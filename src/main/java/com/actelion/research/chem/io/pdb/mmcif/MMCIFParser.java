package com.actelion.research.chem.io.pdb.mmcif;

import com.actelion.research.chem.io.pdb.parser.AtomRecord;
import com.actelion.research.chem.io.pdb.parser.PDBCoordEntryFile;
import com.actelion.research.util.IntArrayComparator;
import com.actelion.research.util.SortedList;

import java.io.*;
import java.net.URI;
import java.net.URLConnection;
import java.nio.charset.StandardCharsets;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

public class MMCIFParser {
	public static PDBCoordEntryFile getFromPDB(String pdbID) throws Exception {
		URLConnection con = new URI("https://files.rcsb.org/download/"+pdbID+".cif.gz").toURL().openConnection();
		return MMCIFParser.parse(new BufferedReader(new InputStreamReader(new GZIPInputStream(con.getInputStream()))));
	}

	public static PDBCoordEntryFile parse(String filename) throws IOException {
		return parse(new File(filename));
	}

	public static PDBCoordEntryFile parse(File file) throws IOException {
		InputStream stream = file.getName().toLowerCase().endsWith(".cif.gz")
						  || file.getName().toLowerCase().endsWith(".mmcif.gz") ?
				new GZIPInputStream(new FileInputStream(file)) : new FileInputStream(file);
		return parse(new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8)));
	}

	public static PDBCoordEntryFile parse(BufferedReader reader) throws IOException {
		PDBCoordEntryFile entryFile = new PDBCoordEntryFile();

		String line;
		while ((line = reader.readLine()) != null) {
			if (line.startsWith("data_")) {
				reader.readLine();	// '#'
			}
			else if (line.startsWith("_")) {
				MMCIFBlock block = new MMCIFBlock(line, reader);

				if (block.is("entry"))
					entryFile.setID(block.get("id"));
				if (block.is("citation"))
					entryFile.setTitle(block.get("title"));
				if (block.is("struct_keywords"))
					entryFile.setKeywords(block.get("text"));
				if (block.is("pdbx_database_status"))
					try { entryFile.setDateDeposition(new SimpleDateFormat("yyyy-MM-dd").parse(block.get("recvd_initial_deposition_date"))); } catch (ParseException pe) {}
				if (block.is("refine")) {
//					entryFile.setExperimentalMethod("pdbx_refine_id");
					entryFile.setResolution(block.get("ls_d_res_high"));
				}
			}
			else if (line.startsWith("loop_")) {
				MMCIFTable table = new MMCIFTable(reader);
				if (table.getName().equals("_atom_site"))
					processAtomTable(table, reader, entryFile);
				else if (table.getName().equals("_struct_conn"))
					processConnectionTable(table, reader, entryFile);
				else
					table.skip(reader);
			}
		}

		return entryFile;
	}

	private static void processAtomTable(MMCIFTable table, BufferedReader reader, PDBCoordEntryFile entryFile) throws IOException {
		TreeSet<AtomRecord> proteinAtoms = new TreeSet<>();
		TreeSet<AtomRecord> hetAtoms = new TreeSet<>();

		String[] row;
		while ((row = table.parseRow(reader)) != null) {
			AtomRecord atomRecord = new AtomRecord(
					Integer.parseInt(row[table.getIndex("id")]),
					row[table.getIndex("label_atom_id")],
					row[table.getIndex("label_alt_id")],
					row[table.getIndex("label_comp_id")],
					row[table.getIndex("label_asym_id")],
					parseInt(row[table.getIndex("label_seq_id")]),
					row[table.getIndex("pdbx_PDB_ins_code")],
					Double.parseDouble(row[table.getIndex("Cartn_x")]),
					Double.parseDouble(row[table.getIndex("Cartn_y")]),
					Double.parseDouble(row[table.getIndex("Cartn_z")]),
					Double.parseDouble(row[table.getIndex("occupancy")]),
					Double.parseDouble(row[table.getIndex("B_iso_or_equiv")]),
					row[table.getIndex("type_symbol")]);
			String group = row[table.getIndex("group_PDB")];
			if (group.equals("ATOM"))
				proteinAtoms.add(atomRecord);
			if (group.equals("HETATM"))
				hetAtoms.add(atomRecord);
		}

		ArrayList<AtomRecord> proteinAtomList = new ArrayList<>(proteinAtoms);
		ArrayList<AtomRecord> hetAtomList = new ArrayList<>(hetAtoms);

		entryFile.setProteinAtoms(proteinAtomList);
		entryFile.setHetAtoms(hetAtomList);
		entryFile.setConnections(new SortedList<>(new IntArrayComparator()));
	}

	private static void processConnectionTable(MMCIFTable table, BufferedReader reader, PDBCoordEntryFile entryFile) throws IOException {
		ArrayList<String[]> connectionList = new ArrayList<>();
		String[] row;
		while ((row = table.parseRow(reader)) != null) {
			String connTypeID = row[table.getIndex("conn_type_id")];
			if ("covale".equals(connTypeID)) {
				String[] atomDescription = new String[2];
				atomDescription[0] = atomDescription(
					row[table.getIndex("ptnr1_label_atom_id")],
					row[table.getIndex("ptnr1_label_comp_id")],
					row[table.getIndex("ptnr1_label_seq_id")],
					row[table.getIndex("ptnr1_label_asym_id")]
				);
				atomDescription[1] = atomDescription(
					row[table.getIndex("ptnr2_label_atom_id")],
					row[table.getIndex("ptnr2_label_comp_id")],
					row[table.getIndex("ptnr2_label_seq_id")],
					row[table.getIndex("ptnr2_label_asym_id")]
				);
				connectionList.add(atomDescription);

//			for (String hn : table.getHeaderNames())
//				System.out.println(hn+":"+row[table.getIndex(hn)]);
//			System.out.println("-----");
			}
		}

		if (!connectionList.isEmpty())
			entryFile.setConnections(connectionList);
	}

	public static String atomDescription(String atomID, String compID, String seqID, String asymID) {
		return atomID
			 + (compID.equals(".") || compID.equals("?") ? "_." : "_"+compID)
			 + (seqID.equals(".") || seqID.equals("?") ? "_0" : "_"+seqID)
			 + (asymID.equals(".") || asymID.equals("?") ? "_." : "_"+asymID);
	}

	private static int parseInt(String s) {
		return s.equals(".") ? 0 : Integer.parseInt(s);
	}
}
