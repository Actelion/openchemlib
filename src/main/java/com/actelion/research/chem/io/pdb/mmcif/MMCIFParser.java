package com.actelion.research.chem.io.pdb.mmcif;

import com.actelion.research.chem.io.pdb.parser.AtomRecord;
import com.actelion.research.chem.io.pdb.parser.PDBCoordEntryFile;
import com.actelion.research.util.SortedList;

import java.io.*;
import java.net.URI;
import java.net.URLConnection;
import java.nio.charset.StandardCharsets;
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
		TreeSet<AtomRecord> proteinAtoms = new TreeSet<>();
		TreeSet<AtomRecord> hetAtoms = new TreeSet<>();

		String line;
		while ((line = reader.readLine()) != null) {
			line = line.trim();
			if (line.equals("loop_")) {
				MMCIFTable table = new MMCIFTable(reader);
				if (table.getName().equals("_atom_site")) {
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
				}
			}
		}

		ArrayList<AtomRecord> protAtomList = new ArrayList<>(proteinAtoms);
		ArrayList<AtomRecord> hetAtomList = new ArrayList<>(hetAtoms);

		PDBCoordEntryFile entryFile = new PDBCoordEntryFile();
		entryFile.setProtAtomRecords(protAtomList);
		entryFile.setHetAtomRecords(hetAtomList);
		entryFile.setLiConnect(new SortedList<>());
		return entryFile;
	}

	private static int parseInt(String s) {
		return s.equals(".") ? 0 : Integer.parseInt(s);
	}
}
