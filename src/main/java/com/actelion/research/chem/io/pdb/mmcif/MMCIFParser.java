package com.actelion.research.chem.io.pdb.mmcif;

import com.actelion.research.chem.io.pdb.parser.AtomRecord;
import com.actelion.research.chem.io.pdb.parser.PDBFileEntry;
import com.actelion.research.util.IntArrayComparator;
import com.actelion.research.util.SortedList;

import java.io.*;
import java.net.URI;
import java.net.URLConnection;
import java.nio.charset.StandardCharsets;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

public class MMCIFParser {
	private ArrayList<AtomRecord> mAtoms;
	private ArrayList<String[]> mTemplateConnections,mNonStandardConnections;
	private final BufferedReader mReader;

	public static PDBFileEntry getFromPDB(String pdbID) throws Exception {
		URLConnection con = new URI("https://files.rcsb.org/download/"+pdbID+".cif.gz").toURL().openConnection();
		return MMCIFParser.parse(new BufferedReader(new InputStreamReader(new GZIPInputStream(con.getInputStream()))));
	}

	public static PDBFileEntry parse(String filename) throws IOException {
		return parse(new File(filename));
	}

	public static PDBFileEntry parse(File file) throws IOException {
		InputStream stream = file.getName().toLowerCase().endsWith(".cif.gz")
						  || file.getName().toLowerCase().endsWith(".mmcif.gz") ?
				new GZIPInputStream(new FileInputStream(file)) : new FileInputStream(file);
		return parse(new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8)));
	}

	public static PDBFileEntry parse(BufferedReader reader) throws IOException {
		return new MMCIFParser(reader).parse();
	}

	private MMCIFParser(BufferedReader reader) {
		mReader = reader;
	}

	private PDBFileEntry parse() throws IOException {
		// ChimeraX Guidelines: https://www.cgl.ucsf.edu/chimerax/docs/devel/modules/mmcif/mmcif_guidelines.html
		PDBFileEntry entryFile = new PDBFileEntry();
		mTemplateConnections = new ArrayList<>();
		mNonStandardConnections = new ArrayList<>();

		String line;
		while ((line = mReader.readLine()) != null) {
			if (line.startsWith("data_")) {
				mReader.readLine();	// '#'
			}
			else if (line.startsWith("_")) {
				MMCIFBlock block = new MMCIFBlock(line, mReader);

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
				MMCIFTable table = new MMCIFTable(mReader);
				if (table.getName().equals("_atom_site"))
					processAtomTable(table);
				else if (table.getName().equals("_chem_comp_bond"))
					processTemplateConnections(table);
				else if (table.getName().equals("_struct_conn"))
					processNonStandardConnections(table);
				else
					table.skip(mReader);
			}
		}

		entryFile.setAtoms(mAtoms);
		entryFile.setNonStandardConnections(mNonStandardConnections);
		entryFile.setTemplateConnections(translateTemplateConnections());

		return entryFile;
	}

	private void processAtomTable(MMCIFTable table) throws IOException {
		TreeSet<AtomRecord> atoms = new TreeSet<>();

		String[] row;
		while ((row = table.parseRow(mReader)) != null)
			// Used in ChimeraX (†:required): id, label_entity_id, label_asym_id†, auth_asym_id, pdbx_PDB_ins_code,
			// label_seq_id†, auth_seq_id, label_alt_id, type_symbol†, label_atom_id†, auth_atom_id, label_comp_id†,
			// auth_comp_id, Cartn_x†, Cartn_y†, Cartn_z†, occupancy, B_iso_or_equiv, pdbx_PDB_model_num
			atoms.add(new AtomRecord(
					row[table.getIndex("group_PDB")].equals("HETATM"),
					Integer.parseInt(row[table.getIndex("id")]),
					row[table.getIndex("label_atom_id")],
					row[table.getIndex("label_alt_id")],
					parseInt(row[table.getIndex("label_seq_id")]),
					row[table.getIndex("auth_comp_id")],
					row[table.getIndex("auth_asym_id")],
					parseInt(row[table.getIndex("auth_seq_id")]),
					row[table.getIndex("pdbx_PDB_ins_code")],
					Double.parseDouble(row[table.getIndex("Cartn_x")]),
					Double.parseDouble(row[table.getIndex("Cartn_y")]),
					Double.parseDouble(row[table.getIndex("Cartn_z")]),
					Double.parseDouble(row[table.getIndex("occupancy")]),
					Double.parseDouble(row[table.getIndex("B_iso_or_equiv")]),
					row[table.getIndex("type_symbol")]));

		mAtoms = new ArrayList<>(atoms);
	}

	private void processTemplateConnections(MMCIFTable table) throws IOException {
		// Used in ChimeraX (†:required): comp_id†, atom_id_1†, atom_id_2†
		String[] row;
		while ((row = table.parseRow(mReader)) != null) {
//			String connTypeID = row[table.getIndex("conn_type_id")];
			String[] connection = new String[3];
			connection[0] = row[table.getIndex("comp_id")];
			connection[1] = row[table.getIndex("atom_id_1")];
			connection[2] = row[table.getIndex("atom_id_2")];
			mTemplateConnections.add(connection);
		}
	}

	private SortedList<int[]> translateTemplateConnections() {
		SortedList<int[]> connections = new SortedList<>(new IntArrayComparator());
		if (!mTemplateConnections.isEmpty()) {
			Comparator<AtomRecord> comparator = (o1, o2) -> {
				int c = o1.getResName().compareTo(o2.getResName());
				if (c != 0) return c;
				c = o1.getLabelAtomName().compareTo(o2.getLabelAtomName());
				if (c != 0) return c;
				c = Integer.compare(o1.getAuthSeqID(), o2.getAuthSeqID());
				if (c != 0) return c;
				return o1.getChainID().compareTo(o2.getChainID());
			};
			SortedList<AtomRecord> sortedAtoms = new SortedList<>(comparator);
			AtomRecord probe = new AtomRecord();
			for (AtomRecord atom : mAtoms)
				sortedAtoms.add(atom);

			for (String[] connection : mTemplateConnections) {
				probe.setAtomAndCompName(connection[1], connection[0]);
				int index1 = sortedAtoms.getIndexAboveEqual(probe);
				probe.setAtomAndCompName(connection[2], connection[0]);
				int index2 = sortedAtoms.getIndexAboveEqual(probe);
				if (index1 != -1 && index2 != -1) {
					for (int i1=index1; i1<sortedAtoms.size()
							&& connection[0].equals(sortedAtoms.get(i1).getResName())
							&& connection[1].equals(sortedAtoms.get(i1).getLabelAtomName()); i1++) {
						AtomRecord atom1 = sortedAtoms.get(i1);
						for (int i2=index2; i2<sortedAtoms.size()
								&& connection[0].equals(sortedAtoms.get(i2).getResName())
								&& connection[2].equals(sortedAtoms.get(i2).getLabelAtomName()); i2++) {
							AtomRecord atom2 = sortedAtoms.get(i2);
							if (atom1.getChainID().equals(atom2.getChainID()) && atom1.getAuthSeqID()==atom2.getAuthSeqID()) {
								int[] bond = new int[2];
								if (atom1.getSerialId() < atom2.getSerialId()) {
									bond[0] = atom1.getSerialId();
									bond[1] = atom2.getSerialId();
								}
								else {
									bond[0] = atom2.getSerialId();
									bond[1] = atom1.getSerialId();
								}
								connections.add(bond);
							}
						}
					}
				}
			}
		}
		return connections;
	}

	private void processNonStandardConnections(MMCIFTable table) throws IOException {
		// Used in ChimeraX (†:required): conn_type_id†, ptnr1_label_asym_id†, pdbx_ptnr1_PDB_ins_code, ptnr1_label_seq_id†,
		// ptnr1_auth_seq_id, pdbx_ptnr1_label_alt_id, ptnr1_label_atom_id†, ptnr1_label_comp_id†, ptnr1_symmetry,
		// ptnr2_label_asym_id†, pdbx_ptnr2_PDB_ins_code, ptnr2_label_seq_id†, ptnr2_auth_seq_id, pdbx_ptnr2_label_alt_id,
		// ptnr2_label_atom_id†, ptnr2_label_comp_id†, ptnr2_symmetry, pdbx_dist_value
		String[] row;
		while ((row = table.parseRow(mReader)) != null) {
			String[] atomDescription = new String[2];
			atomDescription[0] = atomDescription(
				row[table.getIndex("ptnr1_label_atom_id")],
				row[table.getIndex("ptnr1_label_seq_id")],
				row[table.getIndex("ptnr1_auth_comp_id")],
				row[table.getIndex("ptnr1_auth_seq_id")],
				row[table.getIndex("ptnr1_auth_asym_id")]
			);
			atomDescription[1] = atomDescription(
				row[table.getIndex("ptnr2_label_atom_id")],
				row[table.getIndex("ptnr2_label_seq_id")],
				row[table.getIndex("ptnr2_auth_comp_id")],
				row[table.getIndex("ptnr2_auth_seq_id")],
				row[table.getIndex("ptnr2_auth_asym_id")]
			);
			mNonStandardConnections.add(atomDescription);
		}
	}

	public static String atomDescription(String atomID, String labelSeqID, String authCompID, String authSeqID, String authAsymID) {
		return atomID
			 + (labelSeqID.equals(".") || labelSeqID.equals("?") ? "_0" : "_"+labelSeqID)
			 + (authCompID.equals(".") || authCompID.equals("?") ? "_." : "_"+authCompID)
			 + (authSeqID.equals(".") || authSeqID.equals("?") ? "_0" : "_"+authSeqID)
			 + (authAsymID.equals(".") || authAsymID.equals("?") ? "_." : "_"+authAsymID);
	}

	private int parseInt(String s) {
		return s.equals(".") ? 0 : Integer.parseInt(s);
	}
}
