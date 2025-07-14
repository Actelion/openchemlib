package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.io.pdb.converter.BondsCalculator;
import com.actelion.research.chem.io.pdb.converter.BondOrderCalculator;
import com.actelion.research.util.SortedList;

import java.util.*;
import java.util.stream.Collectors;

/**
 * @author JW,TLS
 * December 2019
 * The StructureAssembler class takes a list of AtomRecords and constructs the 3D-Molecules. Protein atoms are grouped 
 * together, HETATM records are grouped according to their connectivity into non-bonded fragments.
 * HETATM molecules that are connected to the protein are merged with it.
 */

public class StructureAssembler {
	public static final boolean USE_NEW_BONDORDER_CALCULATOR = true;

	public static final String PROTEIN_GROUP = "protein";
	public static final String SOLVENT_GROUP = "water";
	public static final String LIGAND_GROUP = "ligand";
	
	private Map<String,List<AtomRecord>> mGroupAtomRecordListMap;
	private final SortedList<int[]> mConnectionList;
	private final List<AtomRecord> mProteinAtomRecordList;
	private final List<AtomRecord> mHetAtomRecordList;
	private final boolean mDetachCovalentLigands;
	private Map<String,List<Molecule3D>> mMolTypeListMap;

	/**
	 * Creates a new StructureAssembler ready to create molecules from the given protein- and het-atoms.
	 * @param connectionBondList sorted bond atom list of disulfide bonds and non-standard-group bonds including all het-atom bonds
	 * @param proteinAtomRecordList
	 * @param hetAtomRecordList
	 * @param detachCovalentLigands	whether covalent ligands shall be detached as ligands or be part of the protein
	 */
	public StructureAssembler(SortedList<int[]> connectionBondList, List<AtomRecord> proteinAtomRecordList,
							  List<AtomRecord> hetAtomRecordList, boolean detachCovalentLigands) {
		mConnectionList = connectionBondList;
		mProteinAtomRecordList = proteinAtomRecordList;
		mHetAtomRecordList = hetAtomRecordList;
		mDetachCovalentLigands = detachCovalentLigands;
	}

	/**
	 * Does the work to build molecules from given protein- and het-atoms.
	 * @return
	 */
	public Map<String,List<Molecule3D>> assemble() {
		// Build map from group names to list(s) of atom records (mGroupAtomRecordListMap).
		buildAtomGroupMap();

		mMolTypeListMap = new HashMap<>();
		mMolTypeListMap.putIfAbsent(SOLVENT_GROUP, new ArrayList<>());
		mMolTypeListMap.putIfAbsent(LIGAND_GROUP, new ArrayList<>());

		// Build protein from atom records and amino acid templates with correct bond orders.
		// For unusual amino acids, where no template exists, PDB files should contain CONECT records for them.
		// If these exist, bonds are created from them, otherwise bondsare created from geometry.
		// Then, bond orders are calculated from geometry and the protein is added to the mMolTypeListMap.
		buildProtein();

		// Ligands and solvents are created from HetAtoms and CONECT records, which result in single bonds
		// unless one of the atoms is a metal atoms, when a zero-oder bond is formed.
		// Then, proper bond orders are calculated from geometry.
		buildHetResidues();

		return mMolTypeListMap;
	}

	private void buildProtein() {
		List<AtomRecord> proteinRecords = mGroupAtomRecordListMap.get(PROTEIN_GROUP);
		Map<String,List<AtomRecord>> residues_ = proteinRecords.stream().collect(Collectors.groupingBy(AtomRecord::getString));
		List<Residue> residues = residues_.values().stream().map(v -> new Residue(v, true, false)).collect(Collectors.toList());
		residues.sort((c1,c2) -> {
			if (!c1.getChainID().equals(c2.getChainID())) //different chains
				return c1.getChainID().compareTo(c2.getChainID());
			else { //same chain
				if (c1.getResnum() != c2.getResnum())
					return Integer.compare(c1.getResnum(), c2.getResnum());
				else { //same chain, same residue number -> check insertion code
					return c1.getInsertionCode().compareTo(c2.getInsertionCode());
				}
			}
			});

		ProteinSynthesizer proteinSynthesizer = new ProteinSynthesizer();
		List<Molecule3D> protMols = new ArrayList<>();
		for (Residue residue : residues) {
			Molecule3D fragment = residue.getMolecule();
			if(fragment.getAtomAmino(0).trim().equals("ACT") || fragment.getAtomAmino(0).trim().equals("LIG")) {
				mMolTypeListMap.get(LIGAND_GROUP).add(fragment);
				continue;
			}
			else if(fragment.getAtomAmino(0).trim().equals("HOH")) {
				mMolTypeListMap.get(SOLVENT_GROUP).add(fragment);
				continue;
			}

			if(proteinSynthesizer.addResidue(fragment)) {
				if(residue.isTerminal()) {
					protMols.add(proteinSynthesizer.getProtein());
					proteinSynthesizer = new ProteinSynthesizer();
				}
			}
			else { //coupling failed
				protMols.add(proteinSynthesizer.getProtein());
				proteinSynthesizer = new ProteinSynthesizer();
				proteinSynthesizer.addResidue(fragment);
			}
		}
		Molecule3D nextMol = proteinSynthesizer.getProtein();
		if (nextMol!=null && !protMols.contains(nextMol))
			protMols.add(nextMol);

		Molecule3D protein = protMols.stream().reduce((mol1, mol2) -> {
			mol1.addMolecule(mol2);
			return mol1;
			}).get();

		addDisulfideBridges(protein);

		protMols = new ArrayList<>();
		protMols.add(protein);
		mMolTypeListMap.put(PROTEIN_GROUP, protMols);
	}

	/**
	 * Based on distance of all pairs of sulfur atoms
	 * @param protein
	 */
	private void addDisulfideBridges(Molecule3D protein) {
		final double MAX_DISULFIDE_BRIDGE_LENGTH = 3.0;
		protein.ensureHelperArrays(Molecule.cHelperNeighbours);
		ArrayList<Integer> sulfurList = new ArrayList<>();
		for (int atom=0; atom<protein.getAtoms(); atom++)
			if (protein.getAtomicNo(atom) == 16 && protein.getConnAtoms(atom) == 1)
				sulfurList.add(atom);

		for (int i=1; i<sulfurList.size(); i++) {
			int atom1 = sulfurList.get(i);
			Coordinates coords1 = protein.getCoordinates(atom1);
			for (int j=0; j<i; j++) {
				int atom2 = sulfurList.get(j);
				Coordinates coords2 = protein.getCoordinates(atom2);
				double dx = Math.abs(coords2.x - coords1.x);
				if (dx < MAX_DISULFIDE_BRIDGE_LENGTH) {
					double dy = Math.abs(coords2.y - coords1.y);
					if (dy < MAX_DISULFIDE_BRIDGE_LENGTH) {
						double dz = Math.abs(coords2.z - coords1.z);
						if (dz < MAX_DISULFIDE_BRIDGE_LENGTH) {
							if (Math.sqrt(dx*dx + dy*dy + dz*dz) < MAX_DISULFIDE_BRIDGE_LENGTH) {
								protein.addBond(atom1, atom2, Molecule.cBondTypeSingle);
								j = i;	// break loop
							}
						}
					}
				}
			}
		}
	}

	private void buildHetResidues() {
		for(String groupName : mGroupAtomRecordListMap.keySet()) {
			if(!groupName.equals(PROTEIN_GROUP)) {
				List<AtomRecord> atomRecords = mGroupAtomRecordListMap.get(groupName);
				Residue residue = new Residue(atomRecords, false, false);
				Molecule3D fragment = residue.getMolecule();
				if(fragment.getAtomAmino(0).equals("HOH")) {
					mMolTypeListMap.putIfAbsent(SOLVENT_GROUP, new ArrayList<>());
					mMolTypeListMap.get(SOLVENT_GROUP).add(fragment);
				}
				else {
					buildBonds(fragment);
					mMolTypeListMap.putIfAbsent(LIGAND_GROUP, new ArrayList<>());
					mMolTypeListMap.get(LIGAND_GROUP).add(fragment);
				}
			}
		}
	}

	private void buildBonds(Molecule3D mol) {
		if (mConnectionList.size() == 0) {
			if (mol.getAllAtoms() > 1) {
				try {
					if (USE_NEW_BONDORDER_CALCULATOR) {
						BondsCalculator.createBonds(mol, true, null);
						new BondOrderCalculator(mol).calculateBondOrders();
					}
					else {
						BondsCalculator.createBonds(mol, true, null);
						BondsCalculator.calculateBondOrders(mol, true);
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
		else {
			TreeMap<Integer,Integer> sequenceToAtomMap = new TreeMap<>();
			for (int atom=0; atom<mol.getAllAtoms(); atom++) {
				int sequence = mol.getAtomSequence(atom);
				if (sequence != -1)
					sequenceToAtomMap.put(sequence, atom);
			}

			for(int i = 0; i<mConnectionList.size(); i++) {
				int[] bond = mConnectionList.get(i);
				Integer atom1 = sequenceToAtomMap.get(bond[0]);
				Integer atom2 = sequenceToAtomMap.get(bond[1]);
				if (atom1 != null && atom2 != null)
					mol.addBond(atom1, atom2);
			}

			try {
				if (USE_NEW_BONDORDER_CALCULATOR) {
					if (mol.getAllBonds() == 0)    // CONECT records didn't cover this molecule
						BondsCalculator.createBonds(mol, true, null);
					new BondOrderCalculator(mol).calculateBondOrders();
				}
				else {
					if (mol.getAllBonds() == 0)	// CONECT records didn't cover this molecule
						BondsCalculator.createBonds(mol, true, null);
					BondsCalculator.calculateBondOrders(mol, true);
				}
			}
			catch (Exception e) {}
		}
	}

	/**
	 * merge atom groups that are connected by a bond
	 */
	private void buildAtomGroupMap() {
		// build map from group names to corresponding atom lists
		mGroupAtomRecordListMap = new HashMap<>();
		mGroupAtomRecordListMap.put(PROTEIN_GROUP, new ArrayList<>(mProteinAtomRecordList));
		mHetAtomRecordList.forEach(hetAtom -> {
			List<AtomRecord> li = mGroupAtomRecordListMap.computeIfAbsent(hetAtom.getString(), k -> new ArrayList<>());
			li.add(hetAtom);
		});

		// Merge atom groups that are connected by a bond
		for (int i=0; i<mConnectionList.size(); i++) {
			int[] bond = mConnectionList.get(i);
			int atom1 = bond[0];
			int atom2 = bond[1];
			String[] grps = new String[2];
			mGroupAtomRecordListMap.forEach((k, v) -> {
				List<Integer> atoms = v.stream().map(e -> e.getSerialId()).collect(Collectors.toList());
				if(atoms.contains(atom1))
					grps[0] = k;
				if(atoms.contains(atom2))
					grps[1] = k;
				});

			if(!grps[0].equals(grps[1])) {
				if(grps[0].equals(PROTEIN_GROUP)) {
					if (!mDetachCovalentLigands) {
						mGroupAtomRecordListMap.get(grps[0]).addAll(mGroupAtomRecordListMap.get(grps[1]));
						mGroupAtomRecordListMap.remove(grps[1]);
					}
				}
				else if(grps[1].equals(PROTEIN_GROUP)) {
					if (!mDetachCovalentLigands) {
						mGroupAtomRecordListMap.get(grps[1]).addAll(mGroupAtomRecordListMap.get(grps[0]));
						mGroupAtomRecordListMap.remove(grps[0]);
					}
				}
				else {
					mGroupAtomRecordListMap.get(grps[0]).addAll(mGroupAtomRecordListMap.get(grps[1]));
					mGroupAtomRecordListMap.remove(grps[1]);
				}
			}
		}
	}
}
