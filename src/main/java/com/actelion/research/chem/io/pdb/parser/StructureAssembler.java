package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.conf.VDWRadii;
import com.actelion.research.chem.io.pdb.converter.BondOrderCalculator;
import com.actelion.research.chem.io.pdb.converter.BondsCalculator;
import com.actelion.research.util.IntArrayComparator;
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
	
	private Map<String,List<AtomRecord>> mGroupNameToAtomListMap;
	private final SortedList<int[]> mTemplateConnectionList,mNonStandardConnectionList;
	private final List<AtomRecord> mProteinAtomList;
	private final List<AtomRecord> mHetAtomList;
	private ArrayList<String> mCovalentLigandGroupList;
	private final boolean mDetachCovalentLigands;
	private Map<String,List<Molecule3D>> mTypeNameToMoleculeListMap;

	/**
	 * Creates a new StructureAssembler ready to create molecules from the given protein- and het-atoms.
	 * @param templateConnections sorted bond atom list of disulfide bonds and non-standard-group bonds plus all or some het-atom bonds
	 * @param nonStandardConnections whether connections include all het-atom bonds (pdb-file) or just some (mmcif-file)
	 * @param proteinAtomRecords
	 * @param hetAtomRecords
	 * @param detachCovalentLigands	whether covalent ligands shall be detached as ligands or be part of the protein
	 */
	public StructureAssembler(SortedList<int[]> templateConnections, SortedList<int[]> nonStandardConnections,
							  List<AtomRecord> proteinAtomRecords,
							  List<AtomRecord> hetAtomRecords, boolean detachCovalentLigands) {
		mTemplateConnectionList = templateConnections;
		mNonStandardConnectionList = nonStandardConnections == null ? new SortedList<>(new IntArrayComparator()) : nonStandardConnections;
		mProteinAtomList = proteinAtomRecords;
		mHetAtomList = hetAtomRecords;
		mDetachCovalentLigands = detachCovalentLigands;
	}

	/**
	 * Does the work to build molecules from given protein- and het-atoms.
	 * @return
	 */
	public Map<String,List<Molecule3D>> assemble() {
		detectMetalLigandBonds();

		// Build map from group names to list(s) of atom records (mGroupAtomRecordListMap).
		buildAtomGroupMap();

		mTypeNameToMoleculeListMap = new HashMap<>();
		mTypeNameToMoleculeListMap.putIfAbsent(SOLVENT_GROUP, new ArrayList<>());
		mTypeNameToMoleculeListMap.putIfAbsent(LIGAND_GROUP, new ArrayList<>());

		// Build protein from atom records and amino acid templates with correct bond orders.
		// For unusual amino acids, where no template exists, PDB/MMCIF files should contain CONECT records for them.
		// If these exist, bonds are created from them, otherwise bonds are created from geometry.
		// Then, bond orders are calculated from geometry and the protein is added to the mTypeNameToMoleculeListMap.
		buildProtein();

		// Ligands and solvents are created from HetAtoms and CONECT records, which result in single bonds
		// unless one of the atoms is a metal atoms, when a zero-oder bond is formed.
		// Then, proper bond orders are calculated from geometry.
		buildHetResidues();

		return mTypeNameToMoleculeListMap;
	}


	private void detectMetalLigandBonds() {
		ArrayList<AtomRecord> metalAtoms = new ArrayList<>();
		ArrayList<AtomRecord> heteroAtoms = new ArrayList<>();
		for (AtomRecord atom : mProteinAtomList) {
			if (Molecule.isAtomicNoMetal(atom.getAtomicNo()))
				metalAtoms.add(atom);
			else if (Molecule.isAtomicNoElectronegative(atom.getAtomicNo()) && !atom.getResName().startsWith("HOH"))
				heteroAtoms.add(atom);
		}
		for (AtomRecord atom : mHetAtomList) {
			if (Molecule.isAtomicNoMetal(atom.getAtomicNo()))
				metalAtoms.add(atom);
			else if (Molecule.isAtomicNoElectronegative(atom.getAtomicNo()))
				heteroAtoms.add(atom);
		}

		for (AtomRecord metal : metalAtoms) {
			for (AtomRecord hetero : heteroAtoms) {
				double limit = 0.5 * VDWRadii.getVDWRadius(metal.getAtomicNo()) + VDWRadii.getVDWRadius(hetero.getAtomicNo());
				double dx = Math.abs(metal.getX() - hetero.getX());
				if (dx < limit) {
					double dy = Math.abs(metal.getY() - hetero.getY());
					if (dy < limit) {
						double dz = Math.abs(metal.getZ() - hetero.getZ());
						if (dz < limit) {
							if (dx*dx+dy*dy+dz*dz < limit*limit) {
								int[] connection = new int[2];
								connection[0] = metal.getSerialId();
								connection[1] = hetero.getSerialId();
								Arrays.sort(connection);
								mNonStandardConnectionList.add(connection);
							}
						}
					}
				}
			}
		}
	}


	private void buildProtein() {
		List<AtomRecord> proteinAtoms = mGroupNameToAtomListMap.get(PROTEIN_GROUP);
		Map<String,List<AtomRecord>> residues_ = proteinAtoms.stream().collect(Collectors.groupingBy(AtomRecord::getString));
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
				mTypeNameToMoleculeListMap.get(LIGAND_GROUP).add(fragment);
				continue;
			}
			else if(fragment.getAtomAmino(0).trim().equals("HOH")) {
				mTypeNameToMoleculeListMap.get(SOLVENT_GROUP).add(fragment);
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
		mTypeNameToMoleculeListMap.put(PROTEIN_GROUP, protMols);
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
			Coordinates coords1 = protein.getAtomCoordinates(atom1);
			for (int j=0; j<i; j++) {
				int atom2 = sulfurList.get(j);
				Coordinates coords2 = protein.getAtomCoordinates(atom2);
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
		for (String groupName : mGroupNameToAtomListMap.keySet()) {
			if (!groupName.equals(PROTEIN_GROUP)) {
				List<AtomRecord> groupOfHetAtoms = mGroupNameToAtomListMap.get(groupName);
				Residue residue = new Residue(groupOfHetAtoms, false, false);
				Molecule3D fragment = residue.getMolecule();
				if (fragment.getAtomAmino(0).equals("HOH")) {
					mTypeNameToMoleculeListMap.putIfAbsent(SOLVENT_GROUP, new ArrayList<>());
					mTypeNameToMoleculeListMap.get(SOLVENT_GROUP).add(fragment);
				} else {
					buildBonds(fragment);
					if (mCovalentLigandGroupList.contains(groupName))
						fragment.setCovalentLigand(true);
					mTypeNameToMoleculeListMap.putIfAbsent(LIGAND_GROUP, new ArrayList<>());
					mTypeNameToMoleculeListMap.get(LIGAND_GROUP).add(fragment);
				}
			}
		}
	}

	private void buildBonds(Molecule3D mol) {
		if (mol.getAllAtoms() > 1) {
			TreeMap<Integer,Integer> sequenceToAtomMap = new TreeMap<>();
			for (int atom=0; atom<mol.getAllAtoms(); atom++) {
				int sequence = mol.getAtomSequence(atom);
				if (sequence != -1)
					sequenceToAtomMap.put(sequence, atom);
			}

			addConnections(mol, mTemplateConnectionList, sequenceToAtomMap);
			try {
				if (mol.getAllBonds() == 0)    // template records didn't cover this molecule
					BondsCalculator.createBonds(mol, true, null);
				addConnections(mol, mNonStandardConnectionList, sequenceToAtomMap);
				new BondOrderCalculator(mol).calculateBondOrders();
			}
			catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	private void addConnections(Molecule3D mol, SortedList<int[]> connections, TreeMap<Integer,Integer> sequenceToAtomMap) {
		if (connections.size() != 0) {
			for(int i = 0; i<connections.size(); i++) {
				int[] bond = connections.get(i);
				Integer atom1 = sequenceToAtomMap.get(bond[0]);
				Integer atom2 = sequenceToAtomMap.get(bond[1]);
				if (atom1 != null && atom2 != null) {
					double dx = mol.getAtomX(atom1)-mol.getAtomX(atom2);
					double dy = mol.getAtomY(atom1)-mol.getAtomY(atom2);
					double dz = mol.getAtomZ(atom1)-mol.getAtomZ(atom2);
					if (dx*dx+dy*dy+dz*dz < 16)
						mol.addBond(atom1, atom2);
				}
			}
		}
	}

	/**
	 * merge atom groups that are connected by a bond
	 */
	private void buildAtomGroupMap() {
		// build map from group names to corresponding atom lists
		mGroupNameToAtomListMap = new HashMap<>();
		mGroupNameToAtomListMap.put(PROTEIN_GROUP, new ArrayList<>(mProteinAtomList));
		mHetAtomList.forEach(hetAtom -> {
			List<AtomRecord> li = mGroupNameToAtomListMap.computeIfAbsent(hetAtom.getString(), k -> new ArrayList<>());
			li.add(hetAtom);
		});
		mCovalentLigandGroupList = new ArrayList<>();

		// Merge atom groups that are connected by a bond
		// unless one of the groups is all water: then move the connected water atom to the other group
		for (int i=0; i<mTemplateConnectionList.size(); i++)
			handleGroupsIfNeeded(mTemplateConnectionList.get(i));
		for (int i=0; i<mNonStandardConnectionList.size(); i++)
			handleGroupsIfNeeded(mNonStandardConnectionList.get(i));
	}

	private void handleGroupsIfNeeded(int[] bond) {
		int atom1ID = bond[0];
		int atom2ID = bond[1];
		AtomRecord atom1 = null;
		AtomRecord atom2 = null;
		String[] grps = new String[2];

		for (String key : mGroupNameToAtomListMap.keySet()) {
			for (AtomRecord atom : mGroupNameToAtomListMap.get(key)) {
				if (atom.getSerialId() == atom1ID) {
					grps[0] = key;
					atom1 = atom;
				}
				if (atom.getSerialId() == atom2ID) {
					grps[1] = key;
					atom2 = atom;
				}
			}
		}

		if(!grps[0].equals(grps[1])) {
			if(grps[0].equals(PROTEIN_GROUP)) {
				if (!mDetachCovalentLigands) {
					mGroupNameToAtomListMap.get(grps[0]).addAll(mGroupNameToAtomListMap.get(grps[1]));
					mGroupNameToAtomListMap.remove(grps[1]);
				}
				else {
					mCovalentLigandGroupList.add(grps[1]);
				}
			}
			else if(grps[1].equals(PROTEIN_GROUP)) {
				if (!mDetachCovalentLigands) {
					mGroupNameToAtomListMap.get(grps[1]).addAll(mGroupNameToAtomListMap.get(grps[0]));
					mGroupNameToAtomListMap.remove(grps[0]);
				}
				else {
					mCovalentLigandGroupList.add(grps[0]);
				}
			}
			else if (grps[0].startsWith("HOH")) {
				mGroupNameToAtomListMap.get(grps[0]).remove(atom1);
				if (mGroupNameToAtomListMap.get(grps[0]).isEmpty())
					mGroupNameToAtomListMap.remove(grps[0]);
				mGroupNameToAtomListMap.get(grps[1]).add(atom1);
			}
			else if (grps[1].startsWith("HOH")) {
				mGroupNameToAtomListMap.get(grps[1]).remove(atom2);
				if (mGroupNameToAtomListMap.get(grps[1]).isEmpty())
					mGroupNameToAtomListMap.remove(grps[1]);
				mGroupNameToAtomListMap.get(grps[0]).add(atom2);
			}
			else {
				mGroupNameToAtomListMap.get(grps[0]).addAll(mGroupNameToAtomListMap.get(grps[1]));
				mGroupNameToAtomListMap.remove(grps[1]);
			}
		}
	}
}
