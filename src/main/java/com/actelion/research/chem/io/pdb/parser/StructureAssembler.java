package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.HydrogenAssembler;
import com.actelion.research.chem.conf.VDWRadii;
import com.actelion.research.chem.io.pdb.calc.BondOrderCalculator;
import com.actelion.research.chem.io.pdb.calc.BondsCalculator;
import com.actelion.research.util.IntArrayComparator;
import com.actelion.research.util.SortedList;

import java.util.*;

/**
 * @author JW,TLS
 * 2019-2025
 * The StructureAssembler class takes a list of AtomRecords and constructs the 3D-Molecules. Protein atoms are grouped 
 * together, HETATM records are grouped according to their connectivity into non-bonded fragments.
 * HETATM molecules that are connected to the protein are merged with it unless detachCovalentLigands==true.
 */

public class StructureAssembler {
	public static final String PROTEIN_GROUP = "protein";
	public static final String SOLVENT_GROUP = "water";
	public static final String LIGAND_GROUP = "ligand";

	private final SortedList<int[]> mTemplateConnectionList,mNonStandardConnectionList;
	private final List<AtomRecord> mAtomList;
	private TreeSet<String> mCovalentLigandGroupSet;
	private final boolean mDetachCovalentLigands;
	private Map<String,List<Molecule3D>> mTypeNameToMoleculeListMap;
	private final Map<Integer,AtomRecord> mSerial2AtomRecordMap;
	private final Map<String,String> mGroupToReassignedGroupMap;

	/**
	 * Creates a new StructureAssembler ready to create molecules from the given protein- and het-atoms.
	 * @param templateConnections sorted bond atom list of disulfide bonds and non-standard-group bonds plus all or some het-atom bonds
	 * @param nonStandardConnections whether connections include all het-atom bonds (pdb-file) or just some (mmcif-file)
	 * @param atomRecords
	 * @param detachCovalentLigands	whether covalent ligands shall be detached as ligands or be part of the protein
	 */
	public StructureAssembler(SortedList<int[]> templateConnections, SortedList<int[]> nonStandardConnections,
							  List<AtomRecord> atomRecords, boolean detachCovalentLigands) {
		mTemplateConnectionList = templateConnections;
		mNonStandardConnectionList = nonStandardConnections == null ? new SortedList<>(new IntArrayComparator()) : nonStandardConnections;
		mAtomList = atomRecords;
		mDetachCovalentLigands = detachCovalentLigands;
		mGroupToReassignedGroupMap = new TreeMap<>();
		mSerial2AtomRecordMap = new TreeMap<>();
		for (AtomRecord atom : atomRecords)
			mSerial2AtomRecordMap.put(atom.getSerialId(), atom);
	}

	/**
	 * Does the work to build molecules from given protein- and het-atoms.
	 * @return
	 */
	public Map<String,List<Molecule3D>> assemble() {
		detectMetalLigandAndDisulfideBonds();

		// Assign atom groups to other groups if we have custom or template bonds connecting them.
		mergeAtomGroupsByBonds();

		mTypeNameToMoleculeListMap = new HashMap<>();
		mTypeNameToMoleculeListMap.putIfAbsent(SOLVENT_GROUP, new ArrayList<>());
		mTypeNameToMoleculeListMap.putIfAbsent(LIGAND_GROUP, new ArrayList<>());

		TreeMap<String,Residue> proteinResidueMap = new TreeMap<>();
		TreeMap<String,Residue> hetAtomResidueMap = new TreeMap<>();
		for (AtomRecord atom : mAtomList) {
			String residueKey = getAssignedGroup(atom);
			if (PROTEIN_GROUP.equals(residueKey)) {
				String residueName = atom.getString();
				Residue residue = proteinResidueMap.get(residueName);
				if (residue == null)
					proteinResidueMap.put(residueName, new Residue(true, true, false, atom));
				else
					residue.addAtom(atom);
			}
			else {
				Residue residue = hetAtomResidueMap.get(residueKey);
				if (residue == null)
					hetAtomResidueMap.put(residueKey, new Residue(false, false, false, atom));
				else
					residue.addAtom(atom);
			}
		}
		for (Residue residue : proteinResidueMap.values())
			residue.build();
		for (Residue residue : hetAtomResidueMap.values())
			residue.build();

		// Build protein from atom records and amino acid templates with correct bond orders.
		// For unusual amino acids, where no template exists, PDB/MMCIF files should contain CONECT records for them.
		// If these exist, bonds are created from them, otherwise bonds are created from geometry.
		// Then, bond orders are calculated from geometry and the protein is added to the mTypeNameToMoleculeListMap.
		buildProtein(proteinResidueMap);

		// Ligands and solvents are created from HetAtoms and CONECT records, which result in single bonds
		// unless one of the atoms is a metal atoms, when a zero-oder bond is formed.
		// Then, proper bond orders are calculated from geometry.
		buildHetResidues(hetAtomResidueMap);

		return mTypeNameToMoleculeListMap;
	}


	private void detectMetalLigandAndDisulfideBonds() {
		ArrayList<AtomRecord> metalAtoms = new ArrayList<>();
		ArrayList<AtomRecord> heteroAtoms = new ArrayList<>();
		ArrayList<AtomRecord> sulfurAtoms = new ArrayList<>();
		for (AtomRecord atom : mAtomList) {
			if (Molecule.isAtomicNoMetal(atom.getAtomicNo()))
				metalAtoms.add(atom);
			else if (Molecule.isAtomicNoElectronegative(atom.getAtomicNo()))
//			 && (atom.isHetAtom() || !atom.getResName().startsWith("HOH")))
				heteroAtoms.add(atom);
			if (atom.getAtomicNo() == 16)
				sulfurAtoms.add(atom);
		}

		for (AtomRecord metal : metalAtoms) {
			for (AtomRecord hetero : heteroAtoms) {
				double limit = 0.5 * VDWRadii.getVDWRadius(metal.getAtomicNo()) + VDWRadii.getVDWRadius(hetero.getAtomicNo());
				addConnectionIfDistanceSmallerThan(metal, hetero, limit);
			}
		}
		for (int i=1; i<sulfurAtoms.size(); i++)
			if (!sulfurAtoms.get(i).isHetAtom())
				for (int j=0; j<i; j++)
					if (!sulfurAtoms.get(j).isHetAtom())
						addConnectionIfDistanceSmallerThan(sulfurAtoms.get(j), sulfurAtoms.get(i), 3.0);
	}

	private void addConnectionIfDistanceSmallerThan(AtomRecord atom1, AtomRecord atom2, double distance) {
		double dx = Math.abs(atom1.getX() - atom2.getX());
		if (dx < distance) {
			double dy = Math.abs(atom1.getY() - atom2.getY());
			if (dy < distance) {
				double dz = Math.abs(atom1.getZ() - atom2.getZ());
				if (dz < distance) {
					if (dx*dx+dy*dy+dz*dz < distance*distance) {
						int[] connection = new int[2];
						connection[0] = atom1.getSerialId();
						connection[1] = atom2.getSerialId();
						Arrays.sort(connection);
						mNonStandardConnectionList.add(connection);
					}
				}
			}
		}
	}

	private void buildProtein(TreeMap<String,Residue> proteinResidueMap) {
		Residue[] residues = proteinResidueMap.values().toArray(new Residue[0]);
		Arrays.sort(residues, (c1,c2) -> {
			int comparison = c1.getChainID().compareTo(c2.getChainID());
			if (comparison != 0) //different chains
				return comparison;
			comparison = Integer.compare(c1.getResnum(), c2.getResnum());
			if (comparison != 0)
				return comparison;
			return c1.getInsertionCode().compareTo(c2.getInsertionCode());
		} );

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

		buildBonds(protein, true);

		protMols = new ArrayList<>();
		protMols.add(protein);
		mTypeNameToMoleculeListMap.put(PROTEIN_GROUP, protMols);
	}

	private void buildHetResidues(TreeMap<String,Residue> hetAtomResidueMap) {
		for (String groupName : hetAtomResidueMap.keySet()) {
			Molecule3D fragment = hetAtomResidueMap.get(groupName).getMolecule();
			if (fragment.getAtomAmino(0).equals("HOH")) {
				mTypeNameToMoleculeListMap.putIfAbsent(SOLVENT_GROUP, new ArrayList<>());
				mTypeNameToMoleculeListMap.get(SOLVENT_GROUP).add(fragment);
			} else {
				buildBonds(fragment, false);
				if (mCovalentLigandGroupSet.contains(groupName))
					fragment.setCovalentLigand(true);
				mTypeNameToMoleculeListMap.putIfAbsent(LIGAND_GROUP, new ArrayList<>());
				mTypeNameToMoleculeListMap.get(LIGAND_GROUP).add(fragment);
			}
		}
	}

	private void buildBonds(Molecule3D mol, boolean isProtein) {
		if (mol.getAllAtoms() > 1) {
			TreeMap<Integer,Integer> sequenceToAtomMap = new TreeMap<>();
			for (int atom=0; atom<mol.getAllAtoms(); atom++) {
				int sequence = mol.getAtomSequence(atom);
				if (sequence != -1)
					sequenceToAtomMap.put(sequence, atom);
			}
			// make sure not to duplicate bonds
			TreeSet<Long> existingBondSet = new TreeSet<>();
			for (int bond=0; bond<mol.getAllBonds(); bond++)
				 existingBondSet.add(mol.getBondAtom(0, bond) < mol.getBondAtom(1, bond) ?
							   ((long)mol.getBondAtom(0, bond) << 32) + mol.getBondAtom(1, bond)
							 : ((long)mol.getBondAtom(1, bond) << 32) + mol.getBondAtom(0, bond));

			if (!isProtein)
				addConnections(mol, mTemplateConnectionList, sequenceToAtomMap, existingBondSet);
			try {
				if (!isProtein && mol.getAllBonds() < mol.getAllAtoms()) {    // template records didn't properly cover this molecule
					BondsCalculator.createBonds(mol, true, null);
					for (int bond=0; bond<mol.getAllBonds(); bond++)
						existingBondSet.add(mol.getBondAtom(0, bond) < mol.getBondAtom(1, bond) ?
								  ((long)mol.getBondAtom(0, bond) << 32) + mol.getBondAtom(1, bond)
								: ((long)mol.getBondAtom(1, bond) << 32) + mol.getBondAtom(0, bond));
				}
				addConnections(mol, mNonStandardConnectionList, sequenceToAtomMap, existingBondSet);
				if (!isProtein)
					new BondOrderCalculator(mol).calculateBondOrders();
			}
			catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	private void addConnections(Molecule3D mol, SortedList<int[]> connections,
								TreeMap<Integer,Integer> sequenceToAtomMap,
								TreeSet<Long> existingBondSet) {
		if (connections.size() != 0) {
			ArrayList<Integer> newBondList = new ArrayList<>();
			for(int i = 0; i<connections.size(); i++) {
				int[] bond = connections.get(i);
				Integer atom1 = sequenceToAtomMap.get(bond[0]);
				Integer atom2 = sequenceToAtomMap.get(bond[1]);
				if (atom1 != null && atom2 != null) {
					long bondID = (atom1 < atom2) ? ((long)atom1 << 32) + atom2 : ((long)atom2 << 32) + atom1;
					if (!existingBondSet.contains(bondID)) {
						existingBondSet.add(bondID);
						double dx = mol.getAtomX(atom1)-mol.getAtomX(atom2);
						if (dx < 4) {
							double dy = mol.getAtomY(atom1)-mol.getAtomY(atom2);
							if (dy < 4) {
								double dz = mol.getAtomZ(atom1)-mol.getAtomZ(atom2);
								if (dz < 4 && dx*dx+dy*dy+dz*dz < 16) {
									int type = mol.isMetalAtom(atom1) || mol.isMetalAtom(atom2) ? Molecule.cBondTypeMetalLigand : Molecule.cBondTypeSingle;
									newBondList.add(mol.addBondNoChecks(atom1, atom2, type));
								}
							}
						}
					}
				}
			}
			ArrayList<Integer> chargeAdaptionList = new ArrayList<>();
			for (int bond : newBondList) {
				// if we add a metal at a position where it would collide with an implicit hydrogen atom,
				// then add a negative charge to prevent the HydrogenAssembler to add a colliding hydrogen!
				int atom1 = mol.getBondAtom(0, bond);
				int atom2 = mol.getBondAtom(1, bond);
				if (mol.isMetalAtom(atom1) && causesHydrogenCollision(mol, atom2, atom1))
					chargeAdaptionList.add(atom2);
				else if (mol.isMetalAtom(atom2) && causesHydrogenCollision(mol, atom1, atom2))
					chargeAdaptionList.add(atom1);
			}
			// we use this delayed method of charge changes to prevent new calculation of helper arrays with every change
			for (int atom : chargeAdaptionList) {
				if (mol.isElectropositive(atom))
					mol.setAtomCharge(atom, mol.getAtomCharge(atom)+1);
				else
					mol.setAtomCharge(atom, mol.getAtomCharge(atom)-1);
			}
		}
	}

	private boolean causesHydrogenCollision(StereoMolecule mol, int rootAtom, int metalAtom) {
		if (mol.getImplicitHydrogens(rootAtom) > 0) {
			mol.ensureHelperArrays(Molecule.cHelperRings);
			ArrayList<Coordinates> coordsList = new HydrogenAssembler(mol).getImplicitHydrogenPositions(rootAtom);
			for (Coordinates ch : coordsList) {
				Coordinates vm = mol.getAtomCoordinates(rootAtom).subC(mol.getAtomCoordinates(metalAtom));
				Coordinates vh = mol.getAtomCoordinates(rootAtom).subC(ch);
				return vm.getAngle(vh) < Math.PI / 5;
			}
		}
		return false;
	}

	/**
	 * merge atom groups that are connected by a bond
	 */
	private void mergeAtomGroupsByBonds() {
		mCovalentLigandGroupSet = new TreeSet<>();

		// Merge atom groups that are connected by a bond
		// unless one of the groups is all water: then move the connected water atom to the other group
		for (int i=0; i<mTemplateConnectionList.size(); i++)
			mergeConnectedGroups(mTemplateConnectionList.get(i));
		for (int i=0; i<mNonStandardConnectionList.size(); i++)
			mergeConnectedGroups(mNonStandardConnectionList.get(i));
	}

	private void mergeConnectedGroups(int[] bond) {
		AtomRecord[] bondAtom = new AtomRecord[2];
		boolean[] isMetal = new boolean[2];
		String[] grps = new String[2];

		for (int i=0; i<2; i++) {
			bondAtom[i] = mSerial2AtomRecordMap.get(bond[i]);
			grps[i] = getAssignedGroup(bondAtom[i]);
			isMetal[i] = Molecule.isAtomicNoMetal(bondAtom[i].getAtomicNo());
		}

		if (!grps[0].equals(grps[1])) {
			for (int i=0; i<2; i++) {
				if(grps[i].equals(PROTEIN_GROUP)) {
					if (!mDetachCovalentLigands || isMetal[1-i]) {
						reassignGroup(grps[1-i], grps[i]);
					}
					else {
						mCovalentLigandGroupSet.add(grps[1-i]);
						if (!isMetal[i]) {
							bondAtom[i].setCovalentBridgeAtom(bondAtom[1-i]);
							bondAtom[1-i].setCovalentBridgeAtom(bondAtom[i]);
						}
					}
					return;
				}
			}

			for (int i=0; i<2; i++) {
				if (grps[i].startsWith("HOH")) {
					reassignGroup(grps[i], grps[1-i]);
					return;
				}
			}

			if (!isMetal[0] && !isMetal[1]) {
				reassignGroup(grps[1], grps[0]);
				if (mCovalentLigandGroupSet.contains(grps[1])) {
					mCovalentLigandGroupSet.remove(grps[1]);
					mCovalentLigandGroupSet.add(grps[0]);
				}
			}
		}
	}

	private void reassignGroup(String oldGroup, String newGroup) {
		for (String key : mGroupToReassignedGroupMap.keySet())
			if (oldGroup.equalsIgnoreCase(mGroupToReassignedGroupMap.get(key)))
				mGroupToReassignedGroupMap.put(key, newGroup);

		mGroupToReassignedGroupMap.put(oldGroup, newGroup);
	}

	private String getAssignedGroup(AtomRecord atom) {
		if (!atom.isHetAtom())
			return PROTEIN_GROUP;

		String group = atom.getString();
		String reassignedGroup = mGroupToReassignedGroupMap.get(group);
		return reassignedGroup == null ? group : reassignedGroup;
	}
}
