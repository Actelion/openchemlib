package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.io.pdb.converter.BondsCalculator;
import com.actelion.research.util.SortedList;

import java.util.*;
import java.util.stream.Collectors;

/**
 * @author JW
 * December 2019
 * The StructureAssembler class takes a list of AtomRecords and constructs the 3D-Molecules. Protein atoms are grouped 
 * together, HETATM records are grouped according to their connectivity into non-bonded fragments, HETATM molecules that are
 * connected to the protein are merged with it. 
 * 
 */


public class StructureAssembler {
	
	public static final String PROTEIN_GROUP = "protein";
	public static final String SOLVENT_GROUP = "water";
	public static final String LIGAND_GROUP = "ligand";
	
	private Map<String,List<AtomRecord>> groups;
	private final SortedList<int[]> bondList;
	private final List<AtomRecord> protAtomRecords;
	private final List<AtomRecord> hetAtomRecords;
	private boolean detachCovalentLigands;
	private Map<String,List<Molecule3D>> mols;


	public StructureAssembler(SortedList<int[]> bondList, List<AtomRecord> protAtomRecords, List<AtomRecord> hetAtomRecords) {
		this.bondList = bondList;
		this.protAtomRecords = protAtomRecords;
		this.hetAtomRecords = hetAtomRecords;
	}

	public void setDetachCovalentLigands(boolean b) {
		detachCovalentLigands = b;
	}

	public Map<String,List<Molecule3D>> assemble() {
		groups = new HashMap<>();
		mols = new HashMap<>();

		group();
		List<Molecule3D> protMols = new ArrayList<>();
		mols.putIfAbsent(SOLVENT_GROUP, new ArrayList<>());
		mols.putIfAbsent(LIGAND_GROUP, new ArrayList<>());
		protMols.add(buildProtein());
		mols.put(PROTEIN_GROUP, protMols);
		buildHetResidues();
		mols.forEach((k,v) -> v.forEach(e -> coupleBonds(e)));
		return mols;
	}
	
	private void group() {
		groups.put(PROTEIN_GROUP, new ArrayList<>(protAtomRecords));
		
		hetAtomRecords.forEach(e -> {
			List<AtomRecord> li = groups.computeIfAbsent(e.getString(), k -> new ArrayList<>());
			li.add(e);
		});

		for(int i=0; i<bondList.size(); i++)
			processBond(bondList.get(i));
	}
	
	private Molecule3D buildProtein() {
		List<AtomRecord> proteinRecords = groups.get(PROTEIN_GROUP);
		Map<String,List<AtomRecord>> residues_ = proteinRecords.stream().collect(Collectors.groupingBy(AtomRecord::getString));
		List<Residue> residues = residues_.values().stream().map(v -> new Residue(v, true, false)).collect(Collectors.toList());
		residues.sort((c1,c2) -> {
				if(!c1.getChainID().equals(c2.getChainID())) //different chains
					return c1.getChainID().compareTo(c2.getChainID());
				else { //same chain
					if(c1.getResnum()!=c2.getResnum())
						return Integer.compare(c1.getResnum(), c2.getResnum());
					else { //same chain, same residue number -> check insertion code
						return c1.getInsertionCode().compareTo(c2.getInsertionCode());
					}
				}
			});

		ProteinSynthesizer proteinSynthesizer = new ProteinSynthesizer();
		List<Molecule3D> protMols = new ArrayList<Molecule3D>();
		for(Residue residue : residues) {
			Molecule3D fragment = residue.getMolecule();
			if(fragment.getAtomAmino(0).trim().equals("ACT") || fragment.getAtomAmino(0).trim().equals("LIG")) {
				mols.get(LIGAND_GROUP).add(fragment);
				continue;
			}
			else if(fragment.getAtomAmino(0).trim().equals("HOH")) {
				mols.get(SOLVENT_GROUP).add(fragment);
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
		if(nextMol!=null && !protMols.contains(nextMol))
				protMols.add(nextMol);
		Molecule3D protein = protMols.stream().reduce((mol1,mol2) ->{
			mol1.addMolecule(mol2);
			return mol1;})
				.get();
		return protein;
		}

	private void buildHetResidues() {
		for(String group : groups.keySet()) {
			if(!group.equals(PROTEIN_GROUP)) {
				List<AtomRecord> records = groups.get(group);
				Residue atomGroup = new Residue(records, false, false);
				Molecule3D fragment = atomGroup.getMolecule();
				if(fragment.getAtomAmino(0).equals("HOH")) {
					mols.putIfAbsent(SOLVENT_GROUP, new ArrayList<Molecule3D>());
					mols.get(SOLVENT_GROUP).add(fragment);
				}
				else {
					mols.putIfAbsent(LIGAND_GROUP, new ArrayList<Molecule3D>());
					mols.get(LIGAND_GROUP).add(fragment);
				}
			}
		}
	}
	
	private void coupleBonds(Molecule3D mol) {
		if (bondList.size() == 0) {
			if (mol.getAllAtoms() > 1) {
				try {
					BondsCalculator.createBonds(mol, true, null);
					BondsCalculator.calculateBondOrders(mol, true);
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

			for(int i=0; i<bondList.size(); i++) {
				int[] bond = bondList.get(i);
				Integer atom1 = sequenceToAtomMap.get(bond[0]);
				Integer atom2 = sequenceToAtomMap.get(bond[1]);
				if (atom1 != null && atom2 != null)
					mol.addBond(atom1, atom2);
			}

			try {
				if (mol.getAllBonds() == 0)	// CONECT records didn't cover this molecule
					BondsCalculator.createBonds(mol, true, null);
				BondsCalculator.calculateBondOrders(mol, true);
			}
			catch (Exception e) {}
		}
	}

	/**
	 * merge atom groups that are connected by a bond
	 * @param bond
	 */
	private void processBond(int[] bond) {
		int atom1 = bond[0];
		int atom2 = bond[1];
		String[] grps = new String[2];
		groups.forEach((k,v) -> {
			List<Integer> atoms = v.stream().map(e -> e.getSerialId()).collect(Collectors.toList());
			if(atoms.contains(atom1))
				grps[0] = k;
			if(atoms.contains(atom2))
				grps[1] = k;
			});

		if(!grps[0].equals(grps[1])) {
			if(grps[0].equals(PROTEIN_GROUP)) {
				if (!detachCovalentLigands) {
					groups.get(grps[0]).addAll(groups.get(grps[1]));
					groups.remove(grps[1]);
				}
			}
			else if(grps[1].equals(PROTEIN_GROUP)) {
				if (!detachCovalentLigands) {
					groups.get(grps[1]).addAll(groups.get(grps[0]));
					groups.remove(grps[0]);
				}
			}
			else {
				groups.get(grps[0]).addAll(groups.get(grps[1]));
				groups.remove(grps[1]);
			}
		}
	}
}
