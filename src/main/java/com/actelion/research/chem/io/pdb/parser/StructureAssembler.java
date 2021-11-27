package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.Molecule3D;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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
	private List<int[]> bondList;
	private List<AtomRecord> atomRecords;
	private List<AtomRecord> hetAtomRecords;
	Map<String,List<Molecule3D>> mols;

	
	public StructureAssembler(List<int[]> bondList, List<AtomRecord> atomRecords, List<AtomRecord> hetAtomRecords) {
		this.bondList = bondList;
		this.atomRecords = atomRecords;
		this.hetAtomRecords = hetAtomRecords;
		groups = new HashMap<String,List<AtomRecord>>();
		mols = new HashMap<String,List<Molecule3D>>();
	}
	
	
	public Map<String,List<Molecule3D>> assemble() {
		group();
		List<Molecule3D> protMols = new ArrayList<>();
		mols.putIfAbsent(SOLVENT_GROUP, new ArrayList<Molecule3D>());
		mols.putIfAbsent(LIGAND_GROUP, new ArrayList<Molecule3D>());
		protMols.add(buildProtein());
		mols.put(PROTEIN_GROUP, protMols);
		buildHetResidues();
		mols.forEach((k,v) -> v.forEach(e -> coupleBonds(e)));
		return mols;
		
	}
	
	private void group() {
		groups.put(PROTEIN_GROUP, atomRecords);
		
		hetAtomRecords.forEach(e -> { 
			String s = e.getString();
			if(groups.get(s)!=null) {
				List<AtomRecord> li = groups.get(s);
				li.add(e);
			}
			else { 
				List<AtomRecord> li = new ArrayList<AtomRecord>();
				li.add(e);
				groups.put(s, li);
			}
		});
		for(int[] bond : bondList) {
			try {
				 processBond(bond);
			}
			catch(Exception e) {
				continue;
			}
		}
	}
	
	private Molecule3D buildProtein() {
		ProteinSynthesizer proteinSynthesizer = new ProteinSynthesizer();
		Map<String,List<AtomRecord>> residues_;
		List<AtomRecord> proteinRecords = groups.get(PROTEIN_GROUP);
		residues_ = proteinRecords.stream().collect(Collectors.groupingBy(AtomRecord::getString));
		List<Residue> residues = residues_.values().stream().map(v -> new Residue(v)).collect(Collectors.toList());
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
				boolean coupled = proteinSynthesizer.addResidue(fragment);
				if(coupled) { 
					if(residue.isTerminal()) {
						protMols.add(proteinSynthesizer.retrieveProtein());
						proteinSynthesizer = new ProteinSynthesizer();
					}
					else 
						continue;
				}
				else { //coupling failed
					protMols.add(proteinSynthesizer.retrieveProtein());
					proteinSynthesizer = new ProteinSynthesizer();
					proteinSynthesizer.addResidue(fragment);
						
				}
			}
		Molecule3D nextMol = proteinSynthesizer.retrieveProtein();
		if(nextMol!=null && !protMols.contains(nextMol))
				protMols.add(nextMol);
		Molecule3D protein = protMols.stream().reduce((mol1,mol2) ->{
			mol1.addMolecule(mol2);
			return mol1;})
				.get();
//		protein.ensureHelperArrays(Molecule.cHelperCIP);    // very expensive. Should not be done here just in case somebody might need parities
		return protein;
		}
	
	private void buildHetResidues() {
		for(String group : groups.keySet()) {
			if(group.equals(PROTEIN_GROUP))
				continue;
			else {
				List<AtomRecord> records = groups.get(group);
				Residue atomGroup = new Residue(records);
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
		for(int[] bond:bondList) {
			int [] bondedAtoms = {-1,-1};
			IntStream.range(0,mol.getAllAtoms()).forEach( e -> {
				int pdbAtomID = mol.getAtomSequence(e);
				if(pdbAtomID==bond[0])
					bondedAtoms[0]=e;
				else if(pdbAtomID==bond[1])
					bondedAtoms[1]=e;
			});
			if(bondedAtoms[0]!=-1 && bondedAtoms[1]!=-1)
				mol.addBond(bondedAtoms[0], bondedAtoms[1]);		
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

		if(grps[0].equals(grps[1]))
			return;
		else {
			if(grps[0].equals(PROTEIN_GROUP)) {
				groups.get(grps[0]).addAll(groups.get(grps[1]));
				groups.remove(grps[1]);
				}
			else if(grps[1].equals(PROTEIN_GROUP)) {
				groups.get(grps[1]).addAll(groups.get(grps[0]));
				groups.remove(grps[0]);
				}
			else {
				groups.get(grps[0]).addAll(groups.get(grps[1]));
				groups.remove(grps[1]);
				}
		}
	}
}
