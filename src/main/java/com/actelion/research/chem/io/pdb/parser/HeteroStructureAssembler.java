package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.Molecule3D;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Cloned from StructureAssembler (Joel Wahl)
 *
 * Considers only hetero atom records
 *
 * Was implemented to extract covalent bond ligands. i.e. Nirmatrelvir in 8h4y.pdb
 *
 * @author MvK
 * June 2024
 */


public class HeteroStructureAssembler {


	private Map<String,List<AtomRecord>> groups;
	private List<int[]> bondList;

	private List<AtomRecord> hetAtomRecords;
	Map<String,List<Molecule3D>> mols;


	public HeteroStructureAssembler(List<int[]> bondList, List<AtomRecord> hetAtomRecords) {
		this.bondList = bondList;
		this.hetAtomRecords = hetAtomRecords;
		groups = new HashMap<String,List<AtomRecord>>();
		mols = new HashMap<String,List<Molecule3D>>();
	}
	
	
	public Map<String,List<Molecule3D>> assemble() {
		group();

		mols.putIfAbsent(StructureAssembler.SOLVENT_GROUP, new ArrayList<Molecule3D>());
		mols.putIfAbsent(StructureAssembler.LIGAND_GROUP, new ArrayList<Molecule3D>());
		buildHetResidues();
		mols.forEach((k,v) -> v.forEach(e -> coupleBonds(e)));
		return mols;
		
	}
	
	private void group() {

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
	

	private void buildHetResidues() {
		for(String group : groups.keySet()) {
			List<AtomRecord> records = groups.get(group);
			Residue atomGroup = new Residue(records);
			Molecule3D fragment = atomGroup.getMolecule();
			if(fragment.getAtomAmino(0).equals("HOH")) {
				mols.putIfAbsent(StructureAssembler.SOLVENT_GROUP, new ArrayList<Molecule3D>());
				mols.get(StructureAssembler.SOLVENT_GROUP).add(fragment);
			}
			else {
				mols.putIfAbsent(StructureAssembler.LIGAND_GROUP, new ArrayList<Molecule3D>());
				mols.get(StructureAssembler.LIGAND_GROUP).add(fragment);
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
			groups.get(grps[0]).addAll(groups.get(grps[1]));
			groups.remove(grps[1]);
		}
	}
}
