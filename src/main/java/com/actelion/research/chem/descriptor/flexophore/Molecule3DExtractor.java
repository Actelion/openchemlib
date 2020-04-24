/*
 * Copyright (c) 2020.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
 */

package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.Molecule3D;

import java.util.*;

public class Molecule3DExtractor {
	
	private Molecule3D molecule3D;
	
	private List<Molecule3D> liFFMoleculeFrag;
	
	private List<int []> liMapFrag2Mol;
	
	public Molecule3DExtractor(Molecule3D molecule3D) {
		
		this.molecule3D = molecule3D;
		
		init();
	}
	
	public Molecule3D getMolecule(int index){
		return liFFMoleculeFrag.get(index);
	}
	
	public int [] getMap(int index){
		return liMapFrag2Mol.get(index);
	}
	
	public int size(){
		return liFFMoleculeFrag.size();
	}

	private void init() {
		List<Integer> liIndexAtom = new ArrayList<Integer>();
		for (int i = 0; i < molecule3D.getAllAtoms(); i++) {
			liIndexAtom.add(i);
		}
		
		List<List<Integer>> liliIndexAtomsFrags = new ArrayList<List<Integer>>(); 
		
		while(!liIndexAtom.isEmpty()) {
			int indexAtomStart = liIndexAtom.remove(0);
			
			HashSet<Integer> hsIndexAtomFrag = extract(molecule3D, indexAtomStart);
			
			liliIndexAtomsFrags.add(new ArrayList<Integer>(hsIndexAtomFrag));
			
			for (int i = liIndexAtom.size()-1; i >= 0; i--) {
				
				if(hsIndexAtomFrag.contains(liIndexAtom.get(i)))
					liIndexAtom.remove(i);
			}
		}
		
		extractWithMapGeneration(molecule3D, liliIndexAtomsFrags);
	}
	
	private static HashSet<Integer> extract(Molecule3D ff, int indexAtomStart) {
		
		HashSet<Integer> hsVisisted = new HashSet<Integer>();
		
		hsVisisted.add(indexAtomStart);
		
		List<Integer> liQueue = new ArrayList<Integer>();
		
		liQueue.add(indexAtomStart);
		
		while(!liQueue.isEmpty()){
			int indexAtom = liQueue.remove(0);
			
			int nConn = ff.getAllConnAtoms(indexAtom);
			
			for (int i = 0; i < nConn; i++) {
				int indexAtomConn = ff.getConnAtom(indexAtom, i);
				
				if(!hsVisisted.contains(indexAtomConn)) {
					hsVisisted.add(indexAtomConn);
					
					liQueue.add(indexAtomConn);
				}
			}
		}
		
		return hsVisisted;
	}
	
	private void extractWithMapGeneration(Molecule3D molecule3D, List<List<Integer>> liliIndexAtomsFrags) {
		
		List<Molecule3D> liFFMolecule = new ArrayList<Molecule3D>();
		
		List<int[]> liMap = new ArrayList<int[]>();
		
		for (int i = 0; i < liliIndexAtomsFrags.size(); i++) {
			
			List<Integer> liIndexAtomsFrag = liliIndexAtomsFrags.get(i);
			
			Molecule3D molFrag = new Molecule3D();
			
			int [] arrMap = new int [liIndexAtomsFrag.size()];
			
			HashMap<Integer, Integer> hmIndexFF_FFFrag = new HashMap<Integer, Integer>();
			
			for (int j = 0; j < liIndexAtomsFrag.size(); j++) {
				
				int indexAtom = liIndexAtomsFrag.get(j);
				
				int indexAtomFrag = molFrag.addAtom(molecule3D, indexAtom);
								
				if(j!=indexAtomFrag){
					throw new RuntimeException("Error in algorithm.");
				}
				arrMap[j]=indexAtom;
				
				hmIndexFF_FFFrag.put(indexAtom, j);
				
			}
			
			for (int j = 0; j < molecule3D.getAllBonds(); j++) {
				int indexAtom1 = molecule3D.getBondAtom(0, j);
				
				if(hmIndexFF_FFFrag.containsKey(indexAtom1)){
					int indexAtom2 = molecule3D.getBondAtom(1, j);
					
					if(hmIndexFF_FFFrag.containsKey(indexAtom2)){
						
						molFrag.addBond(hmIndexFF_FFFrag.get(indexAtom1), hmIndexFF_FFFrag.get(indexAtom2), molecule3D.getBondOrder(j));
					}
				}
			}
			
			liFFMolecule.add(molFrag);
			
			liMap.add(arrMap);
			
			
		}
		
		liFFMoleculeFrag = Collections.unmodifiableList(liFFMolecule);
		
		liMapFrag2Mol = Collections.unmodifiableList(liMap);
		
		
	}
	
	/**
	 * 
	 * @param molecule3D
	 * @param liIndexAtomsFrags the first atom in the index list becomes the first atom in the new molecule and so on.
	 * @return
	 */
	public static Molecule3D extract(Molecule3D molecule3D, List<Integer> liIndexAtomsFrags) {
		
		Molecule3D molFrag = new Molecule3D();
		
		int [] arrMap = new int [liIndexAtomsFrags.size()];
		
		HashMap<Integer, Integer> hmIndexFF_FFFrag = new HashMap<Integer, Integer>();
		
		for (int i = 0; i < liIndexAtomsFrags.size(); i++) {
			
			int indexAtom = liIndexAtomsFrags.get(i);
			
			int indexAtomFrag = molFrag.addAtom(molecule3D, indexAtom);
							
			if(i!=indexAtomFrag){
				throw new RuntimeException("Error in algorithm.");
			}
			arrMap[i]=indexAtom;
			
			hmIndexFF_FFFrag.put(indexAtom, i);
			
		}
		
		for (int i = 0; i < molecule3D.getAllBonds(); i++) {
			int indexAtom1 = molecule3D.getBondAtom(0, i);
			
			if(hmIndexFF_FFFrag.containsKey(indexAtom1)){
				int indexAtom2 = molecule3D.getBondAtom(1, i);
				
				if(hmIndexFF_FFFrag.containsKey(indexAtom2)){
					
					molFrag.addBond(hmIndexFF_FFFrag.get(indexAtom1), hmIndexFF_FFFrag.get(indexAtom2), molecule3D.getBondOrder(i));
				}
			}
		}
		
		return molFrag;
		
	}

}
