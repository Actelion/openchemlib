package com.actelion.research.chem.mcs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import com.actelion.research.chem.StereoMolecule;

/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/
public class SubStructSearchExhaustiveTreeWalker {
	
	private static final boolean VERBOSE = false;
	
	private static final boolean OVERLAP = false;
	
	private static final boolean precheckTotalFormula = false;
	
	private static final int CARDINALITY_RING_ATOM = 16;
	
	private static final int CARDINALITY_HETERO_ATOM = 32;
	
	private static final int ELEMENTS = 83;
	
	private static final int MAX_NUM_ATOMS = 256;
	
	static final int LEN_RING_CLOSURES = 8;
	
	static final int LEN_MAPPING_BLOCK = 8;
	
	private StereoMolecule mol;
	
	private StereoMolecule frag;
	
	private int [] arrTotalFormularMol;
	
	private int [] arrTotalFormularFrag;
	
	private int [] arrCardinalityFragment;
	
	private boolean [] arrVisited;
	
	private int [] arrMapMolecule;
	
	private List<int []> liMatches; 

	private PossibleMappingsFrag2Mol [] arrProcessingArray;
	
	private int [] arrFragIndexConnected;
	
	private int [] arrMapFrag2ProcessingArray;

	
	private boolean maxCapacityMatchListContainerReached;
	
	private LinkedList<PossibleMappingsFrag2Mol> liQueue;
	
	private List<Integer> liFragmentAtomNeighbours;
	
	private final int maxCapacity;
	
	public SubStructSearchExhaustiveTreeWalker() {
		
		arrTotalFormularMol = new int [ELEMENTS];
		
		arrTotalFormularFrag = new int [ELEMENTS];
		
		arrCardinalityFragment = new int [MAX_NUM_ATOMS];
		
		arrFragIndexConnected = new int [MAX_NUM_ATOMS];
		
		arrMapFrag2ProcessingArray = new int [MAX_NUM_ATOMS];
		
		arrMapMolecule = new int [MAX_NUM_ATOMS];
		
		arrVisited = new boolean [MAX_NUM_ATOMS];
		
		arrProcessingArray = new PossibleMappingsFrag2Mol[MAX_NUM_ATOMS];
		
		liQueue = new LinkedList<PossibleMappingsFrag2Mol>();
		
		liFragmentAtomNeighbours = new ArrayList<Integer>();
		
		maxCapacity = MAX_NUM_ATOMS/4;
		
		initProcessingArray();
		
	}
	
	public void setMolecule(StereoMolecule mol){
		
		this.mol=mol;
		
	}
	
	public void setFragment(StereoMolecule frag){
		
		this.frag=frag;
		
		calculateCardinalityFragment();
	}
	
	public int findFragmentInMolecule() {
		
		if(!isFulfillPrechecks()) {
			return 0;
		}
		
		maxCapacityMatchListContainerReached=false;
		
		final int atomsFragment = frag.getAtoms();
		
		final int atomsMolecule = mol.getAtoms();
		
		Arrays.fill(arrMapMolecule, 0, atomsMolecule, -1);
		
		int [] arrMapFragment = new int [atomsFragment];
		
		liMatches = new ArrayList<int[]>();
		
		Arrays.fill(arrMapFragment, -1);
		
		List<Integer> liSortedMoleculeAtomNeighbours = new ArrayList<Integer>();
	
		
		determineProcessingArray();
		
		boolean openEnds=true;
		
		int indexProcessingArray=0;
		
		searchcompleted:
		maxCapacityMatchListContainerReached:
		while(openEnds){
			
			PossibleMappingsFrag2Mol possMapFrag2MolCurrent = arrProcessingArray[indexProcessingArray];
			
			while(possMapFrag2MolCurrent.empty()){
				
				if(VERBOSE) {
					System.err.println("possMapFrag2MolCurrent empty");
					System.out.println();
				}
				
				indexProcessingArray--;
				
				if(indexProcessingArray<0){
					openEnds=false;
					break searchcompleted;
				} else {
					
					possMapFrag2MolCurrent = arrProcessingArray[indexProcessingArray];
					
					final int indexAtomFragment = arrProcessingArray[indexProcessingArray].getIndexAtomFrag();
					
					final int indexMappedAtomMolecule = arrMapFragment[indexAtomFragment];
					
					arrMapFragment[indexAtomFragment]=-1;
					
					arrMapMolecule[indexMappedAtomMolecule]=-1;
					
				}
			}
			
			final int indexMappedAtomMolecule = possMapFrag2MolCurrent.pollIndexMappingAtomMolecule();
				
			final int indexAtomFragment = possMapFrag2MolCurrent.getIndexAtomFrag();
			
			arrMapFragment[indexAtomFragment]=indexMappedAtomMolecule;
				
			arrMapMolecule[indexMappedAtomMolecule]=indexAtomFragment;
				
			
			if(VERBOSE) {
				System.out.println("indexProcessingArray " + indexProcessingArray + " indexAtomFragment " + indexAtomFragment + " indexMappedAtomMolecule " + indexMappedAtomMolecule);
			}
			
			indexProcessingArray++;
			
				
			// One match completed.
			// Each fragment atom is matched to a molecule atom.
			if(indexProcessingArray==atomsFragment){
				
				if(liMatches.size()==maxCapacity){
					
					maxCapacityMatchListContainerReached = true;
					
					break maxCapacityMatchListContainerReached;
				}
				
				liMatches.add(arrMapFragment);
				
				int [] arrMapFragmentNew = new int [atomsFragment];
				
				if(OVERLAP) {
					
					System.arraycopy(arrMapFragment, 0, arrMapFragmentNew, 0, indexProcessingArray);
					
					arrMapFragment=arrMapFragmentNew;
					
					indexProcessingArray--;
					
					// Last entry in map is set to -1.
					arrMapFragment[indexAtomFragment]=-1;
					
					arrMapMolecule[indexMappedAtomMolecule]=-1;
				
				} else {
					indexProcessingArray=0;
					
					arrMapFragmentNew[0]=arrMapFragment[0];
					
					for (int i = 1; i < atomsFragment; i++) {
						
						arrMapFragmentNew[i]=-1;
						
						arrProcessingArray[i].resetMappingMolecules();
					}
					
					arrMapFragment=arrMapFragmentNew;
				}
				
			} else { // Prepare next PossibleMappingsFrag2Mol in arrProcessingArray
				// Find all available neihgbours for the parent molecule atom
				
				final PossibleMappingsFrag2Mol possMapFrag2MolNext = arrProcessingArray[indexProcessingArray];
				
				// System.out.println("Next " + possMapFrag2MolNext.getIndexAtomFrag());
				
				getAvailableMoleculeAtomNeighbours(
						possMapFrag2MolNext, 
						arrMapFragment, 
						arrMapMolecule, 
						liSortedMoleculeAtomNeighbours);
				
				for (int indexAtomMoleculeAvailableNeighbour : liSortedMoleculeAtomNeighbours) {
					
					possMapFrag2MolNext.addIndexMappingAtomMolecule(indexAtomMoleculeAvailableNeighbour);
				}
			}
		}
		
		return liMatches.size();
	}
	
	public List<int []> getMatches(){
		return liMatches;
	}
	
	/**
	 * Generates a tree from the fragment
	 * @param
	 * @return an array with one field for each fragment atom. 
	 */
	private void determineProcessingArray(){
		
		for (int i = 0; i < frag.getAtoms(); i++) {
			arrProcessingArray[i].reset();
			
			arrVisited[i]=false;
		}
		
		determinePossibleMappingsFrag2MolStart();
		
		PossibleMappingsFrag2Mol possMapFrag2MolStart = arrProcessingArray[0];
		
		int indexProcessingArray = 0;
		
		indexProcessingArray++;
		
		arrVisited[possMapFrag2MolStart.getIndexAtomFrag()]=true;
		
		liQueue.clear();
		
		liQueue.add(possMapFrag2MolStart);
		
		liFragmentAtomNeighbours.clear();
		
		while(!liQueue.isEmpty()){
			PossibleMappingsFrag2Mol possibleMappingsFrag2MolParent = liQueue.pollLast();
			
			final int indexAtomFragParent = possibleMappingsFrag2MolParent.getIndexAtomFrag();
			
			getFragmentAtomNeighbours(indexAtomFragParent, liFragmentAtomNeighbours);
			
			for (int indexAtomFragNeighbour : liFragmentAtomNeighbours) {
				
				if(!arrVisited[indexAtomFragNeighbour]){
					
					arrVisited[indexAtomFragNeighbour]=true;
					
					PossibleMappingsFrag2Mol possibleMappingsFrag2Mol = arrProcessingArray[indexProcessingArray];
						
					possibleMappingsFrag2Mol.setIndexAtomFrag(indexAtomFragNeighbour);
						
					possibleMappingsFrag2MolParent.addChild(possibleMappingsFrag2Mol);
					
					indexProcessingArray++;
					
					liQueue.push(possibleMappingsFrag2Mol);
				}
			}
		}
		
		addRingClosures(arrProcessingArray);
		
	}
	
	/**
	 * Ring closures are marked only in the PossibleMappingsFrag2Mol object with the higer index in arrProcessingArray.
	 * @param arrProcessingArray
	 */
	private void addRingClosures(PossibleMappingsFrag2Mol [] arrProcessingArray){
		
		final int atomsFrag = frag.getAtoms();
		
		
		for (int i = 0; i < atomsFrag; i++) {
			arrMapFrag2ProcessingArray[arrProcessingArray[i].getIndexAtomFrag()]=i;
		}
		
		for (int i = atomsFrag-1; i >= 0; i--) {
			final PossibleMappingsFrag2Mol possibleMappingsFrag2Mol = arrProcessingArray[i];
			
			final int indexAtomFragment = possibleMappingsFrag2Mol.getIndexAtomFrag();
			
			int nNeighboursInTree = 0;
			
			if(possibleMappingsFrag2Mol.getParent()!=null){
				nNeighboursInTree++;
			}
			
			final int childs = possibleMappingsFrag2Mol.getChilds(); 
			
			nNeighboursInTree += childs; 
			
			final int nConnected = frag.getConnAtoms(indexAtomFragment);
			
			// Ring closure?
			if(nNeighboursInTree!=nConnected) {
			
				Arrays.fill(arrFragIndexConnected, 0, atomsFrag, 0);
				
				if(possibleMappingsFrag2Mol.getParent()!=null){
					arrFragIndexConnected[possibleMappingsFrag2Mol.getParent().getIndexAtomFrag()]=1;
				}
				
				for (int j = 0; j < childs; j++) {
					
					final PossibleMappingsFrag2Mol possMapFrag2MolChild = possibleMappingsFrag2Mol.getChild(j);
					
					arrFragIndexConnected[possMapFrag2MolChild.getIndexAtomFrag()]=1;
					
				}
				
				for (int j = 0; j < nConnected; j++) {
					final int indexAtomFragConnected = frag.getConnAtom(indexAtomFragment, j);
					
					if(arrFragIndexConnected[indexAtomFragConnected]==0){
						
						final int positionInProcessingArrayConnectedFragmentAtom = arrMapFrag2ProcessingArray[indexAtomFragConnected];
						
						final int positionInProcessingArrayFragmentAtom = possibleMappingsFrag2Mol.getIndexProcessingArray();
						
						if( positionInProcessingArrayConnectedFragmentAtom < positionInProcessingArrayFragmentAtom) {
							
							possibleMappingsFrag2Mol.addIndexCounterAtomFragmentRingClosure(indexAtomFragConnected);
						}
					}
				}
			}
		}
		
	}
	
	/**
	 * Sorting for cardinality needed more time than without. 
	 * @param indexAtomFragment
	 * @param liSortedFragmentAtomNeighbours
	 */
	private void getFragmentAtomNeighbours(int indexAtomFragment, List<Integer> liSortedFragmentAtomNeighbours){
		
		liSortedFragmentAtomNeighbours.clear();
		
		int nConnectedFrag = frag.getConnAtoms(indexAtomFragment);
		
		for (int i = 0; i < nConnectedFrag; i++) {
			
			int indexConnAtomFrag = frag.getConnAtom(indexAtomFragment, i);
			
			liSortedFragmentAtomNeighbours.add(indexConnAtomFrag);
		}
		
	}
	
	/**
	 * Sorting of neighbours according cardinality slows performance down.
	 * @param possMapFrag2MolNext
	 * @param arrMapFragment
	 * @param arrMapMolecule
	 * @param liSortedMoleculeAtomNeighbours
	 */
	private void getAvailableMoleculeAtomNeighbours(
			PossibleMappingsFrag2Mol possMapFrag2MolNext, 
			int [] arrMapFragment, 
			int [] arrMapMolecule, 
			List<Integer> liSortedMoleculeAtomNeighbours){
		
		final PossibleMappingsFrag2Mol possMapFrag2MolParent = possMapFrag2MolNext.getParent();
		
		final int indexAtomMoleculeParent = arrMapFragment[possMapFrag2MolParent.getIndexAtomFrag()];

		final int indexAtomFragment = possMapFrag2MolNext.getIndexAtomFrag();
		
		final int indexAtomFragmentParent = possMapFrag2MolParent.getIndexAtomFrag();
		
		
		liSortedMoleculeAtomNeighbours.clear();
		
		final int nConnectedMol = mol.getConnAtoms(indexAtomMoleculeParent);
		
		final int indexBondFrag = frag.getBond(indexAtomFragmentParent, indexAtomFragment);
		
		for (int i = 0; i < nConnectedMol; i++) {
			
			final int indexConnAtomMol = mol.getConnAtom(indexAtomMoleculeParent, i);
			
			if(arrMapMolecule[indexConnAtomMol]!=-1) {
				continue;
			}
			
			
			if(areAtomsSimilar(indexConnAtomMol, indexAtomFragment)) {
				final int indexBondMol = mol.getBond(indexAtomMoleculeParent, indexConnAtomMol);
				
				if(areBondsSimilar(indexBondMol, indexBondFrag)) {
					boolean ringsOk = true;
					if(possMapFrag2MolNext.getRingClosures()>0) {
						
						if(!isRingClosureBondMatch(possMapFrag2MolNext, indexConnAtomMol, arrMapFragment)){
							ringsOk = false;
						}
					}
					
					if(ringsOk){
						liSortedMoleculeAtomNeighbours.add(indexConnAtomMol);
					}
				}
			}
		}
	}
	
	private boolean isRingClosureBondMatch(PossibleMappingsFrag2Mol possMapFrag2MolNext, int indexConnAtomMol, int [] arrMapFragment){
		
		boolean ringClosuresMatch=true;
		
		final int ringClosures = possMapFrag2MolNext.getRingClosures();
		
		final int indexAtomFragment = possMapFrag2MolNext.getIndexAtomFrag();
		
		for (int i = 0; i < ringClosures; i++) {
			
			final int indexCounterAtomFragmentRingClosure = possMapFrag2MolNext.getIndexCounterAtomFragmentRingClosure(i);
			
			final int indexCounterAtomMoleculeRingClosure = arrMapFragment[indexCounterAtomFragmentRingClosure];
			
			
			final int indexBondFrag = frag.getBond(indexCounterAtomFragmentRingClosure, indexAtomFragment);
			
			final int indexBondMol = mol.getBond(indexCounterAtomMoleculeRingClosure, indexConnAtomMol);
			
			if(indexBondMol==-1){
				ringClosuresMatch=false;
				break;
			}
			
			if(!areBondsSimilar(indexBondMol, indexBondFrag)) {
				ringClosuresMatch=false;
				break;
			}
		}
		
		return ringClosuresMatch;
		
	}
	
	/**
	 * To determine connected atoms and bonds took longer than the saved search time. 
	 * (Asads data set. /home/korffmo1/Projects/Software/Development/SubStructureSearch/asad-Benchmark-9e64792)
	 * @param atomMolecule
	 * @param atomFragment
	 * @return
	 */
	private boolean areAtomsSimilar(int atomMolecule, int atomFragment){
		
		final int moleculeConnAtoms = mol.getConnAtoms(atomMolecule);
		
		final int fragmentConnAtoms = frag.getConnAtoms(atomFragment);
		
		if (fragmentConnAtoms > moleculeConnAtoms)
			return false;

		if (frag.getAtomicNo(atomFragment) == 0)	// attachment points match on any atom
			return true;
		if (mol.getAtomicNo(atomMolecule) == 0)	// real fragment atoms cannot match on attachment points in molecule
			return false;

		
		if(mol.getAtomicNo(atomMolecule)!=frag.getAtomicNo(atomFragment)){
			return false;
		}
		
		return true;
	}
	
	
    private boolean areBondsSimilar(int moleculeBond, int fragmentBond) {
		
    	boolean similar=false;
    	
    	if(mol.isDelocalizedBond(moleculeBond) && frag.isDelocalizedBond(fragmentBond)){
    		
    		similar=true;
    		
    	} else if(mol.getBondType(moleculeBond)==frag.getBondType(fragmentBond)) {
    		
    		similar=true;
    	}
    	
    	if(!similar){
    		return false;
    	}
    	
		return similar;
    }
	
	private boolean isFulfillPrechecks(){
		boolean precheck=true;
		
		if(frag.getAtoms()>mol.getAtoms()){
			return false;
		} else if(frag.getBonds()>mol.getBonds()){
			return false;
		}
		
		if(precheckTotalFormula){
			for (int i = 0; i < arrTotalFormularMol.length; i++) {
				if(arrTotalFormularMol[i]-arrTotalFormularFrag[i]<0) {
					precheck=false;
					break;
				}
			}
		}
		return precheck;
	}

	private void initProcessingArray(){
		
		PossibleMappingsFrag2Mol possMapFrag2MolStart = new PossibleMappingsFrag2Mol(0, MAX_NUM_ATOMS);
		
		arrProcessingArray[0]=possMapFrag2MolStart;
		
		for (int i = 1; i < arrProcessingArray.length; i++) {
			arrProcessingArray[i] = new PossibleMappingsFrag2Mol(i);
		}
	}
	
	/**
	 * 
	 * @return Startobject containing the fragment atom with the highest cardinality and the mapping molecule atoms. 
	 */
	private void determinePossibleMappingsFrag2MolStart() {
		
		int maxCardinality=-1;
		int indexAtomFragMaxCard=-1;
		for (int i = 0; i < frag.getAtoms(); i++) {
			if(arrCardinalityFragment[i]>maxCardinality){
				maxCardinality=arrCardinalityFragment[i];
				indexAtomFragMaxCard=i;
			}
		}
		
		PossibleMappingsFrag2Mol possMapFrag2MolStart = arrProcessingArray[0];
		
		possMapFrag2MolStart.setIndexAtomFrag(indexAtomFragMaxCard);
		
		for (int i = 0; i < mol.getAtoms(); i++) {
			if(areAtomsSimilar(i, indexAtomFragMaxCard)){
			// if(arrAtomMatch[i][indexAtomFragMaxCard]){
				possMapFrag2MolStart.addIndexMappingAtomMolecule(i);
			}
		}
		
		if(VERBOSE){
			if(possMapFrag2MolStart.empty()) {
				System.out.println("No mapping start atom found in molecule.");
			}
		}
		
	}
	
	private void calculateCardinalityFragment(){
		
		for (int i = 0; i < frag.getAtoms(); i++) {
			arrCardinalityFragment[i]=calculateCardinality(frag, i);
		}
	}
	
	private static int calculateCardinality(StereoMolecule mol, int atom){
		
		int cardinality=0;
		
		int atomicNo = mol.getAtomicNo(atom);
		
		if((atomicNo!=1) && (atomicNo != 6)){
			cardinality+=CARDINALITY_HETERO_ATOM;
		}
		
		if(mol.isRingAtom(atom)){
			cardinality+=CARDINALITY_RING_ATOM;
		}
		
		return cardinality;
	}

	/**
	 * The sub structure search algortihm stops if the maximum capacity of the match list container is reached. 
	 * 
	 * @return true if maximum capacity of the match list container is reached. False otherwise.
	 */
	public boolean isMaxCapacityMatchListContainerReached() {
		return maxCapacityMatchListContainerReached;
	}

}
