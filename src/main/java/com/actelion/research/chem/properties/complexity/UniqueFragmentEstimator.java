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

package com.actelion.research.chem.properties.complexity;

import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.mcs.ContainerListWithIntVec;
import com.actelion.research.chem.mcs.ListWithIntVec;
import com.actelion.research.util.SizeOf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

public class UniqueFragmentEstimator {
	
	private static final int SIZE_SUBSET_SOLUTIONS = 10000;
	
	private static final int NUM_BONDS_ESTIMATION_STARTS = 16;
	
	
	private static final int MAX_NUM_BONDS = 256;
	
	private static int CAPACITY_FRAGMENTS = 1000000;
	
	// This list contains the solutions. Each row in the list contains the solutions for the corresponding number of bonds.
	// Consequently the first row is empty.
	private List<ContainerSolutions> liSolutions;
	
	private ExtendedMolecule mol;
	
	private int nBondsMolecule;
	
	private int sizeIntVec;
	
	private int maximumNumberBondsInFrag;
	
	private ListWithIntVec livNeighbours;
	
	private ContainerListWithIntVec containerListWithIntVec;
		
	private boolean fragmentsGenerated;
	
	private int maxCapacityFragments;
	
	private boolean capacityLimitBreakes;
	
	/**
	 * From that number of bonds only a subset of parent solutions will be taken to 
	 * generate new solutions.
	 */
	int numBondsEstimationStarts;
	
	int sizeSubsetSolutions;
	
	
	public UniqueFragmentEstimator() {
		this(CAPACITY_FRAGMENTS);
	}
	
	/**
	 * 
	 * @param totalCapacity maximum number of unique index combinations. 
	 * When the limit is exceeded the variable <code>capacityLimitBreakes</code> is increased by one.
	 */
	public UniqueFragmentEstimator(int totalCapacity) {
		
		maxCapacityFragments = totalCapacity;
				
		int maxSizeIntVec = getSizeArrayLIV();
		
		containerListWithIntVec = new ContainerListWithIntVec(maxSizeIntVec, totalCapacity);
		
		System.out.println("ExhaustiveFragmentGeneratorBonds Used mem " + SizeOf.usedMemoryMB() + "[MB].");
		
		sizeSubsetSolutions = SIZE_SUBSET_SOLUTIONS;
		
		
		numBondsEstimationStarts = NUM_BONDS_ESTIMATION_STARTS;

	}
	
	
	public void set(ExtendedMolecule mol, int nMaximumNumberBonds) {
		
		initHashSet(nMaximumNumberBonds+1);
		
		this.mol = mol;
		
		nBondsMolecule = mol.getBonds();
		
		if(nBondsMolecule>MAX_NUM_BONDS){
			throw new RuntimeException("Maximum number of atoms exceeded.");
		}
		
		sizeIntVec = (nBondsMolecule + Integer.SIZE-1)/ Integer.SIZE;
		
		this.maximumNumberBondsInFrag = nMaximumNumberBonds;
		
		livNeighbours = new ListWithIntVec(sizeIntVec, -1);
		
		fragmentsGenerated = false;
		
		containerListWithIntVec.reset();
		
		capacityLimitBreakes=false;
	}

	private void initHashSet(int nBondsFragsNew){
	
		if(liSolutions==null){
			liSolutions = new ArrayList<ContainerSolutions>();
		}
		
		int sizeOld = liSolutions.size();
		for (int i = 0; i < liSolutions.size(); i++) {
			liSolutions.get(i).clear();
		}
		
		for (int i = sizeOld; i < nBondsFragsNew+1; i++) {
			
			int capacity = 1000;
			if(i > 15) {
				capacity = 1000000;
			}else if(i > 10) {
				capacity = 100000;
			}else if(i > 5) {
				capacity = 10000;
			}
			
			liSolutions.add(new ContainerSolutions(capacity, i));

		}
		
	}
	
	public void generateFragments(){
		
		createInitSolutions();
		
		for (int i = 1; i < maximumNumberBondsInFrag+1; i++) {
			
			getAllPossibleNeigbourCombinations(i);
		}
		
		fragmentsGenerated = true;
		
		System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() getTotalSizeResultList() " + getTotalSizeResultList() + ".");

	}
	
	
	private void createInitSolutions(){
		
		ContainerSolutions csInit = liSolutions.get(1);
		
		for (int i = 0; i < maximumNumberBondsInFrag+1; i++) {
			ListWithIntVec livStart = containerListWithIntVec.get();
			
			// Start solution for this bond.
			livStart.addBit(i);
			
			livStart.calculateHash();
			
			csInit.add(livStart);
		}
		
	}


	/**
	 * Returns list with indices for fragments with <code>size</code> atoms. The indices are coded as bit lists.
	 * @param size
	 * @return
	 */
	public List<ListWithIntVec> get(int size) {
		
		if(!fragmentsGenerated){
			throw new RuntimeException("Fragments have to be generated first. Call generateFragments().");
		}
		
		return liSolutions.get(size).getSolutions();
	}
	
	private void getAllPossibleNeigbourCombinations(int numBondsParent){
		
		
		ContainerSolutions csParent = liSolutions.get(numBondsParent);
		
		List<ListWithIntVec> liSolutionsParent = csParent.getSolutions();
		
		List<ListWithIntVec> liSolutionsParentSubset = null; 
		if(numBondsParent >= numBondsEstimationStarts){
			
			Collections.shuffle(liSolutionsParent);
			
			liSolutionsParentSubset = liSolutionsParent.subList(0, sizeSubsetSolutions);
			
		} else {
			liSolutionsParentSubset = liSolutionsParent;	
		} 
		
		
		ContainerSolutions csChild = liSolutions.get(numBondsParent+1);
		
		if(maximumNumberBondsInFrag==1){
			return;
		}
		
		for (ListWithIntVec livSolutionParent : liSolutionsParentSubset) {
			
			ListWithIntVec livIndexBondsReachable = getAllReachableNeighbourBonds(mol, livSolutionParent);

			for (int i = 0; i < livIndexBondsReachable.size(); i++) {
				
				ListWithIntVec livPathExtended = containerListWithIntVec.getWithCopy(livSolutionParent);
				 
				if(livPathExtended.addBit(livIndexBondsReachable.get(i))){
					 
					livPathExtended.calculateHash();
					
					if(!csChild.add(livPathExtended)){
						
						containerListWithIntVec.back(livPathExtended);
						
						
					} 
				}
			}
		}
			
	}
	
	private int getTotalSizeResultList (){
		int size = 0;
		
		for(ContainerSolutions cs : liSolutions) {
			size += cs.size();
		}
		
		return size;
	}
	
	private void shrinkResultList(){
		
		for (int i = 0; i < liSolutions.size()-1; i++) {
			int size = liSolutions.get(i+1).size();
			
			if(size>0){
				liSolutions.get(i).clear();
			}
		}
	}
	
	private void clearResultList(){
		
		for (int i = 0; i < liSolutions.size(); i++) {
			
			liSolutions.get(i).clear();
			
		}
	}
	
	private final ListWithIntVec getAllReachableNeighbourBonds(ExtendedMolecule mol, ListWithIntVec livIndexBond){
		
		livNeighbours.reset();
		
		for (int i = 0; i < livIndexBond.size(); i++) {
			
			final int indexBond = livIndexBond.get(i);
			
			final int indexAtom1 = mol.getBondAtom(0, indexBond);
			
			final int indexAtom2 = mol.getBondAtom(1, indexBond);
			
			final int nConnected1 = mol.getAllConnAtoms(indexAtom1);
			
			for (int j = 0; j < nConnected1; j++) {
				
				int indexAtConn = mol.getConnAtom(indexAtom1, j);
				
				int indexBondConnected = mol.getBond(indexAtom1, indexAtConn);
				
				if(!livIndexBond.isBitSet(indexBondConnected)) {
					livNeighbours.addBit(indexBondConnected);
				}
			}
			
			final int nConnected2 = mol.getAllConnAtoms(indexAtom2);
			
			for (int j = 0; j < nConnected2; j++) {
				
				int indexAtConn = mol.getConnAtom(indexAtom2, j);
				
				int indexBondConnected = mol.getBond(indexAtom2, indexAtConn);
				
				if(!livIndexBond.isBitSet(indexBondConnected)) {
					livNeighbours.addBit(indexBondConnected);
				}
			}
		}
		
		return livNeighbours;
	}

	/**
	 * If true not all index combinations where generated.
	 * Starts with 0 for each new molecule. 
	 * @return
	 */
	public boolean isCapacityLimitBreakes() {
		return capacityLimitBreakes;
	}
	
	
	public static int getSizeArrayLIV(){
		
		int maxSizeIntVec = (MAX_NUM_BONDS + Integer.SIZE-1)/ Integer.SIZE;
		
		return maxSizeIntVec;

	}
	
	/**
	 * 
	 * ContainerSolutions
	 * @author Modest von Korff
	 * @version 1.0
	 * Mar 1, 2013 MvK Start implementation
	 * Container with unique substructures for a certain number of atoms.
	 */
	private static class ContainerSolutions {
		
		private int bonds;
		
		private HashSet<ListWithIntVec> hsSolution;
				
		public ContainerSolutions(int capacity, int bonds) {
			
			hsSolution = new HashSet<ListWithIntVec>(capacity);
			
			this.bonds = bonds;
		}	
				
		public List<ListWithIntVec> getSolutions(){
			return new ArrayList<ListWithIntVec>(hsSolution);
		}
		
		public boolean add(ListWithIntVec liv){
			return hsSolution.add(liv);
		}
		
		public void clear(){
			hsSolution.clear();
		}
		
		public int size(){
			return hsSolution.size();
		}
		
		public int getBonds(){
			return bonds;
		}
		
		
	}

}
