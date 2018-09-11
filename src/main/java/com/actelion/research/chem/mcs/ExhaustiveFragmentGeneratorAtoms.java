package com.actelion.research.chem.mcs;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.SizeOf;



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
public class ExhaustiveFragmentGeneratorAtoms {
	
	
	private static final boolean DEBUG = false;
	
	
	private static int MAX_NUM_FRAGMENTS = 400000;
	
	private static int CAPACITY_FRAGMENTS = 8000000;
	
	static {
		if(DEBUG){
			MAX_NUM_FRAGMENTS = 20000;
			CAPACITY_FRAGMENTS = 500000;
		}
	}
	
	private static final int MAX_NUM_ATOMS = 256;

	
	private List<HashSet<ListWithIntVec>> liHashSetIntegerList;
	
	private ExtendedMolecule mol;
	
	private int nAtomsMolecule;
	
	private int sizeIntVec;
	
	private int maximumFragmentSite;
	
	private ListWithIntVec livNeighbours;
	
	private ContainerListWithIntVec containerListWithIntVec;
	
	private List<ListWithIntVec> liIntegerListSolution;
	
	private List<ListWithIntVec> liIntegerListTmp;
	
	private boolean fragmentsGenerated;
	
	private int capacityLimitBreakes;
	

	
	public ExhaustiveFragmentGeneratorAtoms() {
		
		liIntegerListSolution = new ArrayList<ListWithIntVec>();
		
		liIntegerListTmp = new ArrayList<ListWithIntVec>();
		
		int maxSizeIntVec = (MAX_NUM_ATOMS + Integer.SIZE-1)/ Integer.SIZE;
		
		containerListWithIntVec = new ContainerListWithIntVec(maxSizeIntVec, CAPACITY_FRAGMENTS);
		
		System.out.println("ExhaustiveFragmentGenerator Used mem " + SizeOf.usedMemoryMB() + "[MB].");
		
	}
	
	public void set(ExtendedMolecule mol, int nMaximumFragmentSize) {
		
		initHashSet(nMaximumFragmentSize);
		
		this.mol = mol;
		
		nAtomsMolecule = mol.getAtoms();
		
		if(nAtomsMolecule>MAX_NUM_ATOMS){
			throw new RuntimeException("Maximum number of atoms exceeded.");
		}
		
		sizeIntVec = (nAtomsMolecule + Integer.SIZE-1)/ Integer.SIZE;
		
		this.maximumFragmentSite = nMaximumFragmentSize;
		
		livNeighbours = new ListWithIntVec(sizeIntVec, -1);
		
		fragmentsGenerated = false;
		
		containerListWithIntVec.reset();
		
		capacityLimitBreakes=0;
	}

	private void initHashSet(int nAtomsFragsNew){
		
		liIntegerListSolution.clear();
		
		liIntegerListTmp.clear();
		
		if(liHashSetIntegerList==null){
			liHashSetIntegerList = new ArrayList<HashSet<ListWithIntVec>>();
		}
		
		int sizeOld = liHashSetIntegerList.size();
		for (int i = 0; i < liHashSetIntegerList.size(); i++) {
			liHashSetIntegerList.get(i).clear();
		}
		
		for (int i = sizeOld; i < nAtomsFragsNew+1; i++) {
			
			int capacity = 1000;
			if(i > 15) {
				capacity = 1000000;
			}else if(i > 10) {
				capacity = 100000;
			}else if(i > 5) {
				capacity = 10000;
			}
			
			liHashSetIntegerList.add(new HashSet<ListWithIntVec>(capacity));

		}
		
	}
	
	public void generateFragments(){
		
		int nAtomsMol = mol.getAtoms();

		for (int i = 0; i < nAtomsMol; i++) {
			
			getAllPossibleNeigbourCombinationsConsidersPrevious(i);
		}
		
		fragmentsGenerated = true;

	}
	
	/**
	 * Returns a bit list. 
	 * @param size
	 * @return
	 */
	public List<ListWithIntVec> get(int size) {
		
		if(!fragmentsGenerated){
			throw new RuntimeException("Fragments have to be generated first. Call generateFragments().");
		}
		
		return new ArrayList<ListWithIntVec>(liHashSetIntegerList.get(size));
	}
	
	public List<ListWithIntVec> getFragmentsForSingleAtom(int indexStartAtom) {
		
		initHashSet(maximumFragmentSite);
		
		return getAllPossibleNeigbourCombinationsConsidersPrevious(indexStartAtom);
	}

	private List<ListWithIntVec> getAllPossibleNeigbourCombinationsConsidersPrevious(int indexStartAtom){
		
		liIntegerListSolution.clear();
		
		ListWithIntVec livStart = containerListWithIntVec.get();
		
		// Start solution for this atom
		livStart.addBit(indexStartAtom);
		
		livStart.calculateHash();
		
		liIntegerListSolution.add(livStart);
		
		if(maximumFragmentSite==1){
			return liIntegerListSolution;
		}
		
		boolean finish=false;
		
		while(!finish){
			
			liIntegerListTmp.clear();
			
			int numAtomsInFragmentParent=liIntegerListSolution.get(0).size();
			
			HashSet<ListWithIntVec> hsIntListGeneratedSolutions = liHashSetIntegerList.get(numAtomsInFragmentParent+1);

			if(hsIntListGeneratedSolutions.size()>=MAX_NUM_FRAGMENTS) {
				
				// System.err.println("ExhaustiveFragmentGenerator getAllPossibleNeigbourCombinationsConsidersPrevious while Break forced");
				
				capacityLimitBreakes++;
				break;
			}
			
			breakLimitFrags:
			for (ListWithIntVec livAtomPath : liIntegerListSolution) {
				 
				ListWithIntVec livIndexAtomsReachable = getAllReachableNeighbours(mol, livAtomPath);
				 
				for (int i = 0; i < livIndexAtomsReachable.size(); i++) {
					
					ListWithIntVec livPathExtended = containerListWithIntVec.getWithCopy(livAtomPath);
					 
					if(livPathExtended.addBit(livIndexAtomsReachable.get(i))){
						 
						livPathExtended.calculateHash();
						
						if(hsIntListGeneratedSolutions.add(livPathExtended)){
							liIntegerListTmp.add(livPathExtended);
							
							if(liIntegerListTmp.size()==MAX_NUM_FRAGMENTS){
								
								// System.err.println("ExhaustiveFragmentGenerator getAllPossibleNeigbourCombinationsConsidersPrevious while Break forced");
								
								capacityLimitBreakes++;

								break breakLimitFrags;
								// System.out.println("liIntegerListTmp " + liIntegerListTmp.size());
							}
						} else {
							containerListWithIntVec.back(livPathExtended);
						}
					} else {
						containerListWithIntVec.back(livPathExtended);
					}
				}
			}
				 
			
			liIntegerListSolution.clear();
			 
			liIntegerListSolution.addAll(liIntegerListTmp);
			
			if(liIntegerListSolution.isEmpty()) {
				finish=true;
			} else {
				int lenFirst=liIntegerListSolution.get(0).size();
				 
				for (ListWithIntVec livAtomPath : liIntegerListSolution) {
					if(livAtomPath.size()==maximumFragmentSite){
						finish=true;
					}
					 
					if(lenFirst!=livAtomPath.size()){
						throw new RuntimeException("Error in algorithm.");
					}
				}
			}
			
			// System.out.println("ExhaustiveFragmentGenerator while Used mem " + SizeOf.usedMemoryMB() + "[MB].");

			
		}
		
		return liIntegerListSolution;
	}
	
	private ListWithIntVec getAllReachableNeighbours(ExtendedMolecule mol, ListWithIntVec livIndexAtom){
		
		livNeighbours.reset();
		
		for (int i = 0; i < livIndexAtom.size(); i++) {
			final int indexAtom = livIndexAtom.get(i);
			
			final int nConnected = mol.getAllConnAtoms(indexAtom);
			
			for (int j = 0; j < nConnected; j++) {
				
				int indexAtConn = mol.getConnAtom(indexAtom, j);
				
				if(!livIndexAtom.isBitSet(indexAtConn)) {
					livNeighbours.addBit(indexAtConn);
				}
			}
		}
		
		return livNeighbours;
	}

	/**
	 * Is increased by one if more than the the allowed number of fragments was generated. 
	 * Starts with 0 for each new molecule. 
	 * @return
	 */
	public int getCapacityLimitBreakes() {
		return capacityLimitBreakes;
	}
	
	public static String getIdCodeFromFragment(StereoMolecule mol, ListWithIntVec livIndexAtom){
		
		boolean [] includeAtom = new boolean [mol.getAtoms()];
		
		for (int i = 0; i < livIndexAtom.size(); i++) {
			includeAtom[livIndexAtom.get(i)]=true;
		}
		
		StereoMolecule frag = new StereoMolecule();
		
		mol.copyMoleculeByAtoms(frag, includeAtom, true, null);
		
		frag.ensureHelperArrays(Molecule.cHelperRings);
		
		Canonizer can = new Canonizer(frag);
		
		return can.getIDCode();
	}


}
