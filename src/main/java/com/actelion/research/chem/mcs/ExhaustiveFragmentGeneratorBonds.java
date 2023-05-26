package com.actelion.research.chem.mcs;

import java.util.Date;
import java.util.List;

import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.properties.complexity.ContainerFragBondsSolutions;
import com.actelion.research.chem.properties.complexity.IBitArray;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.SizeOf;
import com.actelion.research.util.datamodel.IntArray;



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
public class ExhaustiveFragmentGeneratorBonds {
	
	private static boolean ELUSIVE = false;

	public static final long LIMIT_NEIGHBOURS_SINCE_LAST_ADDED = 100 * 1000 * 1000;
		
	// This list contains the solutions. Each row in the list contains the solutions for the corresponding number of atom types.
	// Consequently, the first two rows are empty.
	private ContainerFragBondsSolutions containerDataFragDefByBonds;
	
	private ExtendedMolecule mol;
	
	private int nBondsMolecule;
	
	private int maximumNumberBondsInFrag;
	
	private IntArray arrIndexReachableNeighbours;
	
	
	private boolean fragmentsGenerated;
		
	private boolean capacityLimitBreakes;
		
	private long solutionAdded;

	private int totalMaximumCapacity;
	
	
	/**
	 * @param bits maximum fragment size that can be stored as number of bonds.
	 */
	public ExhaustiveFragmentGeneratorBonds(int bits, int totalMaximumCapacity) {

		ContainerFragBondsSolutions.ELUSIVE = ELUSIVE;

		this.totalMaximumCapacity = totalMaximumCapacity;

		containerDataFragDefByBonds = new ContainerFragBondsSolutions(bits, totalMaximumCapacity);
				
		arrIndexReachableNeighbours = new IntArray();
		
		if(ELUSIVE)
			System.out.println("ExhaustiveFragmentGeneratorBonds constructor finished, used mem " + SizeOf.usedMemoryMB() + "[MB].");
	}

	public void set(ExtendedMolecule mol, int nMaximumNumberBonds) {
		
		this.mol = mol;
		
		nBondsMolecule = mol.getBonds();
		
		if(nBondsMolecule > containerDataFragDefByBonds.getSizeBinaryArray()){
			throw new RuntimeException("Maximum number of bonds exceeded.");
		}
				
		this.maximumNumberBondsInFrag = Math.min(nMaximumNumberBonds, containerDataFragDefByBonds.getMaximumCapacityBondsInFragment());
				
		fragmentsGenerated = false;
				
		capacityLimitBreakes=false;
		
		containerDataFragDefByBonds.setBondsMolecule(nBondsMolecule);
		
	}

	private void initDataContainerAllSingleBonds(){
		
		containerDataFragDefByBonds.reset();
		
		for (int i = 0; i < nBondsMolecule; i++) {
			IBitArray f = containerDataFragDefByBonds.get();
			
			f.setBit(i);
						
			containerDataFragDefByBonds.addFacultative(f);
		}
	}
	
	
	private void initDataContainerOneSingleBond(int indexBond) {
		
		containerDataFragDefByBonds.reset();
				
		IBitArray f = containerDataFragDefByBonds.get();
		
		f.setBit(indexBond);
					
		containerDataFragDefByBonds.addFacultative(f);
	}
	
	public void generateFragmentsAllBonds(){
		
		initDataContainerAllSingleBonds();
		
		generateFragments();
	}
	
	
	public void generateFragmentsForSingleBond(int indexBond) {
		
		initDataContainerOneSingleBond(indexBond);
		
		generateFragments();
		
	}
	
	private void generateFragments(){
		
		int maxNumBondsFrag = Math.min(nBondsMolecule, maximumNumberBondsInFrag);
		
		if(maxNumBondsFrag==1){
			return;
		}
		
		long neighboursTotal = 0;
		long addedTotal = 0;
		if(ELUSIVE) {
			System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() start record capacity " + Formatter.group(containerDataFragDefByBonds.getAvailable()) + ".");
		}
		maximumCapacityBreak:
		for (int i = 1; i < maxNumBondsFrag; i++) {
			
			List<IBitArray> liParent = containerDataFragDefByBonds.getList(i);

			if(ELUSIVE) {
				System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() bonds  " + i + ". Parents " + Formatter.group(liParent.size()) + ".");

			}

			long added = 0;
			
			long neighbours = 0;
			
			long neighboursSinceLastAdded = 0;
			
			for (IBitArray fParent : liParent) {
				
				IntArray arrBondsReachable = getAllReachableNeighbourBonds(mol, fParent);
				 
				for (int j = 0; j < arrBondsReachable.length(); j++) {
										
					IBitArray fChildAddedBond = containerDataFragDefByBonds.getWithCopy(fParent);
					
					fChildAddedBond.setBit(arrBondsReachable.get(j));
					
					// System.out.println(fChildAddedBond.toString());
					
					if(containerDataFragDefByBonds.addFacultative(fChildAddedBond)){
						added++;
						addedTotal++;
						neighboursSinceLastAdded = 0;

						if(addedTotal>totalMaximumCapacity) {
							if(ELUSIVE) {
								System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() maximum capacity (" + totalMaximumCapacity + ") break.");
								log(neighbours, neighboursTotal, added, addedTotal, neighboursSinceLastAdded);
							}
							// Clear all records for this number of bonds.
							containerDataFragDefByBonds.reset(fChildAddedBond.getBitsSet());

							break maximumCapacityBreak;
						}
					}
					neighbours++;
					
					neighboursSinceLastAdded++;
					
					if(neighbours % 50000000 == 0){
						if(ELUSIVE) {
							log(neighbours, neighboursTotal, added, addedTotal, neighboursSinceLastAdded);
						}
					}
				}
				
				// if(neighboursSinceLastAdded > LIMIT_NEIGHBOURS_SINCE_LAST_ADDED) {
				if(neighboursSinceLastAdded > totalMaximumCapacity) {
					System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments(). Break for fragments with " + i + " bonds. Generated  " + totalMaximumCapacity + " neighbours since last add to hash map.");
					break;
				}
			}
			
			// System.out.println("Bonds\t" + i + "\tparents\t" + liParent.size() + "\tneighbours added\t" + added);

			
			neighboursTotal += neighbours;

			
			if(ELUSIVE) {
				log(neighbours, neighboursTotal, added, addedTotal, neighboursSinceLastAdded);
			}
		}
		
		fragmentsGenerated = true;
		
		if(ELUSIVE) {
			System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() getTotalSizeResultList() " + Formatter.group(containerDataFragDefByBonds.getTotalSizeResults()) + ".");
			System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() solutionAdded " + Formatter.group(solutionAdded) + ".");
		}
	}

	private void log(long neighbours, long neighboursTotal, long added, long addedTotal, long neighboursSinceLastAdded){
		System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() " + new Date().toString() + ".");
		System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() neighbours generated " + Formatter.group(neighbours) + ".");
		System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() neighbours generated total " + Formatter.group(neighboursTotal) + ".");
		System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() solutions added " + Formatter.group(added) + ".");
		System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() solutions added total " + Formatter.group(addedTotal) + ".");
		System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() available  " + Formatter.group(containerDataFragDefByBonds.getAvailable()) + ".");
		System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() capacity  " + Formatter.group(containerDataFragDefByBonds.getCapacity()) + ".");
		System.out.println("ExhaustiveFragmentGeneratorBonds generateFragments() neighboursSinceLastAdded  " + Formatter.group(neighboursSinceLastAdded) + ".");
	}


	/**
	 * Returns list with indices for fragments with <code>size</code> bonds. The indices are coded as bit lists.
	 * @param bonds
	 * @return
	 */
	public List<IBitArray> getFragments(int bonds) {
		
		if(!fragmentsGenerated){
			throw new RuntimeException("Fragments have to be generated first. Call generateFragments().");
		}
		
		return containerDataFragDefByBonds.getList(bonds);
	}
	
	/**
	 * get all bonds that can be reached in one step from the input bonds.
	 * @param mol
	 * @param livIndexBond
	 * @return
	 */
	private final IntArray getAllReachableNeighbourBonds(ExtendedMolecule mol, IBitArray livIndexBond){
		
		arrIndexReachableNeighbours.reset();
				
		for (int i = 0; i < nBondsMolecule; i++) {
			
			if(!livIndexBond.isBitSet(i)){
				continue;
			}
			
			final int indexBond = i;
			
			final int indexAtom1 = mol.getBondAtom(0, indexBond);
			
			final int indexAtom2 = mol.getBondAtom(1, indexBond);
			
			final int nConnected1 = mol.getAllConnAtoms(indexAtom1);
			
			for (int j = 0; j < nConnected1; j++) {
				
				int indexAtConn = mol.getConnAtom(indexAtom1, j);
				
				int indexBondConnected = mol.getBond(indexAtom1, indexAtConn);
									
				if(!livIndexBond.isBitSet(indexBondConnected)) {
					
					solutionAdded++;
					
					arrIndexReachableNeighbours.add(indexBondConnected);
				}
			}
			
			final int nConnected2 = mol.getAllConnAtoms(indexAtom2);
			
			for (int j = 0; j < nConnected2; j++) {
				
				int indexAtConn = mol.getConnAtom(indexAtom2, j);
				
				int indexBondConnected = mol.getBond(indexAtom2, indexAtConn);
											
				if(!livIndexBond.isBitSet(indexBondConnected)) {
					
					solutionAdded++;
					
					arrIndexReachableNeighbours.add(indexBondConnected);
				}
			}
		}
		
		return arrIndexReachableNeighbours;
	}

	/**
	 * If true not all index combinations were generated.
	 * Starts with 0 for each new molecule. 
	 * @return
	 */
	public boolean isCapacityLimitBreakes() {
		return capacityLimitBreakes;
	}
	
	
	public int getSizeArrayLIV(){
		
		int maxSizeIntVec = (containerDataFragDefByBonds.getSizeBinaryArray() + Integer.SIZE-1)/ Integer.SIZE;
		
		return maxSizeIntVec;

	}
	
	public int getMaximumCapacityBondsInFragment() {
		return containerDataFragDefByBonds.getMaximumCapacityBondsInFragment();
	}

	public int getMaximumNumberBondsInMolecule(){
		return containerDataFragDefByBonds.getMaximumNumberBondsInMolecule();
	}

	/**
	 * @return the eLUSIVE
	 */
	public static boolean isELUSIVE() {
		return ELUSIVE;
	}


	/**
	 * @param elusive the ELUSIVE to set
	 */
	public static void setELUSIVE(boolean elusive) {
		ELUSIVE = elusive;
		ContainerFragBondsSolutions.ELUSIVE = ELUSIVE;
	}

}
