package com.actelion.research.chem.mcs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.chem.ExtendedMoleculeFunctions;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.datamodel.IntVec;

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
public class MCSFast {
	
	public static final int PAR_CLEAVE_RINGS=0; 
	
	public static final int PAR_KEEP_RINGS=1;
	
	public static final int PAR_KEEP_AROMATIC_RINGS=2; 
	
	private static final boolean DEBUG = false;
	
	private static final int CAPACITY = (60+60)*10; 
	
	private StereoMolecule mol;
	
	private StereoMolecule frag;
	
	private HashSet<IntVec> hsIndexFragCandidates;
	
	private HashSet<IntVec> hsIndexFragGarbage;
	
	private HashSet<IntVec> hsIndexFragSolution;
	
	private SSSearcher sss;
	
	private ComparatorBitsSet comparatorBitsSet;
	
	private StereoMolecule molMCS;
	
	private boolean considerAromaticRings;
	
	private boolean considerRings;
	
	// The largest valid solutions.
	private List<IntVec> liMCSSolutions;
	
	private HashMap<Integer, List<int []>> hmRingBnd_ListRingBnds;
	
	private HashMap<Integer, List<int []>> hmAromaticRingBnd_ListRingBnds;
	
	private RingCollection ringCollection;
	
	public MCSFast() {
		this(PAR_CLEAVE_RINGS);
	}
	
	public MCSFast(int ringStatus) {
		
		considerRings=false;
		
		considerAromaticRings=false;
		
		switch (ringStatus) {
		case PAR_CLEAVE_RINGS:
			break;
		case PAR_KEEP_RINGS:
			considerRings=true;
			break;
		case PAR_KEEP_AROMATIC_RINGS:
			considerAromaticRings=true;
			break;

		default:
			break;
		}
		
		hsIndexFragCandidates = new HashSet<IntVec>(CAPACITY);
		
		hsIndexFragGarbage = new HashSet<IntVec>(CAPACITY);
		
		hsIndexFragSolution = new HashSet<IntVec>(CAPACITY);
		
		liMCSSolutions = new ArrayList<IntVec>(CAPACITY);
		
		sss = new SSSearcher();
		
		comparatorBitsSet = new ComparatorBitsSet();
		
		hmRingBnd_ListRingBnds = new HashMap<Integer, List<int []>>();
		
		hmAromaticRingBnd_ListRingBnds = new HashMap<Integer, List<int []>>();
	}
	
	public void set(StereoMolecule mol, StereoMolecule frag) {
		
		this.mol = mol;
		
		this.frag = frag;
		
		init();
	}
	
	
	private void init(){
		
		hsIndexFragCandidates.clear();
		
		hsIndexFragGarbage.clear();
		
		hsIndexFragSolution.clear();
		
		liMCSSolutions.clear();
		
		hmRingBnd_ListRingBnds.clear();
		
		hmAromaticRingBnd_ListRingBnds.clear();
		
		initCandidates();
		
	}
	
	private void initCandidates(){
		
		ringCollection = frag.getRingSet();
		
		int rings = ringCollection.getSize();
						
		for (int i = 0; i < rings; i++) {
			
			int [] arrIndexBnd = ringCollection.getRingBonds(i);
			
			for (int j = 0; j < arrIndexBnd.length; j++) {
				
				if(!hmRingBnd_ListRingBnds.containsKey(arrIndexBnd[j])){
					
					hmRingBnd_ListRingBnds.put(arrIndexBnd[j], new ArrayList<int []>());
					
				}
				
				List<int []> li = hmRingBnd_ListRingBnds.get(arrIndexBnd[j]);
				
				li.add(arrIndexBnd);
				
			}
			
			if(ringCollection.isAromatic(i)){
				
				for (int j = 0; j < arrIndexBnd.length; j++) {
					
					if(!hmAromaticRingBnd_ListRingBnds.containsKey(arrIndexBnd[j])){
						
						hmAromaticRingBnd_ListRingBnds.put(arrIndexBnd[j], new ArrayList<int []>());
						
					}
					
					List<int []> li = hmAromaticRingBnd_ListRingBnds.get(arrIndexBnd[j]);
					
					li.add(arrIndexBnd);
					
				}
			}
			
			
		}
		
		// Array length for bit list.
		int nInts = (int)(((double)frag.getBonds() / Integer.SIZE) + ((double)(Integer.SIZE-1)/Integer.SIZE));
		
		// The vector is used as bit vector.
		// For each bond one bit.
		for (int i = 0; i < frag.getBonds(); i++) {
			
			IntVec iv = new IntVec(nInts);
			
			setBitAndAddRelatedRingBonds(i, iv);
			
			hsIndexFragCandidates.add(iv);
			
		}
	}
	
	
	private void setBitAndAddRelatedRingBonds(int bit, IntVec iv){
		if(!considerAromaticRings && !considerRings) {
			iv.setBit(bit);
			
		} else if(considerRings) {
			
			iv.setBit(bit);
			
			if(hmRingBnd_ListRingBnds.containsKey(bit)){
				List<int []> li = hmRingBnd_ListRingBnds.get(bit);
				
				for (int i = 0; i < li.size(); i++) {
					int [] a = li.get(i);
					
					for (int j = 0; j < a.length; j++) {
						iv.setBit(a[j]);
					}
				}
				
			}
			
		} else if(considerAromaticRings) {
			
			iv.setBit(bit);
			if(hmAromaticRingBnd_ListRingBnds.containsKey(bit)){
				
				List<int []> li = hmAromaticRingBnd_ListRingBnds.get(bit);
				
				for (int i = 0; i < li.size(); i++) {
					int [] a = li.get(i);
					
					for (int j = 0; j < a.length; j++) {
						iv.setBit(a[j]);
					}
				}
			}
		}
		
		iv.calculateHashCode();
	}
	
	/**
	 * Checks first fragment for being sub structure of molecule. If true, frag is returned.
	 * @return
	 */
	private List<IntVec> getAllSolutionsForCommonSubstructures (){
		
		sss.setMolecule(mol);
		
		frag.setFragment(true);
		
		sss.setFragment(frag);
		
		// The MCS is the complete fragment.
		try {
			if(sss.findFragmentInMolecule()>0){
				molMCS = frag;
				
				List<IntVec> liIndexFragCandidates = new ArrayList<IntVec>(hsIndexFragCandidates);
				
				IntVec iv = liIndexFragCandidates.get(0);
				
				iv.setBits(0, iv.sizeBits());
				
				iv.calculateHashCode();
				
				List<IntVec> li = new ArrayList<IntVec>();
				
				li.add(iv);
				
				return li;
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
				
		int maxSizeCandidates = 0;
		
		while(!hsIndexFragCandidates.isEmpty()){
			
			List<IntVec> liIndexFragCandidates = new ArrayList<IntVec>(hsIndexFragCandidates);
			
			Collections.sort(liIndexFragCandidates, comparatorBitsSet);
			
			IntVec iv = liIndexFragCandidates.get(liIndexFragCandidates.size()-1);
			
			hsIndexFragCandidates.remove(iv);
			
			if(DEBUG) {
				System.out.println("Bits set " + iv.getBitsSet());
				
				if(iv.getBitsSet() == frag.getBonds()){
					System.out.println("Full structure in iv.");
				}
			}
			
			StereoMolecule fragSub = getSubFrag(frag, iv);
			
			sss.setFragment(fragSub);
			
			if(sss.findFragmentInMolecule()>0){

				hsIndexFragSolution.add(iv);

				removeAllSubSolutions(iv);
				
				if(iv.getBitsSet() != frag.getBonds()) {
					List<IntVec> liIV = getAllPlusOneAtomCombinations(iv, frag);
					
					for (IntVec ivPlus : liIV) {
						
						if(DEBUG) {
							if(ivPlus.getBitsSet() == frag.getBonds()){
								System.out.println("Full structure in ivPlus.");
							}
						}
						
						if((!hsIndexFragGarbage.contains(ivPlus)) && (!hsIndexFragSolution.contains(ivPlus))){
							
							hsIndexFragCandidates.add(ivPlus);	
						}
						
					}
					
					if(DEBUG)
						System.out.println("tsIndexFragCandidates " + hsIndexFragCandidates.size());
				}
				
				if(maxSizeCandidates < hsIndexFragCandidates.size())
					maxSizeCandidates = hsIndexFragCandidates.size();
				
			} else {
				hsIndexFragGarbage.add(iv);
			}
		}
		
		if(hsIndexFragSolution.size()==0){
			return null;
		}
		
		return getFinalSolutionSet(hsIndexFragSolution);
	}
	
	/**
	 * All molecules which are sub structures of an other molecule in the list are removed.
	 * @return
	 */
	public List<StereoMolecule> getAllCommonSubstructures (){
		
		List<IntVec> liIndexFragSolution = getAllSolutionsForCommonSubstructures();
		
		if(liIndexFragSolution==null){
			return null;
		}
		
		Collections.sort(liIndexFragSolution, comparatorBitsSet);
		
		IntVec ivMCS = liIndexFragSolution.get(liIndexFragSolution.size()-1);
		
		molMCS = getSubFrag(frag, ivMCS);
		
		List<StereoMolecule> li = new ArrayList<StereoMolecule>();
		
		for (IntVec iv : liIndexFragSolution) {
			li.add(getSubFrag(frag, iv));
		}
		
		return ExtendedMoleculeFunctions.removeSubStructures(li);
	}

	/**
	 * 
	 * @return maximum common substructure at top of list or null of none common MCS was found.
	 */
	public StereoMolecule getMCS(){
		
		List<IntVec> liIndexFragSolution = getAllSolutionsForCommonSubstructures ();
		
		if(liIndexFragSolution==null){
			return null;
		}
		
		Collections.sort(liIndexFragSolution, comparatorBitsSet);
		
		IntVec ivMCS = liIndexFragSolution.get(liIndexFragSolution.size()-1);
		
		molMCS = getSubFrag(frag, ivMCS);
		
		return molMCS;
	}


	private static StereoMolecule getSubFrag(StereoMolecule frag, IntVec iv){
		
		int bonds = frag.getBonds();
		
		HashSet<Integer> hsAtomIndex = new HashSet<Integer>();
		for (int i = 0; i < bonds; i++) {
			if(iv.isBitSet(i)){
				
				int indexAtom1 = frag.getBondAtom(0, i);
				int indexAtom2 = frag.getBondAtom(1, i);
				
				hsAtomIndex.add(indexAtom1);
				
				hsAtomIndex.add(indexAtom2);
			}
		}		
		
		StereoMolecule fragSubBonds = new StereoMolecule(hsAtomIndex.size(), bonds);
		
		fragSubBonds.setFragment(true);
		
		int [] arrMapAtom = new int [frag.getAtoms()];
		
		ArrayUtilsCalc.set(arrMapAtom, -1);
		
		for (int indexAtom : hsAtomIndex) {
			
			int indexAtomNew = fragSubBonds.addAtom(frag.getAtomicNo(indexAtom));
			
			fragSubBonds.setAtomX(indexAtomNew, frag.getAtomX(indexAtom));
			fragSubBonds.setAtomY(indexAtomNew, frag.getAtomY(indexAtom));
			fragSubBonds.setAtomZ(indexAtomNew, frag.getAtomZ(indexAtom));
			
			arrMapAtom[indexAtom]=indexAtomNew;
		}
		
		for (int i = 0; i < bonds; i++) {
			if(iv.isBitSet(i)){
				
				int indexAtom1 = frag.getBondAtom(0, i);
				int indexAtom2 = frag.getBondAtom(1, i);
				
				int indexAtomNew1 = arrMapAtom[indexAtom1];
				int indexAtomNew2 = arrMapAtom[indexAtom2];
				
				int type = frag.getBondType(i);
				
				if(frag.isDelocalizedBond(i)){
					type = Molecule.cBondTypeDelocalized;
					
					// fragSubBonds.setBondQueryFeature(bondIndexNew, Molecule.cBondQFDelocalized, true);
				}
				
				// int bondIndexNew = fragSubBonds.addBond(indexAtomNew1, indexAtomNew2, type);
				fragSubBonds.addBond(indexAtomNew1, indexAtomNew2, type);
				

			}
		}
		
		fragSubBonds.ensureHelperArrays(Molecule.cHelperRings);

		return fragSubBonds;
	}
	

	private List<IntVec> getAllPlusOneAtomCombinations(IntVec iv, StereoMolecule frag){
		
		int bonds = frag.getBonds();
		
		List<IntVec> liIntVec = new ArrayList<IntVec>();
		
		
		HashSet<Integer> hsAtomIndex = new HashSet<Integer>();
		
		for (int i = 0; i < bonds; i++) {
			
			if(iv.isBitSet(i)){
				
				int indexAt1 = frag.getBondAtom(0, i);
				
				int indexAt2 = frag.getBondAtom(1, i);
				
				hsAtomIndex.add(indexAt1);
				hsAtomIndex.add(indexAt2);
			}
		}
		
		for (int i = 0; i < bonds; i++) {
			if(!iv.isBitSet(i)){
				int indexAt1 = frag.getBondAtom(0, i);
				
				int indexAt2 = frag.getBondAtom(1, i);
				
				if(hsAtomIndex.contains(indexAt1) || hsAtomIndex.contains(indexAt2)) {
					IntVec ivPlus = new IntVec(iv.get());
					
					setBitAndAddRelatedRingBonds(i, ivPlus);
					
					liIntVec.add(ivPlus);
				}
			}
		}
		
		return liIntVec;
	}

	/**
	 * Removes all sub-solutions from hsIndexFragCandidates.
	 * @param ivSolution
	 */
	private void removeAllSubSolutions(IntVec ivSolution){
		
		List<IntVec> liIndexFragSolution = new ArrayList<IntVec>(hsIndexFragCandidates);
		
		for (IntVec ivCandidate : liIndexFragSolution) {
			if(isCandidateInSolution(ivSolution, ivCandidate)) {
				
				hsIndexFragCandidates.remove(ivCandidate);
			}
		}
	}
	
	/**
	 * All sub solutions are removed.
	 * @param hsIndexFragSolution
	 * @return
	 */
	private static List<IntVec> getFinalSolutionSet(HashSet<IntVec> hsIndexFragSolution ){
		
		List<IntVec> liIndexFragSolution = new ArrayList<IntVec>(hsIndexFragSolution);
		
		for (int i = liIndexFragSolution.size()-1; i >= 0; i--) {
			IntVec ivCandidate = liIndexFragSolution.get(i);
			
			for (int j = 0; j < liIndexFragSolution.size(); j++) {
				IntVec ivSolution = liIndexFragSolution.get(j);
				
				if(i!=j){
					if(isCandidateInSolution(ivSolution, ivCandidate)) {
						liIndexFragSolution.remove(i);
						break;
					}
				}
			}
		}
		
		return liIndexFragSolution;
	}
	
	private static final boolean isCandidateInSolution(IntVec ivSolution, IntVec ivCandidate){
		
		IntVec iv = IntVec.OR(ivSolution, ivCandidate);
		
		if(iv.equals(ivSolution)){
			return true;
		}
		
		return false;
	}
	
	private static class ComparatorBitsSet implements  Comparator<IntVec> {
		
		public int compare(IntVec iv1, IntVec iv2) {
			
			int bits1 = iv1.getBitsSet();
			
			int bits2 = iv2.getBitsSet();
			
			if(bits1>bits2){
				return 1;
			} else if(bits1<bits2){
				return -1;
			}
			
			return 0;
		}
		
	}
	
	public double getScore(){
		
		double sc = 0;
		
		double nAtmsFrag = frag.getBonds();
		
		double nAtmsMol = mol.getBonds();
		
		double nAtmsMCS = molMCS.getBonds();
		
		sc = nAtmsMCS/Math.max(nAtmsFrag, nAtmsMol);
		
		return sc;
	}

	public boolean isConsiderAromaticRings() {
		return considerAromaticRings;
	}

	public boolean isConsiderRings() {
		return considerRings;
	}
	
//	static class Solution  {
//		
//		IntVec iv;
//		
//		int [] arrMatchMolecule;
//		
//		
//		int hash;
//		
//		public Solution(IntVec iv, int [] arrMatchMolecule) {
//			
//			this.iv = new IntVec(iv);
//			
//			this.arrMatchMolecule = arrMatchMolecule;
//			
//			int [] a = iv.get();
//			
//			int [] arrHash = new int [a.length+arrMatchMolecule.length];
//			
//			System.arraycopy(a, 0, arrHash, 0, a.length);
//			
//			System.arraycopy(arrMatchMolecule, 0, arrHash, a.length, arrMatchMolecule.length);
//			
//			hash = new IntVec(arrHash).hashCode();
//			
//		}
//		
//		
//		public int hashCode() {
//			return hash;
//		}
//		
//		public boolean equals(Object obj) {
//			
//			if(!(obj instanceof Solution)) {
//				return false;
//			}
//			
//			Solution solution = (Solution)obj;
//			
//			if(!iv.equal(solution.iv)){
//				return false;
//			}
//			
//			if(!ArrayUtilsCalc.equals(arrMatchMolecule, solution.arrMatchMolecule)){
//				return false;
//			}
//			
//			return true;
//		}
//		
//
//	}
}
