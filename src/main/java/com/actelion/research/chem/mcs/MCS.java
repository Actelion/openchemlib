package com.actelion.research.chem.mcs;

import com.actelion.research.chem.*;
import com.actelion.research.util.datamodel.IntVec;

import java.util.*;

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
public class MCS {
	public static final int PAR_CLEAVE_RINGS = 0;
	public static final int PAR_KEEP_RINGS = 1;
	public static final int PAR_KEEP_AROMATIC_RINGS = 2;
	private static final boolean DEBUG = false;
	private static final int CAPACITY = (60 + 60) * 10;
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
	private HashMap<Integer, List<int[]>> hmRingBnd_ListRingBnds;
	private HashMap<Integer, List<int[]>> hmAromaticRingBnd_ListRingBnds;
	private RingCollection ringCollection;
	boolean excluded[] = null;
	private int [] arrMatchListFrag2Mol;

	public MCS() {
		this(PAR_CLEAVE_RINGS, null);
	}

	public MCS(int ringStatus) {
		this(ringStatus, null);
	}

	public MCS(int ringStatus, SSSearcher searcher) {
		considerRings = false;
		considerAromaticRings = false;
		switch (ringStatus) {
			case PAR_CLEAVE_RINGS:
				break;
			case PAR_KEEP_RINGS:
				considerRings = true;
				break;
			case PAR_KEEP_AROMATIC_RINGS:
				considerAromaticRings = true;
				break;
			default:
				break;
		}
		hsIndexFragCandidates = new HashSet<IntVec>(CAPACITY);
		hsIndexFragGarbage = new HashSet<IntVec>(CAPACITY);
		hsIndexFragSolution = new HashSet<IntVec>(CAPACITY);
		liMCSSolutions = new ArrayList<IntVec>(CAPACITY);
		if (searcher == null)
			sss = new SSSearcher();
		else
			sss = searcher;
		comparatorBitsSet = new ComparatorBitsSet();
		hmRingBnd_ListRingBnds = new HashMap<Integer, List<int[]>>();
		hmAromaticRingBnd_ListRingBnds = new HashMap<Integer, List<int[]>>();
	}

	public void setSSSearcher(SSSearcher sss) {
		this.sss = sss;
	}

	/**
	 * mol should contain equal or more bonds than frag.
	 *
	 * @param mol
	 * @param frag
	 */
	public void set(StereoMolecule mol, StereoMolecule frag) {
		set(mol, frag, null);
	}

	/**
	 * mol should contain equal or more bonds than frag.
	 * If frag contains more than one molecule only the biggest one is considered.
	 *
	 * @param mol
	 * @param frag
	 * @param excluded
	 */
	public void set(StereoMolecule mol, StereoMolecule frag, boolean excluded[]) {

		StereoMolecule fragBiggestSub = new StereoMolecule(frag);

		fragBiggestSub.ensureHelperArrays(Molecule.cHelperRings);

		fragBiggestSub.stripSmallFragments();

		fragBiggestSub.ensureHelperArrays(Molecule.cHelperRings);

		this.mol = mol;

		this.frag = fragBiggestSub;

		this.excluded = excluded;

		init();
	}

	private void init() {
		hsIndexFragCandidates.clear();
		hsIndexFragGarbage.clear();
		hsIndexFragSolution.clear();
		liMCSSolutions.clear();
		hmRingBnd_ListRingBnds.clear();
		hmAromaticRingBnd_ListRingBnds.clear();
		initCandidates();
	}

	private void initCandidates() {
		ringCollection = frag.getRingSet();
		int rings = ringCollection.getSize();
		for (int i = 0; i < rings; i++) {
			int[] arrIndexBnd = ringCollection.getRingBonds(i);
			for (int j = 0; j < arrIndexBnd.length; j++) {
				if (!hmRingBnd_ListRingBnds.containsKey(arrIndexBnd[j])) {
					hmRingBnd_ListRingBnds.put(arrIndexBnd[j], new ArrayList<int[]>());
				}
				List<int[]> li = hmRingBnd_ListRingBnds.get(arrIndexBnd[j]);
				li.add(arrIndexBnd);
			}
			if (ringCollection.isAromatic(i)) {
				for (int j = 0; j < arrIndexBnd.length; j++) {
					if (!hmAromaticRingBnd_ListRingBnds.containsKey(arrIndexBnd[j])) {
						hmAromaticRingBnd_ListRingBnds.put(arrIndexBnd[j], new ArrayList<int[]>());
					}
					List<int[]> li = hmAromaticRingBnd_ListRingBnds.get(arrIndexBnd[j]);
					li.add(arrIndexBnd);
				}
			}
		}
		// Array length for bit list.
		int nInts = (int) (((double) frag.getBonds() / Integer.SIZE) + ((double) (Integer.SIZE - 1) / Integer.SIZE));
		// The vector is used as bit vector.
		// For each bond one bit.
		for (int i = 0; i < frag.getBonds(); i++) {
			IntVec iv = new IntVec(nInts);
			setBitAndAddRelatedRingBonds(i, iv);
			hsIndexFragCandidates.add(iv);
		}
	}

	private void setBitAndAddRelatedRingBonds(int bit, IntVec iv) {
		if (!considerAromaticRings && !considerRings) {
			iv.setBit(bit);
		} else if (considerRings) {
			iv.setBit(bit);
			if (hmRingBnd_ListRingBnds.containsKey(bit)) {
				List<int[]> li = hmRingBnd_ListRingBnds.get(bit);
				for (int i = 0; i < li.size(); i++) {
					int[] a = li.get(i);
					for (int j = 0; j < a.length; j++) {
						iv.setBit(a[j]);
					}
				}
			}
		} else if (considerAromaticRings) {
			iv.setBit(bit);
			if (hmAromaticRingBnd_ListRingBnds.containsKey(bit)) {
				List<int[]> li = hmAromaticRingBnd_ListRingBnds.get(bit);
				for (int i = 0; i < li.size(); i++) {
					int[] a = li.get(i);
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
	 * The returned IntVec contains the bond information of frag. The IntVec is used bit wise. Each bit corresponds to
	 * one bond in the fragment. The information for the matching substructure in the molecule is not contained in the
	 * IntVec.
	 *
	 * @return
	 */
	private List<IntVec> getAllSolutionsForCommonSubstructures() {
		sss.setMolecule(mol);
		frag.setFragment(true);
		sss.setFragment(frag);
		//
		// The MCS is the complete fragment.
		//
		try {
			if (sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cMatchDBondToDelocalized, excluded) > 0) {
				molMCS = frag;
				List<IntVec> liIndexFragCandidates = new ArrayList<IntVec>(hsIndexFragCandidates);
				if (liIndexFragCandidates != null && !liIndexFragCandidates.isEmpty()) {
					IntVec iv = liIndexFragCandidates.get(0);
					iv.setBits(0, iv.sizeBits());
					iv.calculateHashCode();
					List<IntVec> li = new ArrayList<IntVec>();
					li.add(iv);
					return li;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		int maxSizeCandidates = 0;
		while (!hsIndexFragCandidates.isEmpty()) {
			List<IntVec> liIndexFragCandidates = new ArrayList<IntVec>(hsIndexFragCandidates);
			Collections.sort(liIndexFragCandidates, comparatorBitsSet);
			// Get largest mcs.
			IntVec iv = liIndexFragCandidates.get(liIndexFragCandidates.size() - 1);
			hsIndexFragCandidates.remove(iv);
			if (DEBUG) {
				System.out.println("Bits set " + iv.getBitsSet());
				if (iv.getBitsSet() == frag.getBonds()) {
					System.out.println("Full structure in iv.");
				}
			}
			StereoMolecule fragSub = getSubFrag(frag, iv);
			sss.setFragment(fragSub);
			if (sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode, excluded) > 0) {
				hsIndexFragSolution.add(iv);
				removeAllSubSolutions(iv);
				if (iv.getBitsSet() != frag.getBonds()) {
					List<IntVec> liIV = getAllPlusOneAtomCombinations(iv, frag);
					for (IntVec ivPlus : liIV) {
						if (DEBUG) {
							if (ivPlus.getBitsSet() == frag.getBonds()) {
								System.out.println("Full structure in ivPlus.");
							}
						}
						if ((!hsIndexFragGarbage.contains(ivPlus)) && (!hsIndexFragSolution.contains(ivPlus))) {
							hsIndexFragCandidates.add(ivPlus);
						}
					}
					if (DEBUG) {
						System.out.println("tsIndexFragCandidates " + hsIndexFragCandidates.size());
					}
				}
				if (maxSizeCandidates < hsIndexFragCandidates.size()) {
					maxSizeCandidates = hsIndexFragCandidates.size();
				}
			} else {
				hsIndexFragGarbage.add(iv);
			}
		}
		if (hsIndexFragSolution.size() == 0) {
			return null;
		}
		return getFinalSolutionSet(hsIndexFragSolution);
	}

	/**
	 * All molecules which are sub structures of an other molecule in the list are removed.
	 *
	 * @return
	 */
	public LinkedList<StereoMolecule> getAllCommonSubstructures() {
		List<IntVec> liIndexFragSolution = getAllSolutionsForCommonSubstructures();
		if (liIndexFragSolution == null) {
			return null;
		}
		Collections.sort(liIndexFragSolution, comparatorBitsSet);
		IntVec ivMCS = liIndexFragSolution.get(liIndexFragSolution.size() - 1);
		molMCS = getSubFrag(frag, ivMCS);
		List<StereoMolecule> li = new ArrayList<StereoMolecule>();
		for (IntVec iv : liIndexFragSolution) {
			li.add(getSubFrag(frag, iv));
		}
		return ExtendedMoleculeFunctions.removeSubStructures(li);
	}

	/**
	 * @return maximum common substructure at top of list or null of none common MCS was found.
	 */
	public StereoMolecule getMCS() {
		List<IntVec> liIndexFragSolution = getAllSolutionsForCommonSubstructures();
		if (liIndexFragSolution == null) {
			return null;
		}
		Collections.sort(liIndexFragSolution, comparatorBitsSet);
		IntVec ivMCS = liIndexFragSolution.get(liIndexFragSolution.size() - 1);
		molMCS = getSubFrag(frag, ivMCS);
		return molMCS;
	}

	/**
	 * Calculates the bond arrays for molecule and fragment for their maximum common substructure.
	 *
	 * @param arrBondMCSMol result array. Null or initialized with the length number of bonds in molecule.
	 * @param arrBondFrag   result array. Null or initialized with the length number of bonds in fragment.
	 * @return two arrays, in the first array the bits are set 'true' which corresponds to the bonds for maximum common substructure in the molecule.
	 * In the second are the bits set 'true' which corresponds to the bonds for maximum common substructure in the fragment.
	 */

    public boolean[][] getMCSBondArray(boolean[] arrBondMCSMol, boolean[] arrBondFrag) {
		boolean [][] arrBondMol_Result = new boolean [2][];
		List<IntVec> liIndexFragSolution = getAllSolutionsForCommonSubstructures();
		if(liIndexFragSolution==null){
			return null;
		}
		Collections.sort(liIndexFragSolution, comparatorBitsSet);
		// Get largest mcs.
		IntVec ivMCSLargest = liIndexFragSolution.get(liIndexFragSolution.size()-1);
		int bonds = frag.getBonds();
		if(arrBondFrag == null){
			arrBondFrag = new boolean [bonds];
		} else {
			Arrays.fill(arrBondFrag, false);
		}
		for (int i = 0; i < bonds; i++) {
			if(ivMCSLargest.isBitSet(i)){
				arrBondFrag[i]=true;
			}
		}		
		StereoMolecule fragSub = getSubFrag(frag, ivMCSLargest);
		sss.setFragment(fragSub);
		//
		// The substructure has to be searched in the molecule because the mapping indices are not contained in the
		// IntVec that is the solution from the MCS search.
		arrMatchListFrag2Mol = null;
		if(sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode,excluded)>0){
			ArrayList<int[]> liMatchSubFrag2Mol = sss.getMatchList();
			arrMatchListFrag2Mol = getMappedMatchListFrag2Mol(fragSub, liMatchSubFrag2Mol.get(0));
			int bondsMol = mol.getBonds();
			if(arrBondMCSMol==null) {
				arrBondMCSMol = new boolean [bondsMol];
			} else {
				Arrays.fill(arrBondMCSMol, false);
			}
			getBondArrayMolecule(arrMatchListFrag2Mol, arrBondFrag, arrBondMCSMol);	
		}
		arrBondMol_Result[0]=arrBondMCSMol;
		arrBondMol_Result[1]=arrBondFrag;
		return arrBondMol_Result;
	}

	/**
	 * delivers the match list indices from getMCSBondArray(...)
	 * index is the Fragment atom index, value is the matched atom index in molecule.
	 * @return null if no match was found.
	 */
	public int[] getArrMatchListFrag2Mol() {
		return arrMatchListFrag2Mol;
	}

	private int[] getMappedMatchListFrag2Mol(StereoMolecule fragSub, int[] arrMatchListSubFrag2Mol)
    {
		int [] arrMatchListFrag2Mol = new int [frag.getAtoms()];
//        SSSearcher sss = new SSSearcher();
		sss.setMol(fragSub, frag);
		sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode, null);
		ArrayList<int[]> liMatchSubFrag2Mol = sss.getMatchList();
		int [] arrMatchListSubFrag2Frag = liMatchSubFrag2Mol.get(0);
		HashMap<Integer, Integer> hmIndexSubFrag_IndexFrag = new HashMap<Integer, Integer>();
		for (int i = 0; i < arrMatchListSubFrag2Frag.length; i++) {
			hmIndexSubFrag_IndexFrag.put(i, arrMatchListSubFrag2Frag[i]);
		}
		for (int i = 0; i < arrMatchListSubFrag2Mol.length; i++) {
			int indexAtSubFragment = i;
			int indexAtFragment = hmIndexSubFrag_IndexFrag.get(indexAtSubFragment);
			arrMatchListFrag2Mol[indexAtFragment]=arrMatchListSubFrag2Mol[indexAtSubFragment];
		}
		return arrMatchListFrag2Mol;
	}
	
    private boolean[] getBondArrayMolecule(int[] arrMatchFragment2Mol, boolean[] arrBondMCSFrag, boolean[] arrBondMCSMol)
    {
		for (int i = 0; i < arrMatchFragment2Mol.length; i++) {
			int indexAtFrag = i;
			int indexAtMol = arrMatchFragment2Mol[i];
			int nConn2Frag = frag.getConnAtoms(indexAtFrag);
			for (int j = 0; j < nConn2Frag; j++) {
				int indexAtFragConn = frag.getConnAtom(indexAtFrag, j);
				int indexBondFrag = frag.getBond(indexAtFrag, indexAtFragConn);
				if(arrBondMCSFrag[indexBondFrag]){
					int indexAtMolConn = arrMatchFragment2Mol[indexAtFragConn];
					int indexBondMol = mol.getBond(indexAtMol, indexAtMolConn);
                    if (indexBondMol > -1) {
						arrBondMCSMol[indexBondMol]=true;
                    }
				}
			}
		}
		return arrBondMCSMol;
	}

// The original procedure missed atom charge, mass, etc, which is particularly annoying,
// when providing a custom SSSearcher(), which matches these properties. TLS 11Feb2021
	private static StereoMolecule getSubFrag(StereoMolecule frag, IntVec iv) {
		boolean[] isFragmentAtom = new boolean[frag.getAtoms()];
		int atoms = 0;
		int bonds = frag.getBonds();
		for (int bond=0; bond<bonds; bond++) {
			if (iv.isBitSet(bond)) {
				for (int i=0; i<2; i++) {
					int atom = frag.getBondAtom(i, bond);
					if (!isFragmentAtom[atom]) {
						isFragmentAtom[atom] = true;
						atoms++;
					}
				}
			}
		}

		StereoMolecule fragSubBonds = new StereoMolecule(atoms, bonds);
		fragSubBonds.setFragment(true);
		frag.copyMoleculeByAtoms(fragSubBonds, isFragmentAtom, true, null);
		fragSubBonds.ensureHelperArrays(Molecule.cHelperRings);
		return fragSubBonds;
	}

/*    private static StereoMolecule getSubFrag(StereoMolecule frag, IntVec iv)
    {
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
	}*/
	
    private List<IntVec> getAllPlusOneAtomCombinations(IntVec iv, StereoMolecule frag)
    {
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
     *
	 * @param ivSolution
	 */
    private void removeAllSubSolutions(IntVec ivSolution)
    {
		List<IntVec> liIndexFragSolution = new ArrayList<IntVec>(hsIndexFragCandidates);
		for (IntVec ivCandidate : liIndexFragSolution) {
			if(isCandidateInSolution(ivSolution, ivCandidate)) {
				hsIndexFragCandidates.remove(ivCandidate);
			}
		}
	}
	
	/**
	 * All sub solutions are removed.
     *
	 * @param hsIndexFragSolution
	 * @return
	 */
    private static List<IntVec> getFinalSolutionSet(HashSet<IntVec> hsIndexFragSolution)
    {
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
	
    private static final boolean isCandidateInSolution(IntVec ivSolution, IntVec ivCandidate)
    {
		IntVec iv = IntVec.OR(ivSolution, ivCandidate);
		if(iv.equals(ivSolution)){
			return true;
		}
		return false;
	}
	
    private static class ComparatorBitsSet implements Comparator<IntVec>
    {
        public int compare(IntVec iv1, IntVec iv2)
        {
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
	
	/**
	 * Calculates the score by bonds_mcs/Max(bonds_mol, bonds_frag)
	 * @return 1 if full overlap and 0 if no overlap at all.
	 */
	public double getScore(){
		
		double sc = 0;
		double nBndsFrag = frag.getBonds();
		
		double nBndsMol = mol.getBonds();
		
		double nBndsMCS = molMCS.getBonds();
		
		sc = nBndsMCS/Math.max(nBndsFrag, nBndsMol);
		
		return sc;
	}

    public boolean isConsiderAromaticRings()
    {
		return considerAromaticRings;
	}

    public boolean isConsiderRings()
    {
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
