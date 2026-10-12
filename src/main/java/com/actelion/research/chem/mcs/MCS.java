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
	private final HashSet<IntVec> hsIndexFragCandidates;
	private final HashSet<IntVec> hsIndexFragGarbage;
	private final HashSet<IntVec> hsIndexFragSolution;
	private SSSearcher sss;
	private final ComparatorBitsSet comparatorBitsSet;
	private StereoMolecule molMCS;
	private boolean considerAromaticRings;
	private boolean considerRings;
	// The largest valid solutions.
	private final HashMap<Integer, List<int[]>> hmRingBnd_ListRingBnds;
	private final HashMap<Integer, List<int[]>> hmAromaticRingBnd_ListRingBnds;
	private boolean[] excluded = null;
	private int[] arrMatchListFrag2Mol;

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
	public void set(StereoMolecule mol, StereoMolecule frag, boolean[] excluded) {
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
		hmRingBnd_ListRingBnds.clear();
		hmAromaticRingBnd_ListRingBnds.clear();
		initCandidates();
	}

	private void initCandidates() {
		RingCollection ringCollection = frag.getRingSet();
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
				List<IntVec> liIndexFragCandidates = new ArrayList<>(hsIndexFragCandidates);
				if (!liIndexFragCandidates.isEmpty()) {
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
			StereoMolecule fragSub = getSubFrag(frag, iv, false).fragment;
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
		if (hsIndexFragSolution.isEmpty()) {
			return null;
		}
		return getFinalSolutionSet(hsIndexFragSolution);
	}

	/**
	 * All molecules which are sub structures of another molecule in the list are removed.
	 *
	 * @return
	 */
	public LinkedList<StereoMolecule> getAllCommonSubstructures() {
		List<IntVec> liIndexFragSolution = getAllSolutionsForCommonSubstructures();
		if (liIndexFragSolution == null) {
			return null;
		}
		liIndexFragSolution.sort(comparatorBitsSet);
		IntVec ivMCS = liIndexFragSolution.get(liIndexFragSolution.size() - 1);
		molMCS = getSubFrag(frag, ivMCS, false).fragment;
		List<StereoMolecule> li = new ArrayList<StereoMolecule>();
		for (IntVec iv : liIndexFragSolution) {
			li.add(getSubFrag(frag, iv, false).fragment);
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
		liIndexFragSolution.sort(comparatorBitsSet);
		IntVec ivMCS = liIndexFragSolution.get(liIndexFragSolution.size() - 1);
		molMCS = getSubFrag(frag, ivMCS, false).fragment;
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
		liIndexFragSolution.sort(comparatorBitsSet);
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
		MappedFragment mf = getSubFrag(frag, ivMCSLargest, true);
		sss.setFragment(mf.fragment);
		//
		// The substructure has to be searched in the molecule because the mapping indices are not contained in the
		// IntVec that is the solution from the MCS search.
		arrMatchListFrag2Mol = null;
		if(sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode,excluded)>0){
			ArrayList<int[]> liMatchSubFrag2Mol = sss.getMatchList();
			int [] mapSubFrag2Mol = liMatchSubFrag2Mol.get(0);
			arrMatchListFrag2Mol = new int[frag.getAtoms()];
			Arrays.fill(arrMatchListFrag2Mol, -1);
			for (int indexSubFrag = 0; indexSubFrag < mapSubFrag2Mol.length; indexSubFrag++) {
				for (int indexFrag2SubFrag = 0; indexFrag2SubFrag < mf.map.length; indexFrag2SubFrag++) {
					if(mf.map[indexFrag2SubFrag]==indexSubFrag){
						arrMatchListFrag2Mol[indexFrag2SubFrag]=mapSubFrag2Mol[indexSubFrag];
						break;
					}
				}
			}

			int bondsMol = mol.getBonds();
			if(arrBondMCSMol==null) {
				arrBondMCSMol = new boolean [bondsMol];
			} else {
				Arrays.fill(arrBondMCSMol, false);
			}
			getBondArrayMolecule(arrMatchListFrag2Mol, arrBondMCSMol);
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

	private int[] getMappedMatchListFrag2Mol(MappedFragment mf, int[] arrMatchListSubFrag2Mol) {
		int [] arrMatchListFrag2Mol = new int [frag.getAtoms()];
		Arrays.fill(arrMatchListFrag2Mol, -1);

		sss.setMol(mf.fragment, frag);
		sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode, null);
		ArrayList<int[]> liMatchSubFrag2Mol = sss.getMatchList();
		int [] arrMatchListSubFrag2Frag = liMatchSubFrag2Mol.get(0);
		HashMap<Integer, Integer> hmIndexSubFrag_IndexFrag = new HashMap<Integer, Integer>();
		for (int i = 0; i < arrMatchListSubFrag2Frag.length; i++) {
			hmIndexSubFrag_IndexFrag.put(i, arrMatchListSubFrag2Frag[i]);
		}
		for (int indexAtSubFragment=0; indexAtSubFragment< arrMatchListSubFrag2Mol.length; indexAtSubFragment++) {
			int indexAtFragment = hmIndexSubFrag_IndexFrag.get(indexAtSubFragment);
			arrMatchListFrag2Mol[indexAtFragment]=arrMatchListSubFrag2Mol[indexAtSubFragment];
		}
		return arrMatchListFrag2Mol;
	}


	/**
	 *
	 * @param arrMatchFragment2Mol Index is the Fragment atom index, value is the matched atom index in molecule.
	 * @param arrBondMCSMol The bits are set 'true' which corresponds to the bonds for maximum common substructure in the molecule.
	 * @return
	 */
	private void getBondArrayMolecule(int[] arrMatchFragment2Mol, boolean[] arrBondMCSMol) {
		for (int indexAtFrag1=0; indexAtFrag1<arrMatchFragment2Mol.length; indexAtFrag1++) {
			int indexAtMol1 = arrMatchFragment2Mol[indexAtFrag1];
			if (indexAtMol1 == -1)
				continue;
			for (int indexAtFrag2 = indexAtFrag1+1; indexAtFrag2 < arrMatchFragment2Mol.length; indexAtFrag2++) {

				int indexAtMol2 = arrMatchFragment2Mol[indexAtFrag2];
				if (indexAtMol2 == -1)
					continue;

				int indexBondMol = mol.getBond(indexAtMol1, indexAtMol2);
				if (indexBondMol > -1) {
					arrBondMCSMol[indexBondMol]=true;
				}
			}
		}
	}


// The original procedure missed atom charge, mass, etc., which is particularly annoying,
// when providing a custom SSSearcher(), which matches these properties. TLS 11Feb2021
	private static MappedFragment getSubFrag(StereoMolecule frag, IntVec iv, boolean withMap) {
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

		int [] map = null;
		if(withMap)
			map = new int[frag.getAtoms()];

		StereoMolecule fragSubBonds = new StereoMolecule(atoms, bonds);
		fragSubBonds.setFragment(true);
		frag.copyMoleculeByAtoms(fragSubBonds, isFragmentAtom, true, map);
		fragSubBonds.ensureHelperArrays(Molecule.cHelperRings);

		MappedFragment mf = new MappedFragment(fragSubBonds, map);

		return mf;
	}


    private List<IntVec> getAllPlusOneAtomCombinations(IntVec iv, StereoMolecule frag) {
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
    private void removeAllSubSolutions(IntVec ivSolution) {
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
    private static List<IntVec> getFinalSolutionSet(HashSet<IntVec> hsIndexFragSolution){
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
	
    private static boolean isCandidateInSolution(IntVec ivSolution, IntVec ivCandidate){
		IntVec iv = IntVec.OR(ivSolution, ivCandidate);
		if(iv.equals(ivSolution)){
			return true;
		}
		return false;
	}
	
    private static class ComparatorBitsSet implements Comparator<IntVec> {
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

    public boolean isConsiderRings(){
		return considerRings;
	}

	private static class MappedFragment  {
		StereoMolecule fragment;
		// Index is atom index in fragment, value is atom index in original fragment.
		int [] map;

		public MappedFragment(StereoMolecule fragment, int[] map) {
			this.fragment = fragment;
			this.map = map;
		}
	}
}
