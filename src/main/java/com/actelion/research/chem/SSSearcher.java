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
 * @author Thomas Sander
 */

package com.actelion.research.chem;

import com.actelion.research.util.IntArrayComparator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

public class SSSearcher {
	// CONSTANTS TO DEFINE KIND OF SIMILARITY BETWEEN ATOMS AND BONDS
	public final static int cMatchAtomCharge = 1;
	public final static int cMatchAtomMass = 2;
	public final static int cMatchDBondToDelocalized = 4;
	public final static int cMatchAromDBondToDelocalized = 8;

/*	For index match modes we need to consider the following:
  - If fragment C is SS of fragment B and B is SS of molecule A then C must be SS of A.
	This implies that if X is SS of Y then X is in all respects equal or less exactly
	defined as Y. Consequently, if the SS operator matches any double bonds to aromatic
	double bonds and matches aromatic double bonds to delocaliced bonds, it must match
	any double bond to delocalized bonds (cMatchAromDBondToDelocalized cannot be used for
	index creation). Example: key C=C-N-C=C, query pyrol, molecule indole, key would match
	pyrol but not indol!!!
  - match modes used for the actual atom by atom check must be more or equally restrictive
	than the match mode used for index creation. Otherwise, index keys may filter out
	molecules which would be considered a match with the less strict matching conditions
	of the atom by atom check.
*/

									// match mode to be used for index creation
	public static final int cIndexMatchMode = cMatchDBondToDelocalized;
	public static final int cDefaultMatchMode = cMatchAromDBondToDelocalized;

	public final static int cCountModeExistence     = 1;    // check only, don't create matchList
	public final static int cCountModeFirstMatch	= 2;    // create matchList with just one match
	public final static int cCountModeSeparated		= 3;    // create list of all non-overlapping matches / not optimized for maximum match count
	public final static int cCountModeOverlapping	= 4;    // create list not containing multiple matches sharing exactly the same atoms
	public final static int cCountModeRigorous		= 5;    // create list of all possible matches neglecting any symmetries
	public final static int cCountModeUnique		= 6;    // create list of all distinguishable matches considering symmetries

	// default behaviour for unusual atom masses and atom charges is that
	// - if no atom charge/mass is specified in the query then all charges/masses match
	// - if an atom charge/mass is specified then this charge/mass must match for the atom to match

	// defines details of similarity between atoms and bonds
	private int mDefaultMatchMode;

	// The molecule which is analyzed
	protected StereoMolecule mMolecule;

	// The sub-structure we try to find in the molecule
	protected StereoMolecule mFragment;

	private int[] mMoleculeAtomType;	// atom features required to match
	private int[] mFragmentAtomType;
	private long[] mMoleculeAtomFeatures;	// flags defining given/required atom features
	private long[] mFragmentAtomFeatures;
	private long[] mMoleculeRingFeatures;	// flags defining given/required atom ring size features
	private long[] mFragmentRingFeatures;
	private int[] mMoleculeBondFeatures;	// flags defining given/required bond features
	private int[] mFragmentBondFeatures;

	private int mFragmentExcludeAtoms,mFragmentExcludeBonds;
	private int mFragmentGraphSize;	// the number of wanted atoms & ring closures in fragment graph
	private int mFragmentGraphSizeWithExcludeGroups;	// total number of atoms & ring closures in fragment graph
	private int[] mFragmentGraphAtom;
	private int[] mFragmentGraphParentAtom;
	private int[] mFragmentGraphParentBond;
	private boolean[] mFragmentGraphIsRingClosure;
	private boolean[] mIsExcludeAtom;
	private boolean[] mIsBridgeBondAtom;
	private int[] mFragmentConnAtoms;	// in case of exclude atoms, these are not part of this
	private int[] mMatchTable;
	private int[] mExcludeGroupNo;
	private int[] mExcludeGroupGraphIndex;
	private int[] mFragmentAtomContextRank;

	// depending on the fragment count mode this may contain atom lists
	// of all till now located matching sub-fragments
	private TreeSet<int[]> mSortedMatchSet;
	private ArrayList<int[]> mMatchList;
	private ArrayList<BridgeBond> mBridgeBondList;
	private ArrayList<boolean[]> mBridgeBondAtomList;

	private boolean mMoleculeFeaturesValid;
	private boolean mFragmentFeaturesValid;
	private int mRequiredHelperLevel;
	private int mExcludeGroupCount;

	private volatile boolean mStop;

	/**
	 * Instantiates a SSSearcher object for running sub-structure searches
	 * with one or more sub-structure fragments on one or more molecules.
	 * The search is a pure graph matching algorithm.
	 * For fast sub-structure searches involving an index based pre-screening use the class SSSearcherWithIndex.
	 * For a more high-level structure search supporting multiple cores, sub-structure-, similarity-,
	 * exact-, or tautomer-search use class StructureSearch and related classes.
	 */
	public SSSearcher() {
		mDefaultMatchMode = cDefaultMatchMode;
		mSortedMatchSet = new TreeSet<>(new IntArrayComparator());
		}


	/**
	 * Instantiates a SSSearcher object for running sub-structure searches
	 * with one or more sub-structure fragments on one or more molecules.
	 * The search is a pure graph matching algorithm.
	 * For fast sub-structure searches involving an index based pre-screening use the class SSSearcherWithIndex.
	 * For a more high-level structure search supporting multiple cores, sub-structure-, similarity-,
	 * exact-, or tautomer-search use class StructureSearch and related classes.
	 * @param matchMode combination of cMatchAtomCharge, cMatchAtomMass, cMatchDBondToDelocalized, cMatchAromDBondToDelocalized
	 */
	public SSSearcher(int matchMode) {
		mDefaultMatchMode = matchMode;
		mSortedMatchSet = new TreeSet<>(new IntArrayComparator());
		}


	/**
	 * Defines fragment and molecule before calling isFragmentInMolecule(...)
	 * or findFragmentInMolecule(...).
	 * @param fragment
	 * @param molecule
	 */
	public void setMol(StereoMolecule fragment, StereoMolecule molecule) {
		setMolecule(molecule);
		setFragment(fragment);
		}


	/**
	 * Defines the molecule to be used in isFragmentInMolecule(...)
	 * or findFragmentInMolecule(...).
	 * @param molecule
	 */
	public void setMolecule(StereoMolecule molecule) {
		if (molecule == null || molecule.getAllAtoms() == 0) {
			mMolecule = null;
			return;
			}

		mMolecule = molecule;
		mMoleculeFeaturesValid = false;
		mMolecule.ensureHelperArrays(Molecule.cHelperNeighbours);
		}

	/**
	 * Asks the substructure search to stop without completing as soon as possible.
	 */
	public void stop() {
		mStop = true;
		}

	/**
	 * Defines the fragment to be used in isFragmentInMolecule(...)
	 * or findFragmentInMolecule(...).
	 * @param fragment
	 */
	public void setFragment(StereoMolecule fragment) {
		if (fragment == null || fragment.getAllAtoms() == 0 || !fragment.isFragment()) {
			mFragment = null;
			return;
			}

		mFragment = fragment;
		mFragmentFeaturesValid = false;
		mFragment.ensureHelperArrays(Molecule.cHelperNeighbours);

		mRequiredHelperLevel = Molecule.cHelperRings;
		for (int atom=0; atom<mFragment.getAtoms(); atom++)
			if ((mFragment.getAtomQueryFeatures(atom) & (Molecule.cAtomQFStereoState | Molecule.cAtomQFMatchStereo)) != 0)
				mRequiredHelperLevel = Molecule.cHelperParities;
		for (int bond=0; bond<mFragment.getBonds(); bond++)
			if ((mFragment.getBondQueryFeatures(bond) & Molecule.cBondQFMatchStereo) != 0)
				mRequiredHelperLevel = Molecule.cHelperParities;

		if (mMoleculeFeaturesValid && mRequiredHelperLevel != Molecule.cHelperRings)
			mMolecule.ensureHelperArrays(mRequiredHelperLevel);

		mFragmentExcludeAtoms = 0;
		mFragmentExcludeBonds = 0;

		mIsExcludeAtom = new boolean[mFragment.getAtoms()];
		for (int atom=0; atom<mFragment.getAtoms(); atom++) {
			mIsExcludeAtom[atom] = ((mFragment.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) != 0);
			if (mIsExcludeAtom[atom])
				mFragmentExcludeAtoms++;
			}

		mExcludeGroupCount = 0;
		mExcludeGroupNo = null;

		mFragmentAtomContextRank = null;

		if (mFragmentExcludeAtoms != 0) {
			if (mFragmentExcludeAtoms != 0)
				for (int bond = 0; bond < mFragment.getBonds(); bond++)
					if (mIsExcludeAtom[mFragment.getBondAtom(0, bond)] || mIsExcludeAtom[mFragment.getBondAtom(1, bond)])
						mFragmentExcludeBonds++;

			// find all independent exclude groups
			for (int atom=0; atom<mFragment.getAllAtoms(); atom++)
				mFragment.setAtomMarker(atom, mIsExcludeAtom[atom]);

			mExcludeGroupNo = new int[mFragment.getAllAtoms()];
			mExcludeGroupCount = mFragment.getFragmentNumbers(mExcludeGroupNo, true, false);
			}
		}


	/**
	 * If countMode is cCountModeUnique, then matches are considered distinct, if<br>
	 * - either the list of matching molecule atoms to the query fragment is a different one<br>
	 * - or if the mutual combination of fragment and molecule atom's symmetry rank is different
	 *   (when the list of matched molecule atoms is the same)<br>
	 * For certain situations fragment atoms must be considered different, even if their symmetry rank
	 * is equal, e.g. in the context of a reaction where equivalent reactant atoms end up in different
	 * product environments. This method allows to specify an additional criterion for the uniqueness
	 * comparison to be considered. In case of reactions, these might be the symmetry ranks of the
	 * products atoms mapped to the reactant atoms.
	 * Note: The current implementation only uses the 8 least significant bits of the context rank.
	 */
	public void setFragmentSymmetryConstraints(int[] fragmentContextRank) {
		mFragmentAtomContextRank = fragmentContextRank;
		}

	/**
	 * Build a graph of the query fragment(s) including ring closures as redundant nodes.
	 * If we have exclude groups then these are added at the end of the graph.
	 */
	private void buildFragmentGraph() {
		mFragment.ensureHelperArrays(mRequiredHelperLevel);

		// ringClosures = bonds - atoms + fragments;
		// extreme cases: highly bridged multicycle, e.g. ikosaeder; many exclude atoms (many atoms)
		int graphAllocation = Math.max(mFragment.getAtoms(), mFragment.getBonds()) + 16;

		mFragmentGraphAtom = new int[graphAllocation];
		mFragmentGraphParentAtom = new int[graphAllocation];
		mFragmentGraphParentBond = new int[graphAllocation];
		mFragmentGraphIsRingClosure = new boolean[graphAllocation + 1];

		boolean[] fragmentAtomUsed = new boolean[mFragment.getAtoms()];
		boolean[] fragmentBondUsed = new boolean[mFragment.getBonds()];
		int current = 0;
		for (int atom=0; atom<mFragment.getAtoms(); atom++) {
			if (!mIsExcludeAtom[atom]
			 && !fragmentAtomUsed[atom]) {
				mFragmentGraphAtom[current] = atom;
				mFragmentGraphParentBond[current] = -1;
				mFragmentGraphParentAtom[current] = -1;
				int highest = current;
				while (current <= highest) {
					for (int i=0; i<mFragment.getAllConnAtomsPlusMetalBonds(mFragmentGraphAtom[current]); i++)
						highest = tryAddCandidate(current, highest, i, fragmentAtomUsed, fragmentBondUsed, -1);
					while (mFragmentGraphIsRingClosure[++current]);
					}
				}
			}
		mFragmentGraphSize = current;	// this is the real size of the graph not considering exclude atoms

		if (mFragmentExcludeAtoms != 0) {
			// append all exclude groups including ring closures, one after another, to non-excluded atom tree
			int highest = mFragmentGraphSize - 1;
			for (int excludeGroupNo=0; excludeGroupNo<mExcludeGroupCount; excludeGroupNo++) {
				current = 0;
				while (current <= highest) {
					for (int i=0; i<mFragment.getAllConnAtomsPlusMetalBonds(mFragmentGraphAtom[current]); i++)
						highest = tryAddCandidate(current, highest, i, fragmentAtomUsed, fragmentBondUsed, excludeGroupNo);
					while (mFragmentGraphIsRingClosure[++current]) ;
					}
				}

			// there may still be exclude groups as separated fragments
			for (int atom=0; atom<mFragment.getAtoms(); atom++) {
				if (mIsExcludeAtom[atom] && !fragmentAtomUsed[atom]) {
					mFragmentGraphAtom[current] = atom;
					mFragmentGraphParentBond[current] = -1;
					mFragmentGraphParentAtom[current] = -1;
					highest = current;
					while (current <= highest) {
						for (int i=0; i<mFragment.getAllConnAtomsPlusMetalBonds(mFragmentGraphAtom[current]); i++)
							if (mFragment.getConnAtom(mFragmentGraphAtom[current], i) < mFragment.getAtoms())
								highest = tryAddCandidate(current, highest, i, fragmentAtomUsed, fragmentBondUsed, mExcludeGroupNo[atom]);
						while (mFragmentGraphIsRingClosure[++current]) ;
						}
					}
				}

			mExcludeGroupGraphIndex = new int[mExcludeGroupCount];
			for (int i=0; i<mExcludeGroupCount; i++)
				mExcludeGroupGraphIndex[i] = -1;
			for (int i=mFragmentGraphSize; i<current; i++) {
				int excludeGroupNo = mExcludeGroupNo[mFragmentGraphAtom[i]];
				if (mExcludeGroupGraphIndex[excludeGroupNo] == -1)
					mExcludeGroupGraphIndex[excludeGroupNo] = i;
				}
			}

		mFragmentGraphSizeWithExcludeGroups = current;	// this is the real size of the graph
/*
System.out.print("			"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(mFragmentGraphAtom[i]==-1?"-":Molecule.cAtomLabel[mFragment.getAtomicNo(mFragmentGraphAtom[i])])); System.out.println();
System.out.print("  graphAtom:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(mFragmentGraphAtom[i]==-1?"-":""+mFragmentGraphAtom[i])); System.out.println();
System.out.print(" parentAtom:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(mFragmentGraphParentAtom[i]==-1?"-":""+mFragmentGraphParentAtom[i])); System.out.println();
System.out.print(" parentBond:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(mFragmentGraphParentBond[i]==-1?"-":""+mFragmentGraphParentBond[i])); System.out.println();
System.out.print("ringClosure:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(mFragmentGraphIsRingClosure[i]?"y":"n")); System.out.println();
System.out.println();
*/
		}


	/**
	 * Considers the i-th neighbour of the current graph atom as potential next graph member to add. If<br>
	 * - it is not equal to the parent of the current graph atom and<br>
	 * - if the bond to the candidate is not a bridge bond<br>
	 * - if the bond was not added to the graph already<br>
	 * then attach the candidate to the graph and mark, whether it is a ring closure (and therefore duplicate).
	 * @param current
	 * @param highest
	 * @param i connAtom index of current graph atom
	 * @param fragmentAtomUsed
	 * @param fragmentBondUsed
	 * @return
	 */
	private int tryAddCandidate(int current, int highest, int i, boolean[] fragmentAtomUsed, boolean[] fragmentBondUsed, int excludeGroupNo) {
		int candidate = mFragment.getConnAtom(mFragmentGraphAtom[current], i);
		if ((!mIsExcludeAtom[candidate] || mExcludeGroupNo[candidate] == excludeGroupNo)	// always allow non-exclude atoms, because it may be a ring closure from exclude group to main fragment
		 && candidate != mFragmentGraphParentAtom[current]) {
			int candidateBond = mFragment.getConnBond(mFragmentGraphAtom[current], i);

			if (!fragmentBondUsed[candidateBond]	// if it is a ring closure make sure it is added only once
			 && !mFragment.isBondBridge(candidateBond)) {	// don't consider bridge bonds at this state
				mFragmentGraphAtom[++highest] = candidate;
				mFragmentGraphParentAtom[highest] = mFragmentGraphAtom[current];
				mFragmentGraphParentBond[highest] = candidateBond;
				fragmentBondUsed[candidateBond] = true;
				if (fragmentAtomUsed[candidate])
					mFragmentGraphIsRingClosure[highest] = true;
				else
					fragmentAtomUsed[candidate] = true;
				}
			}
		return highest;
		}

	/**
	 * If the match count mode is one of cCountModeFirstMatch, cCountModeOverlapping,
	 * cCountModeRigorous then this method returns an arraylist of all counted matches,
	 * i.e. int arrays mapping fragment atoms to molecule atoms. Atoms being part of a
	 * matched bridge bond are naturally not covered by the mapping. Atoms being part of a
	 * matching bridge bond are available with getBridgeBondAtomList().<br>
	 * Note: If some query fragment atoms are marked as exclude group, then the respective
	 * matchlist values are -1.
	 * @return list of distinct counted matches.
	 */
	public ArrayList<int[]> getMatchList() {
		return mMatchList;
		}

	/**
	 * getMatchList() doesn't include information about atoms, which are part of a matching bridge bond.
	 * This method returns an atom mask for a given matchNo, where all atoms are flagged that are part of a
	 * matching bridge bond within that match.
	 * If multiple bridges bond matches are possible, for every bridge bond only the shortest bridge is considered.
	 * Multiple bridge matches don't contribute to the multiplicity of match lists, nor are they considered else where.
	 * @param matchNo index of corresponding match from getMatchList()
	 * @return null or atom mask in target atom space
	 */
	public boolean[] getMatchingBridgeBondAtoms(int matchNo) {
		return mBridgeBondAtomList.size() <= matchNo ? null : mBridgeBondAtomList.get(matchNo);
		}

	/**
	 * Fastest check, whether the molecule contains the fragment.
	 * This method uses cCountModeExistance and therefore does not create
	 * any retrievable match list, i.e. mapping from fragment atoms to molecule atoms.
	 * The match mode used is cDefaultMatchMode, unless defined otherwise when
	 * instantiating the SSSearcher.
	 * @return whether fragment was found as sub-structure in molecule
	 */
	public boolean isFragmentInMolecule() {
		return (findFragmentInMolecule(cCountModeExistence, mDefaultMatchMode) > 0);
		}


	/**
	 * Fastest check, whether the molecule contains the fragment.
	 * This method uses cCountModeExistance and therefore does not create
	 * any retrievable match list, i.e. mapping from fragment atoms to molecule atoms.
	 * @param matchMode cDefaultMatchMode or combination of cMatchAtomCharge, cMatchAtomMass, cMatchDBondToDelocalized, cMatchAromDBondToDelocalized
	 * @return whether fragment was found as sub-structure in molecule
	 */
	public boolean isFragmentInMolecule(int matchMode) {
		return (findFragmentInMolecule(cCountModeExistence, matchMode) > 0);
		}

	/**
	 * Locates all matches of the fragment in the molecule that result in distinguishable
	 * sets of molecule atoms. Multiple matches involving the same atoms, e.g. with a benzene ring,
	 * are counted and listed only once. Atom mapping from fragment to molecule
	 * is collected and can be retrieved with getMatchList().
	 * @return count of sub-structure matches of fragment in molecule
	 */
	public int findFragmentInMolecule() {
		return findFragmentInMolecule(cCountModeOverlapping, mDefaultMatchMode);
		}


	/**
	 * Locates all matches of the fragment in the molecule that result in distinguishable
	 * sets of molecule atoms. Multiple matches involving the same atoms, e.g. with a benzene ring,
	 * are counted and listed only once. If count mode is different from cCountModeExistance,
	 * then an atom mapping from fragment to molecule is collected and can be retrieved with getMatchList().
	 * @param countMode one of cCountModeExistance, cCountModeFirstMatch, cCountModeOverlapping, cCountModeRigorous
	 * @param matchMode cDefaultMatchMode or combination of cMatchAtomCharge, cMatchAtomMass, cMatchDBondToDelocalized, cMatchAromDBondToDelocalized
	 * @return count of sub-structure matches of fragment in molecule
	 */
	public int findFragmentInMolecule(int countMode, int matchMode) {
		return findFragmentInMolecule(countMode, matchMode, null);
		}


	/**
	 * Locates all matches of the fragment in the molecule that result in distinguishable
	 * sets of molecule atoms that are not flagged to be excluded from matching.
	 * Multiple matches involving the same atoms, e.g. with a benzene ring,
	 * are counted and listed only once. If count mode is different from cCountModeExistance,
	 * then an atom mapping from fragment to molecule is collected and can be retrieved with getMatchList().
	 * If the query fragment does not contain atoms other than exclude group atoms, then no match is returned.
	 * @param countMode one of cCountModeExistance, cCountModeFirstMatch, cCountModeSeparated, cCountModeOverlapping, cCountModeUnique, cCountModeRigorous
	 * @param matchMode cDefaultMatchMode or combination of cMatchAtomCharge, cMatchAtomMass, cMatchDBondToDelocalized, cMatchAromDBondToDelocalized
	 * @param atomExcluded defines atoms of molecule to be excluded from sub-structure matching
	 * @return count of sub-structure matches of fragment in molecule
	 */
	public int findFragmentInMolecule(int countMode, int matchMode, final boolean[] atomExcluded) {
		mStop = false;
		mMatchList = new ArrayList<>();
		mBridgeBondAtomList = new ArrayList<>();
		mSortedMatchSet.clear();

		if (mMolecule == null
   		 || mFragment == null)
   			return 0;

		if (mFragment.getAtoms() - mFragmentExcludeAtoms > mMolecule.getAtoms()
		 || mFragment.getBonds() - mFragmentExcludeBonds > mMolecule.getBonds())
			return 0;

		if (mFragment.getAtoms() == 0)
			return 0;

/*
System.out.print("  molecule:"); for (int i=0; i<moleculeAtoms; i++) System.out.print(" "+i); System.out.println();
System.out.print("	 label:"); for (int i=0; i<moleculeAtoms; i++) System.out.print(" "+Molecule.cAtomLabel[mMolecule.getAtomicNo(i)]); System.out.println();
for (int j=0; j<moleculeAtoms; j++)
{ System.out.print("conns of "+j+":"); for (int i=0; i<mMolecule.getConnAtoms(j); i++) System.out.print(" "+mMolecule.getConnAtom(j, i)); System.out.println(); }
System.out.println();
*/

		if (countMode == cCountModeUnique)
			mRequiredHelperLevel = Molecule.cHelperSymmetrySimple;

		setupAtomAndBondFeatures(matchMode);

		// atom usage mask in mMolecule
		boolean[] atomUsed = new boolean[mMolecule.getAtoms()];
		if (atomExcluded != null)
			for (int atom=0; atom<mMolecule.getAtoms(); atom++)
				atomUsed[atom] = atomExcluded[atom];

		// mMolecule atom currently matched on mFragment atom
		mMatchTable = new int[mFragment.getAtoms()];
		Arrays.fill(mMatchTable, -1);	// to mark exclude group atoms

		int[] index = new int[mFragmentGraphSizeWithExcludeGroups];
		Arrays.fill(index, -1);
		// contains current molecule atom pointer for graph matching,
		// - in case of sub fragment anchor atom: the current molecule atom index matched to the anchor
		// - otherwise the current connAtom index of the parent atom in the matching graph

		int current = 0;
		while (!mStop) {
/*
System.out.print("  index:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(index[i]==-1?"-":""+index[i])); System.out.println();
System.out.print("		"); for (int i=0; i<current; i++) System.out.print("  "); System.out.println(" ^");
System.out.println();
*/

			if (mFragmentGraphSize != 0) {	// we may have exclude groups only
				int maxIndex = (mFragmentGraphParentAtom[current] == -1) ? mMolecule.getAtoms()
						: mMolecule.getAllConnAtomsPlusMetalBonds(mMatchTable[mFragmentGraphParentAtom[current]]);

				index[current]++;

				if (index[current] == maxIndex) {
					index[current] = -1;
					if (current == 0)
						break;
					current--;
					if (!mFragmentGraphIsRingClosure[current])
						atomUsed[mMatchTable[mFragmentGraphAtom[current]]] = false;
					continue;
					}

				if (mFragmentGraphParentAtom[current] == -1) {	// if current graph atom is sub fragment anchor atom
					if (!atomUsed[index[current]]) {
						if (areAtomsSimilar(index[current], mFragmentGraphAtom[current])) {
							mMatchTable[mFragmentGraphAtom[current]] = index[current];
							atomUsed[index[current]] = true;
							current++;
							}
						}
					}
				else {
					// skip plain hydrogens
					if (mMolecule.getConnAtom(mMatchTable[mFragmentGraphParentAtom[current]], index[current]) >= mMolecule.getAtoms())
						continue;

					int candidate = mMolecule.getConnAtom(mMatchTable[mFragmentGraphParentAtom[current]], index[current]);
					if (!mFragmentGraphIsRingClosure[current]) {	// current graph position is not an anchor
						if (!atomUsed[candidate]) {
							if (areAtomsSimilar(candidate, mFragmentGraphAtom[current])
									&& areBondsSimilar(mMolecule.getConnBond(mMatchTable[mFragmentGraphParentAtom[current]], index[current]), mFragmentGraphParentBond[current])) {
								atomUsed[candidate] = true;
								mMatchTable[mFragmentGraphAtom[current]] = candidate;
								current++;
								}
							}
						}
					else {	// current graph position is ringClosure
						if (candidate == mMatchTable[mFragmentGraphAtom[current]]
								&& areBondsSimilar(mMolecule.getConnBond(mMatchTable[mFragmentGraphParentAtom[current]], index[current]), mFragmentGraphParentBond[current])) {
							current++;
							}
						}
					}
				}

			if (current == mFragmentGraphSize) {
				if (doTHParitiesMatch(-1)
				 && doEZParitiesMatch(-1)
				 && doBridgeBondsMatch(atomUsed, -1)) {
					// we currently have a match not considering exclude groups

					boolean isExcludedMatch = false;
					for (int excludeGroup=0; excludeGroup<mExcludeGroupCount; excludeGroup++) {
						if (isExcludeGroupMatch(atomUsed, index, excludeGroup)) {
							isExcludedMatch = true;
							break;
							}
						}

					if (countMode == cCountModeExistence && !isExcludedMatch)
						return 1;

					if (!isExcludedMatch) {
						addMatchIfQualifies(countMode);

						if (countMode == cCountModeFirstMatch)
							return 1;
						}
					}

				if (current == 0)	// if all fragment atoms are part of an exclude group
					break;

				current--;
				if (!mFragmentGraphIsRingClosure[current])
					atomUsed[mMatchTable[mFragmentGraphAtom[current]]] = false;
				}
			}

		return mMatchList.size();
		}


	private void addMatchIfQualifies(int countMode) {
		if (countMode == cCountModeFirstMatch
		 || countMode == cCountModeRigorous) {
			// count every match (even permutations of same atoms)
			addMatchAtoms();
			}
		else if (countMode == cCountModeOverlapping) {
			int[] sortedMatch = getSortedMatch(copyOf(mMatchTable, mMatchTable.length));
			if (!mSortedMatchSet.contains(sortedMatch)) {
				mSortedMatchSet.add(sortedMatch);
				addMatchAtoms();
				}
			}
		else if (countMode == cCountModeSeparated) {
			int[] sortedMatch = getSortedMatch(copyOf(mMatchTable, mMatchTable.length));
			if (!mSortedMatchSet.contains(sortedMatch)) {
				boolean found = false;
				for (int[] existing:mSortedMatchSet) {
					int existingIndex = 0;
					for (int atom:sortedMatch) {
						while (existingIndex < existing.length && existing[existingIndex] < atom)
							existingIndex++;
						if (existingIndex < existing.length) {
							if (atom == existing[existingIndex]) {
								found = true;
								break;
								}
							}
						}
					if (found)
						break;
					}
				if (!found) {
					mSortedMatchSet.add(sortedMatch);
					addMatchAtoms();
					}
				}
			}
		else if (countMode == cCountModeUnique) {
			int[] sortedMatch = getSortedSymmetryMatch(copyOf(mMatchTable, mMatchTable.length));
			if (!mSortedMatchSet.contains(sortedMatch)) {
				mSortedMatchSet.add(sortedMatch);
				addMatchAtoms();
				}
			}
		}

	private void addMatchAtoms() {
		mMatchList.add(copyOf(mMatchTable, mMatchTable.length));
		if (mBridgeBondList != null)
			mBridgeBondAtomList.add(copyOf(mIsBridgeBondAtom, mIsBridgeBondAtom.length));
		}

	/**
	 * @return sorted match atoms without excluded atoms
	 */
	private int[] getSortedMatch(int[] match) {
		int count = 0;
		for (int atom:match)
			if (atom == -1)
				count++;

		if (count != 0) {
			int[] oldMatch = match;
			match = new int[oldMatch.length - count];
			int index = 0;
			for (int atom:oldMatch)
				if (atom != -1)
					match[index++] = atom;
			}

		Arrays.sort(match);
		return match;
		}

	/**
	 * @return sorted match atoms without excluded atoms
	 */
	private int[] getSortedSymmetryMatch(int[] match) {
		int count = 0;
		for (int atom:match)
			if (atom == -1)
				count++;

		int[] symmetryMatch = new int[match.length - count];
		int index = 0;
		for (int i=0; i<match.length; i++) {
			if (match[i] != -1) {
				symmetryMatch[index] = (mFragment.getSymmetryRank(i) << 16) | mMolecule.getSymmetryRank(match[i]);
				if (mFragmentAtomContextRank != null)
					symmetryMatch[index] |= mFragmentAtomContextRank[i] << 24;
				index++;
				}
			}

		Arrays.sort(symmetryMatch);
		return symmetryMatch;
		}

	public boolean areAtomsSimilar(int moleculeAtom, int fragmentAtom) {
		int moleculeConnAtoms = mMolecule.getConnAtoms(moleculeAtom);
		int fragmentConnAtoms = mFragmentConnAtoms[fragmentAtom];

		if (fragmentConnAtoms > moleculeConnAtoms)
			return false;

		long moleculeQF = mMolecule.getAtomQueryFeatures(moleculeAtom);
		long fragmentQF = mFragment.getAtomQueryFeatures(fragmentAtom);

		int[] fragmentList = mFragment.getAtomList(fragmentAtom);
		int[] moleculeList = mMolecule.getAtomList(moleculeAtom);

		// check atomicNo's considering all 16 combinations of: Any set or not, atom list given or not
		if ((fragmentQF & Molecule.cAtomQFAny) != 0) {
			if (fragmentList != null) {
				if ((moleculeQF & Molecule.cAtomQFAny) != 0) {
					if (moleculeList == null)
						return false;

					if (!isSubListOf(fragmentList, moleculeList))
						return false;
					}
				else {
					if (moleculeList != null) {
						if (listsOverlap(moleculeList, fragmentList))
							return false;
						}
					else {
						if (isListMember(mMolecule.getAtomicNo(moleculeAtom), fragmentList))
							return false;
						}
					}
				}
			}	// 4 cases with Any set in fragment but no fragment exclude list given don't need to be checked
		else {
				// cAtomQFAny not set in fragment but Any set in molecule
			if ((moleculeQF & Molecule.cAtomQFAny) != 0)
				return false;	// regardless of possibly given lists these 4 cases cannot match

				// remaining cases: Any neither set in fragment nor molecule, but check for lists
			if (fragmentList != null) {
				if (moleculeList != null) {
					if (!isSubListOf(moleculeList, fragmentList))
						return false;
					}
				else {
					if (!isListMember(mMolecule.getAtomicNo(moleculeAtom), fragmentList))
						return false;
					}
				}
			else {
				if (moleculeList != null)
					return false;

				if (mMoleculeAtomType[moleculeAtom] != mFragmentAtomType[fragmentAtom])
					return false;
				}
			}	// end atomicNo, cAtomQFAny and AtomList checking

		if ((moleculeQF | fragmentQF) != 0) {
			if ((fragmentQF & Molecule.cAtomQFNoMoreNeighbours) != 0) {
				if (mMolecule.isFragment()
				 && (moleculeQF & Molecule.cAtomQFNoMoreNeighbours) == 0)
					return false;
				else if (fragmentConnAtoms != moleculeConnAtoms)
					return false;
				}
			if ((fragmentQF & Molecule.cAtomQFMoreNeighbours) != 0) {
				if ((fragmentConnAtoms >= moleculeConnAtoms)
				 && (moleculeQF & Molecule.cAtomQFMoreNeighbours) == 0)
					return false;
				}
			}

		// molecule features must be a subset of the fragment features
		if ((mMoleculeAtomFeatures[moleculeAtom] & ~mFragmentAtomFeatures[fragmentAtom]) != 0)
			return false;

		// all ring sizes found for fragment atom must all exist for molecule atom
		if ((mFragmentRingFeatures[fragmentAtom] & ~mMoleculeRingFeatures[moleculeAtom]) != 0)
			return false;

		long fragmentRingQF = fragmentQF & Molecule.cAtomQFNewRingSize;
		if (mMolecule.isFragment()) {
			// For a fragment in fragment search, the query fragment must not be more restrictive than the target.
			// Thus, if we have molecule ring features and no restriction on the query or more allowed features
			// on the query then don't consider the atom a match.
			long moleculeRingQF = fragmentQF & Molecule.cAtomQFNewRingSize;
			if (moleculeRingQF != 0 && (fragmentRingQF == 0 || (fragmentRingQF & ~moleculeRingQF) != 0))
				return false;
			}
		else {
			// at least one of the ring sizes defined in ring query features must match one of the ring sizes found in molecule atom
			if (fragmentRingQF != 0 && (fragmentRingQF & mMoleculeRingFeatures[moleculeAtom]) == 0)
				return false;
			}

		if (mFragment.getAtomCharge(fragmentAtom) != 0
		 && mFragment.getAtomCharge(fragmentAtom) != mMolecule.getAtomCharge(moleculeAtom))
			return false;
		if (mFragment.getAtomMass(fragmentAtom) != 0
		 && mFragment.getAtomMass(fragmentAtom) != mMolecule.getAtomMass(moleculeAtom))
			return false;
		if (mFragment.getAtomRadical(fragmentAtom) != 0
		 && mFragment.getAtomRadical(fragmentAtom) != mMolecule.getAtomRadical(moleculeAtom))
			return false;

		int smallestRingSize = (int)((mFragment.getAtomQueryFeatures(fragmentAtom) & Molecule.cAtomQFSmallRingSize) >> Molecule.cAtomQFSmallRingSizeShift);
		if (smallestRingSize != 0) {
			if (!mMolecule.isFragment()) {
				if (mMolecule.getAtomRingSize(moleculeAtom) != smallestRingSize)
					return false;
				}
			else {
				int targetRingSize = (int)((mMolecule.getAtomQueryFeatures(moleculeAtom) & Molecule.cAtomQFSmallRingSize) >> Molecule.cAtomQFSmallRingSizeShift);
				if (smallestRingSize != targetRingSize)
					return false;
				}
			}

		return true;
		}


	/**
	 * Check whether all stereo centers required to match are compatible
	 * with a match. First all individual stereo centers are checked whether
	 * a molecule atom's ESR setting and parity is in no respect less specific
	 * than the query's respective stereo center. This is based on the MDL
	 * rules for SSS matching of ESR types. In a second step it is checked
	 * whether all stereo match requiring atoms of any ESR group have the same
	 * relative configurations.
	 * @return true if all tetrahedral parities match
	 */
	private boolean doTHParitiesMatch(int excludeGroupNo) {
		int esrGroupAtomCount = 0;
		for (int fragmentAtom=0; fragmentAtom<mFragment.getAtoms(); fragmentAtom++) {
			if ((mExcludeGroupNo == null || mExcludeGroupNo[fragmentAtom] == excludeGroupNo)
			 && (mFragment.getAtomQueryFeatures(fragmentAtom) & Molecule.cAtomQFMatchStereo) != 0) {
				int moleculeAtom = mMatchTable[fragmentAtom];
				int fragmentParity = mFragment.getAtomParity(fragmentAtom);
				int moleculeParity = mMolecule.getAtomParity(moleculeAtom);

				// always consider as match if fragment atom is no stereo center
				if (fragmentParity == Molecule.cAtomParityNone)
			   		continue;

				// unknown fragment centers match everything
				if (fragmentParity == Molecule.cAtomParityUnknown)
			   		continue;

				// Here the fragment center is clearly specified as either 1 or 2.
				// Thus, unknown molecule centers should not be considered a match.
				// A molecule may not have a stereo center here for symmetry reasons,
				// but match 100% anyway. In this case we interpret: the user wants a
				// stereo center. Therefore, we don't consider no stereo centers a match.
				// Into the bargain: parities within idcodes don't include 'unknown',
				// because the information is implicit: if a stereo center within an idcode
				// has no 1 or 2 parity, then it is automatically treated to be unknown.
				// If idcodes are parsed with given or created coordinates, then implicit
				// unknowns are converted to explicit ones. However, if idcodes are intentionally
				// parsed without giving coordinates, then an unknown stereo center looks
				// like a no-stereo-center, because parities were taken from the idcode and
				// never calculated by ensureHelperArrays().
				if (moleculeParity == Molecule.cAtomParityNone
				 || moleculeParity == Molecule.cAtomParityUnknown)
			   		return false;

				// From here both, fragment and molecule, have a defined parity: 1 or 2.

				if (mFragment.getAtomESRType(fragmentAtom) == Molecule.cESRTypeAnd) {
					esrGroupAtomCount++;
					continue;
					}

				if (mMolecule.getAtomESRType(moleculeAtom) == Molecule.cESRTypeAnd)
			   		return false;

				if (mFragment.getAtomESRType(fragmentAtom) == Molecule.cESRTypeOr) {
					esrGroupAtomCount++;
					continue;
					}

				if (mMolecule.getAtomESRType(moleculeAtom) == Molecule.cESRTypeOr)
					return false;

				if (isTHParityInversion(fragmentAtom) == (fragmentParity == moleculeParity))
					return false;   // inverted parity found
				}
			}

		// now checking relative configurations of ESR groups...
		if (esrGroupAtomCount != 0) {
			int[] esrAtom = new int[esrGroupAtomCount];
			int esrAtomIndex = 0;
			for (int fragmentAtom=0; fragmentAtom<mFragment.getAtoms(); fragmentAtom++) {
				if ((mExcludeGroupNo == null || mExcludeGroupNo[fragmentAtom] == excludeGroupNo)
				 && (mFragment.getAtomQueryFeatures(fragmentAtom) & Molecule.cAtomQFMatchStereo) != 0) {
					int fragmentParity = mFragment.getAtomParity(fragmentAtom);
					if (fragmentParity != Molecule.cAtomParityNone
					 && fragmentParity != Molecule.cAtomParityUnknown) {
						esrAtom[esrAtomIndex++] = (mFragment.getAtomESRGroup(fragmentAtom) << 24)
												| (mFragment.getAtomESRType(fragmentAtom) << 22)
												| fragmentAtom;
						}
					}
				}
			Arrays.sort(esrAtom);
			esrAtomIndex = 0;
			while (esrAtomIndex < esrAtom.length) {
				int fragmentBaseAtom = esrAtom[esrAtomIndex] & 0x003FFFFF;
				int moleculeBaseAtom = mMatchTable[fragmentBaseAtom];
				int baseGroupAndType = esrAtom[esrAtomIndex] & 0xFFC00000;
				boolean baseParityComparison = isTHParityInversion(fragmentBaseAtom)
									^ (mFragment.getAtomParity(fragmentBaseAtom)
									== mMolecule.getAtomParity(moleculeBaseAtom));
				for (esrAtomIndex++; esrAtomIndex<esrAtom.length
						&& (esrAtom[esrAtomIndex] & 0xFFC00000)==baseGroupAndType; esrAtomIndex++) {
					int fragmentAtom = esrAtom[esrAtomIndex] & 0x003FFFFF;
					int moleculeAtom = mMatchTable[fragmentAtom];
					if (mMolecule.getAtomESRType(moleculeAtom) != mMolecule.getAtomESRType(moleculeBaseAtom)
					 || mMolecule.getAtomESRGroup(moleculeAtom) != mMolecule.getAtomESRGroup(moleculeBaseAtom))
						return false;
					boolean parityComparison = isTHParityInversion(fragmentAtom)
											 ^ (mFragment.getAtomParity(fragmentAtom)
											 == mMolecule.getAtomParity(moleculeAtom));
					if (parityComparison != baseParityComparison)
						return false;
					}
				}
			}

		return true;
		}


	private boolean isTHParityInversion(int fragmentAtom) {
		boolean inversion = false;
		if (mFragment.getAtomPi(fragmentAtom) == 0) {
			for (int i=1; i<mFragment.getConnAtoms(fragmentAtom); i++) {
				for (int j=0; j<i; j++) {
					int connAtom1 = mFragment.getConnAtom(fragmentAtom,i);
					int connAtom2 = mFragment.getConnAtom(fragmentAtom,j);
					if ((mMatchTable[connAtom1] > mMatchTable[connAtom2])
					  ^ (connAtom1 > connAtom2))
						inversion = !inversion;
					}
				}
			}
		else {  // allene parities
			for (int i=0; i<mFragment.getConnAtoms(fragmentAtom); i++) {
				int connAtom = mFragment.getConnAtom(fragmentAtom,i);
				int neighbours = 0;
				int[] neighbour = new int[3];
				for (int j=0; j<mFragment.getConnAtoms(connAtom); j++) {
					neighbour[neighbours] = mFragment.getConnAtom(connAtom,j);
					if (neighbour[neighbours] != fragmentAtom)
						neighbours++;
					}
				if (neighbours == 2
				 && ((mMatchTable[neighbour[0]] > mMatchTable[neighbour[1]])
					^(neighbour[0] > neighbour[1])))
					inversion = !inversion;
				}
			}
		return inversion;
		}


	/**
	 * Check whether all double bond parities required to match are compatible
	 * with a match.
	 * @return true if all E/Z parities match
	 */
	private boolean doEZParitiesMatch(int excludeGroupNo) {
		for (int fragmentBond=0; fragmentBond<mFragment.getBonds(); fragmentBond++) {
			if ((mFragment.getBondQueryFeatures(fragmentBond) & Molecule.cBondQFMatchStereo) != 0) {
				int fragmentParity = mFragment.getBondParity(fragmentBond);

				// always consider as match if fragment bond parity is none
				if (fragmentParity == Molecule.cBondParityNone)
			   		continue;

				int fragmentAtom1 = mFragment.getBondAtom(0, fragmentBond);
				int fragmentAtom2 = mFragment.getBondAtom(1, fragmentBond);

				if (mExcludeGroupNo == null
				 ||	(excludeGroupNo == -1 && mExcludeGroupNo[fragmentAtom1] == -1 && mExcludeGroupNo[fragmentAtom2] == -1)
				 || (excludeGroupNo != -1 && (mExcludeGroupNo[fragmentAtom1] == excludeGroupNo || mExcludeGroupNo[fragmentAtom2] == excludeGroupNo))) {
//				if ((mIsExcludeAtom[fragmentAtom1] || mIsExcludeAtom[fragmentAtom2]) == isExcludeGroup) {
					int moleculeAtom1 = mMatchTable[fragmentAtom1];
					int moleculeAtom2 = mMatchTable[fragmentAtom2];
					int moleculeBond = mMolecule.getBond(moleculeAtom1, moleculeAtom2);

					// consider as match if an E/Z-bond with defined parity atom matches on one with no parity
					int moleculeParity = mMolecule.getBondParity(moleculeBond);
					if (moleculeParity == Molecule.cBondParityNone) {
						if (mMolecule.isSmallRingBond(moleculeBond))
							moleculeParity = calculateImplicitSmallRingBondParity(moleculeBond);
						if (moleculeParity == Molecule.cBondParityNone)
							continue;
						}

					// unknown molecule E/Z-bonds need to match everything because
					// parities from idcodes don't include them, i.e. depending
					// on the source of the parities, unknown E/Z-bonds may look as
					// non-E/Z-bonds. Both must retrieve the same results.
					if (fragmentParity == Molecule.cBondParityUnknown)
						continue;

					if (moleculeParity == Molecule.cBondParityUnknown)
						continue;

					if (isEZParityInversion(fragmentBond, moleculeBond) == (fragmentParity == moleculeParity))
						return false;   // inverted parity found
					}
				}
			}

		return true;
		}

	/**
	 * IDcodes don't store double bond conformations if they are implicitly given due to small ring membership.
	 * Therefore, if we need to match EZ-configurations and the molecule bond is a small ring bond, we may have
	 * no given configuration, if the molecules was directly created from an idcode.
	 * @param moleculeBond
	 * @return
	 */
	private int calculateImplicitSmallRingBondParity(int moleculeBond) {
		RingCollection ringSet = mMolecule.getRingSet();
		for (int r=0; r<ringSet.getSize(); r++) {
			if (ringSet.isBondMember(r, moleculeBond)) {
				int[] relevantAtom = new int[2];
				for (int i=0; i<2; i++) {
					relevantAtom[i] = Integer.MAX_VALUE;
					int bondAtom = mMolecule.getBondAtom(i, moleculeBond);
					for (int j=0; j<mMolecule.getConnAtoms(bondAtom); j++) {
						int atom = mMolecule.getConnAtom(bondAtom, j);
						if (atom != mMolecule.getBondAtom(1-i, moleculeBond)
								&& relevantAtom[i] > atom)
							relevantAtom[i] = atom;
						}
					}

				int memberCount = 0;
				if (ringSet.isAtomMember(r, relevantAtom[0]))
					memberCount++;
				if (ringSet.isAtomMember(r, relevantAtom[1]))
					memberCount++;
				if (memberCount == 2)
					return Molecule.cBondCIPParityZorM;
				if (memberCount == 1)
					return Molecule.cBondCIPParityEorP;

				return Molecule.cBondCIPParityZorM;
				}
			}

		return Molecule.cBondCIPParityNone;
		}


//	private boolean isEZParityInversion(int fragmentBond, int moleculeBond) {
//		boolean inversion = false;
//		for (int i=0; i<2; i++) {
//			int fragmentAtom = mFragment.getBondAtom(i, fragmentBond);
//			int moleculeAtom = mMatchTable[fragmentAtom];
//			if (mFragment.getConnAtoms(fragmentAtom) == 2) {
//				if (mMolecule.getConnAtoms(moleculeAtom) == 2)
//					continue;
//
//				int fragmentNeighbour = -1;
//				for (int j=0; j<2; j++)
//					if (mFragment.getConnBond(fragmentAtom, j) != fragmentBond)
//						fragmentNeighbour = mFragment.getConnAtom(fragmentAtom, j);
//
//				int moleculeNeighbours = 0;
//				int[] moleculeNeighbour = new int[2];
//				for (int j=0; j<3; j++)
//					if (mMolecule.getConnBond(moleculeAtom, j) != moleculeBond)
//						moleculeNeighbour[moleculeNeighbours++] = mMolecule.getConnAtom(moleculeAtom, j);
//
//				if (mMatchTable[fragmentNeighbour] != moleculeNeighbour[0])
//					inversion = !inversion;
//				}
//			else if (mFragment.getConnAtoms(fragmentAtom) == 3
//	   			  && mMolecule.getConnAtoms(moleculeAtom) == 3) {
//				int[] fragmentNeighbour = new int[2];
//				int fragmentNeighbours = 0;
//				for (int j=0; j<3; j++)
//					if (mFragment.getConnBond(fragmentAtom, j) != fragmentBond)
//						fragmentNeighbour[fragmentNeighbours++] = mFragment.getConnAtom(fragmentAtom, j);
//				if ((mMatchTable[fragmentNeighbour[0]] > mMatchTable[fragmentNeighbour[1]])
//				  ^ (fragmentNeighbour[0] > fragmentNeighbour[1]))
//					inversion = !inversion;
//				}
//			}
//		return inversion;
//		}


	private boolean isEZParityInversion(int fragmentBond, int moleculeBond) {
		boolean inversion = false;
		for (int i=0; i<2; i++) {
			int fragmentAtom = mFragment.getBondAtom(i, fragmentBond);
			int moleculeAtom = mMatchTable[fragmentAtom];
			if (mMolecule.getConnAtoms(moleculeAtom) > 2) {
				int otherFragmentAtom = mFragment.getBondAtom(1-i, fragmentBond);
				int lowFragmentNeighbour = Integer.MAX_VALUE;
				for (int j=0; j<mFragment.getConnAtoms(fragmentAtom); j++) {
					int fragmentNeighbour = mFragment.getConnAtom(fragmentAtom, j);
					if (fragmentNeighbour != otherFragmentAtom
					 && lowFragmentNeighbour > fragmentNeighbour)
						lowFragmentNeighbour = fragmentNeighbour;
					}

				int otherMoleculeAtom = mMatchTable[otherFragmentAtom];
				int lowMoleculeNeighbour = Integer.MAX_VALUE;
				for (int j=0; j<mMolecule.getConnAtoms(moleculeAtom); j++) {
					int moleculeNeighbour = mMolecule.getConnAtom(moleculeAtom, j);
					if (moleculeNeighbour != otherMoleculeAtom
					 && lowMoleculeNeighbour > moleculeNeighbour)
						lowMoleculeNeighbour = moleculeNeighbour;
					}

				if (mMatchTable[lowFragmentNeighbour] != lowMoleculeNeighbour)
					inversion = !inversion;
				}
			}
		return inversion;
		}

	/**
	 * Starting from a full match of the fragment without exclude groups, this method continues
	 * the graph matching to find, whether the specified exclude group can also be matched.
	 * @param atomUsed
	 * @return
	 */
	private boolean isExcludeGroupMatch(boolean[] atomUsed, int[] index, int excludeGroupNo) {
		int excludeGroupGraphBase = mExcludeGroupGraphIndex[excludeGroupNo];
		int excludeGroupGraphMax = excludeGroupGraphBase + 1;
		while (excludeGroupGraphMax < mFragmentGraphSizeWithExcludeGroups
		 && mExcludeGroupNo[mFragmentGraphAtom[excludeGroupGraphMax]] == excludeGroupNo)
			excludeGroupGraphMax++;

		for (int i=excludeGroupGraphBase; i<excludeGroupGraphMax; i++)
			index[i] = -1;

		int current = excludeGroupGraphBase;

		while (true) {
/*
System.out.print("  index:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(index[i]==-1?"-":""+index[i])); System.out.println();
System.out.print("		"); for (int i=0; i<current; i++) System.out.print("  "); System.out.println(" ^");
System.out.println();
*/
			int maxIndex = (mFragmentGraphParentAtom[current] == -1) ? mMolecule.getAtoms()
					: mMolecule.getAllConnAtomsPlusMetalBonds(mMatchTable[mFragmentGraphParentAtom[current]]);

			index[current]++;

			if (index[current] == maxIndex) {
				index[current] = -1;
				if (current == excludeGroupGraphBase)
					break;
				current--;
				if (!mFragmentGraphIsRingClosure[current]) {
					atomUsed[mMatchTable[mFragmentGraphAtom[current]]] = false;
					mMatchTable[mFragmentGraphAtom[current]] = -1;
					}
				continue;
				}

			if (mFragmentGraphParentAtom[current] == -1) {	// if current graph atom is sub fragment anchor atom
				if (!atomUsed[index[current]]) {
					if (areAtomsSimilar(index[current], mFragmentGraphAtom[current])) {
						mMatchTable[mFragmentGraphAtom[current]] = index[current];
						atomUsed[index[current]] = true;
						current++;
					}
				}
			}
			else {
				// skip plain hydrogens
				if (mMolecule.getConnAtom(mMatchTable[mFragmentGraphParentAtom[current]], index[current]) >= mMolecule.getAtoms()) {
					index[current]++;
					continue;
					}

				int candidate = mMolecule.getConnAtom(mMatchTable[mFragmentGraphParentAtom[current]], index[current]);
				if (!mFragmentGraphIsRingClosure[current]) {	// current graph position is not an anchor
					if (!atomUsed[candidate]) {
						if (areAtomsSimilar(candidate, mFragmentGraphAtom[current])
						 && areBondsSimilar(mMolecule.getConnBond(mMatchTable[mFragmentGraphParentAtom[current]], index[current]), mFragmentGraphParentBond[current])) {
							atomUsed[candidate] = true;
							mMatchTable[mFragmentGraphAtom[current]] = candidate;
							current++;
							}
						}
					}
				else {	// current graph position is ringClosure
					if (candidate == mMatchTable[mFragmentGraphAtom[current]]
					 && areBondsSimilar(mMolecule.getConnBond(mMatchTable[mFragmentGraphParentAtom[current]], index[current]), mFragmentGraphParentBond[current])) {
						current++;
						}
					}
				}

			if (current == excludeGroupGraphMax) {
				if (doTHParitiesMatch(excludeGroupNo)
				 && doEZParitiesMatch(excludeGroupNo)
				 && doBridgeBondsMatch(atomUsed, excludeGroupNo)) {

					// remove match table entries for exclude atoms
					for (int i=excludeGroupGraphBase; i<excludeGroupGraphMax; i++) {
						if (!mFragmentGraphIsRingClosure[i]) {
							int atom = mFragmentGraphAtom[i];
							atomUsed[mMatchTable[atom]] = false;
							mMatchTable[atom] = -1;
							}
						}

					return true;
					}

				current--;
				if (!mFragmentGraphIsRingClosure[current]) {
					atomUsed[mMatchTable[mFragmentGraphAtom[current]]] = false;
					mMatchTable[mFragmentGraphAtom[current]] = -1;
					}
				}
			}

		return false;
		}

	/**
	 * Currently bridge bonds are considered as follows:<br>
	 * - Bridge bonds in the molecule are not matched to any fragment bond<br>
	 * - Bridge bonds in the fragment are not considered in the graph-matching and<br>
	 *   checked after an otherwise successful match, whether the shortest path of<br>
	 *   unmatched atoms satisfies the min- and max-atom criteria.<br>
	 * For fragment and(!) molecule containing bridge bonds we would need a more complex
	 * handling: - for every fragment bridge check all unmatched molecule paths and count
	 * for every path min- and max-lengths (considering all bridge bonds within path).
	 * Consider a match, if min- and max-atoms range fits into fragments bridge bond range.
	 * If we have multiple bridge bonds in fragment and partially overlapping bridge matches
	 * in the molecule, it would get very nasty...
	 * @param moleculeAtomUsed
	 * @param excludeGroupNo
	 * @return
	 */
	private boolean doBridgeBondsMatch(boolean[] moleculeAtomUsed, int excludeGroupNo) {
		if (mBridgeBondList != null) {
			mIsBridgeBondAtom = new boolean[moleculeAtomUsed.length];
			for (BridgeBond bb:mBridgeBondList) {
				if (mExcludeGroupNo == null
				 ||	(excludeGroupNo == -1 && mExcludeGroupNo[bb.atom1] == -1 && mExcludeGroupNo[bb.atom2] == -1)
				 || (excludeGroupNo != -1 && (mExcludeGroupNo[bb.atom1] == excludeGroupNo || mExcludeGroupNo[bb.atom2] == excludeGroupNo))) {
					int[] pathAtom = new int[bb.maxBridgeSize+2];
					int bridgeSize = mMolecule.getPath(pathAtom, mMatchTable[bb.atom1], mMatchTable[bb.atom2], bb.maxBridgeSize+1, moleculeAtomUsed, null) - 1;
					if (bridgeSize < bb.minBridgeSize
					 || bridgeSize > bb.maxBridgeSize)
						return false;

					for (int i=1; i<=bridgeSize; i++)
						mIsBridgeBondAtom[pathAtom[i]] = true;
					}
				}
			}

		return true;
		}


	/**
	 * @param moleculeBond flag list of fragment bond features (features present)
	 * @param fragmentBond flag list of molecule bond features (features allowed)
	 * @return true if all molecule bond features are present in fragment bond
	 */
	public boolean areBondsSimilar(int moleculeBond, int fragmentBond) {
		int molDefaults = mMoleculeBondFeatures[moleculeBond];
		int frgDefaults = mFragmentBondFeatures[fragmentBond];

		if ((mFragment.getBondQueryFeatures(fragmentBond) & Molecule.cBondQFMatchFormalOrder) != 0) {
			int molBondType = mMolecule.getBondTypeSimple(moleculeBond);
			int frgBondType = mFragment.getBondTypeSimple(fragmentBond);
			int frgBondTypes = mFragment.getBondQueryFeatures(fragmentBond) & Molecule.cBondQFBondTypes;
			if (molBondType != frgBondType
			 && !(molBondType == Molecule.cBondTypeSingle && (frgBondTypes & Molecule.cBondTypeSingle) != 0)
			 && !(molBondType == Molecule.cBondTypeDouble && (frgBondTypes & Molecule.cBondTypeDouble) != 0)
			 && !(molBondType == Molecule.cBondTypeTriple && (frgBondTypes & Molecule.cBondTypeTriple) != 0)
			 && !(molBondType == Molecule.cBondTypeQuadruple && (frgBondTypes & Molecule.cBondTypeQuadruple) != 0)
			 && !(molBondType == Molecule.cBondTypeQuintuple && (frgBondTypes & Molecule.cBondTypeQuintuple) != 0)
			 && !(molBondType == Molecule.cBondTypeMetalLigand && (frgBondTypes & Molecule.cBondTypeMetalLigand) != 0)
			 && !(molBondType == Molecule.cBondTypeDelocalized && (frgBondTypes & Molecule.cBondTypeDelocalized) != 0))
				return false;

			molDefaults &= ~Molecule.cBondQFBondTypes;
			frgDefaults &= ~Molecule.cBondQFBondTypes;
			}

		if ((molDefaults & ~frgDefaults) != 0)
			return false;

		int ringSize = (mFragment.getBondQueryFeatures(fragmentBond) & Molecule.cBondQFRingSize) >> Molecule.cBondQFRingSizeShift;
		if (ringSize != 0) {
			if (mMolecule.isFragment()
			 && ringSize == (mMolecule.getBondQueryFeatures(fragmentBond) & Molecule.cBondQFRingSize) >> Molecule.cBondQFRingSizeShift)
				return true;

			if (ringSize <= 2) {    // ring size 8-11 is encoded as 1; ring size >=12 is encoded as 2
				int moleculeRingSize = mMolecule.getBondRingSize(moleculeBond);
				if (ringSize == 1)
					return (moleculeRingSize >= 8) && (moleculeRingSize <= 12);
				else
					return moleculeRingSize >= 12;
				}

			boolean found = false;
			RingCollection ringSet = mMolecule.getRingSet();
			for (int i=0; i<ringSet.getSize(); i++) {
				if (ringSet.getRingSize(i) == ringSize) {
					if (ringSet.isBondMember(i, moleculeBond)) {
						found = true;
						break;
						}
					}
				}
			if (!found)
				return false;
			}

		return true;
		}


	private boolean isSubListOf(int[] list1, int[] list2) {
			// returns true if list2 contains at least all members of list1
		int i2 = 0;
		for (int i1=0; i1<list1.length; i1++) {
			int atomicNo1 = list1[i1];
			while (list2[i2] < atomicNo1) {
				i2++;
				if (i2 == list2.length)
					return false;
				}
			if (list2[i2] > atomicNo1)
				return false;
			}
		return true;
		}


	private boolean listsOverlap(int[] list1, int[] list2) {
			// returns true if both lists share at least one common member
		int i1 = 0;
		int i2 = 0;

		while (i1<list1.length && i2<list2.length) {
			int atomicNo1 = list1[i1];
			int atomicNo2 = list2[i2];

			if (atomicNo1 == atomicNo2)
				return true;

			if (atomicNo1 < atomicNo2)
				i1++;
			else
				i2++;
			}

		return false;
		}


	private boolean isListMember(int atomicNo, int[] list) {
		// returns true if list contains atomicNo
		for (int i=0; i<list.length; i++)
			if (list[i] == atomicNo)
				return true;

		return false;
		}


	public void setupAtomAndBondFeatures(int matchMode) {
		if (!mMoleculeFeaturesValid) {
			setupMoleculeFeatures(matchMode);
			mMoleculeFeaturesValid = true;
			}

		if (!mFragmentFeaturesValid) {
			setupFragmentFeatures(matchMode);

			buildFragmentGraph();
			buildBridgeBondList();

			mFragmentFeaturesValid = true;
			}
		}

	private void setupMoleculeFeatures(int matchMode) {
		mMolecule.ensureHelperArrays(mRequiredHelperLevel);
		int nTotalMoleculeAtoms = mMolecule.getAtoms();

		mMoleculeAtomType = new int[nTotalMoleculeAtoms];
		mMoleculeAtomFeatures = new long[nTotalMoleculeAtoms];

		for (int atom=0; atom<nTotalMoleculeAtoms; atom++) {
			mMoleculeAtomFeatures[atom] = ((getAtomQueryDefaults(mMolecule, atom)
				| mMolecule.getAtomQueryFeatures(atom))
					& Molecule.cAtomQFSimpleFeatures)
					^ Molecule.cAtomQFNarrowing;

			mMoleculeAtomType[atom] = mMolecule.getAtomicNo(atom);

			if ((matchMode & cMatchAtomCharge) != 0)
				mMoleculeAtomType[atom] += (mMolecule.getAtomCharge(atom) + 16) << 8;

			if ((matchMode & cMatchAtomMass) != 0)
				mMoleculeAtomType[atom] += mMolecule.getAtomMass(atom) << 16;
			}

		mMoleculeRingFeatures = new long[nTotalMoleculeAtoms];
		RingCollection ringSet = mMolecule.getRingSet();
		for (int i=0; i<ringSet.getSize(); i++) {
			int ringSize = ringSet.getRingSize(i);
			for (int atom:ringSet.getRingAtoms(i)) {
				if (ringSize == 3)
					mMoleculeRingFeatures[atom] |= Molecule.cAtomQFRingSize3;
				else if (ringSize == 4)
					mMoleculeRingFeatures[atom] |= Molecule.cAtomQFRingSize4;
				else if (ringSize == 5)
					mMoleculeRingFeatures[atom] |= Molecule.cAtomQFRingSize5;
				else if (ringSize == 6)
					mMoleculeRingFeatures[atom] |= Molecule.cAtomQFRingSize6;
				else if (ringSize == 7)
					mMoleculeRingFeatures[atom] |= Molecule.cAtomQFRingSize7;
				}
			}
		for (int atom=0; atom<nTotalMoleculeAtoms; atom++) {
			int ringSize = mMolecule.getAtomRingSize(atom);
			if (ringSize == 0)
				mMoleculeRingFeatures[atom] |= Molecule.cAtomQFRingSize0;
			else if (ringSize > 7)
				mMoleculeRingFeatures[atom] |= Molecule.cAtomQFRingSizeLarge;
			}

		int nTotalMoleculeBonds = mMolecule.getBonds();

		mMoleculeBondFeatures = new int[nTotalMoleculeBonds];

		for (int bond=0; bond<nTotalMoleculeBonds; bond++)
			mMoleculeBondFeatures[bond] = (getBondQueryDefaults(mMolecule, bond)
					| mMolecule.getBondQueryFeatures(bond))
					& (Molecule.cBondQFSimpleFeatures | Molecule.cBondQFBridge)
					^ Molecule.cBondQFNarrowing;
			// include cBondQFBridge features to make sure that bridge bonds in molecule are never matched directly
			}

	private void setupFragmentFeatures(int matchMode) {
		long[] atomFeaturesWithoutExcludeAtoms = null;
		int[] bondFeaturesWithoutExcludeAtoms = null;
		int[] atomTypeWithoutExcludeAtoms = null;

		mFragment.ensureHelperArrays(mRequiredHelperLevel);
		mFragmentConnAtoms = new int[mFragment.getAtoms()];
		for (int atom=0; atom<mFragment.getAtoms(); atom++)
			mFragmentConnAtoms[atom] = mFragment.getConnAtoms(atom);

		// If we have exclude groups, we need to determine atom and bond features without the exclude atoms and
		// map them to the original fragment atom and bond indexes.
		if (mFragmentExcludeAtoms != 0) {
			StereoMolecule fragmentWithoutExcludeGroups = new StereoMolecule(mFragment.getAllAtoms(), mFragment.getAllBonds());
			boolean[] isNonExcludeAtom = new boolean[mFragment.getAllAtoms()];
			for (int atom=0; atom<mFragment.getAllAtoms(); atom++)
				isNonExcludeAtom[atom] = !mIsExcludeAtom[atom];
			mFragment.copyMoleculeByAtoms(fragmentWithoutExcludeGroups, isNonExcludeAtom, true, null);

			fragmentWithoutExcludeGroups.ensureHelperArrays(mRequiredHelperLevel);
			setupFragmentFeatures(fragmentWithoutExcludeGroups, matchMode);
			atomFeaturesWithoutExcludeAtoms = mFragmentAtomFeatures;
			bondFeaturesWithoutExcludeAtoms = mFragmentBondFeatures;
			atomTypeWithoutExcludeAtoms = mFragmentAtomType;

			int index = 0;
			for (int atom=0; atom<mFragment.getAtoms(); atom++)
				if (!mIsExcludeAtom[atom])
					mFragmentConnAtoms[atom] = fragmentWithoutExcludeGroups.getConnAtoms(index++);
			}

		// We cannot skip this if we have exclude groups:
		// We need to determine exclude groups atom & bond features for the exclude group matching phase.
		setupFragmentFeatures(mFragment, matchMode);

		if (mFragmentExcludeAtoms != 0) {
			int index = 0;
			for (int atom=0; atom<mFragment.getAllAtoms(); atom++) {
				if (!mIsExcludeAtom[atom]) {
					mFragmentAtomFeatures[atom] = atomFeaturesWithoutExcludeAtoms[index];
					mFragmentAtomType[atom] = atomTypeWithoutExcludeAtoms[index++];
					}
				}
			index = 0;
			for (int bond=0; bond<mFragment.getAllBonds(); bond++) {
				if (!mIsExcludeAtom[mFragment.getBondAtom(0, bond)] && !mIsExcludeAtom[mFragment.getBondAtom(1, bond)]) {
					mFragmentBondFeatures[bond] = bondFeaturesWithoutExcludeAtoms[index++];
					}
				}
			}
		}

	private void setupFragmentFeatures(StereoMolecule fragment, int matchMode) {
		int nTotalFragmentAtoms = fragment.getAtoms();

		mFragmentAtomFeatures = new long[fragment.getAtoms()];
		mFragmentAtomType = new int[fragment.getAtoms()];

		for (int atom=0; atom<nTotalFragmentAtoms; atom++) {
			mFragmentAtomFeatures[atom] = ((getAtomQueryDefaults(fragment, atom)
				| fragment.getAtomQueryFeatures(atom))
					& Molecule.cAtomQFSimpleFeatures)
					^ Molecule.cAtomQFNarrowing;
			mFragmentAtomType[atom] = fragment.getAtomicNo(atom);

			if ((matchMode & cMatchAtomCharge) != 0)
				mFragmentAtomType[atom] += (fragment.getAtomCharge(atom) + 16) << 8;

			if ((matchMode & cMatchAtomMass) != 0)
				mFragmentAtomType[atom] += fragment.getAtomMass(atom) << 16;
			}

		mFragmentRingFeatures = new long[fragment.getAtoms()];
		RingCollection ringSet = fragment.getRingSet();
		for (int i=0; i<ringSet.getSize(); i++) {
			boolean containsBridgeBond = false;
			for (int bond:ringSet.getRingBonds(i)) {
				if (fragment.isBondBridge(bond)) {
					containsBridgeBond = true;
					break;
					}
				}
			if (!containsBridgeBond) {
				int ringSize = ringSet.getRingSize(i);
				for (int atom:ringSet.getRingAtoms(i)) {
					if (ringSize == 3)
						mFragmentRingFeatures[atom] |= Molecule.cAtomQFRingSize3;
					else if (ringSize == 4)
						mFragmentRingFeatures[atom] |= Molecule.cAtomQFRingSize4;
					else if (ringSize == 5)
						mFragmentRingFeatures[atom] |= Molecule.cAtomQFRingSize5;
					else if (ringSize == 6)
						mFragmentRingFeatures[atom] |= Molecule.cAtomQFRingSize6;
					else if (ringSize == 7)
						mFragmentRingFeatures[atom] |= Molecule.cAtomQFRingSize7;
					}
				}
			}
// Cannot require that, because if a molecule atom is also part of a small ring,
// then the large ring membership is not known anymore
//		for (int atom=0; atom<nTotalFragmentAtoms; atom++)
//			if (fragment.getAtomRingSize(atom) > 7)
//				mFragmentRingFeatures[atom] |= Molecule.cAtomQFRingSizeLarge;

		int nTotalFragmentBonds = fragment.getBonds();

		mFragmentBondFeatures = new int[fragment.getBonds()];

		for (int bond=0; bond<nTotalFragmentBonds; bond++) {
			mFragmentBondFeatures[bond] = (getBondQueryDefaults(fragment, bond)
					| fragment.getBondQueryFeatures(bond))
					& Molecule.cBondQFSimpleFeatures
					^ Molecule.cBondQFNarrowing;

			// match fragment's single/double bonds to delocalized molecule bonds also
			if ((matchMode & cMatchDBondToDelocalized) != 0) {
				if ((mFragmentBondFeatures[bond] & Molecule.cBondTypeDouble) != 0)
					mFragmentBondFeatures[bond] |= Molecule.cBondTypeDelocalized;
				}
			else if ((matchMode & cMatchAromDBondToDelocalized) != 0) {
				if ((mFragmentBondFeatures[bond] & Molecule.cBondTypeDouble) != 0
						&& fragment.isAromaticBond(bond))
					mFragmentBondFeatures[bond] |= Molecule.cBondTypeDelocalized;
				}
			}
		}

	/**
	 * Generates inherent feature flags of a given atom.
	 * @param mol molecule or fragment of the SSS
	 * @param atom the atom of which to generate feature flags
	 * @return atom features independent of query features
	 */
	private long getAtomQueryDefaults(StereoMolecule mol, int atom) {
		long queryDefaults = 0;

		if (!mol.isFragment()) {
			if (mol.isHeteroAromaticAtom(atom))
				queryDefaults |= (Molecule.cAtomQFAromatic
								| Molecule.cAtomQFHeteroAromatic);
			else if (mol.isAromaticAtom(atom))
				queryDefaults |= Molecule.cAtomQFAromatic;
			else
				queryDefaults |= Molecule.cAtomQFNotAromatic;

			if (mol.isAtomStereoCenter(atom))
				queryDefaults |= Molecule.cAtomQFIsStereo;
			else
				queryDefaults |= Molecule.cAtomQFIsNotStereo;

			int ringBonds = mol.getAtomRingBondCount(atom);
			if (ringBonds == 0)
				queryDefaults |= (Molecule.cAtomQFNot2RingBonds
								| Molecule.cAtomQFNot3RingBonds
								| Molecule.cAtomQFNot4RingBonds);
			else if (ringBonds == 2)
				queryDefaults |= (Molecule.cAtomQFNotChain
								| Molecule.cAtomQFNot3RingBonds
								| Molecule.cAtomQFNot4RingBonds);
			else if (ringBonds == 3)
				queryDefaults |= (Molecule.cAtomQFNotChain
								| Molecule.cAtomQFNot2RingBonds
								| Molecule.cAtomQFNot4RingBonds);
			else
				queryDefaults |= (Molecule.cAtomQFNotChain
								| Molecule.cAtomQFNot2RingBonds
								| Molecule.cAtomQFNot3RingBonds);

			int charge = mol.getAtomCharge(atom);
			if (charge == 0)
				queryDefaults |= (Molecule.cAtomQFNotChargeNeg | Molecule.cAtomQFNotChargePos);
			else if (charge < 0)
				queryDefaults |= (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargePos);
			else if (charge > 0)
				queryDefaults |= (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargeNeg);

			int hydrogens = mol.getAllHydrogens(atom);
			switch (hydrogens) {
			case 0:
				queryDefaults |= (Molecule.cAtomQFHydrogen & ~Molecule.cAtomQFNot0Hydrogen);
				break;
			case 1:
				queryDefaults |= (Molecule.cAtomQFHydrogen & ~Molecule.cAtomQFNot1Hydrogen);
				break;
			case 2:
				queryDefaults |= (Molecule.cAtomQFHydrogen & ~Molecule.cAtomQFNot2Hydrogen);
				break;
			default:
				queryDefaults |= (Molecule.cAtomQFHydrogen & ~Molecule.cAtomQFNot3Hydrogen);
				break;
				}

			int neighbours = mol.getConnAtoms(atom);
			switch (neighbours) {
			case 0:
				queryDefaults |= (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot0Neighbours);
				break;
			case 1:
				queryDefaults |= (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot1Neighbour);
				break;
			case 2:
				queryDefaults |= (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot2Neighbours);
				break;
			case 3:
				queryDefaults |= (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot3Neighbours);
				break;
			default:
				queryDefaults |= (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot4Neighbours);
				break;
				}

			int eValue = mol.getAtomElectronegativeNeighbours(atom);
			switch (eValue) {
				case 0:
					queryDefaults |= (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot0ENeighbours);
					break;
				case 1:
					queryDefaults |= (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot1ENeighbour);
					break;
				case 2:
					queryDefaults |= (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot2ENeighbours);
					break;
				case 3:
					queryDefaults |= (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot3ENeighbours);
					break;
				default:
					queryDefaults |= (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot4ENeighbours);
					break;
				}

			int piElectrons = mol.getAtomPi(atom);
			switch (piElectrons) {
			case 0:
				queryDefaults |= (Molecule.cAtomQFNot1PiElectron
								| Molecule.cAtomQFNot2PiElectrons);
				break;
			case 1:
				queryDefaults |= (Molecule.cAtomQFNot0PiElectrons
								| Molecule.cAtomQFNot2PiElectrons);
				break;
			default:
				queryDefaults |= (Molecule.cAtomQFNot0PiElectrons
								| Molecule.cAtomQFNot1PiElectron);
				break;
				}
			}
		else {	// The fragments implicit features are not really necessary,
				  // but may speed up the graph matching.
			if (mol.isHeteroAromaticAtom(atom))
				queryDefaults |= (Molecule.cAtomQFAromatic
							  | Molecule.cAtomQFHeteroAromatic);
			else if (mol.isAromaticAtom(atom))
				queryDefaults |= Molecule.cAtomQFAromatic;

			int ringBonds = mol.getAtomRingBondCount(atom);
			if (ringBonds != 0) {
				queryDefaults |= Molecule.cAtomQFNotChain;
				if (ringBonds > 2)
					queryDefaults |= Molecule.cAtomQFNot2RingBonds;
				if (ringBonds > 3)
					queryDefaults |= Molecule.cAtomQFNot3RingBonds;
				}

			int charge = mol.getAtomCharge(atom);
			if (charge < 0)
				queryDefaults |= (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargePos);
			else if (charge > 0)
				queryDefaults |= (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargeNeg);

			// fragment atoms with n neighbours must have
			// at least n neighbours in a matching molecule
			int neighbours = mol.getConnAtoms(atom);
			switch (neighbours) {
			case 0:
				break;
			case 1:
				queryDefaults |= (Molecule.cAtomQFNot0Neighbours);
				break;
			case 2:
				queryDefaults |= (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour);
				break;
			case 3:
				queryDefaults |= (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour | Molecule.cAtomQFNot2Neighbours);
				break;
			default:
				queryDefaults |= (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot4Neighbours);
				break;
				}

			int eValue = mol.getAtomElectronegativeNeighbours(atom);
			switch (eValue) {
				case 0:
					break;
				case 1:
					queryDefaults |= (Molecule.cAtomQFNot0ENeighbours);
					break;
				case 2:
					queryDefaults |= (Molecule.cAtomQFNot0ENeighbours | Molecule.cAtomQFNot1ENeighbour);
					break;
				case 3:
					queryDefaults |= (Molecule.cAtomQFNot0ENeighbours | Molecule.cAtomQFNot1ENeighbour | Molecule.cAtomQFNot2ENeighbours);
					break;
				default:
					queryDefaults |= (Molecule.cAtomQFENeighbours & ~Molecule.cAtomQFNot4ENeighbours);
					break;
				}

			int piElectrons = mol.getAtomPi(atom);
			if (piElectrons > 0)
				queryDefaults |= Molecule.cAtomQFNot0PiElectrons;
			if (piElectrons > 1)
				queryDefaults |= Molecule.cAtomQFNot1PiElectron;
			}

		return queryDefaults;
		}


	/**
	 * Generates inherent feature flags of a given bond.
	 * @param mol molecule or fragment of the SSS
	 * @param bond the bond of which to generate feature flags
	 * @return bond features independent of query features
	 */
	private int getBondQueryDefaults(StereoMolecule mol, int bond) {
		int queryDefaults = 0;

		if (mol.isDelocalizedBond(bond)
		 || mol.getBondType(bond) == Molecule.cBondTypeDelocalized)
			queryDefaults |= Molecule.cBondTypeDelocalized;
		else switch (mol.getBondOrder(bond)) {
			case 0:
				queryDefaults |= Molecule.cBondTypeMetalLigand;
				break;
			case 1:
				queryDefaults |= Molecule.cBondTypeSingle;
				break;
			case 2:
				queryDefaults |= Molecule.cBondTypeDouble;
				break;
			case 3:
				queryDefaults |= Molecule.cBondTypeTriple;
				break;
			case 4:
				queryDefaults |= Molecule.cBondTypeQuadruple;
				break;
			case 5:
				queryDefaults |= Molecule.cBondTypeQuintuple;
				break;
			}

		if (mol.isRingBond(bond))
			queryDefaults |= Molecule.cBondQFRing;
		else if (!mol.isFragment())
			queryDefaults |= Molecule.cBondQFNotRing;

		if (mol.isAromaticBond(bond))
			queryDefaults |= Molecule.cBondQFAromatic;
		else if (!mol.isFragment())
			queryDefaults |= Molecule.cBondQFNotAromatic;

		return queryDefaults;
		}


	private void buildBridgeBondList() {
		mBridgeBondList = null;
		for (int bond=0; bond<mFragment.getBonds(); bond++) {
			if (mFragment.isBondBridge(bond)) {
				if (mBridgeBondList == null)
					mBridgeBondList = new ArrayList<>();
				BridgeBond bridgeBond = new BridgeBond();
				bridgeBond.atom1 = mFragment.getBondAtom(0, bond);
				bridgeBond.atom2 = mFragment.getBondAtom(1, bond);
				bridgeBond.minBridgeSize = mFragment.getBondBridgeMinSize(bond);
				bridgeBond.maxBridgeSize = mFragment.getBondBridgeMaxSize(bond);
				mBridgeBondList.add(bridgeBond);
				}
			}
		}

	private class BridgeBond {
		int atom1,atom2,minBridgeSize,maxBridgeSize;
		}

	private static int[] copyOf(int[] original, int newLength) {
		int[] copy = new int[newLength];
		System.arraycopy(original, 0, copy, 0, Math.min(original.length, newLength));
		return copy;
		}

	private static boolean[] copyOf(boolean[] original, int newLength) {
		boolean[] copy = new boolean[newLength];
		System.arraycopy(original, 0, copy, 0, Math.min(original.length, newLength));
		return copy;
		}
	}
