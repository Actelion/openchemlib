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

package com.actelion.research.chem;

import com.actelion.research.util.IntArrayComparator;

import java.util.*;

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
	than the match mode used for index creation. Otherwise index keys may filter out
	molecules which would be considered a match with the less strict matching consitions
	of the atom by atom check.
*/

									// match mode to be used for index creation
	public static final int cIndexMatchMode = cMatchDBondToDelocalized;
	public static final int cDefaultMatchMode = cMatchAromDBondToDelocalized;

	public final static int cCountModeExistance		= 1;
	public final static int cCountModeFirstMatch	= 2;
	public final static int cCountModeOverlapping	= 3;
	public final static int cCountModeRigorous		= 4;

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
	private int[] mMoleculeAtomFeatures;	// flags defining given/required atom features
	private int[] mFragmentAtomFeatures;
	private int[] mMoleculeBondFeatures;	// flags defining given/required bond features
	private int[] mFragmentBondFeatures;

	private int mFragmentExcludeAtoms,mFragmentExcludeBonds;
	private int mFragmentGraphSize;	// the number of wanted atoms & ring closures in fragment graph
	private int mFragmentGraphSizeWithExcludeGroup;	// total number of atoms & ring closures in fragment graph
	private int[] mFragmentGraphAtom;
	private int[] mFragmentGraphParentAtom;
	private int[] mFragmentGraphParentBond;
	private boolean[] mFragmentGraphIsRingClosure;
	private boolean[] mIsExcludeAtom;
	private int[] mFragmentConnAtoms;	// in case of exclude atoms, these are not part of this
	private int[] mMatchTable;

	// depending on the fragment count mode this may contain atom lists
	// of all till now located matching sub-fragments
	private TreeSet<int[]> mSortedMatchSet,mExcludedMatchSet;
	private ArrayList<int[]> mMatchList;
	private ArrayList<BridgeBond> mBridgeBondList;

	private boolean mMoleculeFeaturesValid;
	private boolean mFragmentFeaturesValid;
	private int mRequiredHelperLevel;

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
		mMatchList = new ArrayList<int[]>();
		mSortedMatchSet = new TreeSet<int[]>(new IntArrayComparator());
		mExcludedMatchSet = new TreeSet<int[]>(new IntArrayComparator());
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
		mMatchList = new ArrayList<int[]>();
		mSortedMatchSet = new TreeSet<int[]>(new IntArrayComparator());
		mExcludedMatchSet = new TreeSet<int[]>(new IntArrayComparator());
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
			if ((mFragment.getAtomQueryFeatures(atom) & Molecule.cAtomQFMatchStereo) != 0)
				mRequiredHelperLevel = Molecule.cHelperParities;
		for (int bond=0; bond<mFragment.getBonds(); bond++)
			if ((mFragment.getBondQueryFeatures(bond) & Molecule.cBondQFMatchStereo) != 0)
				mRequiredHelperLevel = Molecule.cHelperParities;

		if (mMoleculeFeaturesValid && mRequiredHelperLevel != Molecule.cHelperRings)
			mMolecule.ensureHelperArrays(mRequiredHelperLevel);

		mFragmentExcludeAtoms = 0;
		mIsExcludeAtom = new boolean[mFragment.getAtoms()];
		for (int atom=0; atom<mFragment.getAtoms(); atom++) {
			mIsExcludeAtom[atom] = ((mFragment.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) != 0);
			if (mIsExcludeAtom[atom])
				mFragmentExcludeAtoms++;
			}

		mFragmentExcludeBonds = 0;
		if (mFragmentExcludeAtoms != 0)
			for (int bond = 0; bond < mFragment.getBonds(); bond++)
				if (mIsExcludeAtom[mFragment.getBondAtom(0, bond)] || mIsExcludeAtom[mFragment.getBondAtom(1, bond)])
					mFragmentExcludeBonds++;
		}


	/**
	 * Build a graph of the query fragment(s) including ring closures as redundant nodes.
	 * If we have exclude groups then these are added at the end of the graph.
	 */
	private void buildFragmentGraph() {
		mFragment.ensureHelperArrays(mRequiredHelperLevel);

		int graphAllocation = mFragment.getBonds() + 12;    // 12 is max number of separated fragments within mFragment
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
					for (int i=0; i<mFragment.getAllConnAtomsPlusMetalBonds(mFragmentGraphAtom[current]); i++) {
						int candidate = mFragment.getConnAtom(mFragmentGraphAtom[current], i);
						if (candidate < mFragment.getAtoms() && !mIsExcludeAtom[candidate])
							highest = tryAddCandidate(current, highest, i, fragmentAtomUsed, fragmentBondUsed);
						}
					while (mFragmentGraphIsRingClosure[++current]);
					}
				}
			}
		mFragmentGraphSize = current;	// this is the real size of the graph not considering exclude atoms

		if (mFragmentExcludeAtoms != 0) {
			// append all exclude atoms and ring closures to non-exclude atom tree
			int highest = current - 1;
			current = 0;
			while (current <= highest) {
				for (int i = 0; i < mFragment.getAllConnAtomsPlusMetalBonds(mFragmentGraphAtom[current]); i++) {
					int candidate = mFragment.getConnAtom(mFragmentGraphAtom[current], i);
					if (candidate < mFragment.getAtoms() && (mIsExcludeAtom[candidate] || mIsExcludeAtom[mFragmentGraphAtom[current]]))
						highest = tryAddCandidate(current, highest, i, fragmentAtomUsed, fragmentBondUsed);
					}
				while (mFragmentGraphIsRingClosure[++current]) ;
				}

			// there may still be exclude groups as separated fragments
			for (int atom = 0; atom < mFragment.getAtoms(); atom++) {
				if (mIsExcludeAtom[atom]
						&& !fragmentAtomUsed[atom]) {
					mFragmentGraphAtom[current] = atom;
					mFragmentGraphParentBond[current] = -1;
					mFragmentGraphParentAtom[current] = -1;
					highest = current;
					while (current <= highest) {
						for (int i = 0; i < mFragment.getAllConnAtomsPlusMetalBonds(mFragmentGraphAtom[current]); i++)
							if (mFragment.getConnAtom(mFragmentGraphAtom[current], i) < mFragment.getAtoms())
								highest = tryAddCandidate(current, highest, i, fragmentAtomUsed, fragmentBondUsed);

						while (mFragmentGraphIsRingClosure[++current]) ;
						}
					}
				}
			}

		mFragmentGraphSizeWithExcludeGroup = current;	// this is the real size of the graph
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
	private int tryAddCandidate(int current, int highest, int i, boolean[] fragmentAtomUsed, boolean[] fragmentBondUsed) {
		int candidate = mFragment.getConnAtom(mFragmentGraphAtom[current], i);
		if (candidate != mFragmentGraphParentAtom[current]) {
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
	 * matched bridge bond are naturally not covered by the mapping.<br>
	 * Note: If some query fragment atoms are marked as exclude group, then the respective
	 * matchlist values are -1.
	 * @return list of distinct counted matches.
	 */
	public ArrayList<int[]> getMatchList() {
		return mMatchList;
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
		return (findFragmentInMolecule(cCountModeExistance, mDefaultMatchMode) > 0);
		}


	/**
	 * Fastest check, whether the molecule contains the fragment.
	 * This method uses cCountModeExistance and therefore does not create
	 * any retrievable match list, i.e. mapping from fragment atoms to molecule atoms.
	 * @param matchMode cDefaultMatchMode or combination of cMatchAtomCharge, cMatchAtomMass, cMatchDBondToDelocalized, cMatchAromDBondToDelocalized
	 * @return whether fragment was found as sub-structure in molecule
	 */
	public boolean isFragmentInMolecule(int matchMode) {
		return (findFragmentInMolecule(cCountModeExistance, matchMode) > 0);
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
	 * @param countMode one of cCountModeExistance, cCountModeFirstMatch, cCountModeOverlapping, cCountModeRigorous
	 * @param matchMode cDefaultMatchMode or combination of cMatchAtomCharge, cMatchAtomMass, cMatchDBondToDelocalized, cMatchAromDBondToDelocalized
	 * @param atomExcluded defines atoms of molecule to be excluded from sub-structure matching
	 * @return count of sub-structure matches of fragment in molecule
	 */
	public int findFragmentInMolecule(int countMode, int matchMode, final boolean[] atomExcluded) {
		mMatchList = new ArrayList<int[]>();
		mSortedMatchSet.clear();
		mExcludedMatchSet.clear();

		if (mMolecule == null
   		 || mFragment == null)
   			return 0;

		if (mFragment.getAtoms() - mFragmentExcludeAtoms > mMolecule.getAtoms()
		 || mFragment.getBonds() - mFragmentExcludeBonds > mMolecule.getBonds())
			return 0;

		if (mFragment.getAtoms() - mFragmentExcludeAtoms == 0)
			return 0;

/*
System.out.print("  molecule:"); for (int i=0; i<moleculeAtoms; i++) System.out.print(" "+i); System.out.println();
System.out.print("	 label:"); for (int i=0; i<moleculeAtoms; i++) System.out.print(" "+Molecule.cAtomLabel[mMolecule.getAtomicNo(i)]); System.out.println();
for (int j=0; j<moleculeAtoms; j++)
{ System.out.print("conns of "+j+":"); for (int i=0; i<mMolecule.getConnAtoms(j); i++) System.out.print(" "+mMolecule.getConnAtom(j, i)); System.out.println(); }
System.out.println();
*/

		setupAtomAndBondFeatures(matchMode);

		// atom usage mask in mMolecule
		boolean[] atomUsed = new boolean[mMolecule.getAtoms()];
		if (atomExcluded != null)
			for (int atom=0; atom<mMolecule.getAtoms(); atom++)
				atomUsed[atom] = atomExcluded[atom];

		// mMolecule atom currently matched on mFragment atom
		mMatchTable = new int[mFragment.getAtoms()];
		Arrays.fill(mMatchTable, -1);	// to mark exclude group atoms

		int[] index = new int[mFragmentGraphSizeWithExcludeGroup];
		Arrays.fill(index, -1);
		// contains current molecule atom pointer for graph matching,
		// - in case of sub fragment anchor atom: the current molecule atom index matched to the anchor
		// - otherwise the current connAtom index of the parent atom in the matching graph

		int current = 0;
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

			if (current == mFragmentGraphSize) {
				if (doTHParitiesMatch(false)
				 && doEZParitiesMatch(false)
				 && doBridgeBondsMatch(atomUsed, false)) {
					// we currently have a match
					if (countMode == cCountModeExistance && mFragmentExcludeAtoms == 0)
						return 1;

					boolean isExcludedMatch = false;
					if (mFragmentExcludeAtoms != 0) {
						// If we have (an) exclude group(s) then we need to check all permutations
						// of symmetrical matches (those covering the same atoms), whether the graph
						// matching of one of those can be extended to match the exclude group also.
						// If at least one permutations has the exclude group attached, then we need
						// to eliminate all matches on the same atom list off the match list.
						// Therefore we cannot return after the first match is found and must check
						// every match, whether it can be extended to include the exclude group(s).
						// In this case we call it an excluded match.
						int[] sortedMatch = copyOf(mMatchTable, mMatchTable.length);
						Arrays.sort(sortedMatch);
						if (mExcludedMatchSet.contains(sortedMatch)) {
							isExcludedMatch = true;
							}
						else if (doExcludeGroupsMatch(atomUsed, index)) {
							mExcludedMatchSet.add(sortedMatch);
							Comparator<int[]> comparator = new IntArrayComparator();
							int[] tempMatch = new int[sortedMatch.length];
							for (int i=mMatchList.size()-1; i>=0; i--) {
								int[] match = mMatchList.get(i);
								System.arraycopy(match, 0, tempMatch, 0, tempMatch.length);
								Arrays.sort(tempMatch);
								if (comparator.compare(tempMatch, sortedMatch) == 0)
									mMatchList.remove(i);
								}
							isExcludedMatch = true;
							}
						}

					if (!isExcludedMatch) {
						addMatchIfQualifies(countMode);

						if (countMode == cCountModeFirstMatch && mFragmentExcludeAtoms == 0)
							return 1;
						}
					}

				current--;
				if (!mFragmentGraphIsRingClosure[current])
					atomUsed[mMatchTable[mFragmentGraphAtom[current]]] = false;
				}
			}

		return mMatchList.size();
		}


	private void addMatchIfQualifies(int countMode) {
		int[] match = copyOf(mMatchTable, mMatchTable.length);

		if (countMode == cCountModeFirstMatch
		 || countMode == cCountModeRigorous) {
			// count every match (even permutations of same atoms)
			mMatchList.add(match);
			return;
			}

		if (mFragmentExcludeAtoms != 0	// store matches as indication that we have found something
		 || countMode == cCountModeOverlapping) {
			Arrays.sort(match);
			if (!mSortedMatchSet.contains(match)) {
				mSortedMatchSet.add(match);
				mMatchList.add(copyOf(mMatchTable, mMatchTable.length));
				}
			return;
			}

//		if (cCountModeSeparated) {
//	not yet supported
//			}

		return;
		}


	public boolean areAtomsSimilar(int moleculeAtom, int fragmentAtom) {
		int moleculeConnAtoms = mMolecule.getConnAtoms(moleculeAtom);
		int fragmentConnAtoms = mFragmentConnAtoms[fragmentAtom];

		if (fragmentConnAtoms > moleculeConnAtoms)
			return false;

		int moleculeQF = mMolecule.getAtomQueryFeatures(moleculeAtom);
		int fragmentQF = mFragment.getAtomQueryFeatures(fragmentAtom);

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

		if ((mMoleculeAtomFeatures[moleculeAtom] & ~mFragmentAtomFeatures[fragmentAtom]) != 0)
			return false;

		if (mFragment.getAtomCharge(fragmentAtom) != 0
		 && mFragment.getAtomCharge(fragmentAtom) != mMolecule.getAtomCharge(moleculeAtom))
			return false;
		if (mFragment.getAtomMass(fragmentAtom) != 0
		 && mFragment.getAtomMass(fragmentAtom) != mMolecule.getAtomMass(moleculeAtom))
			return false;
		int ringSize = (mFragment.getAtomQueryFeatures(fragmentAtom) & Molecule.cAtomQFRingSize) >> Molecule.cAtomQFRingSizeShift;
		if (ringSize != 0) {
			if (mMolecule.isFragment()
			 && ringSize == (mMolecule.getAtomQueryFeatures(fragmentAtom) & Molecule.cAtomQFRingSize) >> Molecule.cAtomQFRingSizeShift)
				return true;

			boolean found = false;
			RingCollection ringSet = mMolecule.getRingSet();
			for (int i=0; i<ringSet.getSize(); i++) {
				if (ringSet.getRingSize(i) == ringSize) {
					if (ringSet.isAtomMember(i, moleculeAtom)) {
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
	private boolean doTHParitiesMatch(boolean isExcludeGroup) {
		int esrGroupAtomCount = 0;
		for (int fragmentAtom=0; fragmentAtom<mFragment.getAtoms(); fragmentAtom++) {
			if (mIsExcludeAtom[fragmentAtom] == isExcludeGroup
			 && (mFragment.getAtomQueryFeatures(fragmentAtom) & Molecule.cAtomQFMatchStereo) != 0) {
				int moleculeAtom = mMatchTable[fragmentAtom];
				int fragmentParity = mFragment.getAtomParity(fragmentAtom);
				int moleculeParity = mMolecule.getAtomParity(moleculeAtom);

					// always consider as match if fragment atom is no stereo center
				if (fragmentParity == Molecule.cAtomParityNone)
			   		continue;

				// consider as match if assymetric fragment atom matches on non-stereo-center
				if (moleculeParity == Molecule.cAtomParityNone)
			   		continue;

				// unknown molecule centers need to match everything because
				// parities from idcodes don't include them, i.e. depending
				// on the source of the parities, unknown centers may look as
				// no centers. Both must retrieve the same results.
				if (fragmentParity == Molecule.cAtomParityUnknown)
			   		continue;

				if (moleculeParity == Molecule.cAtomParityUnknown)
			   		continue;

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
				if (mIsExcludeAtom[fragmentAtom] == isExcludeGroup
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
	private boolean doEZParitiesMatch(boolean isExcludeGroup) {
		for (int fragmentBond=0; fragmentBond<mFragment.getBonds(); fragmentBond++) {
			if ((mFragment.getBondQueryFeatures(fragmentBond) & Molecule.cBondQFMatchStereo) != 0) {
				int fragmentParity = mFragment.getBondParity(fragmentBond);

				// always consider as match if fragment bond parity is none
				if (fragmentParity == Molecule.cBondParityNone)
			   		continue;

				int fragmentAtom1 = mFragment.getBondAtom(0, fragmentBond);
				int fragmentAtom2 = mFragment.getBondAtom(1, fragmentBond);

				if ((mIsExcludeAtom[fragmentAtom1] || mIsExcludeAtom[fragmentAtom2]) == isExcludeGroup) {
					int moleculeAtom1 = mMatchTable[fragmentAtom1];
					int moleculeAtom2 = mMatchTable[fragmentAtom2];
					int moleculeBond = mMolecule.getBond(moleculeAtom1, moleculeAtom2);

					// consider as match if an E/Z-bond with defined parity atom matches on one with no parity
					int moleculeParity = mMolecule.getBondParity(moleculeBond);
					if (moleculeParity == Molecule.cBondParityNone)
						continue;

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


	private boolean isEZParityInversion(int fragmentBond, int moleculeBond) {
		boolean inversion = false;
		for (int i=0; i<2; i++) {
			int fragmentAtom = mFragment.getBondAtom(i, fragmentBond);
			int moleculeAtom = mMatchTable[fragmentAtom];
			if (mFragment.getConnAtoms(fragmentAtom) == 2) {
				if (mMolecule.getConnAtoms(moleculeAtom) == 2)
					continue;

				int fragmentNeighbour = -1;
				for (int j=0; j<2; j++)
					if (mFragment.getConnBond(fragmentAtom, j) != fragmentBond)
						fragmentNeighbour = mFragment.getConnAtom(fragmentAtom, j);

				int moleculeNeighbours = 0;
				int[] moleculeNeighbour = new int[2];
				for (int j=0; j<3; j++)
					if (mMolecule.getConnBond(moleculeAtom, j) != moleculeBond)
						moleculeNeighbour[moleculeNeighbours++] = mMolecule.getConnAtom(moleculeAtom, j);

				if (mMatchTable[fragmentNeighbour] != moleculeNeighbour[0])
					inversion = !inversion;
				}
			else if (mFragment.getConnAtoms(fragmentAtom) == 3
	   			  && mMolecule.getConnAtoms(moleculeAtom) == 3) {
				int[] fragmentNeighbour = new int[2];
				int fragmentNeighbours = 0;
				for (int j=0; j<3; j++)
					if (mFragment.getConnBond(fragmentAtom, j) != fragmentBond)
						fragmentNeighbour[fragmentNeighbours++] = mFragment.getConnAtom(fragmentAtom, j);
				if ((mMatchTable[fragmentNeighbour[0]] > mMatchTable[fragmentNeighbour[1]])
				  ^ (fragmentNeighbour[0] > fragmentNeighbour[1]))
					inversion = !inversion;
				}
			}
		return inversion;
		}


	/**
	 * Starting from a full match of the fragment without exclude groups, this method
	 * continues the graph matching to find, whether the attached exclude group(s) can
	 * also be matched.
	 * @param atomUsed
	 * @return
	 */
	private boolean doExcludeGroupsMatch(boolean[] atomUsed, int[] index) {
		for (int i=mFragmentGraphSize; i<mFragmentGraphSizeWithExcludeGroup; i++)
			index[i] = -1;

		int current = mFragmentGraphSize;

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
				if (current == mFragmentGraphSize)
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

			if (current == mFragmentGraphSizeWithExcludeGroup) {
				if (doTHParitiesMatch(true)
				 && doEZParitiesMatch(true)
				 && doBridgeBondsMatch(atomUsed, true)) {

					// remove match table entries for exclude atoms
					for (int atom=0; atom<mFragment.getAtoms(); atom++) {
						if (mIsExcludeAtom[atom]) {
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
	 * @return
	 */
	private boolean doBridgeBondsMatch(boolean[] moleculeAtomUsed, boolean isExcludeFragment) {
		if (mBridgeBondList != null) {
			for (BridgeBond bb:mBridgeBondList) {
				if ((mIsExcludeAtom[bb.atom1] || mIsExcludeAtom[bb.atom2]) == isExcludeFragment) {
					int bridgeSize = mMolecule.getPathLength(mMatchTable[bb.atom1], mMatchTable[bb.atom2], bb.maxBridgeSize+1, moleculeAtomUsed) - 1;
					if (bridgeSize < bb.minBridgeSize
					 || bridgeSize > bb.maxBridgeSize)
						return false;
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
		if ((mMoleculeBondFeatures[moleculeBond] & ~mFragmentBondFeatures[fragmentBond]) != 0)
			return false;

		int ringSize = (mFragment.getBondQueryFeatures(fragmentBond) & Molecule.cBondQFRingSize) >> Molecule.cBondQFRingSizeShift;
		if (ringSize != 0) {
			if (mMolecule.isFragment()
			 && ringSize == (mMolecule.getBondQueryFeatures(fragmentBond) & Molecule.cBondQFRingSize) >> Molecule.cBondQFRingSizeShift)
				return true;

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
		mMoleculeAtomFeatures = new int[nTotalMoleculeAtoms];

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
		int[] atomFeaturesWithoutExcludeAtoms = null;
		int[] bondFeaturesWithoutExcludeAtoms = null;
		int[] atomTypeWithoutExcludeAtoms = null;

		mFragment.ensureHelperArrays(mRequiredHelperLevel);
		mFragmentConnAtoms = new int[mFragment.getAtoms()];
		for (int atom=0; atom<mFragment.getAtoms(); atom++)
			mFragmentConnAtoms[atom] = mFragment.getConnAtoms(atom);

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

		mFragmentAtomFeatures = new int[fragment.getAtoms()];
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

		int nTotalFragmentBonds = fragment.getBonds();

		mFragmentBondFeatures = new int[fragment.getBonds()];

		for (int bond=0; bond<nTotalFragmentBonds; bond++) {
			mFragmentBondFeatures[bond] = (getBondQueryDefaults(fragment, bond)
					| fragment.getBondQueryFeatures(bond))
					& Molecule.cBondQFSimpleFeatures
					^ Molecule.cBondQFNarrowing;

			// match fragment's single/double bonds to delocalized molecule bonds also
			if ((matchMode & cMatchDBondToDelocalized) != 0) {
				if ((mFragmentBondFeatures[bond] & Molecule.cBondQFDouble) != 0)
					mFragmentBondFeatures[bond] |= Molecule.cBondQFDelocalized;
				}
			else if ((matchMode & cMatchAromDBondToDelocalized) != 0) {
				if ((mFragmentBondFeatures[bond] & Molecule.cBondQFDouble) != 0
						&& fragment.isAromaticBond(bond))
					mFragmentBondFeatures[bond] |= Molecule.cBondQFDelocalized;
				}
			}
		}

	/**
	 * Generates inherent feature flags of a given atom.
	 * @param mol molecule or fragment of the SSS
	 * @param atom the atom of which to generate feature flags
	 * @return atom features independent of query features
	 */
	private int getAtomQueryDefaults(StereoMolecule mol, int atom) {
		int queryDefaults = 0;

		if (!mol.isFragment()) {
			if (mol.isAromaticAtom(atom))
				queryDefaults |= Molecule.cAtomQFAromatic;
			else
				queryDefaults |= Molecule.cAtomQFNotAromatic;

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
			if (mol.isAromaticAtom(atom))
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
			}

		int piElectrons = mol.getAtomPi(atom);
		if (piElectrons > 0)
			queryDefaults |= Molecule.cAtomQFNot0PiElectrons;
		if (piElectrons > 1)
			queryDefaults |= Molecule.cAtomQFNot1PiElectron;

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
			queryDefaults |= Molecule.cBondQFDelocalized;
		else switch (mol.getBondOrder(bond)) {
			case 0:
				queryDefaults |= Molecule.cBondTypeMetalLigand;
				break;
			case 1:
				queryDefaults |= Molecule.cBondQFSingle;
				break;
			case 2:
				queryDefaults |= Molecule.cBondQFDouble;
				break;
			case 3:
				queryDefaults |= Molecule.cBondQFTriple;
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
					mBridgeBondList = new ArrayList<BridgeBond>();
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
     System.arraycopy(original, 0, copy, 0,
                      Math.min(original.length, newLength));
     return copy;
 }

	}
