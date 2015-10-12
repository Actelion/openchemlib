/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.chem;

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

private int mFragmentGraphSize;	// the number of atoms & ring closures in fragment graph
private int[] mFragmentGraphAtom;
private int[] mFragmentGraphParentAtom;
private int[] mFragmentGraphParentBond;
private boolean[] mFragmentGraphIsRingClosure;
private int[] mMatchTable;

// depending on the fragment count mode this may contain atom lists
// of all till now located matching sub-fragments
private ArrayList<int[]> mMatchList,mSortedMatchList;
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
		if (molecule.getAllAtoms() == 0) {
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
		if (fragment.getAllAtoms() == 0 || !fragment.isFragment()) {
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
		}


	private void buildFragmentGraph() {
		// build a graph of the query fragment(s) including ring closures as redundant nodes
		int graphAllocation = mFragment.getBonds()+12;	// 12 is max number of separated fragments within mFragment
		mFragmentGraphAtom = new int[graphAllocation];
		mFragmentGraphParentAtom = new int[graphAllocation];
		mFragmentGraphParentBond = new int[graphAllocation];
		mFragmentGraphIsRingClosure = new boolean[graphAllocation+1];
		boolean[] fragmentAtomUsed = new boolean[mFragment.getAtoms()];
		int current = 0;
		for (int atom=0; atom<mFragment.getAtoms(); atom++) {
			if (!fragmentAtomUsed[atom]) {
			    mFragmentGraphAtom[current] = atom;
			    mFragmentGraphParentBond[current] = -1;
				mFragmentGraphParentAtom[current] = -1;
				int highest = current;
				while (current <= highest) {
					for (int i=0; i<mFragment.getConnAtoms(mFragmentGraphAtom[current]); i++) {
						int candidate = mFragment.getConnAtom(mFragmentGraphAtom[current], i);
						if (candidate != mFragmentGraphParentAtom[current]
						      // if it is a ring closure make sure it is added only once
						 && (!fragmentAtomUsed[candidate] || candidate>mFragmentGraphAtom[current])) {
                            int candidateBond = mFragment.getConnBond(mFragmentGraphAtom[current], i);
                                // don't consider bridge bonds at this state
                            if (!mFragment.isBondBridge(candidateBond)) {
    						    mFragmentGraphAtom[++highest] = candidate;
    						    mFragmentGraphParentAtom[highest] = mFragmentGraphAtom[current];
    						    mFragmentGraphParentBond[highest] = candidateBond;
    							if (fragmentAtomUsed[candidate])
    							    mFragmentGraphIsRingClosure[highest] = true;
    							else
    								fragmentAtomUsed[candidate] = true;
                                }
							}
						}
					while (mFragmentGraphIsRingClosure[++current]);
					}
				}
			}

		mFragmentGraphSize = current;	// this is the real size of the graph
/*
System.out.print("            "); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(mFragmentGraphAtom[i]==-1?"-":Molecule.cAtomLabel[mFragment.getAtomicNo(mFragmentGraphAtom[i])])); System.out.println();
System.out.print("  graphAtom:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(mFragmentGraphAtom[i]==-1?"-":""+mFragmentGraphAtom[i])); System.out.println();
System.out.print(" parentAtom:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(mFragmentGraphParentAtom[i]==-1?"-":""+mFragmentGraphParentAtom[i])); System.out.println();
System.out.print(" parentBond:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(mFragmentGraphParentBond[i]==-1?"-":""+mFragmentGraphParentBond[i])); System.out.println();
System.out.print("ringClosure:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(mFragmentGraphIsRingClosure[i]?"y":"n")); System.out.println();
System.out.println();
*/
		}


	/**
	 * If the match count mode is one of cCountModeFirstMatch, cCountModeOverlapping,
	 * cCountModeRigorous then this method returns an arraylist of all counted matches,
	 * i.e. an int array mapping fragment atoms to molecule atoms. Atoms being part of a
	 * matched bridge bond are naturally not covered by the mapping.
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
	 * @param countMode one of cCountModeExistance, cCountModeFirstMatch, cCountModeOverlapping, cCountModeRigorous
	 * @param matchMode cDefaultMatchMode or combination of cMatchAtomCharge, cMatchAtomMass, cMatchDBondToDelocalized, cMatchAromDBondToDelocalized
	 * @param atomExcluded defines atoms of molecule to be excluded from sub-structure matching
	 * @return count of sub-structure matches of fragment in molecule
	 */
    public int findFragmentInMolecule(int countMode, int matchMode, final boolean[] atomExcluded) {
        mMatchList = null;
        mSortedMatchList = null;

        if (mMolecule == null
   		 || mFragment == null)
   			return 0;

	    if (mFragment.getAtoms() > mMolecule.getAtoms()
		 || mFragment.getBonds() > mMolecule.getBonds())
			return 0;

/*
System.out.print("  molecule:"); for (int i=0; i<moleculeAtoms; i++) System.out.print(" "+i); System.out.println();
System.out.print("     label:"); for (int i=0; i<moleculeAtoms; i++) System.out.print(" "+Molecule.cAtomLabel[mMolecule.getAtomicNo(i)]); System.out.println();
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

		int[] index = new int[mFragmentGraphSize];
		for (int i=0; i<mFragmentGraphSize; i++)
		    index[i] = -1;
			// contains current molecule atom pointer for graph matching,
			// - in case of sub fragment anchor atom: the current molecule atom index matched to the anchor
			// - otherwise the current connAtom index of the parent atom in the matching graph

		int current = 0;
		while (true) {
/*
System.out.print("  index:"); for (int i=0; i<mFragmentGraphSize; i++) System.out.print(" "+(index[i]==-1?"-":""+index[i])); System.out.println();
System.out.print("        "); for (int i=0; i<current; i++) System.out.print("  "); System.out.println(" ^");
System.out.println();
*/
			index[current]++;
			int maxIndex = (mFragmentGraphParentAtom[current] == -1) ?
			        mMolecule.getAtoms() : mMolecule.getConnAtoms(mMatchTable[mFragmentGraphParentAtom[current]]);
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
		    else if (!mFragmentGraphIsRingClosure[current]) {	// current graph position is not an anchor
	            int candidate = mMolecule.getConnAtom(mMatchTable[mFragmentGraphParentAtom[current]], index[current]);
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
	            int candidate = mMolecule.getConnAtom(mMatchTable[mFragmentGraphParentAtom[current]], index[current]);
	            if (candidate == mMatchTable[mFragmentGraphAtom[current]]
	             && areBondsSimilar(mMolecule.getConnBond(mMatchTable[mFragmentGraphParentAtom[current]], index[current]), mFragmentGraphParentBond[current])) {
	                current++;
		        	}
		    	}

		    if (current == mFragmentGraphSize) {
		        if (doTHParitiesMatch()
		         && doEZParitiesMatch()
                 && doBridgeBondsMatch(atomUsed)) {
			        	// we currently have a match
		            if (countMode == cCountModeExistance)
		                return 1;

		            	// add match list
		            int[] matchList = new int[mFragment.getAtoms()];
		            for (int i=0; i<mFragment.getAtoms(); i++)
	                    matchList[i] = mMatchTable[i];

                    if (mMatchList == null)
		                mMatchList = new ArrayList<int[]>();

                    addMatchIfQualifies(matchList, countMode);

                    if (countMode == cCountModeFirstMatch)
                    	return 1;
		        	}

	            current--;
	            if (!mFragmentGraphIsRingClosure[current])
	                atomUsed[mMatchTable[mFragmentGraphAtom[current]]] = false;
		    	}
			}
		
		return (mMatchList == null) ? 0 : mMatchList.size();
		}


	private void addMatchIfQualifies(int[] matchTable, int countMode) {
		if (countMode == cCountModeFirstMatch || countMode == cCountModeRigorous) {
			// count every match (even assymmetrical permutation of same atoms)
            mMatchList.add(matchTable);
			return;
            }

		if (countMode == cCountModeOverlapping) {
            Arrays.sort(matchTable);

            boolean found = false;
            if (mSortedMatchList == null) {
                mSortedMatchList = new ArrayList<int[]>();
                }
            else {
    			for (int f=0; f<mSortedMatchList.size(); f++) {
                    found = true;
                    int[] fragment = mSortedMatchList.get(f);
    				for (int i=0; i<fragment.length; i++) {
    					if (fragment[i] != matchTable[i]) {
    						found = false;
    						break;
    						}
    					}
                    if (found)
                        break;
                    }
                }
			if (!found) {
                mSortedMatchList.add(matchTable);
                matchTable = new int[mFragment.getAtoms()];
                System.arraycopy(mMatchTable, 0, matchTable, 0, matchTable.length);
                mMatchList.add(matchTable);
                }
			return;
			}

//		if (cCountModeSeparated) {
//	not yet supported
//			}

		return;
		}


	protected boolean areAtomsSimilar(int moleculeAtom, int fragmentAtom) {
		int moleculeConnAtoms = mMolecule.getConnAtoms(moleculeAtom);
		int fragmentConnAtoms = mFragment.getConnAtoms(fragmentAtom);
		if (fragmentConnAtoms > moleculeConnAtoms)
			return false;

/* In order to explicitly search and match attachment points (e.g. used by search substructure and replace)
 * this behaviour was changed. TLS 29 Jan 2015
		if (mFragment.getAtomicNo(fragmentAtom) == 0)	// attachment points match on any atom
			return true;
		if (mMolecule.getAtomicNo(moleculeAtom) == 0)	// real fragment atoms cannot match on attachment points in molecule
			return false;
*/
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
	private boolean doTHParitiesMatch() {
	    int esrGroupAtomCount = 0;
	    for (int fragmentAtom=0; fragmentAtom<mFragment.getAtoms(); fragmentAtom++) {
			if ((mFragment.getAtomQueryFeatures(fragmentAtom) & Molecule.cAtomQFMatchStereo) != 0) {
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
	            if ((mFragment.getAtomQueryFeatures(fragmentAtom) & Molecule.cAtomQFMatchStereo) != 0) {
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
	private boolean doEZParitiesMatch() {
	    for (int fragmentBond=0; fragmentBond<mFragment.getBonds(); fragmentBond++) {
			if ((mFragment.getBondQueryFeatures(fragmentBond) & Molecule.cBondQFMatchStereo) != 0) {
				int fragmentParity = mFragment.getBondParity(fragmentBond);

				// always consider as match if fragment bond parity is none
			    if (fragmentParity == Molecule.cBondParityNone)
			   	    continue;

				int fragmentAtom1 = mFragment.getBondAtom(0, fragmentBond);
				int fragmentAtom2 = mFragment.getBondAtom(1, fragmentBond);

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


	private boolean doBridgeBondsMatch(boolean[] moleculeAtomUsed) {
		// Currently bridge bonds are considered as follows:
		// - Bridge bonds in the molecule are not matched to any fragment bond
		// - Bridge bonds in the fragment are not considered in the graph-matching and
		//   checked after an otherwise successful match, whether the shortest path of
		//   unmatched atoms satisfies the min- and max-atom criteria.
		// For fragment and(!) molecule containing bridge bonds we would need a more complex
		// handling: - for every fragment bridge check all unmatched molecule paths and count
		// for every path min- and max-lengths (considering all bridge bonds within path).
		// Consider a match, if min- and max-atoms range fits into fragments bridge bond range.
		// If we have multiple bridge bonds in fragment and partially overlapping bridge matches
		// in the molecule, it would get very nasty...
        if (mBridgeBondList != null) {
            for (BridgeBond bb:mBridgeBondList) {
                int bridgeSize = mMolecule.getPathLength(mMatchTable[bb.atom1], mMatchTable[bb.atom2], bb.maxBridgeSize, moleculeAtomUsed) - 1;
                if (bridgeSize < bb.minBridgeSize
                 || bridgeSize > bb.maxBridgeSize)
                    return false;
                }
            }

        return true;
        }


    /**
     * @param moleculeBond flag list of fragment bond features (features present)
     * @param fragmentBond flag list of molecule bond features (features allowed)
     * @return true if all molecule bond features are present in fragment bond
     */
    protected boolean areBondsSimilar(int moleculeBond, int fragmentBond) {
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


	private void setupAtomAndBondFeatures(int matchMode) {
	    if (!mMoleculeFeaturesValid) {
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

			mMoleculeFeaturesValid = true;
			}

		if (!mFragmentFeaturesValid) {
		    mFragment.ensureHelperArrays(mRequiredHelperLevel);
			int nTotalFragmentAtoms = mFragment.getAtoms();

			mFragmentAtomFeatures = new int[nTotalFragmentAtoms];
			mFragmentAtomType = new int[nTotalFragmentAtoms];

			for (int atom=0; atom<nTotalFragmentAtoms; atom++) {
				mFragmentAtomFeatures[atom] = ((getAtomQueryDefaults(mFragment, atom)
											  | mFragment.getAtomQueryFeatures(atom))
											 & Molecule.cAtomQFSimpleFeatures)
											^ Molecule.cAtomQFNarrowing;
				mFragmentAtomType[atom] = mFragment.getAtomicNo(atom);

				if ((matchMode & cMatchAtomCharge) != 0)
					mFragmentAtomType[atom] += (mFragment.getAtomCharge(atom) + 16) << 8;

				if ((matchMode & cMatchAtomMass) != 0)
					mFragmentAtomType[atom] += mFragment.getAtomMass(atom) << 16;
				}

			int nTotalFragmentBonds = mFragment.getBonds();

			mFragmentBondFeatures = new int[nTotalFragmentBonds];

			for (int bond=0; bond<nTotalFragmentBonds; bond++) {
				mFragmentBondFeatures[bond] = (getBondQueryDefaults(mFragment, bond)
											 | mFragment.getBondQueryFeatures(bond))
											& Molecule.cBondQFSimpleFeatures
										   ^ Molecule.cBondQFNarrowing;

			// match fragment's single/double bonds to delocalized molecule bonds also
				if ((matchMode & cMatchDBondToDelocalized) != 0) {
					if ((mFragmentBondFeatures[bond] & Molecule.cBondQFDouble) != 0)
						mFragmentBondFeatures[bond] |= Molecule.cBondQFDelocalized;
					}
				else if ((matchMode & cMatchAromDBondToDelocalized) != 0) {
					if ((mFragmentBondFeatures[bond] & Molecule.cBondQFDouble) != 0
					 && mFragment.isAromaticBond(bond))
						mFragmentBondFeatures[bond] |= Molecule.cBondQFDelocalized;
					}
				}

			buildFragmentGraph();
            buildBridgeBondList();

			mFragmentFeaturesValid = true;
			}
		}


    /**
     * Generates inherent feature flags of a given atom. 
     * @param mol molecule or fragment of the SSS
     * @param atm the atom of which to generate feature flags
     * @return atom features independent of query features
     */
	private int getAtomQueryDefaults(StereoMolecule mol, int atm) {
		int queryDefaults = 0;

		if (!mol.isFragment()) {
			if (mol.isAromaticAtom(atm))
				queryDefaults |= Molecule.cAtomQFAromatic;
			else
				queryDefaults |= Molecule.cAtomQFNotAromatic;

			int ringBonds = mol.getAtomRingBondCount(atm);
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

			int charge = mol.getAtomCharge(atm);
			if (charge == 0)
	            queryDefaults |= (Molecule.cAtomQFNotChargeNeg | Molecule.cAtomQFNotChargePos);
			else if (charge < 0)
	            queryDefaults |= (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargePos);
			else if (charge > 0)
	            queryDefaults |= (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargeNeg);

			int hydrogens = mol.getAllHydrogens(atm);
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

			int neighbours = mol.getConnAtoms(atm);
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

			int piElectrons = mol.getAtomPi(atm);
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
		else {    // The fragments implicit features are not really necessary,
		          // but may speed up the graph matching.
            if (mol.isAromaticAtom(atm))
                queryDefaults |= Molecule.cAtomQFAromatic;

            int ringBonds = mol.getAtomRingBondCount(atm);
            if (ringBonds != 0) {
                queryDefaults |= Molecule.cAtomQFNotChain;
                if (ringBonds > 2)
                    queryDefaults |= Molecule.cAtomQFNot2RingBonds;
                if (ringBonds > 3)
                    queryDefaults |= Molecule.cAtomQFNot3RingBonds;
                }

    		int charge = mol.getAtomCharge(atm);
    		if (charge < 0)
                queryDefaults |= (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargePos);
    		else if (charge > 0)
                queryDefaults |= (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargeNeg);

            // fragment atoms with n neighbours must have
            // at least n neighbours in a matching molecule
            int neighbours = mol.getConnAtoms(atm);
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

        int piElectrons = mol.getAtomPi(atm);
        if (piElectrons > 0)
            queryDefaults |= Molecule.cAtomQFNot0PiElectrons;
        if (piElectrons > 1)
            queryDefaults |= Molecule.cAtomQFNot1PiElectron;

		return queryDefaults;
		}


	/**
	 * Generates inherent feature flags of a given bond. 
	 * @param mol molecule or fragment of the SSS
	 * @param bnd the bond of which to generate feature flags
	 * @return bond features independent of query features
	 */
	private int getBondQueryDefaults(StereoMolecule mol, int bnd) {
		int queryDefaults = 0;

		if (mol.isDelocalizedBond(bnd)
		 || mol.getBondType(bnd) == Molecule.cBondTypeDelocalized)
			queryDefaults |= Molecule.cBondQFDelocalized;
		else switch (mol.getBondOrder(bnd)) {
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

		if (mol.isRingBond(bnd))
			queryDefaults |= Molecule.cBondQFRing;
		else if (!mol.isFragment())
			queryDefaults |= Molecule.cBondQFNotRing;

		if (mol.isAromaticBond(bnd))
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
	}
