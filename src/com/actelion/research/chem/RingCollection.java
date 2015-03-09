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

import java.util.ArrayList;

public class RingCollection {
    public static final int MAX_SMALL_RING_SIZE = 7;
//  public static final int MAX_LARGE_RING_SIZE = 24;	// disabled ring size limit, because atom ring flags must be
    													// reliable, e.g. for substituent detection. TLS 20130613

    private static final int MODE_SMALL_RINGS = 1;
    private static final int MODE_LARGE_RINGS = 2;
    private static final int MODE_AROMATICITY = 4;
    public static final int MODE_SMALL_RINGS_ONLY = MODE_SMALL_RINGS;
    public static final int MODE_SMALL_AND_LARGE_RINGS = MODE_SMALL_RINGS
                                                       | MODE_LARGE_RINGS;
    public static final int MODE_SMALL_RINGS_AND_AROMATICITY = MODE_SMALL_RINGS
                                                             | MODE_AROMATICITY;
    public static final int MODE_SMALL_AND_LARGE_RINGS_AND_AROMATICITY = MODE_SMALL_RINGS
                                                                       | MODE_LARGE_RINGS
                                                                       | MODE_AROMATICITY;

    private ExtendedMolecule mMol;
    private ArrayList<int[]> mRingAtomSet;
    private ArrayList<int[]> mRingBondSet;
    private int[] mAtomRingSize;
    private int[] mBondRingSize;
    private int[] mHeteroPosition;
    private boolean[] mAromaticityHandled;
    private boolean[] mIsAromatic;
    private boolean[] mIsDelocalized;

    /**
     * Generates the complete set of small rings, which contains all rings
     * up to 7 members.<br> If mode includes LARGE_RINGS, then it determines
     * for every atom and bond the size of the smallest ring, which they are
     * a member of.<br>If mode includes AROMATICITY then every small ring
     * is checked, whether it is aromatic.
     * @param mol
     * @param mode one of the public MODE_ options
     */
	public RingCollection(ExtendedMolecule mol, int mode) {
		mMol = mol;
		mRingAtomSet = new ArrayList<int[]>();
        mRingBondSet = new ArrayList<int[]>();

        mAtomRingSize = new int[mMol.getAtoms()];
        mBondRingSize = new int[mMol.getBonds()];

        mMol.ensureHelperArrays(ExtendedMolecule.cHelperNeighbours);

        boolean[] isConfirmedChainAtom = new boolean[mMol.getAtoms()];
        boolean[] isConfirmedChainBond = new boolean[mMol.getBonds()];

        boolean found;
        do {    // detect atoms of side chains as non-ring-atoms
            found = false;
            for (int atom=0; atom<mMol.getAtoms(); atom++) {
                if (!isConfirmedChainAtom[atom]) {
                    int potentialRingNeighbours = 0;
                    for (int i=0; i<mMol.getConnAtoms(atom); i++)
                        if (!isConfirmedChainAtom[mMol.getConnAtom(atom, i)])
                            potentialRingNeighbours++;
    
                    if (potentialRingNeighbours < 2) {
                        isConfirmedChainAtom[atom] = true;
                        for (int i=0; i<mMol.getConnAtoms(atom); i++)
                            isConfirmedChainBond[mMol.getConnBond(atom, i)] = true;

                        found = true;
                        }
                    }
                }
            } while (found);

                // generate graph of potential ring atoms to find ring closure bonds
        int startAtom = 0;  // simply take the first potential ring atom as graph base
        while ((startAtom < mMol.getAtoms()) && isConfirmedChainAtom[startAtom])
            startAtom++;
        if (startAtom == mMol.getAtoms())
            return;         // no rings found

                // find all rings with less than 8 members of all closure bonds
        int graphAtom[] = new int[mMol.getAtoms()];
        graphAtom[0] = startAtom;
        int fragmentNo[] = new int[mMol.getAtoms()];
        fragmentNo[startAtom] = 1;
        int current = 0;
        int highest = 0;
        int noOfFragments = 1;
        while (current <= highest) {
            for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
                int candidate = mMol.getConnAtom(graphAtom[current], i);

                if (fragmentNo[candidate] != 0) {   // closure bond
                    addSmallRingsToSet(mMol.getConnBond(graphAtom[current], i), isConfirmedChainAtom);
                    continue;
                    }

                if (!isConfirmedChainAtom[candidate]) {
                    fragmentNo[candidate] = noOfFragments;
                    graphAtom[++highest] = candidate;
                    }
                }
            current++;
            if (current > highest) {
                    // if run out of atoms look for new base atom of other fragment
                for (int atom=0; atom<mMol.getAtoms(); atom++) {
                    if (fragmentNo[atom] == 0 && !isConfirmedChainAtom[atom]) {
                        fragmentNo[atom] = ++noOfFragments;
                        graphAtom[++highest] = atom;
                        break;
                        }
                    }
                }
            }

        if ((mode & MODE_AROMATICITY) != 0) {
            mAromaticityHandled = new boolean[mRingAtomSet.size()];
            mIsAromatic = new boolean[mRingAtomSet.size()];
            mIsDelocalized = new boolean[mRingAtomSet.size()];
            mHeteroPosition = new int[mRingAtomSet.size()];
            determineAromaticity();
            }

        // find large rings by examining every potential ring bond
        // which is not a member of a small ring
        if ((mode & MODE_LARGE_RINGS) != 0) {
            for (int bond=0; bond<mMol.getBonds(); bond++) {
                if (!isConfirmedChainBond[bond]) {
                    int ringAtom[] = findSmallestRing(bond, isConfirmedChainAtom);
                    if (ringAtom != null)
                        updateRingSizes(ringAtom, getRingBonds(ringAtom));
                    }
                }
            }
        }


    private int[] findSmallestRing(int bond, boolean[] isConfirmedChainAtom) {
        // find smallest ring of given bond
    int atom1 = mMol.getBondAtom(0, bond);
    int atom2 = mMol.getBondAtom(1, bond);
    int graphAtom[] = new int[mMol.getAtoms()];
    int graphLevel[] = new int[mMol.getAtoms()];
    int graphParent[] = new int[mMol.getAtoms()];
    graphAtom[0] = atom1;
    graphAtom[1] = atom2;
    graphLevel[atom1] = 1;
    graphLevel[atom2] = 2;
    graphParent[atom1] = -1;
    graphParent[atom2] = atom1;
    int current = 1;
    int highest = 1;
    while (current <= highest) {
//        if (graphLevel[graphAtom[current]] > MAX_LARGE_RING_SIZE)
//            return null;		// disabled ring size limit, TLS 20130613
        for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
            int candidate = mMol.getConnAtom(graphAtom[current], i);
            if ((current > 1) && candidate == atom1) {
                int ringAtom[] = new int[graphLevel[graphAtom[current]]];
                int atom = graphAtom[current];
                for (int j=0; j<ringAtom.length; j++) {
                    ringAtom[j] = atom;
                    atom = graphParent[atom];
                    }
                return ringAtom;
                }
            if (graphLevel[candidate] == 0 && !isConfirmedChainAtom[candidate]) {
                graphAtom[++highest] = candidate;
                graphLevel[candidate] = graphLevel[graphAtom[current]] + 1;
                graphParent[candidate] = graphAtom[current];
                }
            }
        current++;
        }
    return null;
    }


    /**
     * An atom's ring size is the size of the smallest ring the atom is a member of.
     * If the atom doesn't belong to any ring the ring size is 0. If an atom is member
     * of rings larger than 7 members only and if the mode parameter of the constructor
     * didn't include LARGE_RINGS, then the returned ring size is also 0.
     * @param atom
     * @return ring size or 0
     */
    public int getAtomRingSize(int atom) {
        return mAtomRingSize[atom];
        }


    /**
     * A bond's ring size is the size of the smallest ring the bond is a member of.
     * If the bond doesn't belong to any ring the ring size is 0. If a bond is member
     * of rings larger than 7 members only and if the mode parameter of the constructor
     * didn't include LARGE_RINGS, then the returned ring size is also 0.
     * @param bond
     * @return ring size or 0
     */
    public int getBondRingSize(int bond) {
        return mBondRingSize[bond];
        }


    private void addSmallRingsToSet(int closureBond, boolean[] isConfirmedChainAtom) {
        int atom1 = mMol.getBondAtom(0, closureBond);
        int atom2 = mMol.getBondAtom(1, closureBond);
        int graphAtom[] = new int[mMol.getAtoms()];
        int graphLevel[] = new int[mMol.getAtoms()];
        int graphParent[] = new int[mMol.getAtoms()];
        graphAtom[0] = atom1;
        graphAtom[1] = atom2;
        graphLevel[atom1] = 1;
        graphLevel[atom2] = 2;
        graphParent[atom1] = -1;
        graphParent[atom2] = atom1;
        int current = 1;
        int highest = 1;
        while (current <= highest) {
            if (graphLevel[graphAtom[current]] > MAX_SMALL_RING_SIZE)
                return;
            for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
                int candidate = mMol.getConnAtom(graphAtom[current], i);
                if ((current > 1) && candidate == atom1) {
                    int ringAtom[] = new int[graphLevel[graphAtom[current]]];
                    int atom = graphAtom[current];
                    for (int j=0; j<ringAtom.length; j++) {
                        ringAtom[j] = atom;
                        atom = graphParent[atom];
                        }
                    addRingIfNew(ringAtom);
                    continue;
                    }
                if (graphLevel[candidate] == 0 && !isConfirmedChainAtom[candidate]) {
                    graphAtom[++highest] = candidate;
                    graphLevel[candidate] = graphLevel[graphAtom[current]] + 1;
                    graphParent[candidate] = graphAtom[current];
                    }
                }
            current++;
            }
        }


	private void addRingIfNew(int ringAtom[]) {
		int ringAtoms = ringAtom.length;
		int lowAtom = mMol.getMaxAtoms();
		int lowIndex = 0;
		for (int i=0; i<ringAtoms; i++) {
			if (lowAtom > ringAtom[i]) {
				lowAtom = ringAtom[i];
				lowIndex = i;
				}
			}

		int sortedRing[] = new int[ringAtoms];
		int leftIndex = (lowIndex > 0) ? lowIndex - 1 : ringAtoms - 1;
		int rightIndex = (lowIndex < ringAtoms - 1) ? lowIndex + 1 : 0;
		boolean inverse = (ringAtom[leftIndex] < ringAtom[rightIndex]);
		for (int i=0; i<ringAtoms; i++) {
			sortedRing[i] = ringAtom[lowIndex];
			if (inverse) {
				if (--lowIndex < 0)
					lowIndex = ringAtoms - 1;
				}
			else {
				if (++lowIndex == ringAtoms)
					lowIndex = 0;
				}
			}

		for (int i=0; i<mRingAtomSet.size(); i++) {
			int ringOfSet[] = mRingAtomSet.get(i);
			if (ringOfSet.length != ringAtoms)
				continue;
			boolean equal = true;
			for (int j=0; j<ringAtoms; j++) {
				if (ringOfSet[j] != sortedRing[j]) {
					equal = false;
					break;
					}
				}
			if (equal)
				return;
			}

		mRingAtomSet.add(sortedRing);
        int[] ringBond = getRingBonds(sortedRing);
        mRingBondSet.add(ringBond);

        updateRingSizes(sortedRing, ringBond);
	    }


    public int getSize() {
		return mRingAtomSet.size();
		}


	public int[] getRingAtoms(int ringNo) {
		return mRingAtomSet.get(ringNo);
		}


	public int[] getRingBonds(int ringNo) {
		return mRingBondSet.get(ringNo);
		}


	public int getRingSize(int ringNo) {
		return mRingBondSet.get(ringNo).length;
		}


	/**
	 * Return whether the ring is considered aromatic.
	 * If the mode parameter passed to the constructor didn't include AROMATICITY, then
	 * then false is returned.
	 * @param ringNo
	 * @return
	 */
	public boolean isAromatic(int ringNo) {
		return mIsAromatic[ringNo];
		}


	/**
	 * Return whether the ring is considered delocalized, which are 6-membered
	 * aromatic rings with no preference concerning where the double bonds are located.
	 * Pyrrole bonds are not considered delocalized.
	 * If the mode parameter passed to the constructor didn't include AROMATICITY, then
	 * then false is returned.
	 * @param ringNo
	 * @return
	 */
	public boolean isDelocalized(int ringNo) {
		return mIsDelocalized[ringNo];
		}


	public int getAtomIndex(int ringNo, int atom) {
		int[] ringAtom = mRingAtomSet.get(ringNo);
		for (int i=0; i<ringAtom.length; i++)
			if (atom == ringAtom[i])
				return i;

		return -1;
		}

	
	public int getBondIndex(int ringNo, int bond) {
		int[] ringBond = mRingBondSet.get(ringNo);
		for (int i=0; i<ringBond.length; i++)
			if (bond == ringBond[i])
				return i;

		return -1;
		}


	/**
	 * Adds or subtracts the ring size from index to move it
	 * into the valid range from 0 to ringSize-1.
	 * @param ringNo
	 * @param index
	 * @return
	 */
	public int validateMemberIndex(int ringNo, int index) {
		int ringSize = mRingBondSet.get(ringNo).length;
		while (index >= ringSize)
			index -= ringSize;
		while (index < 0)
			index += ringSize;
		return index;
		}

	/**
     * Returns the position of the electron pair providing hetero atom
     * or carbenium atom in case of 5-membered, respective 7-membered
     * aromatic ring.
     * @param ringNo
     * @return position index referring to ringAtom array
     */
    public int getHeteroPosition(int ringNo) {
        return mHeteroPosition[ringNo];
        }


    public boolean isAtomMember(int ringNo, int atom) {
		int[] ringAtom = mRingAtomSet.get(ringNo);
		for (int i=0; i<ringAtom.length; i++)
			if (atom == ringAtom[i])
				return true;

		return false;
		}


    public boolean isBondMember(int ringNo, int bond) {
		int[] ringBond = mRingBondSet.get(ringNo);
		for (int i=0; i<ringBond.length; i++)
			if (bond == ringBond[i])
				return true;

		return false;
		}


    private void updateRingSizes(int[] ringAtom, int[] ringBond) {
        int ringSize = ringAtom.length;
        for (int i=0; i<ringSize; i++)
            if (mAtomRingSize[ringAtom[i]] == 0
             || mAtomRingSize[ringAtom[i]] > ringSize)
                mAtomRingSize[ringAtom[i]] = ringSize;

        for (int i=0; i<ringSize; i++)
            if (mBondRingSize[ringBond[i]] == 0
             || mBondRingSize[ringBond[i]] > ringSize)
                mBondRingSize[ringBond[i]] = ringSize;
        }

    
    private int[] getRingBonds(int[] ringAtom) {
        int ringAtoms = ringAtom.length;
        int ringBond[] = new int[ringAtoms];
        for (int i=0; i<ringAtoms; i++) {
            int atom = (i == ringAtoms - 1) ? ringAtom[0] : ringAtom[i+1];
            for (int j=0; j<mMol.getConnAtoms(ringAtom[i]); j++) {
                if (mMol.getConnAtom(ringAtom[i],j) == atom) {
                    ringBond[i] = mMol.getConnBond(ringAtom[i],j);
                    break;
                    }
                }
            }
        return ringBond;
        }


    private void determineAromaticity() {
        int[][] annelatedRing = new int[mRingAtomSet.size()][];
        for (int i=0; i<mRingAtomSet.size(); i++) {
            annelatedRing[i] = new int[mRingAtomSet.get(i).length];
            for (int j=0; j<mRingAtomSet.get(i).length; j++)
                annelatedRing[i][j] = -1;
            }

        int[] ringMembership = new int[mMol.getBonds()];
        for (int ring=0; ring<mRingBondSet.size(); ring++) {
            int[] ringBond = mRingBondSet.get(ring);
            if (ringBond.length >= 5 && ringBond.length <= 7) {
                for (int i=0; i<ringBond.length; i++) {
                    int bond = ringBond[i];
                    if (mMol.getConnAtoms(mMol.getBondAtom(0, bond)) == 3
                     && mMol.getConnAtoms(mMol.getBondAtom(1, bond)) == 3) {
                        if (ringMembership[bond] > 0) {
                            annelatedRing[ringMembership[bond] >>> 16]
                                         [ringMembership[bond] & 0x7FFF] = ring;
                            annelatedRing[ring][i] = (ringMembership[bond] >>> 16);
                            }
                        else {
                            ringMembership[bond] = (ring << 16) + 0x8000 + i;
                            }
                        }
                    }
                }
            }

        int ringsHandled = 0;
        int lastRingsHandled = -1;
        while (ringsHandled > lastRingsHandled) {
            lastRingsHandled = ringsHandled;
            for (int ring=0; ring<mRingAtomSet.size(); ring++) {
                if (!mAromaticityHandled[ring]) {
                    if (determineAromaticity(ring, annelatedRing)) {
                        mAromaticityHandled[ring] = true;
                        ringsHandled++;
                        }
                    }
                }
            }
        }

    private boolean determineAromaticity(int ringNo, int[][] annelatedRing) {
            // returns true if it can successfully determine and set the ring's aromaticity
        int ringAtom[] = mRingAtomSet.get(ringNo);
        int ringBond[] = mRingBondSet.get(ringNo);
        int ringBonds = ringBond.length;
        int bondSequence = 0;
        int aromaticButNotDelocalizedSequence = 0;
        boolean unhandledAnnelatedRingFound = false;
        for (int i=0; i<ringBonds; i++) {
            bondSequence <<= 1;
            aromaticButNotDelocalizedSequence <<= 1;
            if (mMol.getBondOrder(ringBond[i]) > 1
             || mMol.getBondType(ringBond[i]) == Molecule.cBondTypeDelocalized)
                bondSequence |= 1;
            else {
                int annelated = annelatedRing[ringNo][i];
                if (annelated != -1) {
                    if (mAromaticityHandled[annelated]) {
                        if (mIsAromatic[annelated]) {
                            bondSequence |= 1;
                            if (!mIsDelocalized[annelated])
                                aromaticButNotDelocalizedSequence |= 1;
                            }
                        }
                    else {
                        unhandledAnnelatedRingFound = true;
                        }
                    }
                }
            }

        boolean hasDelocalizationLeak = false;
        switch (ringBonds) {
        case 5:
            final int[] cSequence5Ring = {
               10,    // 01010
                5,    // 00101
               18,    // 10010
                9,    // 01001
               20 };  // 01010
            hasDelocalizationLeak = true;
            for (int heteroPosition=0; heteroPosition<5; heteroPosition++) {
                if ((bondSequence & cSequence5Ring[heteroPosition]) == cSequence5Ring[heteroPosition]) {
                    switch (mMol.getAtomicNo(ringAtom[heteroPosition])) {
                    case 6:
                        if (mMol.getAtomCharge(ringAtom[heteroPosition]) == -1) {
                            mIsAromatic[ringNo] = true;
                            mHeteroPosition[ringNo] = heteroPosition;
                            if ((aromaticButNotDelocalizedSequence & cSequence5Ring[heteroPosition]) == 0)
                                hasDelocalizationLeak = false;
                            }
                        break;
                    case 7:
                        if (mMol.getAtomCharge(ringAtom[heteroPosition]) <= 0) {
                            mIsAromatic[ringNo] = true;
                            mHeteroPosition[ringNo] = heteroPosition;
                            }
                        break;
                    case 8:
                        mIsAromatic[ringNo] = true;
                        mHeteroPosition[ringNo] = heteroPosition;
                        break;
                    case 16:
                        if (mMol.getConnAtoms(ringAtom[heteroPosition]) == 2) {
                            mIsAromatic[ringNo] = true;
                            mHeteroPosition[ringNo] = heteroPosition;
                            }
                        break;
                        }
                    }
                }
            break;
        case 6:
            hasDelocalizationLeak = true;
            if ((bondSequence & 21) == 21) {   // 010101
                mIsAromatic[ringNo] = true;
                if ((aromaticButNotDelocalizedSequence & 21) == 0)
                    hasDelocalizationLeak = false;
                }
            if ((bondSequence & 42) == 42) {   // 101010
                mIsAromatic[ringNo] = true;
                if ((aromaticButNotDelocalizedSequence & 42) == 0)
                    hasDelocalizationLeak = false;
                }
            break;
        case 7:
            final int[] cSequence7Ring = {
                   42,    // 0101010
                   21,    // 0010101
                   74,    // 1001010
                   37,    // 0100101
                   82,    // 1010010
                   41,    // 0101001
                   84 };  // 1010100
            hasDelocalizationLeak = true;
            for (int carbeniumPosition=0; carbeniumPosition<7; carbeniumPosition++) {
                if ((bondSequence & cSequence7Ring[carbeniumPosition]) == cSequence7Ring[carbeniumPosition]) {
                    if (mMol.getAtomicNo(ringAtom[carbeniumPosition]) == 6
                     && mMol.getAtomCharge(ringAtom[carbeniumPosition]) == 1) {
                        mIsAromatic[ringNo] = true;
                        mHeteroPosition[ringNo] = carbeniumPosition;
                        if ((aromaticButNotDelocalizedSequence & cSequence7Ring[carbeniumPosition]) == 0)
                            hasDelocalizationLeak = false;
                        }
                    }
                }
            break;
            }
        
        if (mIsAromatic[ringNo] && !hasDelocalizationLeak)
            mIsDelocalized[ringNo] = true;

        if (mIsAromatic[ringNo])
            return true;

        return !unhandledAnnelatedRingFound;
        }
    }