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

import java.util.ArrayList;

public class RingCollection {
	public static final int MAX_SMALL_RING_SIZE = 7;
	public static final int MAX_SMALL_RING_COUNT = 1024; // to prevent explosions with highly connected metal grids, etc.

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
	public static final int MODE_INCLUDE_TAUTOMERIC_BONDS = 8;

	private static final int FEATURES_RING_SIZE = 0x0000FFFF;
	private static final int FEATURES_AROMATIC = 0x00010000;
	private static final int FEATURES_DELOCALIZED = 0x00020000;
	private static final int FEATURES_HETERO_AROMATIC = 0x00040000;

	private ExtendedMolecule mMol;
	private ArrayList<int[]> mRingAtomSet;
	private ArrayList<int[]> mRingBondSet;
	private int[] mAtomRingFeatures;
	private int[] mBondRingFeatures;
	private int[] mHeteroPosition;
	private boolean[] mIsAromatic;
	private boolean[] mIsDelocalized;
	private int mMaxSmallRingSize;

	/**
	 * Generates the complete set of small rings, which don't contain metal atoms
	 * and have up to 7 members.<br> If mode includes LARGE_RINGS, then it determines
	 * for every atom and bond the size of the smallest ring, which they are
	 * a member of.<br>If mode includes AROMATICITY then every small ring
	 * is checked, whether it is aromatic.
	 * @param mol
	 * @param mode one of the public MODE_ options
	 */
	public RingCollection(ExtendedMolecule mol, int mode) {
		this(mol, mode, MAX_SMALL_RING_SIZE);
		}

	/**
	 * Generates the complete set of small rings, which don't contain metal atoms
	 * and have up to 7 members.<br> If mode includes LARGE_RINGS, then it determines
	 * for every atom and bond the size of the smallest ring, which they are
	 * a member of.<br>If mode includes AROMATICITY then every small ring
	 * is checked, whether it is aromatic.
	 * @param mol
	 * @param mode one of the public MODE_ options
	 * @param maxSmallRingSize largest ring size considered a small ring
	 */
	public RingCollection(ExtendedMolecule mol, int mode, int maxSmallRingSize) {
		mMol = mol;
		mMaxSmallRingSize = maxSmallRingSize;
		mRingAtomSet = new ArrayList<>();
		mRingBondSet = new ArrayList<>();

		mAtomRingFeatures = new int[mMol.getAtoms()];
		mBondRingFeatures = new int[mMol.getBonds()];

		mMol.ensureHelperArrays(ExtendedMolecule.cHelperNeighbours);

		boolean[] isConfirmedChainAtom = new boolean[mMol.getAtoms()];
		boolean[] isConfirmedChainBond = new boolean[mMol.getBonds()];

		boolean found;
		do {	// detect atoms of side chains as non-ring-atoms
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
			return;		 // no rings found

				// find all rings with less than 8 members of all closure bonds
		int graphAtom[] = new int[mMol.getAtoms()];
		graphAtom[0] = startAtom;
		int[] parent = new int[mMol.getAtoms()];
		parent[0] = -1;
		int fragmentNo[] = new int[mMol.getAtoms()];
		fragmentNo[startAtom] = 1;
		int current = 0;
		int highest = 0;
		int noOfFragments = 1;
		while (current <= highest) {
			for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
				int candidate = mMol.getConnAtom(graphAtom[current], i);
				if (candidate == parent[graphAtom[current]])
					continue;

				if (fragmentNo[candidate] != 0) {   // closure bond
					addSmallRingsToSet(mMol.getConnBond(graphAtom[current], i), isConfirmedChainAtom);
					continue;
					}

				if (!isConfirmedChainAtom[candidate]) {
					fragmentNo[candidate] = noOfFragments;
					parent[candidate] = graphAtom[current];
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
						parent[atom] = -1;
						break;
						}
					}
				}
			}

		if ((mode & MODE_AROMATICITY) != 0) {
			mIsAromatic = new boolean[mRingAtomSet.size()];
			mIsDelocalized = new boolean[mRingAtomSet.size()];
			mHeteroPosition = new int[mRingAtomSet.size()];
			determineAromaticity(mIsAromatic, mIsDelocalized, mHeteroPosition, (mode & MODE_INCLUDE_TAUTOMERIC_BONDS) != 0);
			updateAromaticity();
			}

		// find large rings by examining every potential ring bond
		// which is not a member of a small ring
		if ((mode & MODE_LARGE_RINGS) != 0) {
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (!isConfirmedChainBond[bond] && mMol.getBondOrder(bond) != 0) {
					int ringAtom[] = findSmallestRing(bond, isConfirmedChainAtom);
					if (ringAtom != null)
						updateRingSize(ringAtom, getRingBonds(ringAtom));
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
			for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
				int candidate = mMol.getConnAtom(graphAtom[current], i);
				if ((current > 1) && candidate == atom1) {
					int ringAtom[] = new int[graphLevel[graphAtom[current]]];
					int atom = graphAtom[current];
					for (int j = 0; j < ringAtom.length; j++) {
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
	 * of rings larger than 7 members and if the mode parameter of the constructor
	 * didn't include LARGE_RINGS, then the returned ring size is also 0.
	 * @param atom
	 * @return ring size or 0
	 */
	public int getAtomRingSize(int atom) {
		return mAtomRingFeatures[atom] & FEATURES_RING_SIZE;
		}


	/**
	 * A bond's ring size is the size of the smallest ring the bond is a member of.
	 * If the bond doesn't belong to any ring the ring size is 0. If a bond is member
	 * of rings larger than 7 members and if the mode parameter of the constructor
	 * didn't include LARGE_RINGS, then the returned ring size is also 0.
	 * @param bond
	 * @return ring size or 0
	 */
	public int getBondRingSize(int bond) {
		return mBondRingFeatures[bond] & FEATURES_RING_SIZE;
		}


	private void addSmallRingsToSet(int closureBond, boolean[] isConfirmedChainAtom) {
		int[] graphAtom = new int[mMaxSmallRingSize];
		int[] connIndex = new int[mMaxSmallRingSize];
		boolean[] isUsed = new boolean[mMol.getAtoms()];

		int atom1 = mMol.getBondAtom(0, closureBond);
		int atom2 = mMol.getBondAtom(1, closureBond);

		graphAtom[0] = atom1;
		graphAtom[1] = atom2;
		connIndex[1] = -1;
		isUsed[atom2] = true;
		int current = 1;

		while(current >= 1) {
			connIndex[current]++;
			if (connIndex[current] == mMol.getConnAtoms(graphAtom[current])) {
				isUsed[graphAtom[current]] = false;
				current--;
				continue;
				}

			int candidate = mMol.getConnAtom(graphAtom[current], connIndex[current]);
			if (isUsed[candidate] || isConfirmedChainAtom[candidate])
				continue;

			if (candidate == atom1 && current > 1) {
				addRingIfNew(graphAtom, current+1);

				// if we have already such many rings, we only collect the smallest ring to avoid a combinatorial explosion
				if (mRingAtomSet.size() >= MAX_SMALL_RING_COUNT)
					return;

				continue;
				}

			if (current+1 < mMaxSmallRingSize) {
				current++;
				graphAtom[current] = candidate;
				isUsed[candidate] = true;
				connIndex[current] = -1;
				}
			}
		}


	private void addRingIfNew(int ringAtom[], int ringSize) {
		int lowAtom = mMol.getMaxAtoms();
		int lowIndex = 0;
		for (int i=0; i<ringSize; i++) {
			if (lowAtom > ringAtom[i]) {
				lowAtom = ringAtom[i];
				lowIndex = i;
				}
			}

		int sortedRing[] = new int[ringSize];
		int leftIndex = (lowIndex > 0) ? lowIndex - 1 : ringSize - 1;
		int rightIndex = (lowIndex < ringSize - 1) ? lowIndex + 1 : 0;
		boolean inverse = (ringAtom[leftIndex] < ringAtom[rightIndex]);
		for (int i=0; i<ringSize; i++) {
			sortedRing[i] = ringAtom[lowIndex];
			if (inverse) {
				if (--lowIndex < 0)
					lowIndex = ringSize - 1;
				}
			else {
				if (++lowIndex == ringSize)
					lowIndex = 0;
				}
			}

		for (int i=0; i<mRingAtomSet.size(); i++) {
			int ringOfSet[] = mRingAtomSet.get(i);
			if (ringOfSet.length != ringSize)
				continue;
			boolean equal = true;
			for (int j=0; j<ringSize; j++) {
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

		updateRingSize(sortedRing, ringBond);
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
	 * a NullPointerException is raised.
	 * @param ringNo
	 * @return
	 */
	public boolean isAromatic(int ringNo) {
		return mIsAromatic[ringNo];
		}


	/**
	 * Return whether the atom is member if an aromatic ring.
	 * If the mode parameter passed to the constructor didn't include AROMATICITY, then
	 * a NullPointerException is raised.
	 * @param atom
	 * @return
	 */
	public boolean isAromaticAtom(int atom) {
		return (mAtomRingFeatures[atom] & FEATURES_AROMATIC) != 0;
		}


	/**
	 * Return whether the atom is member if a delocalized ring.
	 * If the mode parameter passed to the constructor didn't include AROMATICITY, then
	 * a NullPointerException is raised.
	 * @param atom
	 * @return
	 */
	public boolean isDelocalizedAtom(int atom) {
		return (mAtomRingFeatures[atom] & FEATURES_DELOCALIZED) != 0;
		}


	/**
	 * Return whether the atom is member if a hetero-aromatic ring.
	 * If the mode parameter passed to the constructor didn't include AROMATICITY, then
	 * a NullPointerException is raised.
	 * @param atom
	 * @return
	 */
	public boolean isHeteroAromaticAtom(int atom) {
		return (mAtomRingFeatures[atom] & FEATURES_HETERO_AROMATIC) != 0;
		}


	/**
	 * Return whether the bond is considered aromatic.
	 * If the mode parameter passed to the constructor didn't include AROMATICITY, then
	 * a NullPointerException is raised.
	 * @param bond
	 * @return
	 */
	public boolean isAromaticBond(int bond) {
		return (mBondRingFeatures[bond] & FEATURES_AROMATIC) != 0;
		}


	/**
	 * Return whether the bond is member if a delocalized ring.
	 * If the mode parameter passed to the constructor didn't include AROMATICITY, then
	 * a NullPointerException is raised.
	 * @param bond
	 * @return
	 */
	public boolean isDelocalizedBond(int bond) {
		return (mBondRingFeatures[bond] & FEATURES_DELOCALIZED) != 0;
		}


	/**
	 * Return whether the bond is member if a hetero-aromatic ring.
	 * If the mode parameter passed to the constructor didn't include AROMATICITY, then
	 * a NullPointerException is raised.
	 * @param bond
	 * @return
	 */
	public boolean isHeteroAromaticBond(int bond) {
		return (mBondRingFeatures[bond] & FEATURES_HETERO_AROMATIC) != 0;
		}


	/**
	 * Return whether the ring is considered delocalized, which are 6-membered
	 * aromatic rings with no preference concerning where the double bonds are located.
	 * Pyrrole bonds are not considered delocalized.
	 * If the mode parameter passed to the constructor didn't include AROMATICITY, then
	 * a NullPointerException is raised.
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


	/**
	 * brute force method to check, whether and which ring is shared by two bonds
	 * @param bond1
	 * @param bond2
	 * @return -1 if bond1 and bond2 don't share a common ring
	 */
	public int getSharedRing(int bond1, int bond2) {
		for (int i=0; i<mRingBondSet.size(); i++)
			if (isBondMember(i, bond1) && isBondMember(i, bond2))
				return i;
		return -1;
		}

	private void updateRingSize(int[] ringAtom, int[] ringBond) {
		int ringSize = ringAtom.length;
		for (int i=0; i<ringSize; i++) {
			int currentSize = mAtomRingFeatures[ringAtom[i]] & FEATURES_RING_SIZE;
			if (currentSize == 0 || currentSize>ringSize) {
				mAtomRingFeatures[ringAtom[i]] &= ~FEATURES_RING_SIZE;
				mAtomRingFeatures[ringAtom[i]] |= ringSize;
				}
			}

		for (int i=0; i<ringSize; i++) {
			int currentSize = mBondRingFeatures[ringBond[i]] & FEATURES_RING_SIZE;
			if (currentSize == 0 || currentSize>ringSize) {
				mBondRingFeatures[ringBond[i]] &= ~FEATURES_RING_SIZE;
				mBondRingFeatures[ringBond[i]] |= ringSize;
				}
			}
		}

	private void updateAromaticity() {
		for (int ring=0; ring<mIsAromatic.length; ring++) {
			if (mIsAromatic[ring]) {
				boolean isHeteroAromatic = false;
				for (int atom:mRingAtomSet.get(ring)) {
					mAtomRingFeatures[atom] |= FEATURES_AROMATIC;
					if (qualifiesAsHeteroAtom(atom))
						isHeteroAromatic = true;
					}
				for (int bond:mRingBondSet.get(ring))
					mBondRingFeatures[bond] |= FEATURES_AROMATIC;

				if (mIsDelocalized[ring]) {
					for (int atom:mRingAtomSet.get(ring))
						mAtomRingFeatures[atom] |= FEATURES_DELOCALIZED;
					for (int bond:mRingBondSet.get(ring))
						mBondRingFeatures[bond] |= FEATURES_DELOCALIZED;
					}

				if (isHeteroAromatic) {
					for (int atom:mRingAtomSet.get(ring))
						mAtomRingFeatures[atom] |= FEATURES_HETERO_AROMATIC;
					for (int bond:mRingBondSet.get(ring))
						mBondRingFeatures[bond] |= FEATURES_HETERO_AROMATIC;
					}
				}
			}
		}

	private boolean qualifiesAsHeteroAtom(int atom) {
		if (mMol.isFragment()) {
			if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0)
				return false;
			int[] atomList = mMol.getAtomList(atom);
			if (atomList != null) {
				for (int atomicNo : atomList)
					if (!Molecule.isAtomicNoElectronegative(atomicNo))
						return false;
				return true;
				}
			}
		return Molecule.isAtomicNoElectronegative(mMol.getAtomicNo(atom));
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


	/**
	 *
	 * @param isAromatic empty array sizes as getSize()
	 * @param isDelocalized
	 * @param heteroPosition
	 * @param includeTautomericBonds whether to treat non-methylated amide/thio-amide bonds as pi-bonds
	 */
	public void determineAromaticity(boolean[] isAromatic, boolean[] isDelocalized, int[] heteroPosition,
									  boolean includeTautomericBonds) {
		int[][] annelatedRing = new int[mRingAtomSet.size()][];
		for (int i=0; i<mRingAtomSet.size(); i++) {
			annelatedRing[i] = new int[mRingAtomSet.get(i).length];
			for (int j=0; j<mRingAtomSet.get(i).length; j++)
				annelatedRing[i][j] = -1;
			}

		int[] ringMembership = new int[mMol.getBonds()];
		for (int ring=0; ring<mRingBondSet.size(); ring++) {
			int[] ringBond = mRingBondSet.get(ring);
			if (ringBond.length == 3 || (ringBond.length >= 5 && ringBond.length <= 7)) {
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

		boolean[] aromaticityHandled = new boolean[mRingAtomSet.size()];
		int ringsHandled = 0;
		int lastRingsHandled = -1;
		while (ringsHandled > lastRingsHandled) {
			lastRingsHandled = ringsHandled;
			for (int ring=0; ring<mRingAtomSet.size(); ring++) {
				if (!aromaticityHandled[ring]) {
					if (determineAromaticity(ring, annelatedRing, aromaticityHandled,
							isAromatic, isDelocalized, heteroPosition, includeTautomericBonds)) {
						aromaticityHandled[ring] = true;
						ringsHandled++;
						}
					}
				}
			}
		}

	private boolean determineAromaticity(int ringNo, int[][] annelatedRing, boolean[] aromaticityHandled,
										 boolean []isAromatic, boolean[] isDelocalized, int[] heteroPosition,
										 boolean includeTautomericBonds) {
			// returns true if it can successfully determine and set the ring's aromaticity
		int ringAtom[] = mRingAtomSet.get(ringNo);
		for (int atom:ringAtom)
			if (!qualifiesAsAromaticAtom(atom))
				return true;

		int ringBond[] = mRingBondSet.get(ringNo);
		int ringBonds = ringBond.length;

		int bondSequence = 0;
		int aromaticButNotDelocalizedSequence = 0;
		boolean unhandledAnnelatedRingFound = false;
		for (int i=0; i<ringBonds; i++) {
			bondSequence <<= 1;
			aromaticButNotDelocalizedSequence <<= 1;
			if (qualifiesAsPiBond(ringBond[i])) {
				bondSequence |= 1;
				}
			else if (includeTautomericBonds && qualifiesAsAmideTypeBond(ringBond[i])) {
				bondSequence |= 1;
				aromaticButNotDelocalizedSequence |= 1;
				}
			else {
				int annelated = annelatedRing[ringNo][i];
				if (annelated != -1) {
					if (aromaticityHandled[annelated]) {
						if (isAromatic[annelated]) {
							bondSequence |= 1;
							if (!isDelocalized[annelated])
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
		case 3:
			final int[] cSequence3Ring = {
				2,	 // 010
				1,	 // 001
				4 }; // 100
			hasDelocalizationLeak = true;
			for (int carbeniumPosition=0; carbeniumPosition<3; carbeniumPosition++) {
				if ((bondSequence & cSequence3Ring[carbeniumPosition]) == cSequence3Ring[carbeniumPosition]) {
					if ((mMol.getAtomicNo(ringAtom[carbeniumPosition]) == 6
							&& mMol.getAtomCharge(ringAtom[carbeniumPosition]) == 1)
							|| (mMol.getAtomicNo(ringAtom[carbeniumPosition]) == 5
							&& mMol.getAtomCharge(ringAtom[carbeniumPosition]) == 0)) {
						isAromatic[ringNo] = true;
						heteroPosition[ringNo] = carbeniumPosition;
						if ((aromaticButNotDelocalizedSequence & cSequence3Ring[carbeniumPosition]) == 0)
							hasDelocalizationLeak = false;
						}
					}
				}
			break;
		case 5:
			final int[] cSequence5Ring = {
			   10,	// 01010
				5,	// 00101
			   18,	// 10010
				9,	// 01001
			   20 };// 01010
			hasDelocalizationLeak = true;
			for (int position=0; position<5; position++) {
				if ((bondSequence & cSequence5Ring[position]) == cSequence5Ring[position]) {
					switch (mMol.getAtomicNo(ringAtom[position])) {
					case 6:
						if (mMol.getAtomCharge(ringAtom[position]) == -1) {
							isAromatic[ringNo] = true;
							heteroPosition[ringNo] = position;
							if ((aromaticButNotDelocalizedSequence & cSequence5Ring[position]) == 0)
								hasDelocalizationLeak = false;
							}
						break;
					case 7:
						if (mMol.getAtomCharge(ringAtom[position]) <= 0) {
							isAromatic[ringNo] = true;
							heteroPosition[ringNo] = position;
							}
						break;
					case 8:
						isAromatic[ringNo] = true;
						heteroPosition[ringNo] = position;
						break;
					case 16:
					case 34:
					case 52:
						if (mMol.getConnAtoms(ringAtom[position]) == 2) {
							isAromatic[ringNo] = true;
							heteroPosition[ringNo] = position;
							}
						break;
						}
					}
				}
			break;
		case 6:
			hasDelocalizationLeak = true;
			if ((bondSequence & 21) == 21) {   // 010101
				isAromatic[ringNo] = true;
				if ((aromaticButNotDelocalizedSequence & 21) == 0)
					hasDelocalizationLeak = false;
				}
			if ((bondSequence & 42) == 42) {   // 101010
				isAromatic[ringNo] = true;
				if ((aromaticButNotDelocalizedSequence & 42) == 0)
					hasDelocalizationLeak = false;
				}
			break;
		case 7:
			final int[] cSequence7Ring = {
				   42,	// 0101010
				   21,	// 0010101
				   74,	// 1001010
				   37,	// 0100101
				   82,	// 1010010
				   41,	// 0101001
				   84 };// 1010100
			hasDelocalizationLeak = true;
			for (int carbeniumPosition=0; carbeniumPosition<7; carbeniumPosition++) {
				if ((bondSequence & cSequence7Ring[carbeniumPosition]) == cSequence7Ring[carbeniumPosition]) {
					if ((mMol.getAtomicNo(ringAtom[carbeniumPosition]) == 6
					  && (mMol.getAtomCharge(ringAtom[carbeniumPosition]) == 1
					|| (includeTautomericBonds && hasOxo(ringAtom[carbeniumPosition]))))
					 || (mMol.getAtomicNo(ringAtom[carbeniumPosition]) == 5
					  && mMol.getAtomCharge(ringAtom[carbeniumPosition]) == 0)) {
						isAromatic[ringNo] = true;
						heteroPosition[ringNo] = carbeniumPosition;
						if ((aromaticButNotDelocalizedSequence & cSequence7Ring[carbeniumPosition]) == 0)
							hasDelocalizationLeak = false;
						}
					}
				}
			break;
			}
		
		if (isAromatic[ringNo] && !hasDelocalizationLeak)
			isDelocalized[ringNo] = true;

		if (isAromatic[ringNo])
			return true;

		return !unhandledAnnelatedRingFound;
		}

	private boolean hasOxo(int carbonAtom) {
		for (int i=0; i<mMol.getConnAtoms(carbonAtom); i++)
			if (mMol.getConnBondOrder(carbonAtom, i) == 2
			 && mMol.getAtomicNo(mMol.getConnAtom(carbonAtom, i)) == 8)
				return true;
		return false;
		}

	private boolean qualifiesAsPiBond(int bond) {
		return (mMol.getBondOrder(bond) > 1
			 || mMol.getBondType(bond) == Molecule.cBondTypeDelocalized);
		}

	/**
	 * Checks, whether this bond may contribute pi-electrons from an amide-resonance
	 * to an aromatic ring. According to M J Cook, A R Katritzky, P Linda, R D Tack
	 * J. Chem. Soc., Perkin Trans. 2, 1972, 1295-1301
	 * 2-pyridone and 2-pyridinethione retain most of the aromatic resonance
	 * energy of pyridine unless the nitrogen atom is methylated.
	 * @param bond
	 * @return
	 */
	public boolean qualifiesAsAmideTypeBond(int bond) {
		for (int i=0; i<2; i++) {
			int atom1 = mMol.getBondAtom(i, bond);
			if ((mMol.getAtomicNo(atom1) == 7)
			 && mMol.getConnAtoms(atom1) == 2) {
				int atom2 = mMol.getBondAtom(1-i, bond);
				if (mMol.getAtomicNo(atom2) == 6) {
					for (int j=0; j<mMol.getConnAtoms(atom2); j++) {
						int connAtom = mMol.getConnAtom(atom2, j);
						int connBond = mMol.getConnBond(atom2, j);
						if ((mMol.getAtomicNo(connAtom) == 8 || mMol.getAtomicNo(connAtom) == 16)
						 && mMol.getBondOrder(connBond) == 2
						 && mMol.getConnAtoms(connAtom) == 1)
						return true;
						}
					}
				}
			}

		return false;
		}

	private boolean qualifiesAsAromaticAtom(int atom) {
		// If we have a list or wildcard atom, then the atomicNo is meaningless...
		if (mMol.isFragment()) {
			// We consider wildcard atoms as being compatible with aromaticity
			if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0) {
				// in theory, we must return false, if all atomicNos, which qualify for aromaticity,
				// are part of the exclude list. In reality that is rather unlikely...
				return true;
				}
			else {
				int[] list = mMol.getAtomList(atom);
				if (list != null) {
					for (int atomicNo:list)
						if (qualifiesAsAromaticAtomicNo(atomicNo))
							return true;

					return false;
					}
				}
			}

		return qualifiesAsAromaticAtomicNo(mMol.getAtomicNo(atom));
		}

	public static boolean qualifiesAsAromaticAtomicNo(int atomicNo) {
		return atomicNo == 5
			|| atomicNo == 6
			|| atomicNo == 7
			|| atomicNo == 8
			|| atomicNo == 15   // P
			|| atomicNo == 16   // S
			|| atomicNo == 33   // As
			|| atomicNo == 34;   // Se
		}
	}