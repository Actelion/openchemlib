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

package com.actelion.research.chem.reaction;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.Arrays;

public class MoleculeAutoMapper implements AutoMapper {
	private static final int MASK_ATOM_INDEX = 0x0000FFFF;
	private static final int MASK_ATOM_TYPE = 0xFFFF0000;

	private StereoMolecule		mMol;
    private Canonizer			mCanonizer;
	private int					mCurrentMapNo;
	private int[]				mCounterAtom;
	private boolean[]			mMapNoInUse,mMatchHandled;

	public MoleculeAutoMapper(StereoMolecule mol) {
		mMol = mol;
		}

	public void autoMap() {
		// Removes previously assigned auto-mapping numbers and recursively assigns
		// new mapping numbers to all neighbors of manually assigned seed atoms.
		// Expects existing mapping numbers to occur on two atoms sharing the same atomic number.

		mMol.ensureHelperArrays(Molecule.cHelperRings);
        // Init vars
		mMapNoInUse = new boolean[mMol.getAtoms()+1];
        mMatchHandled = new boolean[mMol.getAtoms()+1];
		mCounterAtom = new int[mMol.getAtoms()];
        mCurrentMapNo = 0;

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			int mapNo = mMol.getAtomMapNo(atom);
			if (mapNo == 0) {
				mCounterAtom[atom] = -1;
				}
			else {
				if (mMol.isAutoMappedAtom(atom)) {
					mMol.setAtomMapNo(atom, 0, false);
					mCounterAtom[atom] = -1;
					}
				else if (!mMapNoInUse[mapNo]) {
					mMapNoInUse[mapNo] = true;

					for (int counterAtom=atom+1; counterAtom<mMol.getAtoms(); counterAtom++) {
						if (mMol.getAtomMapNo(counterAtom) == mapNo) {
							mCounterAtom[atom] = counterAtom;
							mCounterAtom[counterAtom] = atom;
							break;
							}
						}
					if (mCounterAtom[atom] == -1) {
						mMol.setAtomMapNo(atom, 0, false);
						mCounterAtom[atom] = -1;
						}
					}
				}
			}


		matchFragments();

		mCanonizer = new Canonizer(mMol, Canonizer.CREATE_SYMMETRY_RANK);
		boolean found;
		do {
			found = false;
			for (int atom1=0; atom1<mMol.getAtoms(); atom1++) {
				int mapNo = mMol.getAtomMapNo(atom1);
				if (mapNo != 0 && !mMatchHandled[mapNo])
					if (autoMapNeighbors(atom1, mCounterAtom[atom1]))
						found = true;
				}
			} while (found);

		matchFragments();
		}

	private void matchFragments() {
		// locate disconnected unmapped areas as fragments
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			mMol.setAtomMarker(atom, mMol.getAtomMapNo(atom) == 0
					&& (mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) == 0);
		for (int atom=mMol.getAtoms(); atom<mMol.getAllAtoms(); atom++)
			mMol.setAtomMarker(atom, false);
		int[] fragmentNo = new int[mMol.getAllAtoms()];
		int fragmentCount = mMol.getFragmentNumbers(fragmentNo, true, false);	// TODO: check metal bond behaviour
		mMol.removeAtomMarkers();

		if (fragmentCount == 0)
			return;

		StereoMolecule[] fragment = mMol.getFragments(fragmentNo, fragmentCount);

		FragmentSpec[] spec = new FragmentSpec[fragmentCount];
		for (int i=0; i<fragmentCount; i++) {
			fragment[i].ensureHelperArrays(Molecule.cHelperRings);
			spec[i] = new FragmentSpec();
			spec[i].atomMap = new int[fragment[i].getAtoms()];
			spec[i].rank = new int[fragment[i].getAtoms()];
			}

		// create mapping from fragment's atom indices to mMol atom index
		int[] fragmentAtom = new int[fragmentCount];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getAtomMapNo(atom) == 0
			 && (mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) == 0) {
				int fNo = fragmentNo[atom];
				int fAtom = fragmentAtom[fNo]++;

				String label = encodeMappedNeighbors(atom);
				if (label != null) {
					fragment[fNo].setAtomCustomLabel(fAtom, label);

					int highNeighbor = getHighestMappedNeighbor(atom);
					if (spec[fNo].highestMappedNeighbor < highNeighbor)
						spec[fNo].highestMappedNeighbor = highNeighbor;
					}

				spec[fNo].atomMap[fAtom] = atom;
				}
			}

		for (int i=0; i<fragmentCount; i++) {
			Canonizer canonizer = new Canonizer(fragment[i], Canonizer.ENCODE_ATOM_CUSTOM_LABELS);
			spec[i].code = canonizer.getIDCode();
			for (int a=0; a<fragment[i].getAtoms(); a++)
				spec[i].rank[a] = canonizer.getFinalRank()[a];
			}

		for (int i=0; i<fragmentCount-1; i++) {
			if (!spec[i].isMapped) {
				for (int j=i+1; j<fragmentCount; j++) {
					if (!spec[j].isMapped
					 && spec[i].code.equals(spec[j].code)
					 && spec[i].highestMappedNeighbor != spec[j].highestMappedNeighbor) {
						mapFragments(spec[i], spec[j]);
						break;
						}
					}
				}
			}
		}

	private String encodeMappedNeighbors(int atom) {
		int count = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (mMol.getAtomMapNo(mMol.getConnAtom(atom, i)) != 0)
				count++;

		if (count == 0)
			return null;

		if (count == 1) {
			for (int i=0; i<mMol.getConnAtoms(atom); i++) {
				int mapNo = mMol.getAtomMapNo(mMol.getConnAtom(atom, i));
				if (mapNo != 0)
					return encodeBond(mMol.getConnBond(atom, i)) + mapNo;
				}
			}

		String[] code = new String[count];
		count = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int mapNo = mMol.getAtomMapNo(mMol.getConnAtom(atom, i));
			if (mapNo != 0)
				code[count++] = encodeBond(mMol.getConnBond(atom, i)) + mapNo;
			}
		Arrays.sort(code);
		StringBuilder list = new StringBuilder();
		for (int i=0; i<count; i++)
			list.append(code[i]);
		return list.toString();
		}

	private String encodeBond(int bond) {
		if (mMol.isDelocalizedBond(bond))
			return ".";
		int order = mMol.getBondOrder(bond);
		return order == 1 ? "-" : order == 2 ? "=" : "#";
		}

	private int getHighestMappedNeighbor(int atom) {
		int highest = -1;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int connAtom = mMol.getConnAtom(atom, i);
			if (mMol.getAtomMapNo(connAtom) != 0
			 && highest < connAtom)
				highest = connAtom;
			}
		return highest;
		}

	private void mapFragments(FragmentSpec spec1, FragmentSpec spec2) {
		for (int a1=0; a1<spec1.rank.length; a1++) {
			for (int a2=0; a2<spec2.rank.length; a2++) {
				if (spec1.rank[a1] == spec2.rank[a2]) {
					mapAtoms(spec1.atomMap[a1], spec2.atomMap[a2]);
					mMatchHandled[mCurrentMapNo] = true;
					break;
					}
				}
			}
		spec1.isMapped = true;
		spec2.isMapped = true;
		}

	private boolean autoMapNeighbors(int atom1, int atom2) {
		int atomMapNo = Math.abs(mMol.getAtomMapNo(atom1));

		int[] neighborInfo1 = getUnmappedNeighbors(atom1);
		if (neighborInfo1 == null) {
			mMatchHandled[atomMapNo] = true;
			return false;
			}

		int[] neighborInfo2 = getUnmappedNeighbors(atom2);
		if (neighborInfo2 == null) {
			mMatchHandled[atomMapNo] = true;
			return false;
			}

		int index1 = 0;
		int index2 = 0;
		int mappedNeighbors = 0;
		while (index1 < neighborInfo1.length
			&& index2 < neighborInfo2.length) {

			if ((neighborInfo1[index1] & MASK_ATOM_TYPE)
			 == (neighborInfo2[index2] & MASK_ATOM_TYPE)) {
				int count1 = countSimilarNeighbors(neighborInfo1, index1);
				int count2 = countSimilarNeighbors(neighborInfo2, index2);

				mappedNeighbors += tryMapNeighbors(atom1, neighborInfo1, index1, count1, atom2, neighborInfo2, index2, count2);

				index1 += count1;
				index2 += count2;
				continue;
				}

			if (index2 < neighborInfo2.length)
				while ((index1 < neighborInfo1.length)
					&& (neighborInfo1[index1] & MASK_ATOM_TYPE) < (neighborInfo2[index2] & MASK_ATOM_TYPE))
					index1++;
			if (index1 < neighborInfo1.length)
				while ((index2 < neighborInfo2.length)
					&& (neighborInfo2[index2] & MASK_ATOM_TYPE) < (neighborInfo1[index1] & MASK_ATOM_TYPE))
					index2++;
			}

		if (neighborInfo1.length == mappedNeighbors
		 || neighborInfo2.length == mappedNeighbors)
			mMatchHandled[atomMapNo] = true;

		return mappedNeighbors != 0;
		}

	private int countSimilarNeighbors(int[] neighborInfo, int startIndex) {
		int count = 1;
		while ((startIndex+count < neighborInfo.length)
			&& (neighborInfo[startIndex+count] & MASK_ATOM_TYPE) == (neighborInfo[startIndex] & MASK_ATOM_TYPE))
			count++;
		return count;
		}

	private int tryMapNeighbors(int atom1, int[] neighborInfo1, int index1, int count1,
								int atom2, int[] neighborInfo2, int index2, int count2) {
		if (count1 == 1 && count2 == 1) {
			mapAtoms(neighborInfo1[index1] & MASK_ATOM_INDEX, neighborInfo2[index2] & MASK_ATOM_INDEX);
			return 1;
			}

		boolean atom1HasEqualNeighbors = areNeighborsEqual(neighborInfo1, index1, count1);
		boolean atom2HasEqualNeighbors = areNeighborsEqual(neighborInfo2, index2, count2);

		// If at least one atom's neighbors are all symmetrical and the other atom has not more neighbors
		// than this one, then we can match any neighbor of the symmetrical neighbors.
		if ((atom1HasEqualNeighbors && count1 >= count2)
		 || (atom2HasEqualNeighbors && count2 >= count1)) {
			int count = Math.min(count1, count2);
			for (int i=0; i<count; i++)
				mapAtoms(neighborInfo1[index1+i] & MASK_ATOM_INDEX, neighborInfo2[index2+i] & MASK_ATOM_INDEX);

			return count;
			}

		// If we have two matching atoms with no further neighbors and matching bond orders
		for (int i=0; i<count1; i++) {
			int neighbor1 = neighborInfo1[index1+i] & MASK_ATOM_INDEX;
			if (mMol.getConnAtoms(neighbor1) == 1) {
				for (int j=0; j<count2; j++) {
					int neighbor2 = neighborInfo2[index2+j] & MASK_ATOM_INDEX;
					if (mMol.getConnAtoms(neighbor2) == 1) {
						if (mMol.getBondOrder(mMol.getBond(atom1, neighbor1))
						 == mMol.getBondOrder(mMol.getBond(atom2, neighbor2))) {
							mapAtoms(neighbor1, neighbor2);
							return 1;
							// rarely occuring multiple matches are handled in next round
							}
						}
					}
				}
			}

		return 0;
		}

	private boolean areNeighborsEqual(int[] neighborInfo, int startIndex, int count) {
		int symmetryRank = mCanonizer.getSymmetryRank(neighborInfo[startIndex] & MASK_ATOM_INDEX);
		for (int i=1; i<count; i++)
			if (symmetryRank != mCanonizer.getSymmetryRank(neighborInfo[startIndex+i] & MASK_ATOM_INDEX))
				return false;

		return true;
		}

	/**
	 * Creates a sorted list of specifications of all unmapped neighbors.
	 * bit 0-15 contains the neighbor's atom index.
	 * bit 20-27 contains the neighbor's atomic no
	 * @param atom
	 * @return
	 */
	private int[] getUnmappedNeighbors(int atom) {
		int count = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int connAtom = mMol.getConnAtom(atom, i);
			if (mMol.getAtomMapNo(connAtom) == 0
			 && (mMol.getAtomQueryFeatures(connAtom) & Molecule.cAtomQFExcludeGroup) == 0)
				count++;
			}

		if (count == 0)
			return null;

		int[] neighborInfo = new int[count];
		count = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int connAtom = mMol.getConnAtom(atom, i);
			if (mMol.getAtomMapNo(connAtom) == 0
			 && (mMol.getAtomQueryFeatures(connAtom) & Molecule.cAtomQFExcludeGroup) == 0)
				neighborInfo[count++] = (mMol.getAtomicNo(connAtom) << 24)
									  + (mMol.getAtomMass(connAtom) << 16)
									  + connAtom;
			}
		Arrays.sort(neighborInfo);
		return neighborInfo;
		}

	private void mapAtoms(int atom1, int atom2) {
		int newMapNo = getNextFreeMapNo();
		mMol.setAtomMapNo(atom1, newMapNo, true);
		mMol.setAtomMapNo(atom2, newMapNo, true);
		mCounterAtom[atom1] = atom2;
		mCounterAtom[atom2] = atom1;
		}

	private int getNextFreeMapNo() {
		do {
			mCurrentMapNo++;
			} while (mCurrentMapNo<mMapNoInUse.length && mMapNoInUse[mCurrentMapNo]);
		return (mCurrentMapNo<mMapNoInUse.length) ? mCurrentMapNo : -1;
		}
	}

//class SubstituentSpec {
//	public String code;
//	public int[] rank;
//	public int[] atomMap;
//	public boolean isMapped;
//	}

class FragmentSpec {
	public String code;
	public int[] rank;
	public int[] atomMap;
	public int highestMappedNeighbor;
	public boolean isMapped;
	}
