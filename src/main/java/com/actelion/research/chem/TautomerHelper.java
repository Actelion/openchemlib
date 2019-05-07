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

public class TautomerHelper {
	private StereoMolecule mMol;
	private int[] mChainAtom;
	private int[] mChainBond;
	private boolean[] mChainNeedsDoubleBond;
	private boolean[] mChainNeedsDonorAtom;
	private boolean[] mIsTautomerBond;
	private int[] mRegionPiCount;
	private int[] mRegionDCount;
	private int[] mRegionTCount;

	public TautomerHelper(StereoMolecule mol) {
		mMol = mol;
		mMol.ensureHelperArrays(Molecule.cHelperParities);
		mChainAtom = new int[mMol.getAtoms()];
		mChainBond = new int[mMol.getAtoms()];
		mChainNeedsDoubleBond = new boolean[mMol.getAtoms()];
		mChainNeedsDonorAtom = new boolean[mMol.getAtoms()];
		mIsTautomerBond = new boolean[mMol.getBonds()];
		}

	/**
	 * Identifies connected tautomeric regions and assign region numbers to all atoms.
	 * Atoms sharing the same region share the same number.<br>
	 * 0: not member of a tautomer region; 1 and above: region number
	 * @param atomRegionNo int[mol.getAtoms()] filled with 0
	 * @param keepStereoCenters
	 * @return region count
	 */
	public int getAtomRegionNumbers(int[] atomRegionNo, boolean keepStereoCenters) {
		if (!findTautomericBonds(keepStereoCenters))
			return 0;

		int regionCount = assignRegionNumbers(atomRegionNo);
		return regionCount;
		}

	/**
	 * If no tautomers can be formed then the original molecule is returned.
	 * Otherwise the original molecule is copied and normalized to create a
	 * generic tautomer structure without touching the original molecule.
	 * Different tautomers of the same molecule should always result in the same
	 * generic tautomer structure. A generic tautomer contains one or more regions
	 * indicated by bond query features cBondQFSingle & cBondQFDouble. Bond types
	 * of all bonds of any tautomer region are cBondTypeSingle.
	 * The highest ranking atom in every region carries a label defining the number
	 * of double bonds and D and T atoms. The returned molecule has the fragment bit set.
	 * Canonicalizing the returned molecule with Canonizer mode ENCODE_ATOM_CUSTOM_LABELS
	 * produces the same idcode from any tautomer.
	 * If keepStereoCenters is true, then stereo centers with parity 1 or 2, if they are
	 * absolute or if they are part of an AND/OR group with more than one member,
	 * are considered stable (non racemising) and, thus, their proton is not considered
	 * being able to take part in a tautomeric transition.
	 * @param keepStereoCenters if true, then defined (non-racemising) stereo centers cannot be part of tautomeric regions
	 * @return generic tautomer with normalized tautomer regions and custom label to encode pi,D,T counts
	 */
	public StereoMolecule createGenericTautomer(boolean keepStereoCenters) {
		if (!findTautomericBonds(keepStereoCenters))
			return mMol;

		StereoMolecule genericTautomer = mMol.getCompactCopy();
		genericTautomer.setFragment(true);
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mIsTautomerBond[bond]) {
				genericTautomer.setBondType(bond, Molecule.cBondTypeSingle);
				genericTautomer.setBondQueryFeature(bond, Molecule.cBondQFSingle | Molecule.cBondQFDouble, true);
			    }
		    }

		int[] atomRegionNo = new int[mMol.getAtoms()];
		int regionCount = assignRegionNumbers(atomRegionNo);

		// find highest ranking atom in every region
		int[] maxAtom = new int[regionCount];
		int[] maxRank = new int[regionCount];
		int[] atomRank = new Canonizer(genericTautomer).getFinalRank();
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (atomRegionNo[atom] != 0) {
				int regionIndex = atomRegionNo[atom]-1;
				if (maxRank[regionIndex] < atomRank[atom]) {
					maxRank[regionIndex] = atomRank[atom];
					maxAtom[regionIndex] = atom;
					}
				}
			}

		// attach label with region counts to highest ranking atoms
		compileRegionCounts(atomRegionNo, regionCount);
		for (int i=0; i<regionCount; i++) {
			String label = ""+mRegionPiCount[i]+"|"+mRegionDCount[i]+"|"+mRegionTCount[i];
			genericTautomer.setAtomCustomLabel(maxAtom[i], label);
			}

		return genericTautomer;
		}

	/**
	 * Considers connected tautomer bonds to belong to a tautomer region.
	 * All independent tautomer regions are located and member atoms assigned to them.
	 * mAtomRegionNo[] is set accordingly.
	 * 0: not member of a tautomer region; 1 and above: region number
	 * @param atomRegionNo int[mol.getAtoms()] filled with 0
	 * @return number of found tautomer regions
	 */
	private int assignRegionNumbers(int[] atomRegionNo) {
		int[] graphAtom = new int[mMol.getAtoms()];
		boolean[] bondWasSeen = new boolean[mMol.getBonds()];
		int region = 0;

		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (!bondWasSeen[bond] && mIsTautomerBond[bond]) {
				region++;
				atomRegionNo[mMol.getBondAtom(0, bond)] = region;
				atomRegionNo[mMol.getBondAtom(1, bond)] = region;
				bondWasSeen[bond] = true;
				for (int i=0; i<2; i++) {
					int atom = mMol.getBondAtom(i, bond);
					atomRegionNo[atom] = region;
					int current = 0;
					int highest = 0;
					graphAtom[0] = atom;
					while (current <= highest) {
						for (int j=0; j<mMol.getConnAtoms(graphAtom[current]); j++) {
							int connBond = mMol.getConnBond(graphAtom[current], j);
							if (!bondWasSeen[connBond] && mIsTautomerBond[connBond]) {
								bondWasSeen[connBond] = true;
								int connAtom = mMol.getConnAtom(graphAtom[current], j);
								if (atomRegionNo[connAtom] == 0) {
									atomRegionNo[connAtom] = region;
									graphAtom[++highest] = connAtom;
									}
								}
							}
						current++;
						}
					}
				}
			}

		return region;
		}

	/**
	 * Counts for every region: pi-electrons, deuterium atoms, tritium atoms.
	 * Must be called after assignRegionNumbers().
	 * @param atomRegionNo array with valid region numbers
	 */
	private void compileRegionCounts(int[] atomRegionNo, int regionCount) {
		mRegionPiCount = new int[regionCount];
		mRegionDCount = new int[regionCount];
		mRegionTCount = new int[regionCount];
		
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (atomRegionNo[atom] != 0) {
				int regionIndex = atomRegionNo[atom]-1;
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					int connAtom = mMol.getConnAtom(atom, i);
					if (mMol.getAtomicNo(connAtom) == 1) {
						if (mMol.getAtomMass(connAtom) == 2)
							mRegionDCount[regionIndex]++;
						if (mMol.getAtomMass(connAtom) == 3)
							mRegionTCount[regionIndex]++;
						}
					}
				}
			}
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mIsTautomerBond[bond]
			 && mMol.getBondOrder(bond) == 2) {
				mRegionPiCount[atomRegionNo[mMol.getBondAtom(0, bond)]-1] += 2;
				}
			}
		}

	/**
	 * Locates and flags all 
	 * @return true if any tautomeric bonds could be found
	 */
	private boolean findTautomericBonds(boolean keepStereoCenters) {
		boolean[] isProtectedAtom = locateProtectedAtoms(keepStereoCenters);
		mIsTautomerBond = new boolean[mMol.getBonds()];

		// Locate pi-donor fragments as central atom with 1 double bond neighbor
		// and 1 or 2 single bond neighbor(s) with at least one hydrogen attached.
		boolean tautomerFound = false;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (!isProtectedAtom[atom]
			 && isValidHeteroAtom(atom)
			 && (mMol.getAtomPi(atom) != 0 || mMol.getAllHydrogens(atom) != 0)) {
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					int connAtom = mMol.getConnAtom(atom, i);
					int connBondOrder = mMol.getConnBondOrder(atom, i);
					if ((mMol.getAtomPi(atom) == 0 || connBondOrder == 2)
					 && mMol.getAtomPi(connAtom) == 1
					 && !isProtectedAtom[connAtom]) {
						if (findTautomerTreeDepthFirst(atom, connAtom, mMol.getConnBond(atom, i), isProtectedAtom)) {
							tautomerFound = true;
							if (connBondOrder == 2) {
								// in case we have a double bond from the hetero atom, the we need to process other single bonded neighbours
								// of the hereo atom separately in case, the first tree search finds a (vinylog) donor atom.
								for (int j=0; j<mMol.getConnAtoms(atom); j++) {
									if (j != i) {
										int otherConn = mMol.getConnAtom(atom, j);
										int otherBondOrder = mMol.getConnBondOrder(atom, i);
										if (!isProtectedAtom[otherConn]
										 && otherBondOrder == 1
										 && mMol.getAtomPi(otherConn) == 1)
											findTautomerTreeDepthFirst(atom, otherConn, mMol.getConnBond(atom, j), isProtectedAtom);
										}
									}
								}
							}
						}
					}
				}
			}

		return tautomerFound;
		}

	private boolean findTautomerTreeDepthFirst(int heteroAtom, int atom1, int bond1, boolean[] isProtectedAtom) {
		boolean tautomerBondsFound = false;

		int[] connIndex = new int[mMol.getAtoms()];

		boolean[] isChainAtom = new boolean[mMol.getAtoms()];
		isChainAtom[heteroAtom] = true;
		isChainAtom[atom1] = true;

		mChainAtom[0] = atom1;
		mChainBond[0] = bond1;
		mChainNeedsDoubleBond[0] = (mMol.getBondOrder(bond1) == 1);
		mChainNeedsDonorAtom[0] = (mMol.getBondOrder(bond1) == 2);

		int connAtom = -1;
		int connBond = -1;

		int current = 0;
		int highestValid = -1;
		boolean forward = true;
		while (current >= 0) {
			int currentAtom = mChainAtom[current];

			if (forward)
				connIndex[current] = mMol.getConnAtoms(currentAtom);
			connIndex[current]--;

			int desiredBondOrder = mChainNeedsDoubleBond[current] ? 2 : 1;
			while (connIndex[current] >= 0) {
				connAtom = mMol.getConnAtom(currentAtom, connIndex[current]);
				connBond = mMol.getConnBond(currentAtom, connIndex[current]);
				if (!isProtectedAtom[connAtom]
				 && !isChainAtom[connAtom]) {	// don't allow ring closures
					if ((mMol.isDelocalizedBond(connBond) || mMol.getBondOrder(connBond) == desiredBondOrder)
					 && (mMol.getAtomPi(connAtom) != 0
					  || (mChainNeedsDonorAtom[current]) && isValidDonorAtom(connAtom))) {
						// this is a continuation of the delocalized chain or may be an endo-donor as O=C-CHR-C=C
						break;
						}
					}
				connIndex[current]--;
				}
	
			if (connIndex[current] == -1) {
				if (forward && (!mChainNeedsDonorAtom[0] || !mChainNeedsDonorAtom[current])) {
					for (int i=current; i>highestValid; i--) {
						mIsTautomerBond[mChainBond[i]] = true;

						// find and flag exo-donor atoms as CH in NH-C(-CHR2)=C
						if (!mChainNeedsDonorAtom[0] && mChainNeedsDoubleBond[i]) {
							for (int j=0; j<mMol.getConnAtoms(mChainAtom[i]); j++) {
								int donorAtom = mMol.getConnAtom(mChainAtom[i], j);
								if (!isProtectedAtom[donorAtom]
								 && !isChainAtom[donorAtom]
								 && isValidDonorAtom(donorAtom))
									mIsTautomerBond[mMol.getConnBond(mChainAtom[i], j)] = true;
								}
							}
						}
					tautomerBondsFound = true;
					highestValid = current;
					}

				isChainAtom[mChainAtom[current]] = false;
				current--;
				forward = false;
				if (highestValid > current)
					highestValid = current;
				continue;
				}		

			current++;
			forward = true;
			mChainAtom[current] = connAtom;
			mChainBond[current] = connBond;
			mChainNeedsDoubleBond[current] = !mChainNeedsDoubleBond[current-1] && (mMol.getAtomPi(connAtom) != 0);
			mChainNeedsDonorAtom[current] = mChainNeedsDonorAtom[current-1] && (mMol.getAtomPi(connAtom) != 0);
			isChainAtom[connAtom] = true;
			}

		return tautomerBondsFound;
		}

	private boolean[] locateProtectedAtoms(boolean keepStereoCenters) {
		boolean[] isProtectedAtom = new boolean[mMol.getAtoms()];

		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getConnAtoms(atom) > 3)
				isProtectedAtom[atom] = true;

		// protect stereo centers with defined parity
		if (keepStereoCenters)
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (isNonRacemisingStereoCenter(atom))
					isProtectedAtom[atom] = true;

		// protected delocalized 6-membered-ring without hetero-atoms in ring our first shell,
		// that have no ortho or para substitution
		RingCollection ringSet = mMol.getRingSet();
		for (int r=0; r<ringSet.getSize(); r++) {
			if (ringSet.getRingSize(r) == 6
			 && ringSet.isDelocalized(r)) {
				boolean heteroFound = false;
				int oddAndEvenSubstituentMask = 0;
				int[] ringAtom = ringSet.getRingAtoms(r);
				for (int i=0; i<ringAtom.length && !heteroFound; i++) {
					int atom = ringAtom[i];

					int connAtoms = mMol.getConnAtoms(atom);
					if (connAtoms > 2)
						oddAndEvenSubstituentMask |= ((i & 1) != 0) ? 1 : 2;

					for (int j=0; j<connAtoms; j++) {
						if (mMol.getAtomicNo(mMol.getConnAtom(atom, j)) != 6) {
							heteroFound = true;
							break;
							}
						}
					}

				if (!heteroFound && oddAndEvenSubstituentMask != 3)
					for (int atom:ringAtom)
						isProtectedAtom[atom] = true;
				}
			}

		return isProtectedAtom;
		}

	/**
	 * Determines whether this atom is a non-racemizing stereo center,
	 * i.e. it has an assigned parity and it is either absolute or has
	 * a given relative parity towards other atom(s).
	 * @param atom
	 * @return
	 */
	private boolean isNonRacemisingStereoCenter(int atom) {
		if (mMol.getAtomParity(atom) != Molecule.cAtomParity1
		 && mMol.getAtomParity(atom) != Molecule.cAtomParity2)
			return false;

		if (mMol.getAtomESRType(atom) == Molecule.cESRTypeAbs)
			return true;

		int type = mMol.getAtomESRType(atom);
		int group = mMol.getAtomESRGroup(atom);
		for (int i=0; i<mMol.getAtoms(); i++)
			if (mMol.getAtomESRType(i) == type
			 && mMol.getAtomESRGroup(i) == group)
				return true;

		return false;
		}

	private boolean isValidHeteroAtom(int atom) {
		return mMol.getAtomicNo(atom) == 7
			|| mMol.getAtomicNo(atom) == 8
			|| mMol.getAtomicNo(atom) == 16
			|| mMol.getAtomicNo(atom) == 34
			|| mMol.getAtomicNo(atom) == 52;
		}

	private boolean isValidDonorAtom(int atom) {
		return mMol.getAtomPi(atom) == 0
			&& mMol.getAllHydrogens(atom) != 0
			&& (mMol.getAtomicNo(atom) == 6 || isValidHeteroAtom(atom));
		}
	}
