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

public class AromaticityResolver {
	ExtendedMolecule	mMol;
	private boolean		mAllHydrogensAreExplicit;
	private boolean[]	mIsDelocalizedRing,mIsDelocalizedAtom,mIsDelocalizedBond,mIsDelocalizedBridgeHead,
						mIsDelocalizedFiveRingMember,mIsDelocalizedThreeOrSevenRingMember;
    private int mDelocalizedAtoms, mDelocalizedBonds,mPiElectronsAdded;

    /**
     * Creates a new AromaticityResolver for molecule mol.
     * @param mol
     */
    public AromaticityResolver(ExtendedMolecule mol) {
        mMol = mol;
        }

    /**
     * This method promotes all necessary bonds of the defined delocalized part of the molecule
     * from single to double bonds in order to create a valid delocalized system
     * of conjugated single and double bonds.
	 * The delocalized part of the molecule may be defined by passing an array
	 * to isAromaticBond that has all bonds flagged, which are part of a delocalized area.
	 * In this case these bonds are assumed to have bond type cBondTypeSingle.
	 * Alternatively, one may pass null and indicate affected bonds with bond type cBondTypeDelocalized.
     * Non-cyclic atom chains defined to be delocalized are treated depending
     * on whether we have a molecule or a query fragment. For fragments the respective bond
     * types will be set to cBondTypeDelocalized; for molecules the chain will
     * have alternating single and double bonds starting with double at a non-ring end.
     * @return true if all bonds of the delocalized area could be consistently converted. 
     */
	public boolean locateDelocalizedDoubleBonds(boolean[] isAromaticBond) {
		return locateDelocalizedDoubleBonds(isAromaticBond, false, false);
		}

	/**
	 * This method promotes all necessary bonds of the defined delocalized part of the molecule
	 * from single to double bonds in order to create a valid delocalized system
	 * of conjugated single and double bonds.
	 * The delocalized part of the molecule may be defined by passing an array
	 * to isAromaticBond that has all bonds flagged, which are part of a delocalized area.
	 * In this case these bonds are assumed to have bond type cBondTypeSingle.
	 * Alternatively, one may pass null and indicate affected bonds with bond type cBondTypeDelocalized.
	 * Non-cyclic atom chains defined to be delocalized are treated depending
	 * on whether we have a molecule or a query fragment. For fragments the respective bond
	 * types will be set to cBondTypeDelocalized; for molecules the chain will
	 * have alternating single and double bonds starting with double at a non-ring end.
	 * @param isDelocalizedBond if null, then bond type cBondTypeDelocalized is used to indicate delocalized bonds
	 * @param mayChangeAtomCharges true if input molecule doesn't carry atom charges and these may be added to achieve aromaticity
	 * @param allHydrogensAreExplicit true this method can rely on all hydrogens being explicitly present
	 * @return true if all bonds of the delocalized area could be consistently converted.
	 */
	public boolean locateDelocalizedDoubleBonds(boolean[] isDelocalizedBond, boolean mayChangeAtomCharges, boolean allHydrogensAreExplicit) {
		// Helper arrays are calculated once only and stay untouched until entire kekulization is done!
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);

		RingCollection ringSet = new RingCollection(mMol, RingCollection.MODE_SMALL_RINGS_ONLY);
		initialize(isDelocalizedBond, ringSet);
		if (mDelocalizedBonds == 0)
            return true;

		if (mayChangeAtomCharges)
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (mIsDelocalizedAtom[atom] && mMol.getAtomicNo(atom) == 7
				 && (mIsDelocalizedBridgeHead[atom]
				  || (mMol.getConnAtoms(atom) == 3 && !mIsDelocalizedFiveRingMember[atom])))
					mMol.setAtomCharge(atom, 1);

		mAllHydrogensAreExplicit = allHydrogensAreExplicit;

		protectFullValenceAtoms(mayChangeAtomCharges);

        if (mMol.isFragment())
        	promoteDelocalizedChains();

		protectObviousNonRingLeaks(ringSet);

		protectAmideBonds(ringSet);	// using rules for detecting aromatic (thio-)amide bonds
		protectDoubleBondAtoms();

		promoteObviousBonds();

		while (promoteOuterShellDelocalizedRingSystems(ringSet, mayChangeAtomCharges))
			promoteObviousBonds();

		// try to find and promote entirely aromatic 6-rings
		if (promoteSixMemberedAromaticRings(ringSet, mayChangeAtomCharges))
			promoteObviousBonds();

		// try to find and promote other fully aromatic rings
		while (promoteOneDelocalizationLeak(ringSet, mayChangeAtomCharges))
			promoteObviousBonds();

		// Find and promote remaining aromatic bonds
        while (mDelocalizedBonds != 0) {
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (mIsDelocalizedBond[bond]) {
					promoteBond(bond);
					promoteObviousBonds();
					}
				}
            }

		if (mDelocalizedAtoms - mPiElectronsAdded >= 2)
			connectSeparatedSingletons();

		// If finally one or more single delocalized atoms remain, and if these carry an implicit hydrogen,
		// and if it is allowed and possible, then add/remove a charge to make the atom's lone pair contribute
		// to the neighbour atom's delocalization. If modifying charges is not allowed, then remove an implicit
		// hydrogen including one of its binding electrons. An unpaired electron is closer to the intended aromatic
		// state than a wrong hybridisation and one implicit hydrogen too much.
		// In addition, during id-coordinate parsing of implicit hydrogen atoms, the number of implicit hydrogens
		// per non-H atom must exactly match the one given when the idcode with coordinates was encoded.
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mIsDelocalizedAtom[atom] && mMol.getImplicitHydrogens(atom) != 0) {
				if (mayChangeAtomCharges
				 &&	((mMol.getAtomCharge(atom) == 1 && mMol.isElectronegative(atom))
				  || (mMol.getAtomCharge(atom) == -1 && mMol.getAtomicNo(atom) == 5)))
					mMol.setAtomCharge(atom, 0);
				else
					mMol.setAtomRadical(atom, Molecule.cAtomRadicalStateD);
				mPiElectronsAdded++;
				}
			}

		return (mDelocalizedAtoms == mPiElectronsAdded);
		}

	private void initialize(boolean[] isDelocalizedBond, RingCollection ringSet) {
		if (isDelocalizedBond != null) {
			mIsDelocalizedBond = isDelocalizedBond;
		}
		else {
			mIsDelocalizedBond = new boolean[mMol.getBonds()];
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (mMol.getBondType(bond) == Molecule.cBondTypeDelocalized) {
					mIsDelocalizedBond[bond] = true;
					mMol.setBondType(bond, Molecule.cBondTypeSingle);
				}
			}
		}

		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mIsDelocalizedBond[bond])
				mDelocalizedBonds++;

		mIsDelocalizedAtom = new boolean[mMol.getAtoms()];
		mIsDelocalizedBridgeHead = new boolean[mMol.getAtoms()];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			int delocalizedConnBonds = 0;
			for (int i=0; i<mMol.getConnAtoms(atom); i++)
				if (mIsDelocalizedBond[mMol.getConnBond(atom, i)])
					delocalizedConnBonds++;
			if (delocalizedConnBonds > 0) {
				mIsDelocalizedAtom[atom] = true;
				mDelocalizedAtoms++;
				if (delocalizedConnBonds == 3)
					mIsDelocalizedBridgeHead[atom] = true;
			}
		}

		mIsDelocalizedRing = new boolean[ringSet.getSize()];
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			mIsDelocalizedRing[ring] = true;
			for (int bond : ringSet.getRingBonds(ring)) {
				if (!mIsDelocalizedBond[bond]) {
					mIsDelocalizedRing[ring] = false;
					break;
				}
			}
		}

		mIsDelocalizedThreeOrSevenRingMember = new boolean[mMol.getAtoms()];
		mIsDelocalizedFiveRingMember = new boolean[mMol.getAtoms()];
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			if (mIsDelocalizedRing[ring] && ringSet.getRingSize(ring) != 6) {
				for (int atom : ringSet.getRingAtoms(ring)) {
					if (ringSet.getRingSize(ring) == 5)
						mIsDelocalizedFiveRingMember[atom] = true;
					else
						mIsDelocalizedThreeOrSevenRingMember[atom] = true;
				}
			}
		}

		mPiElectronsAdded = 0;
	}

	private boolean promoteSixMemberedAromaticRings(RingCollection ringSet, boolean mayChangeAtomCharges) {
		boolean found = false;
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			if (ringSet.getRingSize(ring) == 6) {
				boolean isQualifyingRing = true;
				int[] ringAtom = ringSet.getRingAtoms(ring);
				int[] ringBond = ringSet.getRingBonds(ring);
				for (int i=0; i<6; i++) {
					if (!checkAtomTypePi1(ringAtom[i], false)
					 || !mIsDelocalizedBond[ringBond[i]]) {
						isQualifyingRing = false;
						break;
					}
				}

				if (isQualifyingRing) {
					if (mayChangeAtomCharges)
						for (int i=0; i<6; i++)
							checkAtomTypePi1(ringAtom[i], true);

					for (int i=0; i<6; i+=2)
						promoteBond(ringBond[i]);

					found = true;
					break;
				}
			}
		}
		return found;
	}

	private boolean promoteOneDelocalizationLeak(RingCollection ringSet, boolean mayChangeAtomCharges) {
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			if (ringSet.getRingSize(ring) != 6 && mIsDelocalizedRing[ring]) {
				boolean isFullyDelocalizedRing = true;
				int[] ringBond = ringSet.getRingBonds(ring);
				for (int i=0; i<ringBond.length; i++) {
					if (!mIsDelocalizedBond[ringBond[i]]) {
						isFullyDelocalizedRing = false;
						break;
					}
				}

				if (isFullyDelocalizedRing) {
					int bestDelocalizationLeakIndex = -1;
					int bestDelocalizationLeakPriority = 0;
					int[] ringAtom = ringSet.getRingAtoms(ring);
					for (int i=0; i<ringAtom.length; i++) {
						int atom = ringAtom[i];
						int priority = mIsDelocalizedFiveRingMember[atom] ?
								checkAtomTypeLeak5(atom, false)
							  : checkAtomTypeLeak7(atom, false);
						if (bestDelocalizationLeakPriority < priority) {
							bestDelocalizationLeakPriority = priority;
							bestDelocalizationLeakIndex = i;
						}
					}
					if (bestDelocalizationLeakIndex != -1) {
						int atom = ringAtom[bestDelocalizationLeakIndex];
						if (mayChangeAtomCharges) {
							if (mIsDelocalizedFiveRingMember[atom])
								checkAtomTypeLeak5(atom, true);
							else
								checkAtomTypeLeak7(atom, true);
						}
						protectAtom(atom);
						return true;
					}
				}
			}
		}
		return false;
	}

	private boolean promoteOuterShellDelocalizedRingSystems(RingCollection ringSet, boolean mayChangeAtomCharges) {
		int[] sharedDelocalizedRingCount = new int[mMol.getBonds()];
		for (int r=0; r<ringSet.getSize(); r++) {
			int[] ringBond = ringSet.getRingBonds(r);
			boolean isDelocalized = true;
			for (int i=0; i<ringBond.length; i++) {
				if (!mIsDelocalizedBond[ringBond[i]]) {
					isDelocalized = false;
					break;
				}
			}
			if (isDelocalized)
				for (int i=0; i<ringBond.length; i++)
					sharedDelocalizedRingCount[ringBond[i]]++;
		}

		int delocalizedBonds = mDelocalizedBonds;

		for (int r=0; r<ringSet.getSize(); r++) {
			boolean found = false;
			int[] ringAtom = ringSet.getRingAtoms(r);
			int[] ringBond = ringSet.getRingBonds(r);
			for (int i=0; i<ringBond.length && !found; i++) {
				if (sharedDelocalizedRingCount[ringBond[i]] > 1) {
					int first = nextIndex(i, ringBond.length);
					if (sharedDelocalizedRingCount[ringBond[first]] == 1) {	// first unshared aromatic bond
						boolean incompatibleAtomFound = false;
						int next = nextIndex(first, ringBond.length);
						while (sharedDelocalizedRingCount[ringBond[next]] == 1) {
							if (!checkAtomTypePi1(ringAtom[next], false)
							 ||	(ringBond.length != 6 && mMol.getAtomicNo(ringAtom[next]) != 6))
								incompatibleAtomFound = true;
							next = nextIndex(next, ringBond.length);
						}
						if (!incompatibleAtomFound) {
							int bridgeBonds = (next > first) ? next - first : next + ringBond.length -first;
							if (bridgeBonds > 2 && (bridgeBonds & 1) == 1) {
								for (int j=1; j<bridgeBonds; j+=2) {
									int index = (first+j < ringBond.length) ? first+j : first+j-ringBond.length;
									if (mayChangeAtomCharges) {
										checkAtomTypePi1(ringAtom[index], true);
										checkAtomTypePi1(ringAtom[index == ringAtom.length-1 ? 0 : index+1], true);
									}
									promoteBond(ringBond[index]);
								}
								found = true;
							}
						}
					}
				}
			}
		}

		return delocalizedBonds != mDelocalizedBonds;
	}

	private int nextIndex(int i, int size) {
		return i == size-1 ? 0 : i+1;
	}

	private void protectFullValenceAtoms(boolean mayChangeAtomCharges) {
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mIsDelocalizedAtom[atom]
			 && mMol.getLowestFreeValence(atom) == 0
			 && (!mayChangeAtomCharges
			  || (mMol.getAtomicNo(atom) == 5 && mMol.getAtomCharge(atom) < 0)
			  || (mMol.getAtomicNo(atom) == 6 || mMol.getAtomicNo(atom) == 14)
			  || (mMol.isElectronegative(atom) && mMol.getAtomCharge(atom) > 0)))
				protectAtom(atom);
		}

	private void protectAtom(int atom) {
		if (mIsDelocalizedAtom[atom]) {
			mIsDelocalizedAtom[atom] = false;
			mDelocalizedAtoms--;
			}
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int connBond = mMol.getConnBond(atom, i);
            if (mIsDelocalizedBond[connBond]) {
                mIsDelocalizedBond[connBond] = false;
                mDelocalizedBonds--;
                }
            }
		}

	private void promoteBond(int bond) {
		if (mMol.getBondType(bond) == Molecule.cBondTypeSingle) {
			mMol.setBondType(bond, Molecule.cBondTypeDouble);
			mPiElectronsAdded += 2;
			}

		for (int i=0; i<2; i++) {
			int bondAtom = mMol.getBondAtom(i, bond);
			mIsDelocalizedAtom[bondAtom] = false;
			for (int j=0; j<mMol.getConnAtoms(bondAtom); j++) {
				int connBond = mMol.getConnBond(bondAtom, j);
                if (mIsDelocalizedBond[connBond]) {
                    mIsDelocalizedBond[connBond] = false;
                    mDelocalizedBonds--;
                    }
                }
			}
		}

	private void promoteObviousBonds() {
			// handle bond orders of aromatic bonds along the chains attached to 5- or 7-membered ring
		boolean terminalAromaticBondFound;
		do {
			terminalAromaticBondFound = false;
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (mIsDelocalizedBond[bond]) {
					boolean isTerminalAromaticBond = false;
					for (int i=0; i<2; i++) {
                        int bondAtom = mMol.getBondAtom(i, bond);
					    boolean aromaticNeighbourFound = false;
						for (int j=0; j<mMol.getConnAtoms(bondAtom); j++) {
							if (bond != mMol.getConnBond(bondAtom, j)
							 && mIsDelocalizedBond[mMol.getConnBond(bondAtom, j)]) {
								aromaticNeighbourFound = true;
								break;
								}
							}
						if (!aromaticNeighbourFound) {
							isTerminalAromaticBond = true;
							break;
							}
						}

					if (isTerminalAromaticBond) {
						terminalAromaticBondFound = true;
						promoteBond(bond);
						}
					}
				}
			} while (terminalAromaticBondFound);
		}

	/**
	 * Protect query features cBondQFDelocalized in open aromatic chains of fragments
	 * with incomplete aromatic rings
	 */
	private void promoteDelocalizedChains() {
        for (int bond=0; bond<mMol.getBonds(); bond++) {
            if (mIsDelocalizedBond[bond]) {
                for (int i=0; i<2; i++) {
                    int terminalAtom = mMol.getBondAtom(i, bond);
                    boolean aromaticNeighbourFound = false;
                    for (int j=0; j<mMol.getConnAtoms(terminalAtom); j++) {
                        if (bond != mMol.getConnBond(terminalAtom, j)
                         && mIsDelocalizedBond[mMol.getConnBond(terminalAtom, j)]) {
                            aromaticNeighbourFound = true;
                            break;
                            }
                        }
                    if (!aromaticNeighbourFound) {
                        int terminalBond = bond;
                        int bridgeAtom = mMol.getBondAtom(1-i, bond);
                        while (terminalBond != -1) {
							mIsDelocalizedAtom[terminalAtom] = false;
                            mIsDelocalizedBond[terminalBond] = false;
                            mDelocalizedBonds--;
                            mMol.setBondType(terminalBond, Molecule.cBondTypeDelocalized);
                            terminalBond = -1;
                            terminalAtom = bridgeAtom;
                            for (int j=0; j<mMol.getConnAtoms(terminalAtom); j++) {
                                if (mIsDelocalizedBond[mMol.getConnBond(terminalAtom, j)]) {
                                    if (terminalBond == -1) {
                                        terminalBond = mMol.getConnBond(terminalAtom, j);
                                        bridgeAtom = mMol.getConnAtom(terminalAtom, j);
                                        }
                                    else {
                                        terminalAtom = -1;	// Stop here! We have hit an aromatic branch.
                                        terminalBond = -1;
                                        break;
                                        }
                                    }
                                }
                            }
						if (terminalAtom != -1)		// Regular end of aromatic chain (no branch).
							mIsDelocalizedAtom[bridgeAtom] = false;
                        break;
                        }
                    }
                }
            }
        }

	private void protectAmideBonds(RingCollection ringSet) {
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mIsDelocalizedBond[bond] && ringSet.qualifiesAsAmideTypeBond(bond)) {
				protectAtom(mMol.getBondAtom(0, bond));
				protectAtom(mMol.getBondAtom(1, bond));
				}
			}
		}

	private void protectDoubleBondAtoms() {
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondOrder(bond) == 2) {
				for (int i=0; i<2; i++) {
					int atom = mMol.getBondAtom(i, bond);
					if (mMol.getAtomicNo(atom) <= 8) {
						for (int j=0; j<mMol.getConnAtoms(atom); j++) {
							int connBond = mMol.getConnBond(atom, j);
							if (mIsDelocalizedBond[connBond]) {
								protectAtom(atom);
								break;
								}
							}
						}
					}
				}
			}
		}

	private void protectObviousNonRingLeaks(RingCollection ringSet) {
		// for all ring atoms add charges and protect preferred delocalization leak atoms
/*		boolean[] isAromaticRingAtom = new boolean[mMol.getAtoms()];
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			int ringSize = ringSet.getRingSize(ring);
			if (ringSize == 3 || ringSize == 5 || ringSize == 6 || ringSize == 7) {
				if (mIsDelocalizedRing[ring]) {
					for (int atom:ringSet.getRingAtoms(ring))
						isAromaticRingAtom[atom] = true;

					boolean possible = true;
					int leakAtom = -1;
					int leakPriority = 0;
					int bridgeHeadCount = 0;

					for (int atom:ringSet.getRingAtoms(ring)) {
						if (mIsDelocalizedBridgeHead[atom])
							bridgeHeadCount++;

						if (ringSize == 6 || mIsDelocalizedBridgeHead[atom]) {	// 6-membered or bridgehead atom
							if (!checkAtomTypePi1(atom, false)) {
								possible = false;
								break;
								}
							}
						else {	// non-bridgehead in 5- or 7-membered ring
							int priority = (ringSize == 5) ?
									checkAtomTypeLeak5(atom, false) : checkAtomTypeLeak7(atom, false);
							if (!checkAtomTypePi1(atom, false)) {
								if (leakPriority == 10) {
									possible = false;
									break;
									}
								leakAtom = atom;
								leakPriority = 20;	// MAX
								}
							else if (leakPriority < priority) {
								leakPriority = priority;
								leakAtom = atom;
								}
							}
						}

					if (possible && (ringSize == 6 || bridgeHeadCount == 0)) {
						for (int atom : ringSet.getRingAtoms(ring)) {
							if (atom == leakAtom) {
								if (ringSize == 5)
									checkAtomTypeLeak5(atom, true);	// 5-membered
								else
									checkAtomTypeLeak7(atom, true);	// 3- or 7-membered

								protectAtom(atom);
								}
							else {
								checkAtomTypePi1(atom, true);
								}
							}
						}
					}
				}
			}*/

		boolean[] isAromaticRingAtom = new boolean[mMol.getAtoms()];
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			int ringSize = ringSet.getRingSize(ring);
			if (ringSize == 3 || ringSize == 5 || ringSize == 6 || ringSize == 7) {
				if (mIsDelocalizedRing[ring]) {
					for (int atom : ringSet.getRingAtoms(ring))
						isAromaticRingAtom[atom] = true;
				}
			}
		}

		// From here locate delocalized strings of atoms, which are not member
		// of an aromatic ring. Protect preferred atoms and add obvious atom charges.

		// count for every atom the number of delocalized bonds attached
		int[] delocalizedNeighbourCount = new int[mMol.getAtoms()];
		boolean[] hasMetalLigandBond = new boolean[mMol.getAtoms()];
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			int atom1 = mMol.getBondAtom(0, bond);
			int atom2 = mMol.getBondAtom(1, bond);
			if (!isAromaticRingAtom[atom1] && !isAromaticRingAtom[atom2]) {
				if (mIsDelocalizedBond[bond]) {
					delocalizedNeighbourCount[atom1]++;
					delocalizedNeighbourCount[atom2]++;
					}
				if (mMol.getBondType(bond) == Molecule.cBondTypeMetalLigand) {
					hasMetalLigandBond[atom1] = true;
					hasMetalLigandBond[atom2] = true;
					}
				}
			}

		// From any delocalized atom with one delocalized neighbor (chain end)
		// locate the path to a branch atom (including) or chain end, whatever comes first.
		// Then mark every second atom from the startAtom as not being capable to assume
		// the role of a delocalization leak (priority:-1).
		int[] priority = new int[mMol.getAtoms()];
		int[] graphAtom = new int[mMol.getAtoms()];
		for (int seedAtom=0; seedAtom<mMol.getAtoms(); seedAtom++) {
			if (delocalizedNeighbourCount[seedAtom] == 1) {
				graphAtom[0] = seedAtom;
				int current = 0;
				int highest = 0;
				while (current <= highest) {
					for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
						if (mIsDelocalizedBond[mMol.getConnBond(graphAtom[current], i)]) {
							int candidate = mMol.getConnAtom(graphAtom[current], i);
							if ((current == 0 || candidate != graphAtom[current-1])
							 && delocalizedNeighbourCount[candidate] != 0) {
								graphAtom[++highest] = candidate;
								if ((delocalizedNeighbourCount[candidate] & 1) != 0) {	// 1 or 3
									for (int j=1; j<highest; j+=2)
										priority[graphAtom[j]] = -1;
									highest = 0;	// to break outer loop
									}
								break;
								}
							}
						}
					current++;
					}
				}
			}

		// For every connected delocalized area not being part of an aromatic ring
		// calculate delocalization leak priorities for all atoms not marked above.
		// Then protect the atom with the highest priority.
		boolean[] atomHandled = new boolean[mMol.getAtoms()];
		for (int seedAtom=0; seedAtom<mMol.getAtoms(); seedAtom++) {
			if (!atomHandled[seedAtom] && delocalizedNeighbourCount[seedAtom] != 0) {
				graphAtom[0] = seedAtom;
				atomHandled[seedAtom] = true;
				int current = 0;
				int highest = 0;
				while (current <= highest) {
					for (int i = 0; i < mMol.getConnAtoms(graphAtom[current]); i++) {
						if (mIsDelocalizedBond[mMol.getConnBond(graphAtom[current], i)]) {
							int candidate = mMol.getConnAtom(graphAtom[current], i);
							if (!atomHandled[candidate]) {
								graphAtom[++highest] = candidate;
								atomHandled[candidate] = true;
								}
							}
						}
					current++;
					}

				// if we have an odd number of delocalized atoms in a region, we need to assign a leak
				if ((highest & 1) == 0) {	// highest is atom count-1

					// check for all potential delocalization leak atoms, whether they are compatible
					for (int i = 0; i <= highest; i++)
						if (priority[graphAtom[i]] == 0)
							priority[graphAtom[i]] = checkAtomTypeLeakNonRing(graphAtom[i], false);

					// check for all atoms, which cannot be the leak, whether they can carry a pi-bond
					boolean isPossible = true;
					for (int i = 0; i <= highest; i++) {
						if (priority[graphAtom[i]] <= 0) {
							if (!checkAtomTypePi1(graphAtom[i], false)) {
								isPossible = false;
								break;
								}
							}
						}

					// find the preferred atom for the leak
					if (isPossible) {
						int maxPriority = 0;
						int maxAtom = -1;
						for (int i = 0; i <= highest; i++) {
							if (maxPriority < priority[graphAtom[i]]) {
								maxPriority = priority[graphAtom[i]];
								maxAtom = graphAtom[i];
								}
							}

						if (maxPriority > 0) {
							checkAtomTypeLeakNonRing(maxAtom, true);
							protectAtom(maxAtom);
							}
						}
					}
				}
			}
		}

	private void connectSeparatedSingletons() {
		for (int atom = 0; atom<mMol.getAtoms(); atom++) {
			mMol.ensureHelperArrays(Molecule.cHelperNeighbours);	// make sure, connBondOrders are
			if (mIsDelocalizedAtom[atom]) {
				boolean found = false;

				int[] graphAtom = new int[mMol.getAtoms()];
				int[] parentAtom = new int[mMol.getAtoms()];
				int[] graphLevel = new int[mMol.getAtoms()];

				graphAtom[0] = atom;
				parentAtom[atom] = -1;
				graphLevel[atom] = 1;

				int current = 0;
				int highest = 0;

				while (current <= highest && !found) {
					int currentAtom = graphAtom[current];
					for (int i=0; i<mMol.getConnAtoms(currentAtom) && !found; i++) {
						// We need alternating bond orders!
						boolean isCompatibleBond = ((graphLevel[currentAtom] & 1) == 1) ^ (mMol.getBondOrder(mMol.getConnBond(currentAtom, i)) > 1);
						int candidate = mMol.getConnAtom(currentAtom, i);
						if (graphLevel[candidate] == 0 && isCompatibleBond) {
							if (mIsDelocalizedAtom[candidate]) {
								if ((graphLevel[currentAtom] & 1) == 1) {	// odd number of bonds
									mIsDelocalizedAtom[atom] = false;
									mIsDelocalizedAtom[candidate] = false;
									mPiElectronsAdded += 2;
									int parent = currentAtom;
									for (int j=0; j<graphLevel[currentAtom]; j++) {
										// trace candidate back to root atom and invert bond orders!
										int bond = mMol.getBond(candidate, parent);
										if (mIsDelocalizedBond[bond]) {
											mIsDelocalizedBond[bond] = false;
											mDelocalizedBonds--;
										}
										mMol.setBondOrder(bond, mMol.getBondOrder(bond) == 1 ? 2 : mMol.getBondOrder(bond)-1);
										candidate = parent;
										parent = parentAtom[candidate];
									}
									found = true;
								}
							}
							else {
								graphAtom[++highest] = candidate;
								parentAtom[candidate] = currentAtom;
								graphLevel[candidate] = graphLevel[currentAtom]+1;
							}
						}
					}
					current++;
				}
			}
		}
	}

	/**
	 * Checks, whether the atom is compatible with an aromatic atom of the type
	 * that carries one half of a delocalized double bond.
	 * @param atom
	 * @param correctCharge if true then may add a charge to make the atom compatible
	 * @return
	 */
	private boolean checkAtomTypePi1(int atom, boolean correctCharge) {
		int atomicNo = mMol.getAtomicNo(atom);
		if ((atomicNo >= 5 && atomicNo <= 8)
		 || atomicNo == 15 || atomicNo == 16 || atomicNo == 33 || atomicNo == 34 || atomicNo == 52) {	// P,S,As,Se,Te

			int freeValence = mMol.getLowestFreeValence(atom);
			if (freeValence != 0)
				return true;

			int charge = mMol.getAtomCharge(atom);
			if (atomicNo == 5 && charge >= 0) {
				if (correctCharge)
					mMol.setAtomCharge(atom, charge-1);
				return true;
				}
			if (atomicNo != 5 && charge <= 0) {
				if (correctCharge)
					mMol.setAtomCharge(atom, charge+1);
				return true;
				}
			}

		return false;
		}

	/**
	 * Checks, whether the atom is compatible with that aromatic atom of
	 * a 5-membered ring that supplies the additional electron pair.
	 * @param atom
	 * @param correctCharge if true then may add a charge to make the atom compatible
	 * @return 0 (not compatible) or priority to be used (higher numbers have higher priority)
	 */
	private int checkAtomTypeLeak5(int atom, boolean correctCharge) {
		if (mIsDelocalizedBridgeHead[atom])
			return 0;

		if (mMol.getAtomicNo(atom) == 7) {
			if (mMol.getAllConnAtoms(atom) == 3)
				return 6;
			else if (mMol.getConnAtoms(atom) == 2)
				return mAllHydrogensAreExplicit ? 0 : 4;
			}
		else if (mMol.getAtomicNo(atom) == 8) {
			return 10;
			}
		else if (mMol.getAtomicNo(atom) == 15 || mMol.getAtomicNo(atom) == 33) {
			if (mMol.getConnAtoms(atom) == 3)
				return 8;
			}
		else if (mMol.getAtomicNo(atom) == 16 || mMol.getAtomicNo(atom) == 34 || mMol.getAtomicNo(atom) == 52) {
			if (mMol.getConnAtoms(atom) == 2)
				return 11;
			if (mMol.getConnAtoms(atom) == 3) {
				if (mMol.getAtomCharge(atom) == 1)
					return 12;
				if (correctCharge)
					mMol.setAtomCharge(atom, 1);
				return 5;
				}
			}
		else if (mMol.getAtomicNo(atom) == 6) {
			if (mMol.getAtomCharge(atom) == -1)
				return mMol.getAllConnAtoms(atom) == 3 ? 16
					 : mMol.getAllConnAtomsPlusMetalBonds(atom) == 3 ? 15 : 14;;
			if (correctCharge)
				mMol.setAtomCharge(atom, -1);
			return (mMol.getAllConnAtoms(atom) != mMol.getAllConnAtomsPlusMetalBonds(atom)) ? 2 : 3;
			}

		return 0;
		}


	/**
	 * Checks, whether the atom is compatible with that aromatic atom of
	 * a 3- or 7-membered ring that supplies the empty orbital.
	 * @param atom
	 * @param correctCharge if true then may add a charge to make the atom compatible
	 * @return 0 (not compatible) or priority to be used (higher numbers have higher priority)
	 */
	private int checkAtomTypeLeak7(int atom, boolean correctCharge) {
		if (mIsDelocalizedBridgeHead[atom])
			return 0;

		if (mAllHydrogensAreExplicit) {
			if (mMol.getAllConnAtoms(atom) != 3)
				return 0;
			}
		else {
			if (mMol.getAllConnAtoms(atom) > 3)
				return 0;
			}

		if (mMol.getAtomicNo(atom) == 6) {
			if (correctCharge)
				mMol.setAtomCharge(atom, 1);
			return 2;
			}
		if (mMol.getAtomicNo(atom) == 5) {
			return 4;
			}

		return 0;
		}

	/**
	 * Checks, whether the atom is compatible with the (typically charged) atom
	 * in a delocalized chain of an odd number of atoms that does not carry a pi bond.
	 * @param atom
	 * @param correctCharge if true then may add a charge to make the atom compatible
	 * @return 0 (not compatible) or priority to be used (higher numbers have higher priority)
	 */
	private int checkAtomTypeLeakNonRing(int atom, boolean correctCharge) {
		if (mMol.getAtomCharge(atom) != 0)
			return 0;

		if (mAllHydrogensAreExplicit) {
			if (mMol.getAtomicNo(atom) == 5) {
				if (mMol.getOccupiedValence(atom) != 2)
					return 0;
				if (correctCharge)
					mMol.setAtomCharge(atom, 1);
				return 1;
				}

			if (mMol.getAtomicNo(atom) == 7) {
				if (mMol.getOccupiedValence(atom) != 2)
					return 0;
				if (correctCharge)
					mMol.setAtomCharge(atom, -1);
				return hasMetalNeighbour(atom) ? 6 : 3;
				}

			if (mMol.getAtomicNo(atom) == 8) {
				if (mMol.getOccupiedValence(atom) != 1)
					return 0;
				if (correctCharge)
					mMol.setAtomCharge(atom, -1);
				return hasMetalNeighbour(atom) ? 7 : 4;
				}

			if (mMol.getAtomicNo(atom) == 16) {
				if (mMol.getOccupiedValence(atom) != 1)
					return 0;
				if (correctCharge)
					mMol.setAtomCharge(atom, -1);
				return hasMetalNeighbour(atom) ? 5 : 2;
				}

			if (mMol.getAtomicNo(atom) == 34) {
				if (mMol.getOccupiedValence(atom) != 1)
					return 0;
				if (correctCharge)
					mMol.setAtomCharge(atom, -1);
				return hasMetalNeighbour(atom) ? 4 : 1;
				}
			}
		else {
			if (mMol.getAtomicNo(atom) == 5) {
				if (mMol.getOccupiedValence(atom) > 2)
					return 0;
				if (correctCharge)
					mMol.setAtomCharge(atom, 1);
				return 1;
				}

			if (mMol.getAtomicNo(atom) == 7) {
				if (mMol.getOccupiedValence(atom) > 2)
					return 0;
				if (correctCharge)
					mMol.setAtomCharge(atom, -1);
				return hasMetalNeighbour(atom) ? 5 : 3;
				}

			if (mMol.getAtomicNo(atom) == 8) {
				if (mMol.getOccupiedValence(atom) > 1)
					return 0;
				if (correctCharge)
					mMol.setAtomCharge(atom, -1);
				return hasMetalNeighbour(atom) ? 7 : 4;
				}

			if (mMol.getAtomicNo(atom) == 16) {
				if (mMol.getOccupiedValence(atom) > 1)
					return 0;
				if (correctCharge)
					mMol.setAtomCharge(atom, -1);
				return hasMetalNeighbour(atom) ? 5 : 2;
				}
			}

		return 0;
		}

	private boolean hasMetalNeighbour(int atom) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (mMol.isMetalAtom(mMol.getConnAtom(atom, i)))
				return true;

		return false;
		}
	}
