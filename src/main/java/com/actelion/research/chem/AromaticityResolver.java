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
	private boolean[]	mIsDelocalizedAtom,mIsDelocalizedBond;
    private int         mAromaticAtoms,mAromaticBonds,mPiElectronsAdded;

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
	 * @param isAromaticBond if null, then bond type cBondTypeDelocalized is used to indicate delocalized bonds
	 * @param mayChangeAtomCharges true if input molecule doesn't carry atom charges and these may be added to achieve aromaticity
	 * @param allHydrogensAreExplicit true this method can rely on all hydrogens being explicitly present
	 * @return true if all bonds of the delocalized area could be consistently converted.
	 */
	public boolean locateDelocalizedDoubleBonds(boolean[] isAromaticBond, boolean mayChangeAtomCharges, boolean allHydrogensAreExplicit) {
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);

		if (isAromaticBond != null) {
			mIsDelocalizedBond = isAromaticBond;
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

		mPiElectronsAdded = 0;

		mIsDelocalizedAtom = new boolean[mMol.getAtoms()];
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mIsDelocalizedBond[bond]) {
				mAromaticBonds++;
				for (int i=0; i<2; i++) {
					if (!mIsDelocalizedAtom[mMol.getBondAtom(i, bond)]) {
						mIsDelocalizedAtom[mMol.getBondAtom(i, bond)] = true;
						mAromaticAtoms++;
						}
					}
				}
			}

		if (mAromaticBonds == 0)
            return true;

		mAllHydrogensAreExplicit = allHydrogensAreExplicit;

		protectFullValenceAtoms(mayChangeAtomCharges);

        if (mMol.isFragment())
        	promoteDelocalizedChains();

		// create small-ring set without aromaticity information
		RingCollection ringSet = new RingCollection(mMol, RingCollection.MODE_SMALL_RINGS_ONLY);

		if (mayChangeAtomCharges)
			addObviousAtomCharges(ringSet);

        // find mandatory conjugation breaking atoms in 3-, 5-, and 7-membered rings,
		// i.e. atoms whose neighbour bonds must be single bonds
		protectObviousDelocalizationLeaks(ringSet);

		protectAmideBonds(ringSet);	// using rules for detecting aromatic (thio-)amide bonds
		protectDoubleBondAtoms();

		promoteObviousBonds();

		while (promoteOuterShellDelocalizedRingSystems(ringSet))
			promoteObviousBonds();

        while (mAromaticBonds != 0) {
            boolean bondsPromoted = false;

            if (!bondsPromoted) {
                // try to find and promote one entire aromatic 6-ring
                for (int ring=0; ring<ringSet.getSize(); ring++) {
                    if (ringSet.getRingSize(ring) == 6) {
                        boolean isAromaticRing = true;
                        int[] ringBond = ringSet.getRingBonds(ring);
                        for (int i=0; i<6; i++) {
                            if (!mIsDelocalizedBond[ringBond[i]]) {
                                isAromaticRing = false;
                                break;
                                }
                            }
    
                        if (isAromaticRing) {
                            for (int i=0; i<6; i+=2)
                                promoteBond(ringBond[i]);
                            bondsPromoted = true;
                            break;
                            }
                        }
                    }
                }


			if (bondsPromoted) {
				promoteObviousBonds();
				continue;
				}

			if (!bondsPromoted) {
                // find and promote one aromatic bond
                // (should never happen, but to prevent an endless loop nonetheless)
                for (int bond=0; bond<mMol.getBonds(); bond++) {
                    if (mIsDelocalizedBond[bond]) {
                        promoteBond(bond);
                        promoteObviousBonds();
                        bondsPromoted = true;
                        break;
                        }
                    }
                }
            }

/*		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mIsDelocalizedAtom[atom] && mMol.getImplicitHydrogens(atom) != 0) {
				mMol.setAtomRadical(atom, Molecule.cAtomRadicalStateD);
				mPiElectronsAdded++;
				}
			}*/

		return (mAromaticAtoms == mPiElectronsAdded);
		}


	private void protectObviousDelocalizationLeaks(RingCollection ringSet) {
		for (int r=0; r<ringSet.getSize(); r++) {
			int ringSize = ringSet.getRingSize(r);
			if (ringSize == 3 || ringSize == 5 || ringSize == 7) {
				int[] ringAtom = ringSet.getRingAtoms(r);
				for (int i=0; i<ringSize; i++) {
					int atom = ringAtom[i];
					if (isAromaticAtom(atom)) {
						if (ringSize == 5) {
							// C(-),N,O,S in cyclopentadienyl, furan, pyrrol, etc.
							if ((mMol.getAtomicNo(atom) == 6 && mMol.getAtomCharge(atom) == -1 && mMol.getAllConnAtoms(atom) == 3)
							 || (mMol.getAtomicNo(atom) == 7 && mMol.getAtomCharge(atom) == 0 && mMol.getAllConnAtoms(atom) == 3)
							 || (mMol.getAtomicNo(atom) == 8 && mMol.getAtomCharge(atom) == 0 && mMol.getConnAtoms(atom) == 2)
							 || (mMol.getAtomicNo(atom) == 16 && mMol.getAtomCharge(atom) == 0 && mMol.getConnAtoms(atom) == 2)
							 || (mMol.getAtomicNo(atom) == 34 && mMol.getAtomCharge(atom) == 0 && mMol.getConnAtoms(atom) == 2))
								protectAtom(atom);
							}
						else {
							// B,C+ in tropylium
							if ((mMol.getAtomicNo(atom) == 5 && mMol.getAtomCharge(atom) == 0 && mMol.getAllConnAtoms(atom) == 3)
							 || (mMol.getAtomicNo(atom) == 6 && mMol.getAtomCharge(atom) == 1))
								protectAtom(atom);
							}
						}
					}
				}
			}

		// in 5-membered rings with 5 delocalized bonds and more than one negative carbon, we choose the most obvious as leak
		for (int r=0; r<ringSet.getSize(); r++) {
			if (ringSet.getRingSize(r) == 5) {
				int[] ringBond = ringSet.getRingBonds(r);
				boolean isDelocalized = true;
				for (int i=0; i<ringBond.length; i++) {
					if (!mIsDelocalizedBond[ringBond[i]]) {
						isDelocalized = false;
						break;
						}
					}
				if (isDelocalized) {
					int[] ringAtom = ringSet.getRingAtoms(r);
					int negativeCarbonPriority = 0;
					int negativeCarbon = -1;
					for (int i=0; i<ringBond.length; i++) {
						if (mMol.getAtomCharge(ringAtom[i]) == -1 && mMol.getAtomicNo(ringAtom[i]) == 6) {
							int priority = mMol.getAllConnAtoms(ringAtom[i]) == 3 ? 3
										 : mMol.getAllConnAtomsPlusMetalBonds(ringAtom[i]) == 3 ? 2 : 1;
							if (negativeCarbonPriority < priority) {
								negativeCarbonPriority = priority;
								negativeCarbon = ringAtom[i];
								}
							}
						}
					if (negativeCarbon != -1)
						protectAtom(negativeCarbon);
					}
				}
			}
		}


	private boolean promoteOuterShellDelocalizedRingSystems(RingCollection ringSet) {
		int[] sharedDelocalizedRingCount = new int[mMol.getBonds()];
		for (int r=0; r<ringSet.getSize(); r++) {
			int[] ringBond = ringSet.getRingBonds(r);
			boolean isDelocalized = true;
			for (int i = 0; i < ringBond.length; i++) {
				if (!mIsDelocalizedBond[ringBond[i]]) {
					isDelocalized = false;
					break;
					}
				}
			if (isDelocalized)
				for (int i=0; i<ringBond.length; i++)
					sharedDelocalizedRingCount[ringBond[i]]++;
			}

		int delocalizedBonds = mAromaticBonds;

		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (sharedDelocalizedRingCount[bond]==1) {
				for (int i=0; i<2 && mIsDelocalizedBond[bond]; i++) {
					int atom1 = mMol.getBondAtom(i, bond);
					int atom2 = mMol.getBondAtom(1-i, bond);
					if (hasSharedDelocalizedBond(atom1, sharedDelocalizedRingCount)
					 && !hasSharedDelocalizedBond(atom2, sharedDelocalizedRingCount)) {
						int connIndex;
						while (-1 != (connIndex = getNextOuterDelocalizedConnIndex(atom2, atom1, sharedDelocalizedRingCount))) {
							int atom3 = mMol.getConnAtom(atom2, connIndex);
							int bond2to3 = mMol.getConnBond(atom2, connIndex);
							if (!mIsDelocalizedBond[bond2to3])
								break;

							promoteBond(bond2to3);
							connIndex = getNextOuterDelocalizedConnIndex(atom3, atom2, sharedDelocalizedRingCount);
							if (connIndex == -1)
								break;

							atom1 = atom3;
							atom2 = mMol.getConnAtom(atom3, connIndex);
							}
						}
					}
				}
			}

		return delocalizedBonds != mAromaticBonds;
		}


	private boolean hasSharedDelocalizedBond(int atom, int[] sharedDelocalizedRingCount) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (sharedDelocalizedRingCount[mMol.getConnBond(atom, i)] > 1)
				return true;
		return false;
		}


	private int getNextOuterDelocalizedConnIndex(int atom, int previousAtom, int[] sharedDelocalizedRingCount) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (sharedDelocalizedRingCount[mMol.getConnBond(atom, i)] == 1
			 && mMol.getConnAtom(atom, i) != previousAtom)
				return i;
		return -1;
		}


	private boolean isAromaticAtom(int atom) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (mIsDelocalizedBond[mMol.getConnBond(atom, i)])
				return true;
		return false;
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
			mAromaticAtoms--;
			}
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int connBond = mMol.getConnBond(atom, i);
            if (mIsDelocalizedBond[connBond]) {
                mIsDelocalizedBond[connBond] = false;
                mAromaticBonds--;
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
                    mAromaticBonds--;
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


	private void promoteDelocalizedChains() {
        // protect query features cBondQFDelocalized in open aromatic chains of fragments
	    // with incomplete aromatic rings
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
                            mIsDelocalizedBond[terminalBond] = false;
                            mAromaticBonds--;
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
                                        terminalAtom = -1;
                                        terminalBond = -1;
                                        break;
                                        }
                                    }
                                }
                            }
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

	private void addObviousAtomCharges(RingCollection ringSet) {
		// count for every atom of how many delocalized rings it is a member
		boolean[] isDelocalized = new boolean[ringSet.getSize()];
		int[] delocalizedRingCount = new int[mMol.getAtoms()];
		for (int r=0; r<ringSet.getSize(); r++) {
			isDelocalized[r] = true;
			for (int bond:ringSet.getRingBonds(r)) {
				if (!mIsDelocalizedBond[bond]) {
					isDelocalized[r] = false;
					break;
					}
				}
			if (isDelocalized[r])
				for (int atom:ringSet.getRingAtoms(r))
					delocalizedRingCount[atom]++;
			}

		// for all ring atoms add charges and protect preferred delocalization leak atoms
		boolean[] isAromaticRingAtom = new boolean[mMol.getAtoms()];
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			int ringSize = ringSet.getRingSize(ring);
			if (ringSize == 3 || ringSize == 5 || ringSize == 6 || ringSize == 7) {
				if (isDelocalized[ring]) {
					for (int atom:ringSet.getRingAtoms(ring))
						isAromaticRingAtom[atom] = true;

					boolean possible = true;
					int leakAtom = -1;
					int leakPriority = 0;

					for (int atom:ringSet.getRingAtoms(ring)) {
						if (ringSize == 6 || delocalizedRingCount[atom] > 1) {	// bridgehead atom
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

					if (possible) {
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
		int graphAtom[] = new int[mMol.getAtoms()];
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
		// Then protect the atom with highest priority.
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
		if (mMol.getAtomicNo(atom) == 7) {
			if (mMol.getAllConnAtoms(atom) == 3)
				return 6;
			else if (mMol.getConnAtoms(atom) == 2)
				return 4;
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
				return 12;
			}
		else if (mMol.getAtomicNo(atom) == 6) {
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
