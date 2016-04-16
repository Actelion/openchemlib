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

public class AromaticityResolver {
	ExtendedMolecule	mMol;
	private boolean[]	mIsAromaticBond;
    private int         mAromaticAtoms,mAromaticBonds,mPiElectronsAdded;

    /**
     * Creates a new AromaticityResolver for molecule mol that assumes that
     * all delocalized bonds of mol have bondType Molecule.cBondTypeDelocalized.
     * Internally this constructor creates a bond masks encoding the
     * delocalization from the bond type, then changes these bond type to
     * Molecule.cBondTypeSingle and calls the other constructor.
     * @param mol
     */
    public AromaticityResolver(ExtendedMolecule mol) {
        mMol = mol;

        mol.ensureHelperArrays(Molecule.cHelperNeighbours);

        mIsAromaticBond = new boolean[mol.getBonds()];
        for (int bond=0; bond<mol.getBonds(); bond++) {
            if (mol.getBondType(bond) == Molecule.cBondTypeDelocalized) {
                mIsAromaticBond[bond] = true;
                mol.setBondType(bond, Molecule.cBondTypeSingle);
                }
            }

        initialize();
    	}

    /**
     * Creates a new AromaticityResolver for molecule mol that assumes that
     * all delocalized bonds are flagged properly in isAromaticBond.
     * BondTypes of these bonds are assumed to be Molecule.cBondTypeSingle.
     * @param mol
     * @param isAromaticBond
     */
    public AromaticityResolver(ExtendedMolecule mol, boolean[] isAromaticBond) {
        mMol = mol;
		mIsAromaticBond = isAromaticBond;

        mMol.ensureHelperArrays(Molecule.cHelperNeighbours);

        initialize();
        }

    private void initialize() {
    	mPiElectronsAdded = 0;

    	boolean[] isAromaticAtom = new boolean[mMol.getAtoms()];
        for (int bond=0; bond<mMol.getBonds(); bond++) {
            if (mIsAromaticBond[bond]) {
            	mAromaticBonds++;
            	for (int i=0; i<2; i++) {
            		if (!isAromaticAtom[mMol.getBondAtom(i, bond)]) {
            			isAromaticAtom[mMol.getBondAtom(i, bond)] = true;
        				mAromaticAtoms++;
            			}
            		}
            	}
            }
    	}

    /**
     * This method promotes all necessary bonds of the defined delocalized part of the molecule
     * from single to double bonds in order to create a valid delocalized system
     * of conjugated single and double bonds.
     * Non-cyclic atom chains defined to be delocalized are treated depending
     * on whether we have a molecule or a query fragment. For fragments the respective bond
     * types will be set to cBondTypeDelocalized; for molecules the chain will
     * have alternating single and double bonds starting with double at a non-ring end.
     * @return true if all bonds of the delocalized area could be consistently converted. 
     */
	public boolean locateDelocalizedDoubleBonds() {
        if (mAromaticBonds == 0)
            return true;

        RingCollection ringSet = new RingCollection(mMol, RingCollection.MODE_SMALL_RINGS_ONLY);

        if (mMol.isFragment())	
        	promoteDelocalizedChains();

        // find mandatory conjugation breaking atoms, i.e. atoms whose neighbour bonds must (!) be single bonds
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (isAromaticAtom(atom)
			 && mMol.getAtomicNo(atom) ==  6) {
	            // C+ in tropylium and C- in cyclopentadienyl
				if ((mMol.getAtomCharge(atom) == -1 && ringSet.getAtomRingSize(atom) == 5)
				 || (mMol.getAtomCharge(atom) == 1 && ringSet.getAtomRingSize(atom) == 7))
					protectAtom(atom);
				}
			}

		promoteObviousBonds();

        while (mAromaticBonds != 0) {
            boolean bondsPromoted = false;
            for (int bond=0; bond<mMol.getBonds(); bond++) {
                if (mIsAromaticBond[bond]) {
                    int aromaticConnBonds = 0;
                    for (int j=0; j<2; j++) {
                        int bondAtom = mMol.getBondAtom(j, bond);
                        for (int k=0; k<mMol.getConnAtoms(bondAtom); k++)
                            if (mIsAromaticBond[mMol.getConnBond(bondAtom, k)])
                                aromaticConnBonds++;
                        }
    
                    if (aromaticConnBonds == 4) {
                        promoteBond(bond);
                        promoteObviousBonds();
                        bondsPromoted = true;
                        break;
                        }
                    }
                }

            if (!bondsPromoted) {
                // try to find and promote one entire aromatic 6-ring
                for (int ring=0; ring<ringSet.getSize(); ring++) {
                    if (ringSet.getRingSize(ring) == 6) {
                        boolean isAromaticRing = true;
                        int[] ringBond = ringSet.getRingBonds(ring);
                        for (int i=0; i<6; i++) {
                            if (!mIsAromaticBond[ringBond[i]]) {
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

            if (!bondsPromoted) {
                // find and promote one aromatic bond
                // (should never happen, but to prevent an endless loop nonetheless)
                for (int bond=0; bond<mMol.getBonds(); bond++) {
                    if (mIsAromaticBond[bond]) {
                        promoteBond(bond);
                        promoteObviousBonds();
                        bondsPromoted = true;
                        break;
                        }
                    }
                }
            }

        return (mAromaticAtoms == mPiElectronsAdded);
		}


	private boolean isAromaticAtom(int atom) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (mIsAromaticBond[mMol.getConnBond(atom, i)])
				return true;
		return false;
		}


	private void protectAtom(int atom) {
        mAromaticAtoms--;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int connBond = mMol.getConnBond(atom, i);
            if (mIsAromaticBond[connBond]) {
                mIsAromaticBond[connBond] = false;
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
			for (int j=0; j<mMol.getConnAtoms(bondAtom); j++) {
				int connBond = mMol.getConnBond(bondAtom, j);
                if (mIsAromaticBond[connBond]) {
                    mIsAromaticBond[connBond] = false;
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
				if (mIsAromaticBond[bond]) {
					boolean isTerminalAromaticBond = false;
					for (int i=0; i<2; i++) {
                        int bondAtom = mMol.getBondAtom(i, bond);
					    boolean aromaticNeighbourFound = false;
						for (int j=0; j<mMol.getConnAtoms(bondAtom); j++) {
							if (bond != mMol.getConnBond(bondAtom, j)
							 && mIsAromaticBond[mMol.getConnBond(bondAtom, j)]) {
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
            if (mIsAromaticBond[bond]) {
                for (int i=0; i<2; i++) {
                    int terminalAtom = mMol.getBondAtom(i, bond);
                    boolean aromaticNeighbourFound = false;
                    for (int j=0; j<mMol.getConnAtoms(terminalAtom); j++) {
                        if (bond != mMol.getConnBond(terminalAtom, j)
                         && mIsAromaticBond[mMol.getConnBond(terminalAtom, j)]) {
                            aromaticNeighbourFound = true;
                            break;
                            }
                        }
                    if (!aromaticNeighbourFound) {
                        int terminalBond = bond;
                        int bridgeAtom = mMol.getBondAtom(1-i, bond);
                        while (terminalBond != -1) {
                            mIsAromaticBond[terminalBond] = false;
                            mAromaticBonds--;
                            mMol.setBondType(terminalBond, Molecule.cBondTypeDelocalized);
                            terminalBond = -1;
                            terminalAtom = bridgeAtom;
                            for (int j=0; j<mMol.getConnAtoms(terminalAtom); j++) {
                                if (mIsAromaticBond[mMol.getConnBond(terminalAtom, j)]) {
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
    }
