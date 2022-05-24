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

@Deprecated // use Canonizer with mode NEGLECT_ANY_STEREO_INFORMATION instead
public class SimpleCanonizer {
    private static final int cIDCodeVersion3 = 9;
    // productive version since May 2006 based on the molfile version 3
    // being compatible with MDL's "Enhanced Stereo Representation"

    private ExtendedMolecule mMol;
	private int mCanRank[];
	private long mCanBaseValue[];

	private boolean mGraphGenerated;
	private int mGraphRings;
	private int mGraphAtom[];
	private int mGraphBond[];
	private int mGraphFrom[];
	private int mGraphClosure[];

    public static final int MAX_ATOM_BITS = 8;

	private String         	mIDCode,mCoordinates;
    private StringBuffer    mEncodingBuffer;
    private int             mEncodingBitsAvail,mEncodingTempData;
    private boolean			mZCoordinatesAvailable;

	public SimpleCanonizer(ExtendedMolecule mol) {
		mMol = mol;
		mMol.ensureHelperArrays(Molecule.cHelperRings);

        for (int atom=0; atom<mMol.getAtoms(); atom++) {
            if (mMol.getAtomZ(atom) != 0.0) {
                mZCoordinatesAvailable = true;
                break;
                }
            }

		mCanRank = new int[mMol.getAllAtoms()];
		mCanBaseValue = new long[mMol.getAtoms()];

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0)
				mCanBaseValue[atom] = 6;
			else
				mCanBaseValue[atom] = mMol.getAtomicNo(atom);
			mCanBaseValue[atom] <<= 8;
			mCanBaseValue[atom] += mMol.getAtomMass(atom);
			mCanBaseValue[atom] <<= 2;
			mCanBaseValue[atom] += mMol.getAtomPi(atom);
			mCanBaseValue[atom] <<= 3;
			mCanBaseValue[atom] += mMol.getConnAtoms(atom);
			mCanBaseValue[atom] <<= 4;
			if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0)
				mCanBaseValue[atom] += 8;
			else
				mCanBaseValue[atom] += (8 + mMol.getAtomCharge(atom));
			mCanBaseValue[atom] <<= 5;
			mCanBaseValue[atom] += mMol.getAtomRingSize(atom);
            mCanBaseValue[atom] <<= 4;
            mCanBaseValue[atom] += (mMol.getAtomAbnormalValence(atom)+1);
			mCanBaseValue[atom] <<= Molecule.cAtomQFNoOfBits;
			mCanBaseValue[atom] += mMol.getAtomQueryFeatures(atom);
			mCanBaseValue[atom] <<= 1;
			if (mMol.getAtomList(atom) != null)
				mCanBaseValue[atom]++;
			}
		int noOfRanks = canPerformFullRanking();


			// ############### begin tie breaking ##############

		while (noOfRanks < mMol.getAtoms()) {
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				mCanBaseValue[atom] = 2 * mCanRank[atom];

			int rank;
			for (rank=1; rank<=noOfRanks; rank++) {
				int rankAtoms = 0;
				for (int atom=0; atom<mMol.getAtoms(); atom++)
					if (mCanRank[atom] == rank)
						rankAtoms++;
				if (rankAtoms > 1)
					break;
				}
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mCanRank[atom] == rank) {
					mCanBaseValue[atom]++;   // select 1st atom of lowest rank for tie breaking
					break;
					}
				}
			noOfRanks = canPerformFullRanking();
			}
		}


	private int canPerformFullRanking() {
		// Does ranking based on sorted rank list of neighbour atoms and an attached
		// flag if neighbour is connected via a non-aromatic double bond.
		// This needs to be called at least once before tie-breaking in order to
		// distinguish atoms connected in a symmetrical way to to anti-aromatic rings
		// with different distribution of pi-bonds.
		int oldNoOfRanks,newNoOfRanks;

		newNoOfRanks = canConsolidate();
		do {
			oldNoOfRanks = newNoOfRanks;
			canCalcNextBaseValues();
			newNoOfRanks = canConsolidate();
			} while (oldNoOfRanks != newNoOfRanks);

		return newNoOfRanks;
		}


	private void canCalcNextBaseValues() {
		int	connRank[] = new int[ExtendedMolecule.cMaxConnAtoms];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
								// generate sorted list of ranks of neighbours
			for (int i=0; i<mMol.getConnAtoms(atom); i++) {
				int rank = 2 * mCanRank[mMol.getConnAtom(atom,i)];
				int connBond = mMol.getConnBond(atom,i);
				if (mMol.getBondOrder(connBond) == 2)
					if (!mMol.isAromaticBond(connBond))
						rank++;		// set a flag for non-aromatic double bond
				int j;
				for (j=0; j<i; j++)
					if (rank < connRank[j])
						break;
				for (int k=i; k>j; k--)
					connRank[k] = connRank[k-1];
				connRank[j] = rank;
				}
			
			mCanBaseValue[atom] = 0;
			for (int i=0; i<Math.min(6, mMol.getConnAtoms(atom)); i++) {
                mCanBaseValue[atom] <<= MAX_ATOM_BITS + 1;
				mCanBaseValue[atom] += connRank[i];
				}
			mCanBaseValue[atom] += (mCanRank[atom] << 55);	// should be 63-8
			}
		}


	private int canConsolidate() {
		int canRank = 0;	// all hydrogens have mCanRank[] = 0
		long lowest;

		while (true) {
			lowest = 0x7fffffffffffffffL;
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (lowest > mCanBaseValue[atom])
					lowest = mCanBaseValue[atom];

			if (lowest != 0x7fffffffffffffffL) {
				canRank++;
				for (int atom=0; atom<mMol.getAtoms(); atom++)
					if (mCanBaseValue[atom] == lowest) {
						mCanBaseValue[atom] = 0x7fffffffffffffffL;
						mCanRank[atom] = canRank;
						}
				continue;
				}
			break;
			}

		return canRank;
		}


	private void generateGraph() {
		if (mMol.getAtoms() == 0)
			return;
		if (mGraphGenerated)
			return;

		int startAtom = 0;
		for (int atom=1; atom<mMol.getAtoms(); atom++)
			if (mCanRank[atom] > mCanRank[startAtom])
				startAtom = atom;

		boolean atomHandled[] = new boolean[mMol.getAtoms()];
		boolean bondHandled[] = new boolean[mMol.getBonds()];
		int newAtomNo[] = new int[mMol.getAtoms()];
		mGraphAtom = new int[mMol.getAtoms()];
		mGraphFrom = new int[mMol.getAtoms()];
		mGraphBond = new int[mMol.getBonds()];
		mGraphAtom[0] = startAtom;
		atomHandled[startAtom] = true;

		int atomsWithoutParents = 1;	// the startatom has no parent
		int firstUnhandled = 0;
		int firstUnused = 1;
		int graphBonds = 0;
		while (firstUnhandled < mMol.getAtoms()) {
			if (firstUnhandled < firstUnused) {	// attach neighbours in rank order to unhandled
				while (true) {
					int highestRankingConnAtom = 0;
					int highestRankingConnBond = 0;
					int highestRank = -1;
					for (int i=0; i<mMol.getConnAtoms(mGraphAtom[firstUnhandled]); i++) {
						int connAtom = mMol.getConnAtom(mGraphAtom[firstUnhandled],i);
						if (!atomHandled[connAtom] && mCanRank[connAtom] > highestRank) {
							highestRankingConnAtom = connAtom;
							highestRankingConnBond = mMol.getConnBond(mGraphAtom[firstUnhandled],i);
							highestRank = mCanRank[connAtom];
							}
						}

					if (highestRank == -1)
						break;

					newAtomNo[highestRankingConnAtom] = firstUnused;
					mGraphFrom[firstUnused] = firstUnhandled;
					mGraphAtom[firstUnused++] = highestRankingConnAtom;
					mGraphBond[graphBonds++] = highestRankingConnBond;
					atomHandled[highestRankingConnAtom] = true;
					bondHandled[highestRankingConnBond] = true;
					}
				firstUnhandled++;
				}
			else {
				int highestRankingAtom = 0;
				int highestRank = -1;
				for (int atom=0; atom<mMol.getAtoms(); atom++) {
					if (!atomHandled[atom] && mCanRank[atom] > highestRank) {
						highestRankingAtom = atom;
						highestRank = mCanRank[atom];
						}
					}
				atomsWithoutParents++;
				newAtomNo[highestRankingAtom] = firstUnused;
				mGraphFrom[firstUnused] = -1;	// no parent atom in graph tree
				mGraphAtom[firstUnused++] = highestRankingAtom;
				atomHandled[highestRankingAtom] = true;
				}
			}

		mGraphClosure = new int[2 * (mMol.getBonds() - graphBonds)];
		mGraphRings = 0;
		while (true) {	// add ring closure bonds (those with lowest new atom numbers first)
			int lowAtomNo1 = mMol.getMaxAtoms();
			int lowAtomNo2 = mMol.getMaxAtoms();
			int lowBond = -1;
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				int loAtom,hiAtom;
				if (!bondHandled[bond]) {
					if (newAtomNo[mMol.getBondAtom(0,bond)]
					  < newAtomNo[mMol.getBondAtom(1,bond)]) {
						loAtom = newAtomNo[mMol.getBondAtom(0,bond)];
						hiAtom = newAtomNo[mMol.getBondAtom(1,bond)];
						}
					else {
						loAtom = newAtomNo[mMol.getBondAtom(1,bond)];
						hiAtom = newAtomNo[mMol.getBondAtom(0,bond)];
						}
					if (loAtom < lowAtomNo1
					 || (loAtom == lowAtomNo1 && hiAtom < lowAtomNo2)) {
						lowAtomNo1 = loAtom;
						lowAtomNo2 = hiAtom;
						lowBond = bond;
						}
					}
				}

			if (lowBond == -1)
				break;

			bondHandled[lowBond] = true;
			mGraphBond[graphBonds++] = lowBond;
			mGraphClosure[2*mGraphRings] = lowAtomNo1;
			mGraphClosure[2*mGraphRings+1] = lowAtomNo2;
			mGraphRings++;
			}

		mGraphGenerated = true;
		}


/*	protected Molecule getCanMolecule() {
		generateGraph();
		create new molecule and return it...
		return mol;
		}*/


	public String getIDCode() {
		generateGraph();
		idCodeCreate();
		return mIDCode.toString();
		}


	private void idCodeCreate() {
        encodeBitsStart();
        encodeBits(cIDCodeVersion3, 4);
        int nbits = Math.max(idGetNeededBits(mMol.getAtoms()),
                             idGetNeededBits(mMol.getBonds()));
        encodeBits(nbits, 4);

		if (nbits == 0) {
			encodeBits(mMol.isFragment() ? 1 : 0, 1);	// query fragment ?
			encodeBits(0, 1);
	        mIDCode = encodeBitsEnd();
	        return;
			}

		int nitrogens,oxygens,otherAtoms,chargedAtoms;
		nitrogens = oxygens = otherAtoms = chargedAtoms = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) == 0) {
				switch (mMol.getAtomicNo(atom)) {
				case 6:
					break;
				case 7:
					nitrogens++;
					break;
				case 8:
					oxygens++;
					break;
				default:
					otherAtoms++;
					break;
					}
				if (mMol.getAtomCharge(atom) != 0)
					chargedAtoms++;
				}
			}


		encodeBits(mMol.getAtoms(), nbits);
		encodeBits(mMol.getBonds(), nbits);
		encodeBits(nitrogens, nbits);
		encodeBits(oxygens, nbits);
		encodeBits(otherAtoms, nbits);
		encodeBits(chargedAtoms, nbits);

		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomicNo(mGraphAtom[atom]) == 7
			 && (mMol.getAtomQueryFeatures(mGraphAtom[atom]) & Molecule.cAtomQFAny) == 0)
				encodeBits(atom, nbits);
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomicNo(mGraphAtom[atom]) == 8
			 && (mMol.getAtomQueryFeatures(mGraphAtom[atom]) & Molecule.cAtomQFAny) == 0)
				encodeBits(atom, nbits);
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomicNo(mGraphAtom[atom]) != 6
			 && mMol.getAtomicNo(mGraphAtom[atom]) != 7
			 && mMol.getAtomicNo(mGraphAtom[atom]) != 8
			 && (mMol.getAtomQueryFeatures(mGraphAtom[atom]) & Molecule.cAtomQFAny) == 0) {
				encodeBits(atom, nbits);
				encodeBits(mMol.getAtomicNo(mGraphAtom[atom]), 8);
				}
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomCharge(mGraphAtom[atom]) != 0
			 && (mMol.getAtomQueryFeatures(mGraphAtom[atom]) & Molecule.cAtomQFAny) == 0) {
				encodeBits(atom, nbits);
				encodeBits(8 + mMol.getAtomCharge(mGraphAtom[atom]), 4);
				}

		int maxdif = 0;
		int base = 0;
		for (int atom=1; atom<mMol.getAtoms(); atom++) {
			int dif;
			if (mGraphFrom[atom] == -1) {
				dif = 0;
				}
			else {
				dif = 1 + mGraphFrom[atom] - base;
				base = mGraphFrom[atom];
				}
			if (maxdif < dif)
				maxdif = dif;
			}

		int dbits = idGetNeededBits(maxdif);
		encodeBits(dbits, 4);
		base = 0;
		for (int atom=1; atom<mMol.getAtoms(); atom++) {
			int dif;
			if (mGraphFrom[atom] == -1) {
				dif = 0;
				}
			else {
				dif = 1 + mGraphFrom[atom] - base;
				base = mGraphFrom[atom];
				}
			encodeBits(dif, dbits);
			}

		for (int i=0; i<2*mGraphRings; i++)
			encodeBits(mGraphClosure[i], nbits);

		for (int bond=0; bond<mMol.getBonds(); bond++) {
			int bondOrder = ((mMol.getBondQueryFeatures(mGraphBond[bond]) & Molecule.cBondQFBridge) != 0
					|| mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeMetalLigand) ?
					1 : (mMol.isDelocalizedBond(mGraphBond[bond])) ?
					0 : mMol.getBondOrder(mGraphBond[bond]);
			encodeBits(bondOrder, 2);
			}

		encodeBits(0, nbits);	// THCount
		encodeBits(0, 1);		// ChiralFlag
		encodeBits(0, nbits);	// EZCount

		encodeBits(mMol.isFragment() ? 1 : 0, 1);	// query fragment ?

		int count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomMass(mGraphAtom[atom]) != 0)
				count++;
		if (count != 0) {
			encodeBits(1, 1);	//	more data to come
			encodeBits(1, 4);	//	1 = datatype 'isotope'
			encodeBits(count, nbits);
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mMol.getAtomMass(mGraphAtom[atom]) != 0) {
					encodeBits(atom, nbits);
					encodeBits(mMol.getAtomMass(mGraphAtom[atom]), 8);
					}
				}
			}

        boolean isSecondFeatureBlock = false;

        if (mMol.isFragment()) {    // QueryFeatures and fragment specific properties
            addAtomQueryFeatures(0, false, nbits, Molecule.cAtomQFNoMoreNeighbours, 1, -1);

            addAtomQueryFeatures(3, false, nbits, Molecule.cAtomQFMoreNeighbours, 1, -1);

            addAtomQueryFeatures(4, false, nbits,
                                 Molecule.cAtomQFRingState,
                                 Molecule.cAtomQFRingStateBits,
                                 Molecule.cAtomQFRingStateShift);

            addAtomQueryFeatures(5, false, nbits,
                                 Molecule.cAtomQFAromState,
                                 Molecule.cAtomQFAromStateBits,
                                 Molecule.cAtomQFAromStateShift);

            addAtomQueryFeatures(6, false, nbits, Molecule.cAtomQFAny, 1, -1);

            addAtomQueryFeatures(7, false, nbits,
                                 Molecule.cAtomQFHydrogen,
                                 Molecule.cAtomQFHydrogenBits,
                                 Molecule.cAtomQFHydrogenShift);

            count = 0;
            for (int atom=0; atom<mMol.getAtoms(); atom++)
                if (mMol.getAtomList(mGraphAtom[atom]) != null)
                    count++;
            if (count > 0) {
                encodeBits(1, 1);   //  more data to come
                encodeBits(8, 4);   //  8 = datatype 'AtomList'
                encodeBits(count, nbits);
                for (int atom=0; atom<mMol.getAtoms(); atom++) {
                    int[] atomList = mMol.getAtomList(mGraphAtom[atom]);
                    if (atomList != null) {
                        encodeBits(atom, nbits);
                        encodeBits(atomList.length, 4);
                        for (int i=0; i<atomList.length; i++)
                            encodeBits(atomList[i], 8);
                        }
                    }
                }

            addBondQueryFeatures(9, false, nbits,
                                 Molecule.cBondQFRingState,
                                 Molecule.cBondQFRingStateBits,
                                 Molecule.cBondQFRingStateShift);

            addBondQueryFeatures(10, false, nbits,
                                 Molecule.cBondQFBondTypes,
                                 Molecule.cBondQFBondTypesBits,
                                 Molecule.cBondQFBondTypesShift);

            addAtomQueryFeatures(11, false, nbits, Molecule.cAtomQFMatchStereo, 1, -1);

            addBondQueryFeatures(12, false, nbits,
                                 Molecule.cBondQFBridge,
                                 Molecule.cBondQFBridgeBits,
                                 Molecule.cBondQFBridgeShift);

            addAtomQueryFeatures(13, false, nbits,
                                 Molecule.cAtomQFPiElectrons,
                                 Molecule.cAtomQFPiElectronBits,
                                 Molecule.cAtomQFPiElectronShift);

            addAtomQueryFeatures(14, false, nbits,
                                 Molecule.cAtomQFNeighbours,
                                 Molecule.cAtomQFNeighbourBits,
                                 Molecule.cAtomQFNeighbourShift);

            isSecondFeatureBlock |= addAtomQueryFeatures(16, isSecondFeatureBlock, nbits,
                                                         Molecule.cAtomQFSmallRingSize,
                                                         Molecule.cAtomQFSmallRingSizeBits,
                                                         Molecule.cAtomQFSmallRingSizeShift);
        	}

        count = 0;
        for (int atom=0; atom<mMol.getAtoms(); atom++)
            if (mMol.getAtomAbnormalValence(mGraphAtom[atom]) != -1)
                count++;
            if (count != 0) {
                isSecondFeatureBlock = ensureSecondFeatureBlock(isSecondFeatureBlock);
                encodeBits(1, 1);   //  more data to come
                encodeBits(1, 4);   //  (17-offset) 17 = datatype 'AtomAbnormalValence'
                encodeBits(count, nbits);
                for (int atom=0; atom<mMol.getAtoms(); atom++) {
                    if (mMol.getAtomAbnormalValence(mGraphAtom[atom]) != -1) {
                        encodeBits(atom, nbits);
                        encodeBits(mMol.getAtomAbnormalValence(mGraphAtom[atom]), 4);
                        }
                    }
                }

		if (mMol.isFragment()) {	// more QueryFeatures and fragment specific properties
			isSecondFeatureBlock |= addAtomQueryFeatures(19, isSecondFeatureBlock, nbits,
														 Molecule.cAtomQFCharge,
														 Molecule.cAtomQFChargeBits,
														 Molecule.cAtomQFChargeShift);

            isSecondFeatureBlock |= addBondQueryFeatures(20, isSecondFeatureBlock, nbits,
            											 Molecule.cBondQFRingSize,
            											 Molecule.cBondQFRingSizeBits,
            											 Molecule.cBondQFRingSizeShift);
			}

		count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
		    if (mMol.getAtomRadical(mGraphAtom[atom]) != 0)
		        count++;
        if (count != 0) {
            isSecondFeatureBlock = ensureSecondFeatureBlock(isSecondFeatureBlock);
            encodeBits(1, 1);   //  more data to come
            encodeBits(5, 4);   //  (21-offset) 21 = datatype 'AtomRadicalState'
            encodeBits(count, nbits);
            for (int atom=0; atom<mMol.getAtoms(); atom++) {
                if (mMol.getAtomRadical(mGraphAtom[atom]) != 0) {
                    encodeBits(atom, nbits);
                    encodeBits(mMol.getAtomRadical(mGraphAtom[atom]) >> Molecule.cAtomRadicalStateShift, 2);
                    }
                }
            }

		if (mMol.isFragment()) {	// more QueryFeatures and fragment specific properties
            isSecondFeatureBlock |= addAtomQueryFeatures(22, isSecondFeatureBlock, nbits, Molecule.cAtomQFFlatNitrogen, 1, -1);
            isSecondFeatureBlock |= addBondQueryFeatures(23, isSecondFeatureBlock, nbits, Molecule.cBondQFMatchStereo, 1, -1);
            isSecondFeatureBlock |= addBondQueryFeatures(24, isSecondFeatureBlock, nbits,
            											 Molecule.cBondQFAromState,
            											 Molecule.cBondQFAromStateBits,
            											 Molecule.cBondQFAromStateShift);
			}

		boolean[] isAromaticSPBond = getAromaticSPBonds();
		if (isAromaticSPBond != null) {
			count = 0;
			for (int bond=0; bond<mMol.getBonds(); bond++)
				if (isAromaticSPBond[mGraphBond[bond]])
					count++;

			isSecondFeatureBlock = ensureSecondFeatureBlock(isSecondFeatureBlock);
			encodeBits(1, 1);   //  more data to come
			encodeBits(10, 4);   //  (26-offset) 26 = datatype 'delocalized high order bond'
			encodeBits(count, nbits);
			for (int bond=0; bond<mMol.getBonds(); bond++)
				if (isAromaticSPBond[mGraphBond[bond]])
					encodeBits(bond, nbits);
			}

		if (mMol.isFragment())	// 27 = datatype 'part of an exclude-group'
			isSecondFeatureBlock |= addAtomQueryFeatures(27, isSecondFeatureBlock, nbits, Molecule.cAtomQFExcludeGroup, 1, -1);

		count = 0;
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeMetalLigand)
				count++;
		isSecondFeatureBlock = ensureSecondFeatureBlock(isSecondFeatureBlock);
		encodeBits(1, 1);   //  more data to come
		encodeBits(12, 4);   //  (28-offset) 28 = datatype 'coordinate bond'
		encodeBits(count, nbits);
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeMetalLigand)
				encodeBits(bond, nbits);

		encodeBits(0, 1);
        mIDCode = encodeBitsEnd();
		}


    private boolean ensureSecondFeatureBlock(boolean isSecondFeatureBlock) {
        if (!isSecondFeatureBlock) {
            encodeBits(1, 1);   //  more data to come
            encodeBits(15, 4);   //  15 = datatype 'start second query feature set'
            }
        return true;
        }


	private boolean addAtomQueryFeatures(int codeNo, boolean isSecondFeatureBlock, int nbits,
                                         long qfMask, int qfBits, int qfShift) {
        int count = 0;
        for (int atom=0; atom<mMol.getAtoms(); atom++)
            if ((mMol.getAtomQueryFeatures(mGraphAtom[atom]) & qfMask) != 0)
                count++;
        
        if (count == 0)
            return false;
        
        if (codeNo > 15) {
            ensureSecondFeatureBlock(isSecondFeatureBlock);
            codeNo -= 16;
            }
        
        encodeBits(1, 1);           //  more data to come
        encodeBits(codeNo, 4);      //  datatype
        encodeBits(count, nbits);
        for (int atom=0; atom<mMol.getAtoms(); atom++) {
        long feature = mMol.getAtomQueryFeatures(mGraphAtom[atom]) & qfMask;
        if (feature != 0) {
            encodeBits(atom, nbits);
            if (qfBits != 1)
                encodeBits(feature >> qfShift, qfBits);
            }
        }
    
    return true;
    }


	private boolean addBondQueryFeatures(int codeNo, boolean isSecondFeatureBlock, int nbits, int qfMask, int qfBits, int qfShift) {
	    int count = 0;
	    for (int bond=0; bond<mMol.getBonds(); bond++)
	        if ((mMol.getBondQueryFeatures(mGraphBond[bond]) & qfMask) != 0)
	            count++;

	    if (count == 0)
	        return false;

        if (codeNo > 15) {
            ensureSecondFeatureBlock(isSecondFeatureBlock);
            codeNo -= 16;
            }
        
	    encodeBits(1, 1);           //  more data to come
	    encodeBits(codeNo, 4);      //  datatype
	    encodeBits(count, nbits);
	    for (int bond=0; bond<mMol.getBonds(); bond++) {
	        int feature = mMol.getBondQueryFeatures(mGraphBond[bond]) & qfMask;
	        if (feature != 0) {
	            encodeBits(bond, nbits);
	            if (qfBits != 1)
	                encodeBits(feature >> qfShift, qfBits);
	            }
	        }

	    return true;
	    }


	private boolean[] getAromaticSPBonds() {
		boolean[] isAromaticSPBond = null;
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getAtomPi(mMol.getBondAtom(0, bond)) == 2
			 && mMol.getAtomPi(mMol.getBondAtom(1, bond)) == 2) {
				if (isAromaticSPBond == null)
					isAromaticSPBond = new boolean[mMol.getBonds()];
				isAromaticSPBond[bond] = true;
				}
			}
		return isAromaticSPBond;
		}


	public String getEncodedCoordinates() {
		if (mCoordinates == null) {
			generateGraph();
			encodeCoordinates(mZCoordinatesAvailable);
			}

		return mCoordinates;
		}


	private void encodeCoordinates(boolean keepAbsoluteValues) {
		if (mMol.getAtoms() == 0) {
            mCoordinates = "";
			return;
		    }

		int resolutionBits = mZCoordinatesAvailable ? 16 : 8;
		encodeBitsStart();
		mEncodingBuffer.append('!');
		encodeBits(mZCoordinatesAvailable ? 1 : 0, 1);
		encodeBits(keepAbsoluteValues ? 1 : 0, 1);
		encodeBits(resolutionBits/2, 4);	// resolution bits devided by 2

		double maxDelta = 0.0;
		for (int i=1; i<mMol.getAtoms(); i++) {
			int atom = mGraphAtom[i];
			int from = (mGraphFrom[i] == -1) ? -1 : mGraphAtom[mGraphFrom[i]];

			double deltaX = (from == -1) ?
							Math.abs(mMol.getAtomX(atom) - mMol.getAtomX(mGraphAtom[0])) / 8.0
						  : Math.abs(mMol.getAtomX(atom) - mMol.getAtomX(from));
			if (maxDelta < deltaX)
				maxDelta = deltaX;

			double deltaY = (from == -1) ?
							Math.abs(mMol.getAtomY(atom) - mMol.getAtomY(mGraphAtom[0])) / 8.0
						  : Math.abs(mMol.getAtomY(atom) - mMol.getAtomY(from));
			if (maxDelta < deltaY)
				maxDelta = deltaY;

			if (mZCoordinatesAvailable) {
				double deltaZ = (from == -1) ?
								Math.abs(mMol.getAtomZ(atom) - mMol.getAtomZ(mGraphAtom[0])) / 8.0
							  : Math.abs(mMol.getAtomZ(atom) - mMol.getAtomZ(from));
				if (maxDelta < deltaZ)
					maxDelta = deltaZ;
				}
			}

        if (maxDelta == 0.0) {
            mCoordinates = "";
            return;
            }

        int binCount = (1 << resolutionBits);
		double increment = maxDelta / (binCount / 2.0 - 1);
		double halfIncrement = increment / 2.0;

		for (int i=1; i<mMol.getAtoms(); i++) {
			int atom = mGraphAtom[i];
			int from = (mGraphFrom[i] == -1) ? -1 : mGraphAtom[mGraphFrom[i]];

			double deltaX = (from == -1) ?
							(mMol.getAtomX(atom) - mMol.getAtomX(mGraphAtom[0])) / 8.0
						   : mMol.getAtomX(atom) - mMol.getAtomX(from);

			double deltaY = (from == -1) ?
							(mMol.getAtomY(atom) - mMol.getAtomY(mGraphAtom[0])) / 8.0
						   : mMol.getAtomY(atom) - mMol.getAtomY(from);

			encodeBits((int)((maxDelta + deltaX + halfIncrement) / increment), resolutionBits);
			encodeBits((int)((maxDelta + deltaY + halfIncrement) / increment), resolutionBits);

			if (mZCoordinatesAvailable) {
				double deltaZ = (from == -1) ?
								(mMol.getAtomZ(atom) - mMol.getAtomZ(mGraphAtom[0])) / 8.0
							   : mMol.getAtomZ(atom) - mMol.getAtomZ(from);

				encodeBits((int)((maxDelta + deltaZ + halfIncrement) / increment), resolutionBits);
				}
			}

		if (keepAbsoluteValues) {
			encodeBits(encodeABVL(mMol.getAverageBondLength(true), binCount), resolutionBits);

			encodeBits(encodeShift(mMol.getAtomX(mGraphAtom[0]), binCount), resolutionBits);
			encodeBits(encodeShift(mMol.getAtomY(mGraphAtom[0]), binCount), resolutionBits);

            if (mZCoordinatesAvailable)
    			encodeBits(encodeShift(mMol.getAtomZ(mGraphAtom[0]), binCount), resolutionBits);
            }

        mCoordinates = encodeBitsEnd();
		}


	/**
	 * Encode a floating point value into an integer with precision proportional to the value itself.
	 * @param value
	 * @return
	 */
	private int encodeABVL(double value, int binCount) {
		return Math.min(binCount-1, Math.max(0, (int)(0.5 + Math.log10(value/0.1) / Math.log10(200/0.1) * (binCount-1))));
		}


	private int encodeShift(double value, int binCount) {
		int halfBinCount = binCount / 2;
		boolean isNegative =  (value < 0);
		value = Math.abs(value);
		float steepness = (float)binCount/100f;
		int intValue = (int)(0.5 + value * (halfBinCount-1) / (value + steepness));
		return isNegative ? halfBinCount + intValue : intValue;
		}


	private void encodeBitsStart() {
        mEncodingBuffer = new StringBuffer();
        mEncodingBitsAvail = 6;
        mEncodingTempData = 0;
        }


    private void encodeBits(long data, int bits) {
//System.out.println(bits+" bits:"+data+"  mode="+mode);
        while (bits != 0) {
            if (mEncodingBitsAvail == 0) {
                mEncodingBuffer.append((char)(mEncodingTempData + 64));
                mEncodingBitsAvail = 6;
                mEncodingTempData = 0;
                }
            mEncodingTempData <<= 1;
            mEncodingTempData |= (data & 1);
            data >>= 1;
            bits--;
            mEncodingBitsAvail--;
            }
        }


    private String encodeBitsEnd() {
        mEncodingTempData <<= mEncodingBitsAvail;
        mEncodingBuffer.append((char)(mEncodingTempData + 64));
        return mEncodingBuffer.toString();
        }


	private int idGetNeededBits(int no) {
		int bits = 0;
		while (no > 0) {
			no >>= 1;
			bits++;
			}
		return bits;
		}
	}
