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

package com.actelion.research.chem.conf;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;

public class TorsionDetail {
    public static final int SYMMETRY_C1C1_OR_C1D1 = 0;  // 0 -> 359 degrees
    public static final int SYMMETRY_C1D2 = 1;  // 0 -> 179 (equal to -180 -> -1)
    public static final int SYMMETRY_D1D1 = 2;  // 0, 1 -> 179 (equals -1 -> -179), 180
    public static final int SYMMETRY_D1D2_OR_D2D2 = 3;  // 0 -> 90 match 180 -> 90, 0 -> -90, -180 -> -90

    private static final int HALF_SYMMETRY_C1 = 0;  // three distinct terminal neighbors
    private static final int HALF_SYMMETRY_D1 = 1;  // e.g. single terminal neighbor or two equal sp3 neighbors
    private static final int HALF_SYMMETRY_D2 = 2;  // two equal sp2 neighbors
//  private static final int HALF_SYMMETRY_D3 = 3;  // for simplicity reasons this is covered by D1

    private static final int[][] SYMMETRY =
        { { SYMMETRY_C1C1_OR_C1D1, SYMMETRY_C1C1_OR_C1D1, SYMMETRY_C1D2 },
          { SYMMETRY_C1C1_OR_C1D1, SYMMETRY_D1D1, SYMMETRY_D1D2_OR_D2D2 },
          { SYMMETRY_C1D2, SYMMETRY_D1D2_OR_D2D2, SYMMETRY_D1D2_OR_D2D2 } };
    private static final int MAX_TRIPLE_BONDS = 8;
    private static final int FRAGMENT_ATOMS = 8+2*MAX_TRIPLE_BONDS;
	private static final int FRAGMENT_BONDS = 13+2*MAX_TRIPLE_BONDS;

	private StereoMolecule mFragment;
	private int[] mCentralAtom,mRearAtom,mRefAtom,mToFragmentAtom,mToMoleculeAtom;
	private int mAlkyneAtomCount;
	private String mID;

	/**
	 * This creates an empty torsion classification detail, which multiply can be used
	 * to classify the environment of a rotatable bond. The result is an identifier
	 * for this torsion situation that may serve to lookup probable torsion angles
	 * for this rotatable bond form a torsion library. Along with the identifier a
	 * 4-atom-chain is uniquely determined, to which torsion angles should refer to.
	 */
	public TorsionDetail() {
		mFragment = new StereoMolecule(FRAGMENT_ATOMS, FRAGMENT_BONDS);
		mToMoleculeAtom = new int[FRAGMENT_ATOMS];
		mCentralAtom = new int[2];
        mRearAtom = new int[2];
		mRefAtom = new int[2];
		}

	/**
	 * @return whether the previous classification created a valid result and ID
	 */
	public boolean isValid() {
		return mID != null;
		}

	/**
	 * Returns the torsion identifier that resulted from the most previous classification.
	 * @return a valid ID or null, if the classification failed
	 */
	public String getID() {
		return mID;
		}

	/**
	 * Returns one of the atoms of the torsion fragment's central rotatable bond.
	 * Two central and two terminal atoms define the strand that torsion angles refer to.
	 * If the rotatable bond is extended by (a) triple bond(s) then the central atom
	 * is one of the first non-sp atoms at the end of the linear atom strand.
	 * @param no 0 or 1; 0 refers to the higher ranking central atom
	 * @return
	 */
	public int getCentralAtom(int no) {
		return mCentralAtom[no];
		}

	/**
	 * Returns that neighbor atom of a central atom that lies in the rotatable bond axis.
	 * Usually this is the other atom of the rotatable bond. However, if the rotatable
	 * bond is extended by (a) triple bond(s) then getRearAtom(no) returns that sp-atom,
	 * which is attached to the atom returned by getCentralAtom(no).
	 * @param no 0 or 1; 0 refers to that atom being connected to central atom 0
	 * @return
	 */
	public int getRearAtom(int no) {
		return mRearAtom[no];
		}

	/**
	 * Returns the reference atom of one part of the torsion fragment. Two central
	 * and two terminal reference atoms define the strand that torsion angles refer to.
	 * @param no 0 or 1; 0 refers to the side with higher ranking central bond atom
	 * @return
	 */
	public int getReferenceAtom(int no) {
		return mRefAtom[no];
		}

	/**
	 * @return count of linear strand of sp-hybridized atoms between torsion atoms (usually 0)
	 */
	public int getAlkyneAtomCount() {
		return mAlkyneAtomCount;
		}

	/**
	 * Returns the direct surrounding of the previously classified rotatable bond.
	 * This fragment may contain query features to describe characteristics, which
	 * have an influence on torsion angles, but are beyond its atom and bond types
	 * (e.g. aromaticity and ring membership). Stereo centers and ESR attributes
	 * are normalized. Thus, atom parities may be inverted to the source molecule.
	 * @return
	 */
	public StereoMolecule getFragment() {
		return mFragment;
		}

	/**
     * Determines uniquely an identifying name for the rotatable bond and its vicinity.
     * If the bond is not a single bond, is aromatic, is a member of a <=5-membered ring,
     * if one of its atoms has 0,3 or more neighbours other than the other bond atom,
     * or if one of its atoms is a stereo center with unknown parity, then the classification fails.
     * The four atoms, that define the torsion angle, are determined and are available through
     * getReferenceAtom(0 or 1) and getCentralAtom(0 or 1). A fragment being unique for the
     * this particular torsion situation is determined.
     * If one of the bond's atoms is a stereo center, this fragment may be inverted for
     * normalization, which would be represented with a trailing '<'. A trailing '>' indicates
     * a non-inverted stereo center. If bond is part of a consecutive sp-sp atom chain, then
     * the classifying fragment covers all linear atoms plus two end atoms not being member
     * of the linear sp-atom strand.
     * If a TorsionDetail is passed, then this will be filled.
	 * @param mol
	 * @param bond the rotatable bond 
	 * @return true if a valid torsion identifier could be determined
	 */
	public boolean classify(StereoMolecule mol, int bond) {
		mFragment.clear();
		mID = null;

        mol.ensureHelperArrays(Molecule.cHelperSymmetrySimple);

        if (mol.getBondOrder(bond) != 1
         || mol.isAromaticBond(bond))
        	return false;

        if (mol.getAtomicNo(mol.getBondAtom(0, bond)) == 1
         || mol.getAtomicNo(mol.getBondAtom(1, bond)) == 1)
        	return false;

        boolean isSmallRingBond = mol.isSmallRingBond(bond);
        if (isSmallRingBond && mol.getBondRingSize(bond) < 6)  // 3- to 5-membered rings
            return false;

        // create fragment with all properties encoded as query features
        boolean[] atomMask = new boolean[mol.getAtoms()];

    	mAlkyneAtomCount = 0;
        for (int i=0; i<2; i++) {
        	mCentralAtom[i] = mol.getBondAtom(i, bond);
        	mRearAtom[i] = mol.getBondAtom(1-i, bond);

        	// walk along sp-chains to first sp2 or sp3 atom
        	while (mol.getAtomPi(mCentralAtom[i]) == 2
        		&& mol.getNonHydrogenNeighbourCount(mCentralAtom[i]) == 2
        		&& mol.getAtomicNo(mCentralAtom[i]) < 10) {
        		for (int j=0; j<mol.getConnAtoms(mCentralAtom[i]); j++) {
        			int connAtom = mol.getConnAtom(mCentralAtom[i], j);
        			if (connAtom != mRearAtom[i] && mol.getAtomicNo(connAtom) != 1) {
        				if (mol.getConnAtoms(connAtom) == 1
        				 || mAlkyneAtomCount == 2*MAX_TRIPLE_BONDS)
        					return false;

        				atomMask[mCentralAtom[i]] = true;
        				mRearAtom[i] = mCentralAtom[i];
        				mCentralAtom[i] = connAtom;
        				mAlkyneAtomCount++;
        				break;
        				}
        			}
        		}

	        int nonHNeighbours = mol.getNonHydrogenNeighbourCount(mCentralAtom[i]);
	        if (nonHNeighbours > 4
             || nonHNeighbours == 1)
        		return false;

        	atomMask[mCentralAtom[i]] = true;
        	}

        for (int i=0; i<2; i++) {
	        for (int j=0; j<mol.getConnAtoms(mCentralAtom[i]); j++) {
		        int connAtom = mol.getConnAtom(mCentralAtom[i], j);
		        if (mol.getAtomicNo(connAtom) != 1)
			        atomMask[connAtom] = true;
	            }
            }

        mToFragmentAtom = new int[mol.getAtoms()];
        mol.copyMoleculeByAtoms(mFragment, atomMask, true, mToFragmentAtom);
        for (int i=0; i<mToFragmentAtom.length; i++)
        	if (mToFragmentAtom[i] != -1)
        		mToMoleculeAtom[mToFragmentAtom[i]] = i;

		mFragment.setFragment(true);

        if (isSmallRingBond) {    // set ringBond flag for all 6-or 7-membered rings
            int bondInFragment = mFragment.getBond(mToFragmentAtom[mCentralAtom[0]], mToFragmentAtom[mCentralAtom[1]]);
            if (bondInFragment != -1) {	// -1 happens, if we have a triple bond in the ring
	            mFragment.setBondQueryFeature(bondInFragment, Molecule.cBondQFRing, true);
	
	            // flag terminal bonds that share the same ring with the central bond
	            RingCollection ringSet = mol.getRingSet();
	            for (int ringNo=0; ringNo<ringSet.getSize(); ringNo++) {
	            	if (ringSet.isBondMember(ringNo, bond)) {
	                    for (int i=0; i<2; i++) {
	                    	for (int j=0; j<mol.getConnAtoms(mCentralAtom[i]); j++) {
	                    		int connAtom = mol.getConnAtom(mCentralAtom[i], j);
	            	            if (connAtom != mRearAtom[i]) {
	                            	if (ringSet.isAtomMember(ringNo, connAtom) && mol.getAtomicNo(connAtom) != 1) {
			                            mFragment.setBondQueryFeature(mFragment.getBond(mToFragmentAtom[mCentralAtom[i]],
							                            mToFragmentAtom[connAtom]), Molecule.cBondQFRing, true);
	                            		break;
	                            		}
	            	            	}
	                    		}
	                    	}
	            		}
	            	}
            	}
            }

        for (int i=0; i<2; i++) {
            if (mol.isFlatNitrogen(mCentralAtom[i]))
                mFragment.setAtomQueryFeature(mToFragmentAtom[mCentralAtom[i]], Molecule.cAtomQFFlatNitrogen, true);

            boolean delocalizedBondFound = false;
            for (int j=0; j<mol.getConnAtoms(mCentralAtom[i]); j++) {
	            int connAtom = mol.getConnAtom(mCentralAtom[i], j);
	            if (connAtom != mRearAtom[i] && mol.getAtomicNo(connAtom) != 1) {
	            	int fragBond = mFragment.getBond(mToFragmentAtom[mCentralAtom[i]], mToFragmentAtom[connAtom]);
	            	if (mFragment.getBondType(fragBond) == Molecule.cBondTypeDelocalized) {
	            		delocalizedBondFound = true;
	            		}
	            	else if (mol.getAtomicNo(connAtom) == 6
	            		  && !mol.isAromaticAtom(mCentralAtom[i])) {	// if not encoded on central atom
	                    long feature = mol.isAromaticAtom(connAtom) ? Molecule.cAtomQFAromatic
	                                                           		: Molecule.cAtomQFNotAromatic;
	                    mFragment.setAtomQueryFeature(mToFragmentAtom[connAtom], feature, true);
	                    }
/*	                if (mol.isElectronegative(connAtom) && mFragment.getFreeValence(mToFragmentAtom[connAtom]) != 0) {
	                    int feature = (mol.getAllHydrogens(connAtom) != 0) ?
	                            Molecule.cAtomQFNot0Hydrogen
	                          : Molecule.cAtomQFNot1Hydrogen + Molecule.cAtomQFNot2Hydrogen + Molecule.cAtomQFNot3Hydrogen;
	                    mFragment.setAtomQueryFeature(mToFragmentAtom[connAtom], feature, true);
	                    }*/
		            int connBond = mol.getConnBond(mCentralAtom[i], j);
		            int ringSize = mol.getBondRingSize(connBond);
		            if (ringSize == 3 || ringSize == 4)
	                    mFragment.setBondQueryFeature(fragBond, ringSize << Molecule.cBondQFRingSizeShift, true);

		            // account for allyl-1,3-strain
		            if (mol.isAromaticBond(connBond)
		             || mol.getConnBondOrder(mCentralAtom[i], j) == 2) {
			            int nonHNeighbours = mol.getNonHydrogenNeighbourCount(connAtom);
		            	boolean hasZ = (nonHNeighbours == 3);
		            	if (!hasZ && nonHNeighbours == 2 && !mol.isRingAtom(connAtom))
		            		hasZ = (mol.getZNeighbour(mCentralAtom[1-i], connBond) != -1);
		            	if (hasZ) {
		            		// there is no query feature 'has-Z-neighbor', thus use 'has3neighbors'
		                    long feature = (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot3Neighbours);
		                    mFragment.setAtomQueryFeature(mToFragmentAtom[connAtom], feature, true);
		            		}
		            	else if (mol.isAromaticBond(connBond)) {
		            		// we show that there is no ortho substituent with 'has2neighbors'
				            long feature = (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot2Neighbours);
		                    mFragment.setAtomQueryFeature(mToFragmentAtom[connAtom], feature, true);
		            		}
		            	}

		            // account for gauche-pentane situations
		            if (mol.getConnBondOrder(mCentralAtom[i], j) == 1) {
		            	if (mol.getNonHydrogenNeighbourCount(connAtom) == 4) {
				            long feature = (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot4Neighbours);
		                    mFragment.setAtomQueryFeature(mToFragmentAtom[connAtom], feature, true);
			            	}
		            	else if (mol.getAtomicNo(connAtom) == 6) {
		                    mFragment.setAtomQueryFeature(mToFragmentAtom[connAtom], Molecule.cAtomQFNot4Neighbours, true);
		            		}
		            	}
	            	}
	            }

            if (!delocalizedBondFound) {
            	if (mol.isAromaticAtom(mCentralAtom[i]))
            		mFragment.setAtomQueryFeature(mToFragmentAtom[mCentralAtom[i]], Molecule.cAtomQFAromatic, true);
            	else
            		mFragment.setAtomQueryFeature(mToFragmentAtom[mCentralAtom[i]], Molecule.cAtomQFNotAromatic, true);
            	}
        	}

    	mFragment.ensureHelperArrays(Molecule.cHelperSymmetrySimple | Molecule.cHelperBitIncludeNitrogenParities);
//System.out.println("idcode:"+new Canonizer(mol).getIDCode()+" mol:"+mol.getAtoms()+" frag:"+mFragment.getAllAtoms());
    	// If we have a stereo center with unknown parity in the fragment(, which should only happen with 2D-molecules),
    	// and if due to symmetry that atom is not a stereo center in the molecule, then we choose parity1.
        for (int i=0; i<2; i++) {
        	int fragmentAtom = mToFragmentAtom[mCentralAtom[i]];
        	if (mFragment.getAtomParity(fragmentAtom) == Molecule.cAtomParityUnknown) {
        		if (mol.getAtomParity(mCentralAtom[i]) == Molecule.cAtomParityUnknown) {
        			return false;
        			}
        		else {
        			int preferredBond = mFragment.getAtomPreferredStereoBond(fragmentAtom);
        			mFragment.setBondType(preferredBond, Molecule.cBondTypeUp);
        			if (mFragment.getBondAtom(0, preferredBond) != fragmentAtom) {
        				mFragment.setBondAtom(1, preferredBond, mFragment.getBondAtom(0, preferredBond));
        				mFragment.setBondAtom(0, preferredBond, fragmentAtom);
        				}
                	mFragment.ensureHelperArrays(Molecule.cHelperSymmetrySimple | Molecule.cHelperBitIncludeNitrogenParities);
        			}
        		}
        	}

        int fatom1 = mToFragmentAtom[mCentralAtom[0]];
        int fatom2 = mToFragmentAtom[mCentralAtom[1]];
        int ratom1 = mToFragmentAtom[mRearAtom[0]];
        int ratom2 = mToFragmentAtom[mRearAtom[1]];

        int esrType1 = mFragment.getAtomESRType(fatom1);
        int esrType2 = mFragment.getAtomESRType(fatom2);

        if (mFragment.isAtomStereoCenter(fatom1) && mFragment.isAtomStereoCenter(fatom2)) {
            if ((esrType1 != Molecule.cESRTypeAbs || esrType2 != Molecule.cESRTypeAbs)
             && (esrType1 != esrType2 || mFragment.getAtomESRGroup(fatom1) != mFragment.getAtomESRGroup(fatom2)))
            	return false;
            }

        boolean esrTypeChanged = false;
        if (mFragment.isAtomStereoCenter(fatom1) && esrType1 != Molecule.cESRTypeAbs) {
        	mFragment.setAtomESR(fatom1, Molecule.cESRTypeAbs, -1);
        	esrTypeChanged = true;
        	}
        if (mFragment.isAtomStereoCenter(fatom2) && esrType2 != Molecule.cESRTypeAbs) {
        	mFragment.setAtomESR(fatom2, Molecule.cESRTypeAbs, -1);
        	esrTypeChanged = true;
        	}
        if (esrTypeChanged)
        	mFragment.ensureHelperArrays(Molecule.cHelperSymmetrySimple | Molecule.cHelperBitIncludeNitrogenParities);

        int rank1 = mFragment.getSymmetryRank(fatom1);
        int rank2 = mFragment.getSymmetryRank(fatom2);
        if (rank1 < rank2) {
            int temp = mCentralAtom[0];
            mCentralAtom[0] = mCentralAtom[1];
            mCentralAtom[1] = temp;
            temp = mRearAtom[0];
            mRearAtom[0] = mRearAtom[1];
            mRearAtom[1] = temp;
            temp = fatom1;
            fatom1 = fatom2;
            fatom2 = temp;
            temp = ratom1;
            ratom1 = ratom2;
            ratom2 = temp;
            }

        // normalize stereo configuration
        boolean isInverted = false;
        if (mFragment.isAtomStereoCenter(fatom1)
         || mFragment.isAtomStereoCenter(fatom2)) {
            if (mFragment.isAtomStereoCenter(fatom1)) {
            	isInverted = (mFragment.getAbsoluteAtomParity(fatom1) == Molecule.cAtomParity1);
                }
            else if (mFragment.isAtomStereoCenter(fatom2)) {
            	isInverted = (mFragment.getAbsoluteAtomParity(fatom2) == Molecule.cAtomParity1);
                }

            if (isInverted) {
                for (int atom=0; atom<mFragment.getAllAtoms(); atom++)
                	mFragment.setAtomX(atom, -mFragment.getAtomX(atom));

                mFragment.ensureHelperArrays(Molecule.cHelperSymmetrySimple | Molecule.cHelperBitIncludeNitrogenParities);
                }
            }

        int fconn1 = getReferenceNeighbor(fatom1, ratom1);
        int fconn2 = getReferenceNeighbor(fatom2, ratom2);

        mRefAtom[0] = (fconn1 == -1) ? -1 : mToMoleculeAtom[fconn1];
        mRefAtom[1] = (fconn2 == -1) ? -1 : mToMoleculeAtom[fconn2];

        String idcode = mFragment.getIDCode();
        if (idcode == null)
            return false;

        int halfSymmetry1 = getHalfSymmetry(fatom1, ratom1);
        int halfSymmetry2 = getHalfSymmetry(fatom2, ratom2);

        int symmetryType;
        if (halfSymmetry1 == HALF_SYMMETRY_C1
         && halfSymmetry2 == HALF_SYMMETRY_C1
         && (mFragment.getChirality() & ~Molecule.cChiralityIsomerCountMask) == Molecule.cChiralityMeso)
            symmetryType = SYMMETRY_D1D1;
        else
            symmetryType = SYMMETRY[halfSymmetry1][halfSymmetry2];

        String symmetryID = (symmetryType == SYMMETRY_C1C1_OR_C1D1) ?  (isInverted ? "<" : ">")
                          : (symmetryType == SYMMETRY_C1D2) ?          (isInverted ? "-" : "+")
                          : (symmetryType == SYMMETRY_D1D2_OR_D2D2) ?                 "=" : "";

        mID = idcode + symmetryID;
        return true;
		}

	/**
     * Tries to uniquely determine one of the terminal neighbor's of atom in mFragment
     * to serve as reference atom that any torsion angles are assigned to.
     * The logic is as follows: If we have one terminal neighbor, this is selected.
     * If there is a neighbor with a unique symmetry rank, then the one with the highest
     * rank is selected. If three neighbors share the same rank, then one of them is
     * selected. If we have two neighbors that share the same rank, then if atom is sp2
     * then one of them is selected. If atom is sp3 then we have to refer to a virtual
     * neighbor (-1) that is assumed to be at the third sp3 position.
     * @param atom one of the bond atoms of the rotatable bond
     * @param remoteBondAtom the remote atom of the rotatable bond
     * @return a unique neighbour atom or -1 if we have a virtual neighbor
     */
    private int getReferenceNeighbor(int atom, int remoteBondAtom) {
        int maxConn = -1;
        int maxRank = -1;
        int symConn = -1;
        boolean[] connHandled = new boolean[mFragment.getConnAtoms(atom)];
        for (int i=0; i<mFragment.getConnAtoms(atom); i++) {
            if (!connHandled[i]) {
                int conn = mFragment.getConnAtom(atom, i);
                if (conn != remoteBondAtom) {
                    int connRank = mFragment.getSymmetryRank(conn);
                    if (maxRank < connRank) {
                        boolean equalRankFound = false;
                        for (int j=i+1; j<mFragment.getConnAtoms(atom); j++) {
                            int candidate = mFragment.getConnAtom(atom, j);
                            if (candidate != remoteBondAtom
                             && mFragment.getSymmetryRank(candidate) == connRank) {
                                connHandled[j] = true;
                                if (equalRankFound)
                                	return conn;	// 3 equal ranking neighbors -> just take one of them
                                equalRankFound = true;
                                }
                            }
                        if (!equalRankFound) {
                            maxRank = connRank;
                            maxConn = conn;
                            }
                        else {
                            symConn = conn;
                        	}
                        }
                    }
                }
            }

        if (maxConn == -1) // no outer neighbor with unique rank found
        	if (isFlatAtom(atom))
        		return symConn;

        return maxConn;	// may be -1 if we have two outer neighbors with same rank and atom is SP3
        }

    /**
     * Checks whether atom is sp2 hybridized.
     * Amide nitrogens are also considered to be sp2.
     * @param atom
     * @return
     */
    private boolean isFlatAtom(int atom) {
	    return (mFragment.getAtomPi(atom) == 1 && mFragment.getAtomicNo(atom)<10)
			    || mFragment.isAromaticAtom(atom)
			    || mFragment.isFlatNitrogen(atom);
        }

    /**
     * Determines the symmetry of one end of the 4-atom sequence,
     * which may be one of:
     * HALF_SYMMETRY_C1: not symmetric due to stereo center or tetrahedral nitrogen.
     * HALF_SYMMETRY_D1: mirror plane due to one terminal atom only
     * or at least 2 symmetrical atoms at sp3 center.
     * HALF_SYMMETRY_D2: two symmetrical atoms at sp2 center.
     * The fragment's helper array level should be cHelperSymmetrySimple.
     * @param atom one of the bond atoms of the rotatable bond
     * @param rearAtom the remote atom of the rotatable bond
     * @return
     */
    private int getHalfSymmetry(int atom, int rearAtom) {
        if (mFragment.getConnAtoms(atom) == 2)
            return HALF_SYMMETRY_D1;

        int[] connAtom = getTerminalAtoms(atom, rearAtom);

        if (mFragment.getConnAtoms(atom) == 3) {
            if (mFragment.getSymmetryRank(connAtom[0]) == mFragment.getSymmetryRank(connAtom[1]))
                return isFlatAtom(atom) ? HALF_SYMMETRY_D2 : HALF_SYMMETRY_D1;
            else
                return isFlatAtom(atom) ? HALF_SYMMETRY_D1 : HALF_SYMMETRY_C1;
            }

        if (mFragment.getConnAtoms(atom) == 4) {
            // two equal ranks with additional neighbor that will serve as reference atom
            for (int i=0; i<connAtom.length; i++) {
                int rank = mFragment.getSymmetryRank(connAtom[i]);
                for (int j=i+1; j<connAtom.length; j++)
                    if (rank == mFragment.getSymmetryRank(connAtom[j]))
                            return HALF_SYMMETRY_D1;
                }
            }

        return HALF_SYMMETRY_C1;
        }

    private int[] getTerminalAtoms(int atom, int rearAtom) {
        int index = 0;
        int[] connAtom = new int[mFragment.getConnAtoms(atom)-1];
        for (int i=0; i<mFragment.getConnAtoms(atom); i++)
            if (mFragment.getConnAtom(atom, i) != rearAtom)
                connAtom[index++] = mFragment.getConnAtom(atom, i);
        return connAtom;
        }
    }
