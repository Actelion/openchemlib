/*
 * Copyright 2013-2020 Thomas Sander, openmolecules.org
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.conf.TorsionDetail;
import com.actelion.research.chem.conf.TorsionPrediction;

/**
 * A RotatableBond knows the two rigid fragments within a molecule
 * that are connected by this bond. It also knows about possible torsion
 * states with associated likelyhoods, which are taken from COD statistics
 * and modified to account for collisions due to bulky groups in the molecule.
 * It knows the smaller half of the molecule and rotate the smaller half to
 * a given torsion angle.
 */
public class RotatableBond {
	private static final short[] SIXTY_DEGREE_TORSION = { 0, 60, 120, 180, 240, 300};
	private static final short[] SIXTY_DEGREE_FREQUENCY = { 17, 17, 17, 17, 17, 17};
	private static final short[][] SIXTY_DEGREE_RANGE = { {-20,20},{40,80},{100,140},{160,200},{220,260},{280,320}};

	private RigidFragment mFragment1,mFragment2;
	private String mTorsionID;
	private int mRotationCenter,mBond,mFragmentNo1,mFragmentNo2;
	private boolean mBondAtomsInFragmentOrder;
	private float mBondRelevance;
	private short[] mDefaultTorsion;
	private short[] mDefaultFrequency;
	private short[][] mDefaultTorsionRange;
	private int[] mTorsionAtom,mRearAtom,mSmallerSideAtomList;

	public RotatableBond(StereoMolecule mol, int bond, int[] fragmentNo, int[] disconnectedFragmentNo,
	                     int disconnectedFragmentSize, RigidFragment[] fragment) {
		this(mol, bond, fragmentNo, disconnectedFragmentNo, disconnectedFragmentSize, fragment, false);
		}

	public RotatableBond(StereoMolecule mol, int bond, int[] fragmentNo, int[] disconnectedFragmentNo,
	                     int disconnectedFragmentSize, RigidFragment[] fragment, boolean use60degreeSteps) {
if (TorsionDB.getTorsionFrequencies("gGP`@dfyjidNcGI[WQCP`<")[0]==-1)
	System.out.println("ERROR");
		mBond = bond;
		mTorsionAtom = new int[4];
		mRearAtom = new int[2];
		TorsionDetail detail = new TorsionDetail();
		if (TorsionDB.getTorsionID(mol, bond, mTorsionAtom, detail) != null) {
			mRearAtom[0] = detail.getRearAtom(0);
			mRearAtom[1] = detail.getRearAtom(1);
			}
		else {
			predictAtomSequence(mol);
			}

		mFragmentNo1 = fragmentNo[mTorsionAtom[1]];
		mFragmentNo2 = fragmentNo[mTorsionAtom[2]];
		mFragment1 = fragment[mFragmentNo1];
		mFragment2 = fragment[mFragmentNo2];

		mBondAtomsInFragmentOrder = (fragmentNo[mol.getBondAtom(0, bond)] == mFragmentNo1);

if (TorsionDB.getTorsionFrequencies("gGP`@dfyjidNcGI[WQCP`<")[0]==-1)
	System.out.println("ERROR");


		if (use60degreeSteps) {
			mDefaultTorsion = SIXTY_DEGREE_TORSION;
			mDefaultFrequency = SIXTY_DEGREE_FREQUENCY;
			mDefaultTorsionRange = SIXTY_DEGREE_RANGE;
			}
		else {
			mTorsionID = detail.getID();
			mDefaultTorsion = TorsionDB.getTorsions(detail.getID());
			if (mDefaultTorsion == null) {
				TorsionPrediction prediction = new TorsionPrediction(mol, mTorsionAtom);
				mDefaultTorsion = prediction.getTorsions();
				mDefaultFrequency = prediction.getTorsionFrequencies();
				mDefaultTorsionRange = prediction.getTorsionRanges();
			} else {
				mDefaultFrequency = TorsionDB.getTorsionFrequencies(detail.getID());
				mDefaultTorsionRange = TorsionDB.getTorsionRanges(detail.getID());
			}
		}

		removeIllegalTorsions(mol);
		removeEquivalentTorsions(mol);

		findSmallerSideAtomList(mol, disconnectedFragmentNo, disconnectedFragmentSize);
		}

	public RigidFragment getFragment(int i) {
		return (i == 0) ? mFragment1 : mFragment2;
		}

	public int getFragmentNo(int i) {
		return (i == 0) ? mFragmentNo1 : mFragmentNo2;
		}

	public boolean areBondAtomsInFragmentOrder() {
		return mBondAtomsInFragmentOrder;
		}

	private void predictAtomSequence(StereoMolecule mol) {
        for (int i=0; i<2; i++) {
    		int centralAtom = mol.getBondAtom(i, mBond);
        	int rearAtom = mol.getBondAtom(1-i, mBond);

        	// walk along sp-chains to first sp2 or sp3 atom
        	while (mol.getAtomPi(centralAtom) == 2
        		&& mol.getConnAtoms(centralAtom) == 2
        		&& mol.getAtomicNo(centralAtom) < 10) {
        		for (int j=0; j<2; j++) {
        			int connAtom = mol.getConnAtom(centralAtom, j);
        			if (connAtom != rearAtom) {
        				rearAtom = centralAtom;
        				centralAtom = connAtom;
        				break;
        				}
        			}
        		}

        	mTorsionAtom[i+1] = centralAtom;
           	mRearAtom[i] = rearAtom;
        	}

    	// A TorsionPrediction does not distinguish hetero atoms from carbons a positions 0 and 3.
        // Therefore we can treat two sp2 neighbors as equivalent when predicting torsions.
        if (mol.getAtomPi(mTorsionAtom[1]) == 0 && mol.getConnAtoms(mTorsionAtom[1]) == 3) {
			mTorsionAtom[0] = -1;
        	}
        else {
			for (int i=0; i<mol.getConnAtoms(mTorsionAtom[1]); i++) {
				int connAtom = mol.getConnAtom(mTorsionAtom[1], i);
				if (connAtom != mTorsionAtom[2]) {
					mTorsionAtom[0] = connAtom;
					break;
					}
				}
        	}
        if (mol.getAtomPi(mTorsionAtom[2]) == 0 && mol.getConnAtoms(mTorsionAtom[2]) == 3) {
			mTorsionAtom[3] = -1;
        	}
        else {
			for (int i=0; i<mol.getConnAtoms(mTorsionAtom[2]); i++) {
				int connAtom = mol.getConnAtom(mTorsionAtom[2], i);
				if (connAtom != mTorsionAtom[1]) {
					mTorsionAtom[3] = connAtom;
					break;
					}
				}
        	}
		}

	private void findSmallerSideAtomList(StereoMolecule mol, int[] disconnectedFragmentNo, int disconnectedFragmentSize) {
		boolean[] isMember = new boolean[mol.getAllAtoms()];
		int memberCount = mol.getSubstituent(mRearAtom[0], mTorsionAtom[1], isMember, null, null);

		int alkyneAtoms = 0;	// if we have an extended linear sp-atom strain
		if (mRearAtom[0] != mTorsionAtom[2])
			alkyneAtoms = mol.getPathLength(mRearAtom[0], mTorsionAtom[2]);

		boolean invert = false;
		if (memberCount > disconnectedFragmentSize-alkyneAtoms-memberCount) {
			memberCount = disconnectedFragmentSize-alkyneAtoms-memberCount;
			invert = true;
			}

		// if invert, then flag all linear alkyne atoms to be avoided
		if (invert && alkyneAtoms != 0) {
			int spAtom = mRearAtom[0];
			int backAtom = mTorsionAtom[1];
        	while (mol.getAtomPi(spAtom) == 2
           		&& mol.getConnAtoms(spAtom) == 2
           		&& mol.getAtomicNo(spAtom) < 10) {
        		isMember[spAtom] = true;
           		for (int j=0; j<2; j++) {
           			int connAtom = mol.getConnAtom(spAtom, j);
           			if (connAtom != backAtom) {
           				backAtom = spAtom;
           				spAtom = connAtom;
           				break;
           				}
           			}
           		}
			}

		int memberNo = 0;
		int fragmentNo = disconnectedFragmentNo[mTorsionAtom[1]];
		mSmallerSideAtomList = new int[memberCount];
		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			if (disconnectedFragmentNo[atom] == fragmentNo && (isMember[atom] ^ invert))
				mSmallerSideAtomList[memberNo++] = atom;

		mBondRelevance = (float)((memberCount == 1) ? 1 : 2*memberCount) / mol.getAtoms();

		mRotationCenter = mTorsionAtom[invert ? 2 : 1];
		}

	/**
	 * @return bond index in molecule
	 */
	public int getBond() {
		return mBond;
		}

	public int[] getTorsionAtoms() {
		return mTorsionAtom;
		}

	public int getTorsionCount() {
		return mDefaultTorsion.length;
		}

	public int getRotationCenter() {
		return mRotationCenter;
		}

	public String getTorsionID() {
		return mTorsionID;
		}

	/**
	 * @return all default torsion angle in degrees
	 */
	public short[] getDefaultTorsions() {
		return mDefaultTorsion;
		}

	/**
	 * @return all default frequencies of the respective torsion anges
	 */
	public short[] getDefaultFrequencies() {
		return mDefaultFrequency;
	}

	/**
	 * @return all default frequencies of the respective torsion anges
	 */
	public short[][] getDefaultTorsionRanges() {
		return mDefaultTorsionRange;
	}

	/**
	 * @return the likelyhood of torsion i among all torsions of this bond
	 **
	public double getTorsionLikelyhood(int t) {
		return mLikelyhood[t];
		}   */

	/**
	 * @return atoms of the smaller half of the molecule excluding anchor atom
	 */
	public int[] getSmallerSideAtoms() {
		return mSmallerSideAtomList;
		}

	/**
	 * The relevance of a rotatable bond and its torsion angle for creating substantially different conformers
	 * depends on how close the bond is to the center of the molecule. Bond relevance values range from
	 * 1.0/atomCount (e.g. bond to methyl group) to 1.0 (bond dividing molecule into two equally large parts).
	 * @return relevance of this bond in the molecule to contribute to conformation change
	 */
	public float getRelevance() {
		return mBondRelevance;
		}

	/**
	 * If we have a BINAP stereo constraint, we have to remove colliding torsions
	 * @param mol
	 */
	private void removeIllegalTorsions(StereoMolecule mol) {
		if (mol.getBondOrder(mBond) == 1
		 && (mol.getBondParity(mBond) == Molecule.cBondParityEor1 || mol.getBondParity(mBond) == Molecule.cBondParityZor2)) {
			boolean inverse = false;
			for (int i=0; i<2; i++) {
				int conn = mTorsionAtom[3*i];
				int atom = mTorsionAtom[1+i];
				int rear = mTorsionAtom[2-i];
				for (int j=0; j<mol.getConnAtoms(atom); j++) {
					int other = mol.getConnAtom(atom, j);
					if (other != rear && other != conn) {
						if (other < conn)
							inverse = !inverse;
						break;
						}
					}
				}
			if (mol.getBondParity(mBond) == Molecule.cBondParityEor1)
				inverse = !inverse;

			// parityEor1 requires torsions values from 0...pi considering lowest atom indexes for mTorsionAtom[0 and 3]
			int count = 0;
			int frequencySum = 0;
			for (int i=0; i<mDefaultTorsion.length; i++) {
				if (mDefaultTorsion[i]<180 ^ inverse) {
					frequencySum += mDefaultFrequency[i];
					count++;
					}
				}

			if (count < mDefaultTorsion.length) {
				short[] newTorsion = new short[count];
				short[] newFrequency = new short[count];
				short[][] newRange = new short[count][];
				count = 0;
				for (int i=0; i<mDefaultTorsion.length; i++) {
					if (mDefaultTorsion[i]<180 ^ inverse) {
						newTorsion[count] = mDefaultTorsion[i];
						newFrequency[count] = (short)(mDefaultFrequency[i] * 100 / frequencySum);
						newRange[count] = mDefaultTorsionRange[i];
						count++;
						}
					}
				mDefaultTorsion = newTorsion;
				mDefaultFrequency = newFrequency;
				mDefaultTorsionRange = newRange;
				}
			}
		}

	/**
	 * For terminal fragments with D2 or D3 symmetry we may remove parts
	 * of the torsion list, because we would get equivalent conformers.
	 */
	private void removeEquivalentTorsions(StereoMolecule mol) {
		final int[][] SYMMETRY_COUNT = {{1,2,3},{2,4,12},{3,12,6}};
		int symmetryCount1 = (mFragment1.getConnectionPointCount() != 1) ? 1
				: countSymmetricalTerminalNeighbors(mol, mTorsionAtom[1], mRearAtom[0]);
		int symmetryCount2 = (mFragment2.getConnectionPointCount() != 1) ? 1
				: countSymmetricalTerminalNeighbors(mol, mTorsionAtom[2], mRearAtom[1]);

		if (symmetryCount1 == 1 && symmetryCount2 == 1)
			return;

		// we assume that we have only 1,2,3 as individual symmetryCounts
		int maxAngle = 360 / Math.max(symmetryCount1, symmetryCount2);

		int count = 0;
		int frequencySum = 0;
		for (int i=0; i<mDefaultTorsion.length && mDefaultTorsion[i] < maxAngle; i++) {
			frequencySum += mDefaultFrequency[i];
			count++;
			}

		if (count == 0)	// should not happen
			return;

		short[] newTorsion = new short[count];
		short[] newFrequency = new short[count];
		short[][] newRange = new short[count][];
        for (int i=0; i<count; i++) {
        	newTorsion[i] = mDefaultTorsion[i];
        	newFrequency[i] = (short)(mDefaultFrequency[i] * 100 / frequencySum);
        	newRange[i] = mDefaultTorsionRange[i];
        	}

		mDefaultTorsion = newTorsion;
		mDefaultFrequency = newFrequency;
		mDefaultTorsionRange = newRange;
		}

	private void sortTorsionsByFrequency() {
		int count = mDefaultFrequency.length;

		short[] newTorsion = new short[count];
		short[] newFrequency = new short[count];
		short[][] newRange = new short[count][];
		boolean[] used = new boolean[count];

		for (int i=0; i<count; i++) {
			int maxFrequency = -1;
			int maxIndex = -1;
			for (int j=0; j<count; j++) {
				if (!used[j] && maxFrequency < mDefaultFrequency[j]) {
					maxFrequency = mDefaultFrequency[j];
					maxIndex = j;
					used[j] = true;
				}
			}

			newTorsion[i] = mDefaultTorsion[maxIndex];
			newFrequency[i] = mDefaultFrequency[maxIndex];
			newRange[i] = mDefaultTorsionRange[maxIndex];
		}

		mDefaultTorsion = newTorsion;
		mDefaultFrequency = newFrequency;
		mDefaultTorsionRange = newRange;
	}

	/**
     * Checks whether all neighbors of atom (not considering rearAtom)
     * have the same symmetry rank. Implicit hydrogens are considered.
     * For sp2 atoms this requires 2 equally ranked neighbors, for sp3
     * atoms there must be three.
     * @param mol
     * @param atom
     * @param rearAtom connected to atom and not considered
     * @return 1,2, or 3
     */
    private int countSymmetricalTerminalNeighbors(StereoMolecule mol, int atom, int rearAtom) {
		if (mol.getAtomPi(atom) == 2)
			return 1;
		if ((mol.getAtomPi(atom) == 1 || mol.isFlatNitrogen(atom)) && mol.getConnAtoms(atom) != 3)
			return 1;
		if (mol.getAtomPi(atom) == 0 && mol.getConnAtoms(atom) != 4)
			return 1;

		int rank = -2;
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int connAtom = mol.getConnAtom(atom, i);
			if (connAtom != rearAtom) {
				if (rank == -2)
					rank = mol.getSymmetryRank(connAtom);
				else if (rank != mol.getSymmetryRank(connAtom))
					return 1;
				}
			}

		return mol.getConnAtoms(atom)-1;
    	}
	}