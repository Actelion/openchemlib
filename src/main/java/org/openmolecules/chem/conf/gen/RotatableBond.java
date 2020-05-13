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

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.*;

import java.util.Random;

/**
 * A RotatableBond knows the two rigid fragments within a molecule
 * that are connected by this bond. It also knows about possible torsion
 * states with associated likelyhoods, which are taken from COD statistics
 * and modified to account for collisions due to bulky groups in the molecule.
 * It knows the smaller half of the molecule and rotate the smaller half to
 * a given torsion angle.
 */
public class RotatableBond {
	private static final double ANGLE_TOLERANCE = 0.001f;	// limit for considering bonds as parallel

	// For every optimum torsion we check for collisions and in case we try, whether left or right
	// range ends are substantially better.
	private static final double ACCEPTABLE_CENTER_STRAIN = 0.05f;	// if we have less strain than this we don't check whether range edges improve strain
	private static final double NECESSARY_EDGE_STRAIN_IMPROVEMENT = 0.02f;	// necessary strain improvement to use edge torsion rather than center torsion
	private static final double MAXIMUM_CENTER_STRAIN = 0.2f;		// if center torsion is above this limit we refuse that torsion

	// If no acceptable torsion remains after checking and attempting to use range edge values,
	// we try to modify the least bad torsion stepwise until it is acceptable.
	private static final int ESCAPE_ANGLE = 8;  // degrees to rotate the rotatable bond to escape collisions
	private static final int ESCAPE_STEPS = 4;	// how often we apply this rotation trying to solve the collision
	private static final double MIN_ESCAPE_GAIN_PER_STEP = 0.05;

	private static final short[] SIXTY_DEGREE_TORSION = { 0, 60, 120, 180, 240, 300};
	private static final short[] SIXTY_DEGREE_FREQUENCY = { 17, 17, 17, 17, 17, 17};
	private static final short[][] SIXTY_DEGREE_RANGE = { {-20,20},{40,80},{100,140},{160,200},{220,260},{280,320}};

	private RigidFragment mFragment1,mFragment2;
	private Random mRandom;
	private int mRotationCenter,mBond,mFragmentNo1,mFragmentNo2;
	private boolean mBondAtomsInFragmentOrder;
	private float mBondRelevance;
	private short[] mTorsion;
	private short[] mFrequency;
	private short[][] mTorsionRange;
	private double[] mLikelyhood; // considering directly connected rigid fragments (frequency and collision strain)
	private int[] mTorsionAtom,mRearAtom,mSmallerSideAtomList;

	public RotatableBond(StereoMolecule mol, int bond, int[] fragmentNo, int[] disconnectedFragmentNo,
	                     int disconnectedFragmentSize, RigidFragment[] fragment, Random random) {
		this(mol, bond, fragmentNo, disconnectedFragmentNo, disconnectedFragmentSize, fragment, random, false);
		}

	public RotatableBond(StereoMolecule mol, int bond, int[] fragmentNo, int[] disconnectedFragmentNo,
	                     int disconnectedFragmentSize, RigidFragment[] fragment, Random random, boolean use60degreeSteps) {
		mBond = bond;
		mRandom = random;
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

		if (use60degreeSteps) {
			mTorsion = SIXTY_DEGREE_TORSION;
			mFrequency = SIXTY_DEGREE_FREQUENCY;
			mTorsionRange = SIXTY_DEGREE_RANGE;
			}
		else {
			mTorsion = TorsionDB.getTorsions(detail.getID());
			if (mTorsion == null) {
				TorsionPrediction prediction = new TorsionPrediction(mol, mTorsionAtom);
				mTorsion = prediction.getTorsions();
				mFrequency = prediction.getTorsionFrequencies();
				mTorsionRange = prediction.getTorsionRanges();
			} else {
				mFrequency = TorsionDB.getTorsionFrequencies(detail.getID());
				mTorsionRange = TorsionDB.getTorsionRanges(detail.getID());
			}
		}

		removeIllegalTorsions(mol);
		removeEquivalentTorsions(mol);
		mLikelyhood = new double[mTorsion.length];

		findSmallerSideAtomList(mol, disconnectedFragmentNo, disconnectedFragmentSize);
		}

	public RigidFragment getFragment(int i) {
		return (i == 0) ? mFragment1 : mFragment2;
		}

	public int getFragmentNo(int i) {
		return (i == 0) ? mFragmentNo1 : mFragmentNo2;
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

	public int getTorsionCount() {
		return mTorsion.length;
		}

	/**
	 * Calculates a random torsion index giving torsions with higher likelyhoods
	 * (i.e. frequencies and collision strains) a higher chance to be selected.
	 * With a progress value of 0.0 selection likelyhoods are proportional to
	 * the torsion frequencies. With increasing progress value the higher frequent
	 * torsions are less and less preferred until 1.0 without any preference.
	 * @param random
	 * @param progress 0...1 
	 */
	public int getLikelyRandomTorsionIndex(double random, double progress) {
		double sum = 0;
		for (int t=0; t<mTorsion.length; t++) {
			double contribution = (1f-progress)*mLikelyhood[t] + progress/mTorsion.length;
			sum += contribution;
			if (random <= sum)
				return t;
			}
		return mTorsion.length-1;  // should never reach this
		}

	/**
	 * @return the i'th torsion angle in degrees
	 */
	public short getTorsion(int t) {
		return mTorsion[t];
		}

	/**
	 * @return the likelyhood of torsion i among all torsions of this bond
	 */
	public double getTorsionLikelyhood(int t) {
		return mLikelyhood[t];
		}

	/**
	 * @return count of atoms of the smaller half of the molecule excluding anchor atom
	 */
	public int getSmallerSideAtomCount() {
		return mSmallerSideAtomList.length;
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
	 * Checks both rigid fragments that are connected by this bond, whether they have
	 * been attached to the conformer yet, i.e. whether their local coordinates have been
	 * copied to conformer and transformed according to the connecting torsion.
	 * If one was already attached, then the other's coordinates are transformed according
	 * to the torsion and copied to the conformer. A likelihood is calculated for every torsion
	 * value based on its frequency and the collision strain of the two fragments' atoms.
	 * If both fragments were not attached yet, then the larger one's coordinates are
	 * copied and the smaller one's coordinates are translated and then copied.
	 * Unlikely torsions, i.e. where collisions strain outweighs frequency, are removed from torsion list.
	 */
	public void connectFragments(Conformer conformer, boolean[] isAttached, int[] fragmentPermutation) {
		if (!isAttached[mFragmentNo1] && !isAttached[mFragmentNo2]) {
			RigidFragment largerFragment = (mFragment1.getCoreSize() > mFragment2.getCoreSize()) ? mFragment1 : mFragment2;
            int largerFragmentNo = (mFragment1.getCoreSize() > mFragment2.getCoreSize()) ? mFragmentNo1 : mFragmentNo2;
			isAttached[largerFragmentNo] = true;
			int fragmentConformer = (fragmentPermutation == null) ? 0 : fragmentPermutation[largerFragmentNo];
			for (int i=0; i<largerFragment.getExtendedSize(); i++) {
				int atom = largerFragment.extendedToOriginalAtom(i);
				conformer.setCoordinates(atom, largerFragment.getExtendedCoordinates(fragmentConformer, i));
				}
			}

		assert(isAttached[mFragmentNo1] ^ isAttached[mFragmentNo2]);

		int rootAtom,rearAtom,fragmentNo,bondAtomIndex;
		RigidFragment fragment;
		if (isAttached[mFragmentNo1]) {
            fragmentNo = mFragmentNo2;
			fragment = mFragment2;
			bondAtomIndex = mBondAtomsInFragmentOrder ? 1 : 0;
			}
		else {
            fragmentNo = mFragmentNo1;
			fragment = mFragment1;
			bondAtomIndex = mBondAtomsInFragmentOrder ? 0 : 1;
			}

		rootAtom = conformer.getMolecule().getBondAtom(bondAtomIndex, mBond);
		rearAtom = conformer.getMolecule().getBondAtom(1-bondAtomIndex, mBond);

		int fragmentConformer = (fragmentPermutation == null) ? 0 : fragmentPermutation[fragmentNo];

		int fRootAtom = fragment.originalToExtendedAtom(rootAtom);
		int fRearAtom = fragment.originalToExtendedAtom(rearAtom);

		Coordinates froot = fragment.getExtendedCoordinates(fragmentConformer, fRootAtom);
		Coordinates root = conformer.getCoordinates(rootAtom);
		Coordinates fuv = froot.subC(fragment.getExtendedCoordinates(fragmentConformer, fRearAtom)).unit();
		Coordinates uv = root.subC(conformer.getCoordinates(rearAtom)).unit();
		double alpha = fuv.getAngle(uv);

		Coordinates[] fcoords = new Coordinates[fragment.getExtendedSize()];

		if (alpha < ANGLE_TOLERANCE) { // special case both axes parallel: we only need to translate
			for (int i=0; i<fragment.getExtendedSize(); i++) {
				int atom = fragment.extendedToOriginalAtom(i);
				if (atom != rootAtom && atom != rearAtom)
					fcoords[i] = (alpha > Math.PI / 2) ?
							froot.subC(fragment.getExtendedCoordinates(fragmentConformer, i))
						  : fragment.getExtendedCoordinates(fragmentConformer, i).subC(froot);
				}
			}
		else {
			Coordinates rotationAxis;
			if (alpha < Math.PI - ANGLE_TOLERANCE) {    // normal case, just rotate around the plane orthogonal
				rotationAxis = uv.cross(fuv);
				}
			else {    // special case both axes anti-parallel: for cross-product we need any vector being different from uv
				if (Math.abs(uv.x) >= Math.abs(uv.y) && Math.abs(uv.x) >= Math.abs(uv.z))
					rotationAxis = new Coordinates(-(uv.y + uv.z) / uv.x, 1.0, 1.0);
				else if (Math.abs(uv.y) >= Math.abs(uv.x) && Math.abs(uv.y) >= Math.abs(uv.z))
					rotationAxis = new Coordinates(1.0, -(uv.x + uv.z) / uv.y, 1.0);
				else
					rotationAxis = new Coordinates(1.0, 1.0, -(uv.x + uv.y) / uv.z);
				}

			double[][] m = getRotationMatrix(rotationAxis.unit(), alpha);
	
			for (int i=0; i<fragment.getExtendedSize(); i++) {
				int atom = fragment.extendedToOriginalAtom(i);
				if (atom != rootAtom && atom != rearAtom)
					fcoords[i] = fragment.getExtendedCoordinates(fragmentConformer, i).subC(froot).rotate(m);
				}
			}
		isAttached[fragmentNo] = true;

		// we need to restore valid coordinates for rootAtom and neighbors
		// to correctly calculate torsion from atom sequence
		for (int i=0; i<conformer.getMolecule().getConnAtoms(rootAtom); i++) {
			int connAtom = conformer.getMolecule().getConnAtom(rootAtom, i);
			if (connAtom != rearAtom) {
				int fAtom = fragment.originalToExtendedAtom(connAtom);
				conformer.setCoordinatesReplace(connAtom, fcoords[fAtom].addC(root));
				}
			}

		double startTorsion = TorsionDB.calculateTorsionExtended(conformer, mTorsionAtom);

		boolean multipleFragments = (mFragment1.getConformerCount() > 1)
								 || (mFragment2.getConformerCount() > 1);

//System.out.print("connectFragments() original torsions:"); for (int t=0; t<mTorsion.length; t++) System.out.print(mTorsion[t]+" "); System.out.println();
		short currentTorsion = -1;
		int bestTorsionIndex = 0;
		int bestTorsionEdgeUsed = 0;	// 0:peakCenter, 1:leftEdge, 2:rightEdge
		double bestTorsionStrain = 0;
		for (int t=0; t<mTorsion.length; t++) {
			currentTorsion = mTorsion[t];
			double[][] m = getRotationMatrix(uv, Math.PI * currentTorsion / 180 - startTorsion);
			for (int i=0; i<fragment.getExtendedSize(); i++) {
				int atom = fragment.extendedToOriginalAtom(i);
				if (atom != rootAtom && atom != rearAtom) {
					conformer.setCoordinates(atom, fcoords[i]);
					conformer.getCoordinates(atom).rotate(m).add(root);
					}
				}

			double centerStrain = calculateCollisionStrain(conformer);
			double usedStrain = centerStrain;	// default
			int torsionEdgeUsed = 0;	// default
			// if the strain is above a certain limit, we investigate whether we should use the
			// limits of the torsion range rather than the central torsion value, which has the
			// highest frequency.
			// If we use torsion range limits, then we need to consider half of the frequency
			// for each of the lower and higher limits.
			if (multipleFragments) {
				mLikelyhood[t] = mFrequency[t];
				}
			else if (centerStrain < ACCEPTABLE_CENTER_STRAIN) {
				double relativeStrain = centerStrain / MAXIMUM_CENTER_STRAIN;
				mLikelyhood[t] = mFrequency[t] * (1.0 - relativeStrain * relativeStrain);
				}
			else {
				boolean foundAlternative = false;
				boolean isFirstAlternative = true;
				for (int r=0; r<2; r++) {
					currentTorsion = mTorsionRange[t][r];
					m = getRotationMatrix(uv, Math.PI * currentTorsion / 180 - startTorsion);
					for (int i=0; i<fragment.getExtendedSize(); i++) {
						int atom = fragment.extendedToOriginalAtom(i);
						if (atom != rootAtom && atom != rearAtom) {
							conformer.setCoordinates(atom, fcoords[i]);
							conformer.getCoordinates(atom).rotate(m).add(root);
							}
						}

					double rangeStrain = calculateCollisionStrain(conformer);
					if (centerStrain - rangeStrain > NECESSARY_EDGE_STRAIN_IMPROVEMENT) {
						if (isFirstAlternative) {
							mTorsion[t] = mTorsionRange[t][r];
							mFrequency[t] = (short)((mFrequency[t]+1) / 2);
							double relativeStrain = rangeStrain / MAXIMUM_CENTER_STRAIN;
							mLikelyhood[t] = mFrequency[t] * (1f - relativeStrain * relativeStrain);
							usedStrain = rangeStrain;
							torsionEdgeUsed = r+1;
							isFirstAlternative = false;
							}
						else {
							double relativeStrain = rangeStrain / MAXIMUM_CENTER_STRAIN;
							insertTorsion(t, r, relativeStrain);
							if (mLikelyhood[t+1] > mLikelyhood[t]) {
								usedStrain = rangeStrain;
								torsionEdgeUsed = r+1;
								}
							t++;
							}
						foundAlternative = true;
						}
					}
				if (!foundAlternative /* && strain < MAXIMUM_CENTER_STRAIN really bad ones should get a negative likelyhood */) {
					double relativeStrain = centerStrain / MAXIMUM_CENTER_STRAIN;
					mLikelyhood[t] = mFrequency[t] * (1f - relativeStrain * relativeStrain);
					}
				}
			if (mLikelyhood[bestTorsionIndex] < mLikelyhood[t]) {
				bestTorsionIndex = t;
				bestTorsionEdgeUsed = torsionEdgeUsed;
				bestTorsionStrain = usedStrain;	// this is the strain with the highest likelyhood (not necessarily the lowest strain)
				}
			}

		double totalLikelyhood = 0f;
		for (int t=0; t<mTorsion.length; t++)
			if (mLikelyhood[t] > 0f)
				totalLikelyhood += mLikelyhood[t];

		// make sure, we have at least one torsion with positive likelyhood, because only those are considered later
		if (mLikelyhood[bestTorsionIndex] <= 0f) {
			mLikelyhood[bestTorsionIndex] = 1.0f;
			int angle = bestTorsionEdgeUsed == 1 ? -ESCAPE_ANGLE
					  : bestTorsionEdgeUsed == 2 ? ESCAPE_ANGLE
					  : (mRandom.nextDouble() < 0.5) ? -ESCAPE_ANGLE : ESCAPE_ANGLE;
			for (int step=1; step<=ESCAPE_STEPS; step++) {
				currentTorsion = (short)(mTorsion[bestTorsionIndex]+angle*step);
				double[][] m = getRotationMatrix(uv, Math.PI * currentTorsion / 180 - startTorsion);
				for (int i=0; i<fragment.getExtendedSize(); i++) {
					int atom = fragment.extendedToOriginalAtom(i);
					if (atom != rootAtom && atom != rearAtom) {
						conformer.setCoordinates(atom, fcoords[i]);
						conformer.getCoordinates(atom).rotate(m).add(root);
						}
					}

				double escapeStrain = calculateCollisionStrain(conformer);
				if (bestTorsionStrain - escapeStrain < MIN_ESCAPE_GAIN_PER_STEP)
					break;

				mTorsion[bestTorsionIndex] = currentTorsion;
				}
			}
		else {
			for (int t = 0; t < mTorsion.length; t++)
				mLikelyhood[t] /= totalLikelyhood;
			}

		bestTorsionIndex = removeUnlikelyTorsions(bestTorsionIndex);

//System.out.print("connectFragments() applied torsions:"); for (int t=0; t<mTorsion.length; t++) System.out.print(mTorsion[t]+" "); System.out.println();
//System.out.print("connectFragments() applied likelyhoods:"); for (int t=0; t<mTorsion.length; t++) System.out.print(mLikelyhood[t]+" "); System.out.println();
//System.out.println();

		if (currentTorsion != mTorsion[bestTorsionIndex]) {
			double[][] m = getRotationMatrix(uv, Math.PI * mTorsion[bestTorsionIndex] / 180 - startTorsion);
			for (int i=0; i<fragment.getExtendedSize(); i++) {
				int atom = fragment.extendedToOriginalAtom(i);
				if (atom != rootAtom && atom != rearAtom) {
					conformer.setCoordinates(atom, fcoords[i]);
					conformer.getCoordinates(atom).rotate(m).add(root);
					}
				}
			}
		conformer.setBondTorsion(mBond, mTorsion[bestTorsionIndex]);
		}

	private double calculateCollisionStrain(Conformer conformer) {
		double panalty = 0f;
		int bondAtom1 = conformer.getMolecule().getBondAtom(0, mBond);
		int bondAtom2 = conformer.getMolecule().getBondAtom(1, mBond);
		for (int i=0; i<mFragment1.getCoreSize(); i++) {
			int atom1 = mFragment1.coreToOriginalAtom(i);
			if (atom1 != bondAtom1 && atom1 != bondAtom2) {
				double vdwr1 = ConformerGenerator.getToleratedVDWRadius(conformer.getMolecule().getAtomicNo(atom1));
				for (int j=0; j<mFragment2.getCoreSize(); j++) {
					int atom2 = mFragment2.coreToOriginalAtom(j);
					if (atom2 != bondAtom1 && atom2 != bondAtom2) {
						double minDistance = vdwr1 + ConformerGenerator.getToleratedVDWRadius(
								conformer.getMolecule().getAtomicNo(atom2));
						double dx = Math.abs(conformer.getX(atom1) - conformer.getX(atom2));
						if (dx < minDistance) {
							double dy = Math.abs(conformer.getY(atom1) - conformer.getY(atom2));
							if (dy < minDistance) {
								double dz = Math.abs(conformer.getZ(atom1) - conformer.getZ(atom2));
								if (dz < minDistance) {
									double distance = Math.sqrt(dx*dx+dy*dy+dz*dz);
									if (distance < minDistance) {
										double p = (minDistance - distance) / minDistance;
										panalty += p*p;
										}
									}
								}
							}
						}
					}
				}
			}
		return panalty;
		}

	/**
	 * Removes all torsions with non-positive likelyhoods.
	 * Sorts torsions, frequencies and likelyhoods by descending likelyhoods.
	 * Translates originalTorsionIndex to new sorted and resized torsion values.
	 * @param originalTorsionIndex
	 * @return adapted torsion index
	 */
	private int removeUnlikelyTorsions(int originalTorsionIndex) {
		int newTorsionIndex = -1;
		int count = 0;
		for (int t=0; t<mTorsion.length; t++)
			if (mLikelyhood[t] > 0f)
				count++;

		short[] newTorsion = new short[count];
		short[] newFrequency = new short[count];
		double[] newLikelyhood = new double[count];
		for (int i=0; i<count; i++) {
			double maxLikelyhood = 0f;
			int maxIndex = -1;
			for (int t=0; t<mTorsion.length; t++) {
				if (maxLikelyhood < mLikelyhood[t]) {
					maxLikelyhood = mLikelyhood[t];
					maxIndex = t;
					}
				}
			newTorsion[i] = mTorsion[maxIndex];
			newFrequency[i] = mFrequency[maxIndex];
			newLikelyhood[i] = mLikelyhood[maxIndex];

			if (maxIndex == originalTorsionIndex)
				newTorsionIndex = i;

			mLikelyhood[maxIndex] = 0f;
			}
		mTorsion = newTorsion;
		mFrequency = newFrequency;
		mLikelyhood = newLikelyhood;

		return newTorsionIndex;
		}

	private void insertTorsion(int t, int r, double relativeStrain) {
		short[] newTorsion = new short[mTorsion.length+1];
		short[][] newRange = new short[mTorsion.length+1][];
		short[] newFrequency = new short[mTorsion.length+1];
		double[] newLikelyhood = new double[mTorsion.length+1];
		int oldIndex = 0;

		short torsion = mTorsionRange[t][r];
		short frequency = (short)((mFrequency[t]+1) / 2);
		double likelyhood = mFrequency[t] * (1.0 - relativeStrain * relativeStrain);

		for (int i=0; i<=mTorsion.length; i++) {
			if (i == t+1) {
				newTorsion[i] = torsion;
				newRange[i] = new short[2];
				newRange[i][0] = torsion;
				newRange[i][1] = torsion;
				newFrequency[i] = frequency;
				newLikelyhood[i] = likelyhood;
				}
			else {
				newTorsion[i] = mTorsion[oldIndex];
				newRange[i] = mTorsionRange[oldIndex];
				newFrequency[i] = mFrequency[oldIndex];
				newLikelyhood[i] = mLikelyhood[oldIndex];
				oldIndex++;
				}
			}

		mTorsion = newTorsion;
		mTorsionRange = newRange;
		mFrequency = newFrequency;
		mLikelyhood = newLikelyhood;
		}

	private static double[][] getRotationMatrix(Coordinates n, double alpha) {
		double sinAlpha = Math.sin(alpha);
		double cosAlpha = Math.cos(alpha);
		double invcosAlpha = 1.0-cosAlpha;

		// rotation matrix is:  m11 m12 m13
		//					    m21 m22 m23
		//					    m31 m32 m33
		double[][] m = new double[3][3];
		m[0][0] = n.x*n.x*invcosAlpha+cosAlpha;
		m[1][1] = n.y*n.y*invcosAlpha+cosAlpha;
		m[2][2] = n.z*n.z*invcosAlpha+cosAlpha;
		m[0][1] = n.x*n.y*invcosAlpha-n.z*sinAlpha;
		m[1][2] = n.y*n.z*invcosAlpha-n.x*sinAlpha;
		m[2][0] = n.z*n.x*invcosAlpha-n.y*sinAlpha;
		m[0][2] = n.x*n.z*invcosAlpha+n.y*sinAlpha;
		m[1][0] = n.y*n.x*invcosAlpha+n.z*sinAlpha;
		m[2][1] = n.z*n.y*invcosAlpha+n.x*sinAlpha;
		return m;
		}

	/**
	 * Rotate the smaller side of the molecule around this bond
	 * to reach the torsion angle defined by torsionIndex.
	 * @param conformer
	 * @param torsionIndex
	 */
	public void rotateToIndex(Conformer conformer, int torsionIndex) {
		rotateTo(conformer, mTorsion[torsionIndex]);
		}

	/**
	 * Rotate the smaller side of the molecule around this bond
	 * to reach the defined torsion angle.
	 * @param conformer
	 * @param torsion in degrees (0 ... 359)
	 */
	public void rotateTo(Conformer conformer, short torsion) {
		if (torsion != conformer.getBondTorsion(mBond)) {
			int deltaTorsion = torsion - conformer.getBondTorsion(mBond);
			rotateSmallerSide(conformer, Math.PI * deltaTorsion / 180.0);
			conformer.setBondTorsion(mBond, torsion);
			}
		}

	/**
	 * If we have a BINAP stereo contraint, we have to remove colliding torsions
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
			for (int i=0; i<mTorsion.length; i++) {
				if (mTorsion[i]<180 ^ inverse) {
					frequencySum += mFrequency[i];
					count++;
					}
				}

			if (count < mTorsion.length) {
				short[] newTorsion = new short[count];
				short[] newFrequency = new short[count];
				short[][] newRange = new short[count][];
				count = 0;
				for (int i=0; i<mTorsion.length; i++) {
					if (mTorsion[i]<180 ^ inverse) {
						newTorsion[count] = mTorsion[i];
						newFrequency[count] = (short)(mFrequency[i] * 100 / frequencySum);
						newRange[count] = mTorsionRange[i];
						count++;
						}
					}
				mTorsion = newTorsion;
				mFrequency = newFrequency;
				mTorsionRange = newRange;
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
		for (int i=0; i<mTorsion.length && mTorsion[i] < maxAngle; i++) {
			frequencySum += mFrequency[i];
			count++;
			}

		if (count == 0)	// should not happen
			return;

		short[] newTorsion = new short[count];
		short[] newFrequency = new short[count];
		short[][] newRange = new short[count][];
        for (int i=0; i<count; i++) {
        	newTorsion[i] = mTorsion[i];
        	newFrequency[i] = (short)(mFrequency[i] * 100 / frequencySum);
        	newRange[i] = mTorsionRange[i];
        	}

        mTorsion = newTorsion;
        mFrequency = newFrequency;
        mTorsionRange = newRange;
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

	/**
	 * Rotates all atoms in atomList around the axis leading from atom1 through atom2
	 * by angle. The coordinates of atom1 and atom2 are not touched.
	 * @param conformer
	 * @param alpha
	 */
	private void rotateSmallerSide(Conformer conformer, double alpha) {
		Coordinates t2 = conformer.getCoordinates(mTorsionAtom[2]);
		Coordinates unit = t2.subC(conformer.getCoordinates(mTorsionAtom[1])).unit();

		double[][] m = getRotationMatrix(unit, (mRotationCenter == mTorsionAtom[1]) ? alpha : -alpha);

		for (int atom:mSmallerSideAtomList) {
			if (atom != mRotationCenter) {
				double x = conformer.getX(atom) - t2.x;
				double y = conformer.getY(atom) - t2.y;
				double z = conformer.getZ(atom) - t2.z;
				conformer.setX(atom, x*m[0][0]+y*m[0][1]+z*m[0][2] + t2.x);
				conformer.setY(atom, x*m[1][0]+y*m[1][1]+z*m[1][2] + t2.y);
				conformer.setZ(atom, x*m[2][0]+y*m[2][1]+z*m[2][2] + t2.z);
				}
			}
		}
	}