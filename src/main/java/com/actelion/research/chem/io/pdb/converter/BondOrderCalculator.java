/*
 * Copyright 2025 Thomas Sander, openmolecules.org
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

package com.actelion.research.chem.io.pdb.converter;

import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.BondLengthSet;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.util.DoubleFormat;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * This class calculates bond orders from atom coordinates only.
 */
public class BondOrderCalculator {
	private static final boolean DEBUG_OUTPUT = false;

	private static final double HALF_PI = 0.5 * Math.PI;
	private static final double SP_BOND_ANGLE_LIMIT = 0.85 * Math.PI;
	private static final double SP2_ANGLE_SUM = 2.000 * Math.PI;
	private static final double SP3_ANGLE_SUM = 1.833 * Math.PI;
	private static final double OUT_OF_PLANE_TORSION_LIMIT_DBOND = 0.25 * Math.PI;	// relaxed setting, because PDB errors are sometimes substantial
	private static final double OUT_OF_PLANE_TORSION_LIMIT_SBOND = 0.10 * Math.PI;	// tolerance increases with decreasing bond length, i.e. higher double bond likelyhood
	private static final double OUT_OF_PLANE_TORSION_FLATNESS_LIMIT = 0.20 * Math.PI;	// value considered to have flatness 0.0
	private static final double METAL_HYDROGEN_COLLISION_ANGLE = 0.25 * Math.PI;
	private static final double SINGLE_DOUBLE_BOND_LENGTH_FACTOR = 0.92;
	private static final double DOUBLE_TRIPLE_BOND_LENGTH_FACTOR = 0.88;
	private static final double MIN_DIF_SINGLE_TO_DOUBLE_BOND_LENGTH = 0.08;
	private static final double MIN_DIF_SINGLE_TO_AROM_BOND_LENGTH = 0.05;
	private static final double MIN_DIF_DOUBLE_TO_TRIPLE_BOND_LENGTH = 0.08;
	private static final double AROMATIC_5RING_BOND_LENGTH_TOLERANCE = 0.07;	// larger values -> more aromatic 5-rings
	private static final double RING_BOND_AROMATICITY_LIMIT = 0.25;
	private static final double SP_BOND_LENGTH_SUM_TOLERANCE = 0.1;
	private static final double SP_BOND_LENGTH_SUM_CONTRIBUTION_FACTOR = 10;
	private static final double PROPARGYL_LIKELYHOOD_INCREASE = 0.6;
	private static final double PATH_START_AND_END_TOLERANCE = 0.1;
	private static final double CHINONE_CONVERSION_MINIMUM_SCORE = 0.02;	// was 0.02
	private static final double HETERO_ATOM_CHARGE_PENALTY = 0.1;
	private static final double DELOCALIZED_ZERO_PI_NITROGEN_MALUS = 0.05;
	private static final double ENOL_TAUTOMER_MALUS = 0.2;
	private static final double PYRIDINONE_BONUS = 0.1;
	private static final double METAL_LIGAND_BONUS = 0.02;
	// TODO do we really want to put a preference on endo-cyclic just because guanine is usually drawn that way?
	private static final double GUANIDINE_AROM_NEIGHBOUR_BONUS = 0.01;	// allowed additional deviation to prefer endo-cyclic guanine type double bonds in resonance to aromatic ring

	private static final int MAX_ASSIGNMENT_ROUNDS = 5;

	private final Molecule3D mMol;
	private double[] mBondLength,mFractionalBondOrder;
	private final double[] mBondFlatness;
	private final int[] mHybridisation;
	private final boolean[] mIsOutOfPlaneBond,mBondIsFinal;
	private boolean[] mIsDelocalizedBond,mIsAromaticBond,mIsAromaticAtom;
	private final Random mRandom;
	private final ArrayList<Neighbour> mNeighbourList;

// TODO remove this
private static String sDebugLigandID;
public static void setDebugLigandID(String id) { sDebugLigandID = id; }

	/**
	 * This class calculates bond orders from atom coordinates only.
	 * Neither bonds nor hydrogen atoms are required, but may be present.
	 * If the molecule contains bonds, then these are expected to be metal and single bonds only.
	 * For all atoms the molecule must contain 3D-coordinates, e.g. taken from a PDB or MMCIF file.
	 * The BondOrderCalculator is not thread-safe, which means every separate thread needs its own
	 * BondOrderCalculator.
	 * @param mol
	 */
	public BondOrderCalculator(Molecule3D mol) throws Exception {
		mMol = mol;
		mRandom = new Random(1234);
		mNeighbourList = new ArrayList<>();

		if (mMol.getAllBonds() == 0)
			BondsCalculator.createBonds(mMol, true, null);

		mMol.ensureHelperArrays(Molecule.cHelperRings);
		mHybridisation = new int[mMol.getAtoms()];
		mBondIsFinal = new boolean[mMol.getBonds()];
		mIsOutOfPlaneBond = new boolean[mMol.getBonds()];
		mBondFlatness = new double[mMol.getBonds()];
	}

	public void calculateBondOrders() {
if (sDebugLigandID != null && !mMol.getName().contains(sDebugLigandID)) return;

		if (mMol.getBonds() == 0)
			return;

		// Don't touch metal bonds!
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mMol.getBondType(bond) == Molecule.cBondTypeMetalLigand)
				mBondIsFinal[bond] = true;

		// Calculate lengths of all bonds
		mBondLength = new double[mMol.getBonds()];
		for (int bond=0; bond<mMol.getBonds(); bond++)
			mBondLength[bond] = mMol.getAtomCoordinates(mMol.getBondAtom(0, bond)).subC(mMol.getAtomCoordinates(mMol.getBondAtom(1, bond))).getLength();

		calculateFractionalBondOrders();

		determineObviousHybridisation();

		// Determine all but terminal SP1 by bond angles
		// Determine out-of-plane torsions for all bonds
		determineOutOfPlaneTorsions();
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mIsOutOfPlaneBond[bond])
				assignBond(bond, 1, true);

		assignObviousBonds();

		determineAromaticBonds();
		boolean success = new AromaticityResolver(mMol, mBondLength).locateDelocalizedDoubleBonds(mIsAromaticBond.clone(), true, false);
		if (!success)
			System.out.println("$$$ WARNING: Assignment of aromatic ring bonds failed.");

		// declare all assigned aromatic bonds as final
		mIsAromaticAtom = new boolean[mMol.getAtoms()];
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mIsAromaticBond[bond]) {
				mIsAromaticAtom[mMol.getBondAtom(0, bond)] = true;
				mIsAromaticAtom[mMol.getBondAtom(1, bond)] = true;
			}
		}

		correctPyridinonesAndAlike();

		// Oxidative bond change of type -=(-=)n- to =-(=-)n= if length statistics suggests it
		// First and last bond must not be aromatic, all other must be!
		correctChinonesAndAlike();

		// aromatic bonds only
		correctTautomers(true);

		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mIsAromaticBond[bond])
				mBondIsFinal[bond] = true;

		assignKnownGroups();

		// Process bonds multiple times in random order and change bond type according to length statistics.
		// Also flip conjugated bonds to opposite tautomer, if length statistics suggests it.
		assignNonAromaticBondOrders();

		convertObviousHydroxyIminesToAmides();

		// doesn't care about final bonds
		correctNonLinearTripleBonds();

		// doesn't care about final bonds
		correctExceededValences();

		// doesn't care about final bonds
		correctForgottenSP2();

		// non-final bonds only
		correctTautomers(false);

		correctAtomCharges();

		correctMetalLigandCharges();

		try { mMol.canonizeCharge(true, false); } catch(Exception e) {}
	}

	private void determineObviousHybridisation() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getConnAtoms(atom) == 2) {
				double angle = GeometryCalculator.getAngle(mMol, mMol.getConnAtom(atom, 0), atom, mMol.getConnAtom(atom, 1));
				if (angle > SP_BOND_ANGLE_LIMIT) {
					double expectedLengthSum = 0;
					double lengthSum = 0;
					for (int i=0; i<2; i++) {
						lengthSum += mBondLength[mMol.getConnBond(atom, i)];
						int connAtom = mMol.getConnAtom(atom, i);
						int bondIndex = BondLengthSet.getBondIndex(2, false, false,
								mMol.getAtomicNo(atom), mMol.getAtomicNo(connAtom), 2, 1,
								mMol.getConnAtoms(atom), mMol.getConnAtoms(connAtom), false);
						if (bondIndex != -1)
							expectedLengthSum += BondLengthSet.getBondLength(bondIndex);
					}
					double angleContribution = (angle - SP_BOND_ANGLE_LIMIT) / (Math.PI - SP_BOND_ANGLE_LIMIT);
					double lengthContribution = SP_BOND_LENGTH_SUM_CONTRIBUTION_FACTOR * (expectedLengthSum - lengthSum + SP_BOND_LENGTH_SUM_TOLERANCE);
					if (2.0 * angleContribution + lengthContribution > 0.0)
						assignAtom(atom, 1);
				}
			}
			else if (mMol.getAllConnAtoms(atom) == 3) {
				double angleSum = GeometryCalculator.getAngle(mMol, mMol.getConnAtom(atom, 0), atom, mMol.getConnAtom(atom, 1))
								+ GeometryCalculator.getAngle(mMol, mMol.getConnAtom(atom, 0), atom, mMol.getConnAtom(atom, 2))
								+ GeometryCalculator.getAngle(mMol, mMol.getConnAtom(atom, 1), atom, mMol.getConnAtom(atom, 2));
				double planarity = (angleSum - SP3_ANGLE_SUM) / (SP2_ANGLE_SUM - SP3_ANGLE_SUM);

				// for carbon atoms we also consider relative bond length of a potential double bond as 20% of the criterion
				double adjustSP2 = 0.0;
				double adjustSP3 = 0.0;
				if (mMol.getAtomicNo(atom) == 6) {
					double fractionalBondOrder = 1.0;
					int conns = mMol.getConnAtoms(atom);
					for (int i=0; i<conns; i++) {
						int connAtom = mMol.getConnAtom(atom, i);
						int connBond = mMol.getConnBond(atom, i);
						if (!mIsOutOfPlaneBond[connBond]
						 && mHybridisation[connAtom] != 3
						 && fractionalBondOrder < mFractionalBondOrder[connBond])
							fractionalBondOrder = Math.min(2.0, mFractionalBondOrder[connBond]);
					}

					adjustSP2 = 0.4 * (fractionalBondOrder - 1.5);
					adjustSP3 = 0.4 * (fractionalBondOrder - 1.5);
//					planarity = 0.2 * (4.0 * planarity + (fractionalBondOrder - 1.0));
				}

				if (planarity + adjustSP2 > 0.90)
					assignAtom(atom, 2);
				else if (planarity + adjustSP3 < 0.5)
					assignAtom(atom, 3);
			}
			else if (mMol.getAllConnAtoms(atom) == 4) {
				assignAtom(atom, 3);
			}
		}
	}

	/**
	 * Determine all but terminal SP1 by bond angles and assign 'flat' and 'out-of-plane' to bonds with non-uncertain torsions.
	 */
	private void determineOutOfPlaneTorsions() {
		int[] tAtom = new int[4];
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			tAtom[1] = mMol.getBondAtom(0, bond);
			tAtom[2] = mMol.getBondAtom(1, bond);
			int conns1 = mMol.getAllConnAtoms(tAtom[1]);
			int conns2 = mMol.getAllConnAtoms(tAtom[2]);

			// If one of the bond's atoms has a 180-degree
			if (mHybridisation[tAtom[1]] == 1 || mHybridisation[tAtom[2]] == 1)
				continue;

			// Both atoms have no more neighbour
			if (conns1 == 1 && conns2 == 1)
				continue;

			double adaptedOutOfPlaneLimit = -1.0;	// adapt limit depending on double bond character of central bond

			if (conns1 > 1 && conns2 > 1) {
				for (int i=0; i<conns1; i++) {
					tAtom[0] = mMol.getConnAtom(tAtom[1], i);
					if (tAtom[0] != tAtom[2]) {
						for (int j=0; j<conns2; j++) {
							tAtom[3] = mMol.getConnAtom(tAtom[2], j);
							if (tAtom[3] != tAtom[1]) {
								double outOfPlaneTorsion = HALF_PI - Math.abs(Math.abs(TorsionDB.calculateTorsionExtended(mMol, tAtom)) - HALF_PI);
								if (outOfPlaneTorsion > OUT_OF_PLANE_TORSION_LIMIT_SBOND) {
									if (adaptedOutOfPlaneLimit == -1.0) {
										double doubleBondLikelyhood = 0.8;	// relaxed default
										int doubleBondIndex = BondLengthSet.getBondIndex(2, false, false,
												mMol.getAtomicNo(tAtom[1]), mMol.getAtomicNo(tAtom[2]), 1, 1, conns1, conns2, false);
										int singleBondIndex = BondLengthSet.getBondIndex(1, false, false,
												mMol.getAtomicNo(tAtom[1]), mMol.getAtomicNo(tAtom[2]), 0, 0, conns1, conns2, false);
										if (singleBondIndex != -1 && doubleBondIndex != -1) {
											double singleBondLength = BondLengthSet.getBondLength(singleBondIndex);
											double doubleBondLength = BondLengthSet.getBondLength(doubleBondIndex);
											doubleBondLikelyhood = (singleBondLength > doubleBondLength + MIN_DIF_SINGLE_TO_DOUBLE_BOND_LENGTH) ?
												(singleBondLength - mBondLength[bond]) / (singleBondLength - doubleBondLength) : 0.0;
										}
										adaptedOutOfPlaneLimit = OUT_OF_PLANE_TORSION_LIMIT_SBOND + doubleBondLikelyhood *
												(OUT_OF_PLANE_TORSION_LIMIT_DBOND - OUT_OF_PLANE_TORSION_LIMIT_SBOND);
									}
									if (outOfPlaneTorsion > adaptedOutOfPlaneLimit)
										mIsOutOfPlaneBond[bond] = true;
								}
								if (mMol.getConnAtoms(tAtom[1]) + mMol.getConnAtoms(tAtom[2]) == 6) {
									double flatness = Math.max(0.01, 1.0 - outOfPlaneTorsion / OUT_OF_PLANE_TORSION_FLATNESS_LIMIT);
									if (mBondFlatness[bond] == 0.0 || mBondFlatness[bond] > flatness)
										mBondFlatness[bond] = flatness;
								}
							}
						}
					}
				}
			}
		}
	}

	private void assignObviousBonds() {
		// All bonds connected to small atom SP3 are single
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mHybridisation[atom] == 3 && mMol.getAtomicNo(atom) <= 9)
				for (int i=0; i<mMol.getConnAtoms(atom); i++)
					assignBond(mMol.getConnBond(atom, i), 1, true);

		// If there is only one qualifying neighbour to SP2-carbon, then make it a double bond and follow-up neighbours
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mHybridisation[atom] == 2 && mMol.getAtomicNo(atom) == 6 && getAtomPi(atom) == 0) {
				int onlyQualifyingConnAtom = -1;
				int onlyQualifyingConnBond = -1;
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					int connAtom = mMol.getConnAtom(atom, i);
					int connBond = mMol.getConnBond(atom, i);
					if (!mBondIsFinal[connBond]
					 && getFreeAtomValence(mMol.getConnAtom(atom, i), -1, false) > 0) {
						if (onlyQualifyingConnBond != -1) {
							onlyQualifyingConnBond = -2;
							break;
						}
						onlyQualifyingConnAtom = connAtom;
						onlyQualifyingConnBond = connBond;
					}
				}
				if (onlyQualifyingConnBond >= 0) {
					finalizeDoubleBondChain(onlyQualifyingConnAtom, onlyQualifyingConnBond);
					for (int i=0; i<mMol.getConnAtoms(atom); i++) {
						int connBond = mMol.getConnBond(atom, i);
						if (connBond != onlyQualifyingConnBond && !mBondIsFinal[connBond])
							assignBond(connBond, 1, true);
					}
				}
			}
		}

		// for sp-atoms with two neighbours of which at least one is not sp,
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			int atomPi = getAtomPi(atom);
			if (mHybridisation[atom] == 1
			 && mMol.getConnAtoms(atom) == 2	// should always be the case
			 && (mHybridisation[mMol.getConnAtom(atom, 0)] != 1 || mHybridisation[mMol.getConnAtom(atom, 1)] != 1)
			 && atomPi < 2
			 && (!mBondIsFinal[mMol.getConnBond(atom, 0)] || !mBondIsFinal[mMol.getConnBond(atom, 1)])) {
				int[] freeValence = new int[2];
				for (int i=0; i<2; i++)
					if (!mBondIsFinal[mMol.getConnBond(atom, i)])
						freeValence[i] = getFreeAtomValence(mMol.getConnAtom(atom, i), -1, false);

				if (freeValence[0]+freeValence[1] < 2-atomPi) {
					System.out.println("$$$ ERROR: Couldn't increase order at sp-hybridized atom.");
					continue;
				}

				boolean done = false;
				for (int i=0; i<2; i++) {
					int connBond = mMol.getConnBond(atom, i);
					if (freeValence[i] <= 0) {
						assignBond(mMol.getConnBond(atom, 1 - i), 4-mMol.getBondOrder(connBond), true);
						done = true;
						break;
					}
				}
				if (!done && atomPi == 0 && freeValence[0] == 1 && freeValence[1] == 1) {
					for (int i=0; i<2; i++)
						assignBond(mMol.getConnBond(atom, i), 2, true);
					done = true;
				}
				if (!done) {
					int centralAtomicNo = mMol.getAtomicNo(atom);
					int[] preferred = new int[2];
					double[][] likelyhood = new double[2][3];	// likelyhoods for single,double,triple on both sides
					for (int i=0; i<2; i++) {
						int connAtom = mMol.getConnAtom(atom, i);
						int connBond = mMol.getConnBond(atom, i);
						int atomicNo = mMol.getAtomicNo(connAtom);
						int connConns = mMol.getConnAtoms(connAtom);
						for (int j=0; j<3; j++) {
							if (freeValence[i] >= j) {
								int connPi = (j==0 && mHybridisation[connAtom] == 2) ? 1 : getAtomPi(connAtom)+j;
								int index = BondLengthSet.getBondIndex(1+j, false, false,
										centralAtomicNo, atomicNo, 2, connPi, 2, connConns, false);
								if (index != -1) {
									double bondLength = BondLengthSet.getBondLength(index);
									likelyhood[i][j] = 1.0 - Math.abs(mBondLength[connBond] - bondLength);
								}
							}
						}
						if (likelyhood[i][0] != 0 && mHybridisation[connAtom] == 2)	// We cheat here: if in resonance with a pi-system,
							likelyhood[i][0] += PROPARGYL_LIKELYHOOD_INCREASE;		// propargyl single bonds are much shorter. We adapt...
						double bestLikelyhood = 0;
						for (int j=0; j<3; j++) {
							if (bestLikelyhood < likelyhood[i][j]) {
								bestLikelyhood = likelyhood[i][j];
								preferred[i] = j;
							}
						}
					}
					while (preferred[0] + preferred[1] > 2) {
						int betterIndex = -1;
						double bestPenalty = Double.MAX_VALUE;
						for (int i=0; i<2; i++) {
							if (preferred[i] > 0 && likelyhood[i][preferred[i]-1] != 0) {
								double penalty = likelyhood[i][preferred[i]] - likelyhood[i][preferred[i]-1];
								if (bestPenalty > penalty) {
									bestPenalty = penalty;
									betterIndex = i;
								}
							}
						}
						if (betterIndex == -1)
							break;	// should never happen
						preferred[betterIndex]--;
					}
					while (preferred[0] + preferred[1] < 2) {
						int betterIndex = -1;
						double bestPenalty = Double.MAX_VALUE;
						for (int i=0; i<2; i++) {
							if (preferred[i] < 2 && likelyhood[i][preferred[i]+1] != 0) {
								double penalty = likelyhood[i][preferred[i]] - likelyhood[i][preferred[i]+1];
								if (bestPenalty > penalty) {
									bestPenalty = penalty;
									betterIndex = i;
								}
							}
						}
						if (betterIndex == -1)
							break;	// should never happen
						preferred[betterIndex]++;
					}
					for (int i=0; i<2; i++)
						assignBond(mMol.getConnBond(atom, i), 1+preferred[i], true);
				}
			}
		}
	}

	/**
	 * Calculates for all non-aromatic bonds a fractional bond order between 1.0 and 3.0 representing
	 * bond order probability values when comparing real bond lengths with the ones from the crystallographic database.
	 * This method assumes that all neighbour bonds are single bonds, if they are not already defined otherwise.
	 */
	private void calculateFractionalBondOrders() {
		mFractionalBondOrder = new double[mMol.getBonds()];

		int[] atom = new int[2];
		int[] piCount = new int[2];
		int[] conns = new int[2];
		int[] freeValence = new int[2];

		for (int bond=0; bond<mMol.getBonds(); bond++) {
			double metalLigandBonus = 0.0;	// if we have a coordinating metal atom, we increase tendency to higher bond order

			for (int i=0; i<2; i++) {
				atom[i] = mMol.getBondAtom(i, bond);
				if (mHybridisation[atom[i]] == 3 && mMol.getAtomicNo(atom[i]) < 14) {
					mFractionalBondOrder[bond] = 1.0;
					break;
				}

				freeValence[i] = getFreeAtomValence(atom[i], bond, true);
				if (freeValence[i] == 0) {
					mFractionalBondOrder[bond] = 1.0;
					break;
				}

				conns[i] = mMol.getConnAtoms(atom[i]);
				for (int k=0; k<conns[i]; k++) {
					int connBond = mMol.getConnBond(atom[i], k);
					if (connBond != bond)
						piCount[i] += mMol.getBondOrder(connBond) - 1;
				}

				if (mMol.getAllConnAtomsPlusMetalBonds(atom[i]) > mMol.getAllConnAtoms(atom[i]))
					metalLigandBonus = METAL_LIGAND_BONUS;
			}
			if (mFractionalBondOrder[bond] != 0.0)
				continue;

			int doubleBondIndex = BondLengthSet.getBondIndex(2, false, false,
					mMol.getAtomicNo(atom[0]), mMol.getAtomicNo(atom[1]), piCount[0]+1, piCount[1]+1, conns[0], conns[1], false);
			int singleBondIndex = BondLengthSet.getBondIndex(1, false, false,
					mMol.getAtomicNo(atom[0]), mMol.getAtomicNo(atom[1]), piCount[0], piCount[1], conns[0], conns[1], false);

			if (singleBondIndex == -1 && doubleBondIndex == -1) {
				mFractionalBondOrder[bond] = 1.5;	// we don't know
				continue;
			}

			double singleBondLength = BondLengthSet.getBondLength(singleBondIndex);
			double doubleBondLength = BondLengthSet.getBondLength(doubleBondIndex);

			if (singleBondIndex == -1)
				singleBondLength = doubleBondLength / SINGLE_DOUBLE_BOND_LENGTH_FACTOR;
			if (doubleBondIndex == -1)
				doubleBondLength = SINGLE_DOUBLE_BOND_LENGTH_FACTOR * singleBondLength;

			double realBondLength = mBondLength[bond] - metalLigandBonus;

			if (freeValence[0] >= 2 && freeValence[1] >= 2
					&& (conns[0] == 1 || mHybridisation[atom[0]] == 1)
					&& (conns[1] == 1 || mHybridisation[atom[1]] == 1)) {
				int tripleBondIndex = BondLengthSet.getBondIndex(3, false, false,
						mMol.getAtomicNo(atom[0]), mMol.getAtomicNo(atom[1]), 2, 2, conns[0], conns[1], false);
				double tripleBondLength = (tripleBondIndex == -1) ?
						DOUBLE_TRIPLE_BOND_LENGTH_FACTOR * doubleBondLength : BondLengthSet.getBondLength(tripleBondIndex);
				if (doubleBondLength - tripleBondLength < MIN_DIF_DOUBLE_TO_TRIPLE_BOND_LENGTH)
					tripleBondLength = doubleBondLength - MIN_DIF_DOUBLE_TO_TRIPLE_BOND_LENGTH;
				if (realBondLength < doubleBondLength) {
					mFractionalBondOrder[bond] = 2.0 + (doubleBondLength - realBondLength) / (doubleBondLength - tripleBondLength);
					continue;
				}
			}

			mFractionalBondOrder[bond] = 1.0 + (singleBondLength - realBondLength) / (singleBondLength - doubleBondLength);
		}
	}

	/**
	 * Call this method to set bond to a double bond and to continue
	 * setting obvious bond orders from the far atom of this bond.
	 * Continuation depends on knowledge of the hybridisation of atom
	 * and may end there or recursively crawl along a sp1 chain.
	 * @param atom
	 * @param bond
	 */
	private void finalizeDoubleBondChain(int atom, int bond) {
		boolean again;
		do {
			again = false;
			assignBond(bond, 2, true);
			if (mHybridisation[atom] == 2) {
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					int connBond = mMol.getConnBond(atom, i);
					if (connBond != bond && !mBondIsFinal[connBond])
						assignBond(connBond, 1, true);
				}
			}
			else if (mHybridisation[atom] == 1) {
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					int connBond = mMol.getConnBond(atom, i);
					if (connBond != bond && !mBondIsFinal[connBond]) {
						atom = mMol.getConnAtom(atom, i);
						bond = connBond;
						again = true;
						break;
					}
				}
			}
		} while (again);
	}

	private void assignKnownGroups() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mIsAromaticAtom[atom])
				continue;

			// N-C(=N)-N
			if (mMol.getAtomicNo(atom) == 6 && mMol.getConnAtoms(atom) == 3 && mHybridisation[atom] == 2) {
				Neighbour[] nitrogens = getNeighbourAtoms(atom, 7, 0, 0, 0, true);
				if (nitrogens.length == 3) {
					// if one of the guanidine nitrogens is connected to an aromatic ring (usually in guanine),
					// then the double bond is in resonance to that even if its bond length is slightly longer than the exocyclic one
					Neighbour preferredNitrogen = null;
					double preferredDeviation = 100;
					for (Neighbour nitrogen : nitrogens) {
						if (mBondIsFinal[nitrogen.getBond()]) {
							if (mMol.getBondOrder(nitrogen.getBond()) >= 2) {
								preferredNitrogen = nitrogen;
								break;
							}
							continue;
						}
						double deviation = nitrogen.getBondLengthDeviation();
						if (hasAromaticNeighbour(nitrogen.getAtom()))
							deviation -= GUANIDINE_AROM_NEIGHBOUR_BONUS;
						if (deviation < preferredDeviation) {
							preferredNitrogen = nitrogen;
							preferredDeviation = deviation;
						}
					}
					for (Neighbour nitrogen : nitrogens)
						assignBond(nitrogen.getBond(), nitrogen == preferredNitrogen ? 2 : 1, true);
					continue;
				}
			}

			// any-C(=O)-O
			if (mMol.getAtomicNo(atom) == 6 && mMol.getConnAtoms(atom) == 3 && mHybridisation[atom] == 2) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 2) {
					assignDoubleToFirstNeighbours(oxygen, 1);
					continue;
				}
			}

			// any-C(=O)-X
			if (mMol.getAtomicNo(atom) == 6 && mMol.getConnAtoms(atom) == 3 && mHybridisation[atom] == 2) {
				Neighbour[] electronegative = getNeighbourAtoms(atom, -1, 1, 0, 0, true);
				if (electronegative.length >= 2 && mMol.getAtomicNo(electronegative[0].getAtom()) == 8) {
					assignDoubleToFirstNeighbours(electronegative, 1);
					continue;
				}
			}

			// any-N(=O)-O
			if (mMol.getAtomicNo(atom) == 7 && mMol.getConnAtoms(atom) == 3) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 2) {
					mMol.setAtomCharge(atom, 1);
					assignDoubleToFirstNeighbours(oxygen, 1);
					for (int i=0; i<oxygen.length; i++) {
						if (mMol.getBondOrder(oxygen[i].getBond()) == 1) {
							mMol.setAtomCharge(oxygen[i].getAtom(), -1);
							break;
						}
					}
					continue;
				}
			}

			// any-P(=O)-any
			if (mMol.getAtomicNo(atom) == 15 && mMol.getConnAtoms(atom) == 3) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 1) {
					assignDoubleToFirstNeighbours(oxygen, 1);
					continue;
				}
			}

			// any-P(=O)(-any)-any
			if (mMol.getAtomicNo(atom) == 15 && mMol.getConnAtoms(atom) == 4) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 1) {
					assignDoubleToFirstNeighbours(oxygen, 1);
					continue;
				}
			}

			// any-S(=O)-any
			if (mMol.getAtomicNo(atom) == 16 && mMol.getConnAtoms(atom) == 3) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 1) {
					assignDoubleToFirstNeighbours(oxygen, 1);
					continue;
				}
			}

			// any-S(=O)(=O)-any
			if (mMol.getAtomicNo(atom) == 16 && mMol.getConnAtoms(atom) == 4) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 2) {
					assignDoubleToFirstNeighbours(oxygen, 2);
					continue;
				}
			}
		}
	}

	private void assignDoubleToFirstNeighbours(Neighbour[] neighbour, int doubleBondCount) {
		int doubleBondsAssigned = 0;
		for (int i=0; i<neighbour.length; i++) {
			if (mBondIsFinal[neighbour[i].getBond()]) {
				if (mMol.getBondOrder(neighbour[i].getBond()) >= 2)
					doubleBondsAssigned++;
				continue;
			}
			if (doubleBondsAssigned<doubleBondCount && getFreeAtomValence(neighbour[i].getAtom(), -1, false) > 0) {
				assignBond(neighbour[i].getBond(), 2, true);
				doubleBondsAssigned++;
				continue;
			}
			assignBond(neighbour[i].getBond(), 1, true);
		}
	}

	private boolean hasAromaticNeighbour(int atom) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (mIsAromaticAtom[mMol.getConnAtom(atom, i)])
				return true;

		return false;
	}

	private void determineAromaticBonds() {
		RingCollection ringSet = mMol.getRingSet();
		boolean[] isAromaticRing = new boolean[ringSet.getSize()];
		Arrays.fill(isAromaticRing, true);	// default

		boolean[] isDelocalizedBond = new boolean[mMol.getBonds()];
		boolean[] hasBulkySubstituent = new boolean[8];

		for (int r=0; r<ringSet.getSize(); r++) {
			int[] ringAtom = ringSet.getRingAtoms(r);
			int[] ringBond = ringSet.getRingBonds(r);

			if (ringBond.length <= 4) {
				isAromaticRing[r] = false;
				continue;
			}

			for (int bond : ringBond) {
				if (mBondIsFinal[bond]
				 || mIsOutOfPlaneBond[bond]) {
					isAromaticRing[r] = false;
					break;
				}
			}
			if (!isAromaticRing[r])
				continue;

			double angleSum = 0.0;
			for (int i = 0; i<ringAtom.length; i++)
				angleSum += GeometryCalculator.getAngle(mMol, ringAtom[(i + 1) % ringAtom.length], ringAtom[i], ringAtom[(i + ringAtom.length - 1) % ringAtom.length]);

			double angleSumDeviation = Math.PI * (ringAtom.length - 2) - angleSum;
			double flatness = Math.max(0.0, 1.0 - angleSumDeviation /
				 ((ringAtom.length == 4) ? 0.10 * Math.PI
				: (ringAtom.length == 5) ? 0.067 * Math.PI    // tiny margin: 108 degrees is flat 5-ring; 109.5 degrees is perfect sp3
				: (ringAtom.length == 6) ? 0.125 * Math.PI : 0.133 * Math.PI));

			double[] bondAromaticity = new double[ringAtom.length];
			for (int i=0; i<ringAtom.length; i++)
				bondAromaticity[i] = getBondAromaticityFromLength(ringBond[i], ringBond.length);
			Arrays.sort(bondAromaticity);

			// we use the second-worst bond aromaticity to tolerate one too long bond
			if (flatness * bondAromaticity[1] < RING_BOND_AROMATICITY_LIMIT) {
				isAromaticRing[r] = false;
				continue;
			}

			int oxaCount = 0;
			int carbonCount = 0;
			int metalNeighbourCount = 0;
			for (int atom : ringAtom) {
				int atomicNo = mMol.getAtomicNo(atom);
				if (atomicNo != 5 && atomicNo != 6 && atomicNo != 7 && atomicNo != 8 && atomicNo != 14 && atomicNo != 15 && atomicNo != 16) {
					isAromaticRing[r] = false;
					break;
				}
				if (!isSPSeTe(atom) && (mHybridisation[atom] == 1 || mHybridisation[atom] == 3)) {
					isAromaticRing[r] = false;
					break;
				}
				if (atomicNo == 6) {
					carbonCount++;
				}
				if (atomicNo == 8 || atomicNo == 16) {
					oxaCount++;
				}
				if (mMol.getAllConnAtoms(atom) != mMol.getAllConnAtomsPlusMetalBonds(atom)) {
					metalNeighbourCount++;
				}
			}
			if (ringBond.length == 5 && carbonCount == 5 && metalNeighbourCount == 0)
				isAromaticRing[r] = false;
			if (oxaCount > (ringBond.length == 6 ? 0 : 1))
				isAromaticRing[r] = false;
			if (!isAromaticRing[r])
				continue;

			if (ringSet.getRingSize(r) == 6)
				for (int bond : ringSet.getRingBonds(r))
					isDelocalizedBond[bond] = true;
		}

		// flat 5-rings may not be aromatic; we check bond lengths and set ring as non-aromatic in case of too long bond
		for (int r=0; r<ringSet.getSize(); r++) {
			int[] ringBond = ringSet.getRingBonds(r);
			if (ringBond.length == 5 && isAromaticRing[r]) {
				for (int bond : ringBond) {
					if (!isDelocalizedBond[bond]) {
						int atom1 = mMol.getBondAtom(0, bond);
						int atom2 = mMol.getBondAtom(1, bond);
						int atomicNo1 = mMol.getAtomicNo(atom1);
						int atomicNo2 = mMol.getAtomicNo(atom2);
						int conns1 = mMol.getConnAtoms(atom1);
						int conns2 = mMol.getConnAtoms(atom2);
						boolean isAtom1SOrO = (mMol.getAtomicNo(atom1) == 8 || mMol.getAtomicNo(atom1) == 16) && conns1 == 2;
						boolean isAtom2SOrO = (mMol.getAtomicNo(atom2) == 8 || mMol.getAtomicNo(atom2) == 16) && conns2 == 2;

						// check against length of aromatic single bond (don't check for double bond, because that would be shorter anyway)
						int bondIndex = BondLengthSet.getBondIndex(1, true, false, atomicNo1, atomicNo2, isAtom1SOrO ? 0 : 1, isAtom2SOrO ? 0 : 1, conns1, conns2, false);
						if (bondIndex != -1) {
							double bondLength = BondLengthSet.getBondLength(bondIndex);
							if (mBondLength[bond] - bondLength >AROMATIC_5RING_BOND_LENGTH_TOLERANCE) {
								isAromaticRing[r] = false;
								break;
							}
						}
					}
				}
			}
		}

		mIsAromaticBond = new boolean[mMol.getBonds()];
		mIsDelocalizedBond = new boolean[mMol.getBonds()];
		for (int r=0; r<ringSet.getSize(); r++) {
			if (isAromaticRing[r]) {
				for (int atom : ringSet.getRingAtoms(r))
					if ((!isSPSeTe(atom) || mMol.getConnAtoms(atom) == 2))
						assignAtom(atom, 2);

				for (int bond : mMol.getRingSet().getRingBonds(r)) {
					mIsAromaticBond[bond] = true;
					if (ringSet.getRingSize(r) == 6)
						mIsDelocalizedBond[bond] = true;
				}
			}
		}

		markHemeBridgesAsAromatic(isAromaticRing);
	}

	private double getBondAromaticityFromLength(int bond, int ringSize) {
		int atom1 = mMol.getBondAtom(0, bond);
		int atom2 = mMol.getBondAtom(1, bond);
		int pi1 = mMol.isElectronegative(atom1) ? 0 : 1;
		int pi2 = mMol.isElectronegative(atom2) ? 0 : 1;
		if (pi1 + pi2 == 0) {
			if (mMol.getAtomicNo(atom1) == 7)
				pi1 = 1;
			else if (mMol.getAtomicNo(atom2) == 7)
				pi2 = 1;
		}
		int conns1 = mMol.getConnAtoms(atom1);
		int conns2 = mMol.getConnAtoms(atom2);
		int atomicNo1 = mMol.getAtomicNo(atom1);
		int atomicNo2 = mMol.getAtomicNo(atom2);
		int index1 = BondLengthSet.getBondIndex(1, true,
				ringSize == 6, atomicNo1, atomicNo2, pi1, pi2, conns1, conns2, false);
		if (index1 == -1)
			return 1.0;
		int index2 = BondLengthSet.getBondIndex(1, false,
				false, atomicNo1, atomicNo2, 0, 0, conns1, conns2, false);
		if (index2 == -1)
			return 1.0;
		float aromLength = BondLengthSet.getBondLength(index1);
		float nonAromLength = BondLengthSet.getBondLength(index2);

		// No show-stopper, if difference between aromatic and non-aromatic bond lengths is insignificant; typically N-N, N-O
		if (nonAromLength - aromLength < MIN_DIF_SINGLE_TO_AROM_BOND_LENGTH)
			return 1.0;

		// if we have bond flatness information, then we include it as 50% of the result value
		double aromaticity = Math.max(0.0, (nonAromLength - mBondLength[bond]) / (nonAromLength - aromLength));
		if (mBondFlatness[bond] != 0)
			aromaticity = 0.5 * (mBondFlatness[bond] + aromaticity);

		return aromaticity;
	}

	private int getExocyclicIndex(int[] ringAtom, int i) {
		if (mMol.getConnAtoms(ringAtom[i]) > 2) {
			int[] ringConn = ringNeighbourAtoms(ringAtom, i);
			for (int j=0; j<mMol.getConnAtoms(ringAtom[i]); j++) {
				int connAtom = mMol.getConnAtom(ringAtom[i], j);
				if (connAtom != ringConn[0] && connAtom != ringConn[1])
					return j;
			}
		}
		return -1;
	}

	private void markHemeBridgesAsAromatic(boolean[] isAromaticRing) {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (qualifiesForHemeBridge(atom)) {
				int leftNextRoot = findNextHemeBridgeRootAtom(mMol.getConnAtom(atom, 0));
				if (leftNextRoot != -1) {
					int rightNextRoot = findNextHemeBridgeRootAtom(mMol.getConnAtom(atom, 1));
					if (rightNextRoot != -1) {
						int leftBridgeAtom = findNextHemeBridgeAtom(leftNextRoot);
						if (leftBridgeAtom != -1) {
							int rightBridgeAtom = findNextHemeBridgeAtom(rightNextRoot);
							if (rightBridgeAtom != -1) {
								int leftRemoteRoot = findNextHemeBridgeRootAtom(mMol.getConnAtom(
										leftBridgeAtom,mMol.getConnAtom(leftBridgeAtom, 0) == leftNextRoot ? 1 : 0));
								if (leftRemoteRoot != -1) {
									int rightRemoteRoot = findNextHemeBridgeRootAtom(mMol.getConnAtom(
											rightBridgeAtom,mMol.getConnAtom(rightBridgeAtom, 0) == rightNextRoot ? 1 : 0));
									if (rightRemoteRoot != -1) {
										int leftRemoteBridgeAtom = findNextHemeBridgeAtom(leftRemoteRoot);
										if (leftRemoteBridgeAtom != -1) {
											int rightRemoteBridgeAtom = findNextHemeBridgeAtom(rightRemoteRoot);
											if (rightRemoteBridgeAtom == leftRemoteBridgeAtom) {
												mIsAromaticBond[mMol.getConnBond(atom, 0)] = true;
												mIsAromaticBond[mMol.getConnBond(atom, 1)] = true;
												mIsAromaticBond[mMol.getConnBond(leftBridgeAtom, 0)] = true;
												mIsAromaticBond[mMol.getConnBond(leftBridgeAtom, 1)] = true;
												mIsAromaticBond[mMol.getConnBond(rightBridgeAtom, 0)] = true;
												mIsAromaticBond[mMol.getConnBond(rightBridgeAtom, 1)] = true;
												mIsAromaticBond[mMol.getConnBond(leftRemoteBridgeAtom, 0)] = true;
												mIsAromaticBond[mMol.getConnBond(leftRemoteBridgeAtom, 1)] = true;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	private boolean qualifiesForHemeBridge(int atom) {
		if (mMol.getConnAtoms(atom) == 2 && mHybridisation[atom] != 3) {
			int connAtom1 = mMol.getConnAtom(atom, 0);
			int connAtom2 = mMol.getConnAtom(atom, 1);
			if (mMol.getAtomRingSize(connAtom1) == 5
			 && (mIsAromaticBond[mMol.getConnBond(connAtom1, 0)] || mIsAromaticBond[mMol.getConnBond(connAtom1, 1)])
			 && mMol.getAtomRingSize(connAtom2) == 5
			 && (mIsAromaticBond[mMol.getConnBond(connAtom2, 0)] || mIsAromaticBond[mMol.getConnBond(connAtom2, 1)])) {
				int connBond1 = mMol.getConnBond(atom, 0);
				int connBond2 = mMol.getConnBond(atom, 1);
				if (!mIsAromaticBond[connBond1] && !mIsOutOfPlaneBond[connBond1]
				 && !mIsAromaticBond[connBond2] && !mIsOutOfPlaneBond[connBond2])
					return true;
			}
		}
		return false;
	}

	private int findNextHemeBridgeRootAtom(int atom) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int nitrogen = mMol.getConnAtom(atom, i);
			if (mMol.getAtomicNo(nitrogen) == 7
			 && mMol.getAtomRingSize(nitrogen) == 5
			 && mIsAromaticBond[mMol.getConnBond(atom, i)]) {
				for (int j=0; j<mMol.getConnAtoms(nitrogen); j++) {
					int connAtom = mMol.getConnAtom(nitrogen, j);
					if (connAtom != atom
					 && mMol.getConnAtoms(connAtom) == 3
					 && mMol.getAtomRingSize(connAtom) == 5
					 && mIsAromaticBond[mMol.getConnBond(nitrogen, j)]) {
						return connAtom;
					}
				}
			}
		}
		return -1;
	}

	private int findNextHemeBridgeAtom(int atom) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int connAtom = mMol.getConnAtom(atom, i);
			if (qualifiesForHemeBridge(connAtom))
				return connAtom;
			}
		return -1;
	}

	private int[] ringNeighbourAtoms(int[] ringAtom, int i) {
		int[] atom = new int[2];
		atom[0] = ringAtom[i==0 ? ringAtom.length-1 : i-1];
		atom[1] = ringAtom[i==ringAtom.length-1 ? 0 : i+1];
		return atom;
	}

	private boolean isSPSeTe(int atom) {
		return mMol.getAtomicNo(atom) == 15
			|| mMol.getAtomicNo(atom) == 16
			|| mMol.getAtomicNo(atom) == 34
			|| mMol.getAtomicNo(atom) == 52;
	}

	private int getAtomPi(int atom) {
		int pi = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			pi += mMol.getBondOrder(mMol.getConnBond(atom, i)) - 1;
		return pi;
	}

	/**
	 * If we have sp2 carbons without pi-electrons after assigning all double bonds to match bond lengths,
	 * then we probably have an issue with bonds lengths. We collect all neighbours that are not sp3, are not connected
	 * with non-flat torsions, that have a free valence. Then we convert that neighbour with the highest bond order deviation
	 * from the expected single-bond length into a double bond.
	 */
	private void correctForgottenSP2() {
		GraphWalker walker = new GraphWalker(mMol, 12) {
			@Override
			public boolean qualifiesAsNext(int parentAtom, int atom, int bond, int size) {
				return (((size & 1) == 0) ^ (mMol.getBondOrder(bond) > 1))
					&& !mIsOutOfPlaneBond[bond]
					&& mMol.getConnAtoms(atom) <= 3;
			}

			@Override
			public boolean qualifiesAsFinal(int parentAtom, int atom, int bond, int size) {
				return (((size & 1) == 0)	// even number of atoms: single bonded end and oxidation on both ends
					 && mHybridisation[atom] != 3
					 && getAtomPi(atom) == 0
					 && getFreeAtomValence(atom, -1, false) > 0)
					|| ((size & 1) == 1    // odd number of atoms: double bonded path end
					 && (mHybridisation[atom] != 2 || (mIsAromaticAtom[atom] && mMol.getAtomicNo(atom) != 6)));
			}

			@Override
			public double calculatePathScore(int finalAtom, int[] graphParent) {
				double score = 0.0;
				int atom2 = finalAtom;
				int atom1 = graphParent[atom2];
				while (atom1 != -1) {
					int bond = mMol.getBond(atom1, atom2);
					score += calculateBondChangeScore(bond,
							mMol.getBondOrder(bond) == 1 ? 2 : mMol.getBondOrder(bond)-1,
							0,
							graphParent[atom1] != -1 ? atom1 : -1,
							atom2 != finalAtom ? atom2 : -1);
					atom2 = atom1;
					atom1 = graphParent[atom1];
				}
				return score;
			}
		};

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mHybridisation[atom] == 2 && mMol.getAtomicNo(atom) == 6 && getAtomPi(atom) == 0) {
				walker.find(atom, -1);
				int[] bestPath = walker.getBestScoringPath();
				if (bestPath != null) {
					for (int i=1; i<bestPath.length; i++) {
						// force to override potentially final bonds
						int bond = mMol.getBond(bestPath[i-1], bestPath[i]);
						mMol.setBondOrder(bond, mMol.getBondOrder(bond) == 1 ? 2 : mMol.getBondOrder(bond)-1);
//TLS						mBondIsFinal[bond] = true;
					}
				}
			}
		}
	}

	private void correctTautomers(boolean aromatic) {
		GraphWalker walker = createTautomerWalker(aromatic, 0.0, true, true);

		boolean changed=true;
		for (int i=0; changed && i<MAX_ASSIGNMENT_ROUNDS; i++) {
			changed = false;
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				changed |= tryChangeTautomer(mMol.getBondAtom(0, bond), mMol.getBondAtom(1, bond), bond, walker)
						|| tryChangeTautomer(mMol.getBondAtom(1, bond), mMol.getBondAtom(0, bond), bond, walker);
			}
		}
	}

	private void correctPyridinonesAndAlike() {
		GraphWalker walker = createTautomerWalker(true, PYRIDINONE_BONUS, false, false);

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getConnAtoms(atom) == 1) {
				int atomicNo = mMol.getAtomicNo(atom);
				int connAtom = mMol.getConnAtom(atom, 0);
				int connBond = mMol.getConnBond(atom, 0);
				if ((atomicNo == 8 || atomicNo == 16 || atomicNo == 34)
				 && mMol.getBondOrder(connBond) == 1
				 && mIsAromaticAtom[connAtom]) {
					if (tryChangeTautomer(atom, connAtom, connBond, walker))
						mBondIsFinal[connBond] = true;
				}
			}
		}
	}

	/**
	 * Chinones have been treated as aromatic ring and end up as 'hydro-chinones'.
	 * We correct here by detecting too short HO-C(atom) and oxydizing along the path
	 * of conjugated single/double bonds to a qualifying partner atom.
	 */
	private void correctChinonesAndAlike() {
		GraphWalker walker = new GraphWalker(mMol, 12) {
			@Override
			public boolean qualifiesAsNext(int parentAtom, int atom, int bond, int size) {
				return (((size & 1) == 0) ^ (mMol.getBondOrder(bond) > 1))
					&& !mIsOutOfPlaneBond[bond]
					&& mMol.getConnAtoms(atom) <= 3;
			}

			@Override
			public boolean qualifiesAsFinal(int parentAtom, int atom, int bond, int size) {
				return (size & 1) == 0	// even number of atoms (converts hydrochinone to chinones)
					&& mHybridisation[atom] != 3
					&& getAtomPi(atom) == 0
					&& !mBondIsFinal[bond]
					&& !mIsAromaticBond[bond]
					&& calculateBondChangeScore(bond, (size & 1) == 0 ? 2 : 1, 0, parentAtom, -1) > 0.0;
			}

			@Override
			public double calculatePathScore(int finalAtom, int[] graphParent) {
				double score = 0.0;
				int atom2 = finalAtom;
				int atom1 = graphParent[atom2];
				while (atom1 != -1) {
					int bond = mMol.getBond(atom1, atom2);
					score += calculateBondChangeScore(bond,
							mMol.getBondOrder(bond) == 1 ? 2 : mMol.getBondOrder(bond)-1,
							2,
							graphParent[atom1] != -1 ? atom1 : -1,
							atom2 != finalAtom ? atom2 : -1);
					atom2 = atom1;
					atom1 = graphParent[atom1];
				}
				return score;
			}
		};

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mHybridisation[atom] != 3
			 && getAtomPi(atom) == 0) {
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					int bond = mMol.getConnBond(atom, i);
					if (mMol.getBondOrder(bond) == 1
					 && !mBondIsFinal[bond]
					 && !mIsAromaticBond[bond]) {
						int connAtom = mMol.getConnAtom(atom, i);
						if (getAtomPi(connAtom) != 0
						 && calculateBondChangeScore(bond, 2, 0, connAtom, -1) > 0.0) {
							walker.find(atom, connAtom);
							if (!Double.isNaN(walker.getBestScore()) && walker.getBestScore() > CHINONE_CONVERSION_MINIMUM_SCORE) {
								int[] bestPath = walker.getBestScoringPath();
								if (bestPath != null) {
									for (int j=1; j<bestPath.length; j++) {
										int pathBond = mMol.getBond(bestPath[j-1], bestPath[j]);
										assignBond(pathBond, mMol.getBondOrder(pathBond) == 1 ? 2 : mMol.getBondOrder(pathBond)-1, false);
									}
								}
							}
						}
					}
				}
			}
		}

//		for (int bond=0; bond<mMol.getBonds(); bond++) {
//			if (mMol.getBondOrder(bond) == 1
//			 && !mBondIsFinal[bond]
//			 && !mIsAromaticBond[bond]) {
//				for (int i=0; i<2; i++) {
//					int atom = mMol.getBondAtom(i, bond);
//					int connAtom = mMol.getBondAtom(1-i, bond);
//					if (mHybridisation[atom] != 3
//					 && getAtomPi(atom) == 0
//					 && getAtomPi(connAtom) != 0
//					 && calculateBondChangeScore(bond, 2, 0, connAtom, -1) > 0.0) {
//						walker.find(atom, connAtom);
//						if (!Double.isNaN(walker.getBestScore()) && walker.getBestScore() > CHINONE_CONVERSION_MINIMUM_SCORE) {
//							int[] bestPath = walker.getBestScoringPath();
//							if (bestPath != null) {
//								for (int j=1; j<bestPath.length; j++) {
//									int pathBond = mMol.getBond(bestPath[j-1], bestPath[j]);
//									assignBond(pathBond, mMol.getBondOrder(pathBond) == 1 ? 2 : mMol.getBondOrder(pathBond)-1, false);
//								}
//							}
//						}
//					}
//				}
//			}
//		}
	}

	/**
	 * Creates a walker that finds chains of conjugated double bonds, which start with a single bond,
	 * ends with a double bond, and can be potentially flipped into the inverse tautomer.
	 */
	private GraphWalker createTautomerWalker(final boolean aromatic, final double startScore,
											 final boolean applyZeroPiNitrogenMalus,
											 final boolean firstBondIsAromatic) {
		return new GraphWalker(mMol, 12) {
			@Override
			public boolean qualifiesAsFirst(int parentAtom, int atom, int bond) {
				return (!aromatic || firstBondIsAromatic == mIsAromaticBond[bond])
					  && !mBondIsFinal[bond]
					  && mMol.getBondOrder(bond) == 1
					  && !mIsOutOfPlaneBond[bond]
					  && mHybridisation[atom] != 3
					  && getAtomPi(atom) != 0
					  && mMol.getConnAtoms(atom)<=3;
			}

			@Override
			public boolean qualifiesAsNext(int parentAtom, int atom, int bond, int size) {
				return aromatic == mIsAromaticBond[bond]
					&& !mBondIsFinal[bond]
					&& (((size & 1) == 0) ^ (mMol.getBondOrder(bond)>1))
					&& !mIsOutOfPlaneBond[bond]
					&& mHybridisation[atom] != 3
					&& getAtomPi(atom) != 0
					&& mMol.getConnAtoms(atom)<=3;
			}

			@Override
			public boolean qualifiesAsFinal(int parentAtom, int atom, int bond, int size) {
				return (size & 1) == 1    // enoles, imines, etc. and vinylogic variants
					&& (mHybridisation[atom] != 2 || mMol.getAtomicNo(atom) != 6)
					&& calculateBondChangeScore(bond, (size & 1) == 0 ? 2 : 1, 0, parentAtom, -1) > -PATH_START_AND_END_TOLERANCE;
			}

			@Override
			public double calculatePathScore(int finalAtom, int[] graphParent) {
				double score = startScore;

				// for tautomers within the aromatic part we give a malus for nitrogen in 6-membered rings if we remove the pi-bond
				if (aromatic && applyZeroPiNitrogenMalus && mMol.getAtomicNo(finalAtom) == 7 && mMol.getAtomRingSize(finalAtom) == 6)
					score -= DELOCALIZED_ZERO_PI_NITROGEN_MALUS;

				// for non-aromatic tautomers with an enol at the end, we give a malus
				if (!aromatic && mMol.getAtomicNo(finalAtom) == 8 && mMol.getConnAtoms(finalAtom) == 1)
					score -= ENOL_TAUTOMER_MALUS;

				int atom2 = finalAtom;
				int atom1 = graphParent[atom2];
				while (atom1 != -1) {
					int bond = mMol.getBond(atom1, atom2);
					score += calculateBondChangeScore(bond,
							mMol.getBondOrder(bond) == 1 ? 2 : mMol.getBondOrder(bond) - 1,
							firstBondIsAromatic ? 0 : 2,
							graphParent[atom1] != -1 ? atom1 : -1,
							atom2 != finalAtom ? atom2 : -1);
					atom2 = atom1;
					atom1 = graphParent[atom1];
				}

				// for tautomers within the aromatic part we give a bonus for nitrogen in 6-membered rings if we re-establish a pi-bond
				if (aromatic && applyZeroPiNitrogenMalus && mMol.getAtomicNo(atom2) == 7 && mMol.getAtomRingSize(atom2) == 6)
					score += DELOCALIZED_ZERO_PI_NITROGEN_MALUS;

				// for non-aromatic tautomers with an enol at the end, we give a malus
				if (!aromatic && mMol.getAtomicNo(atom2) == 8 && mMol.getConnAtoms(atom2) == 1)
					score += ENOL_TAUTOMER_MALUS;

				if (DEBUG_OUTPUT) {
 System.out.print("tauto path:");
 while (finalAtom != -1) {
  System.out.print(" " + finalAtom);
  finalAtom = graphParent[finalAtom];
 }
 System.out.println(" score:" + score);
}
				return score;
			}
		};
	}

	private boolean tryChangeBondOrder(int bond) {
// using pre-calculated value is faster, but creates slightly worse results than using live up-to-date neighbour bond orders
//		int order = (mFractionalBondOrder[bond] < 1.5) ? 1 : (mFractionalBondOrder[bond] < 2.5) ? 2 : 3;
		int order = guessBondOrderFromBondLength(bond);		// this uses live up-to-date neighbour bond orders
		if (order == mMol.getBondOrder(bond))
			return false;

		assignBond(bond, order, false);
		return true;
	}

	/**
	 * Uses the given tautomer walker to detect potential tautomer situations at atom,
	 * scores all potential tautomers, which can be created by shifting double bonds one position,
	 * and, if a score suggests it, shifts double bonds to generate the respective tautomer.
	 */
	private boolean tryChangeTautomer(int atom, int nextAtom, int bond, GraphWalker walker) {
		if (walker.qualifiesAsFirst(atom, nextAtom, bond)
		 && getAtomPi(atom) == 0
		 && getFreeAtomValence(atom, -1, false) != 0
		 && calculateBondChangeScore(bond, 2, 0, nextAtom, -1) > -PATH_START_AND_END_TOLERANCE) {
			walker.find(atom, nextAtom);
			if (!Double.isNaN(walker.getBestScore()) && walker.getBestScore() > 0.0) {
				int[] bestPath = walker.getBestScoringPath();
				if (bestPath != null) {
if (DEBUG_OUTPUT) {
 System.out.print("tauto flip:");
 for (int pa : bestPath)
 System.out.print(" " + pa);
 System.out.println(" score:" + walker.getBestScore());
}
					for (int k=1; k<bestPath.length; k++) {
						int pathBond = mMol.getBond(bestPath[k-1], bestPath[k]);
						assignBond(pathBond, mMol.getBondOrder(pathBond) == 1 ? 2 : mMol.getBondOrder(pathBond)-1, false);
						if (k == bestPath.length-1) {
							int  finalAtom = bestPath[k];
							if (mMol.getAtomicNo(finalAtom) == 6)
								mHybridisation[finalAtom] = 0;    // remove sp2 indication
							if (mMol.getAtomCharge(finalAtom) == 1 && mMol.isElectronegative(finalAtom))
								mMol.setAtomCharge(finalAtom, 0);
						}
					}
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Calculate the improvement in Angstrom regarding difference of bond length to expected bond length.
	 * This method automatically considers the correctly updated pi-electron count at both bond atoms
	 * after the bond order change for the lookup in the statistics bond length table.
	 * If we consider flipping double bonds along a chain of conjugated bonds, we may disable the
	 * automatic pi-change because atom's pi-electron count is not changed when shifting double bonds by
	 * one position.
	 * @param bond
	 * @param newOrder
	 * @param aromaticityChange 0:none; 1:aromatic; 2: non-aromatic
	 * @param noPiChangeAtom1 -1 or atom where PI change should be neglected
	 * @param noPiChangeAtom2 -1 or second atom where PI change should be neglected
	 * @return improvement regarding how much the bond length meets the new bond order's expected length
	 */
	private double calculateBondChangeScore(int bond, int newOrder, int aromaticityChange, int noPiChangeAtom1, int noPiChangeAtom2) {
		int oldOrder = mMol.getBondOrder(bond);
		if (oldOrder == newOrder)
			return 0.0;

		int atom1 = mMol.getBondAtom(0, bond);
		int atom2 = mMol.getBondAtom(1, bond);
		int atomicNo1 = mMol.getAtomicNo(atom1);
		int atomicNo2 = mMol.getAtomicNo(atom2);
		int conns1 = mMol.getConnAtoms(atom1);
		int conns2 = mMol.getConnAtoms(atom2);
		int pi1 = getAtomPi(atom1);
		int pi2 = getAtomPi(atom2);
		int deltaPi1 = (atom1 == noPiChangeAtom1 || atom1 == noPiChangeAtom2) ? 0 : newOrder - oldOrder;
		int deltaPi2 = (atom2 == noPiChangeAtom1 || atom2 == noPiChangeAtom2) ? 0 : newOrder - oldOrder;
		boolean oldDelocalized = mIsDelocalizedBond[bond] && pi1 != 0 && pi2 != 0;
		boolean newAromaticity = (aromaticityChange == 0) ? mIsAromaticBond[bond] : aromaticityChange == 1;
		boolean newDelocalized = (aromaticityChange == 0) ? mIsDelocalizedBond[bond] && pi1+deltaPi1 != 0 && pi2+deltaPi2 != 0 : aromaticityChange == 1;
		int oldBondIndex = BondLengthSet.getBondIndex(oldOrder, mIsAromaticBond[bond], oldDelocalized,
				atomicNo1, atomicNo2, pi1, pi2, conns1, conns2, false);
		int newBondIndex = BondLengthSet.getBondIndex(newOrder, newAromaticity, newDelocalized,
				atomicNo1, atomicNo2, pi1+deltaPi1, pi2+deltaPi2, conns1, conns2, false);
		if (oldBondIndex == -1 || newBondIndex == -1)
			return 0.0;

		double nitrogenPenalty = 0.0;
		if (Molecule.isAtomicNoElectronegative(atomicNo1) && conns1 == 3 && atom1 != noPiChangeAtom1 && atom1 != noPiChangeAtom2 && !hasNegativeNeighbour(atom1))
			nitrogenPenalty += (newOrder > oldOrder) ? HETERO_ATOM_CHARGE_PENALTY : -HETERO_ATOM_CHARGE_PENALTY;
		if (Molecule.isAtomicNoElectronegative(atomicNo2) && conns2 == 3 && atom2 != noPiChangeAtom1 && atom2 != noPiChangeAtom2 && !hasNegativeNeighbour(atom2))
			nitrogenPenalty += (newOrder > oldOrder) ? HETERO_ATOM_CHARGE_PENALTY : -HETERO_ATOM_CHARGE_PENALTY;

		double newExpectedLength = BondLengthSet.getBondLength(newBondIndex);
		double oldExpectedLength = BondLengthSet.getBondLength(oldBondIndex);
		double newDistance = (newOrder > oldOrder) ? mBondLength[bond] - newExpectedLength : newExpectedLength - mBondLength[bond];
		double oldDistance = (oldOrder > newOrder) ? mBondLength[bond] - oldExpectedLength : oldExpectedLength - mBondLength[bond];
		if ((newOrder > oldOrder) == (newExpectedLength > oldExpectedLength)) {
			newDistance = 0.0;
			oldDistance = 0.0;
		}

		// we multiply with a value close to 1.0 to slightly prefer bonds that are far off, if they go in the right direction
		double newPenalty = Math.max(0, newDistance) * (1.0 + newDistance);
		double oldPenalty = Math.max(0, oldDistance) * (1.0 + oldDistance);

if (DEBUG_OUTPUT) {
//if (bond==13 || bond==14) {
 String oldCrit = "oldCrit:"+oldOrder+" arom:"+mIsAromaticBond[bond]+" delo:"+mIsDelocalizedBond[bond]+" atomicNo:"+atomicNo1+","+atomicNo2
	+" pi:"+pi1+","+pi2+" conns:"+conns1+","+conns2;
 String newCrit = "newCrit:"+newOrder+" arom:"+newAromaticity+" delo:"+newDelocalized+" atomicNo:"+atomicNo1+","+atomicNo2
	+" pi:"+(pi1+deltaPi1)+","+(pi2+deltaPi2)+" conns:"+conns1+","+conns2;
 System.out.println("calcBondScore bond:" + bond + " len:" + DoubleFormat.toString(mBondLength[bond]) + " oldOrder:" + oldOrder + " newOrder:" + newOrder
	+ " oldLen:" + BondLengthSet.getBondLength(oldBondIndex) + " newLen:" + BondLengthSet.getBondLength(newBondIndex)
	+ " oldPen:" + DoubleFormat.toString(oldPenalty) + " newPen:" + DoubleFormat.toString(newPenalty) + " nitPen:" + nitrogenPenalty
	+ " score:" + DoubleFormat.toString(oldPenalty - newPenalty - nitrogenPenalty));
 System.out.println(" "+oldCrit+"; "+newCrit);
//}
}

		return oldPenalty - newPenalty - nitrogenPenalty;
	}

	private boolean hasNegativeNeighbour(int atom) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (mMol.getAtomCharge(mMol.getConnAtom(atom, i)) < 0)
				return true;
		return false;
	}

	/**
	 * In random order go through all non-aromatic, hitherto not assigned bonds and
	 * assign that bond order, which best explains the bond length considering for both
	 * bond atoms their neighbour counts and pi-electron counts.
	 * Since atom pi counts are changed by this procedure, and pi count has an influence
	 * on the predicted bond orders, we repeat the procedure until we have a steady state
	 * of predicted bond orders.
	 */
	private void assignNonAromaticBondOrders() {
		GraphWalker tautomerWalker = createTautomerWalker(false, 0.0, false, false);

		int[] randomAndBondIndex = new int[mMol.getBonds()];
		for (int i=0; i<MAX_ASSIGNMENT_ROUNDS; i++) {
			for (int bond=0; bond<mMol.getBonds(); bond++)
				randomAndBondIndex[bond] = bond + (mRandom.nextInt(0x00007FFF) << 16);
			Arrays.sort(randomAndBondIndex);

			boolean changedAny = false;
			for (int rabi : randomAndBondIndex) {
				int bond = rabi & 0x0000FFFF;
				if (!mBondIsFinal[bond]) {
					boolean changed = tryChangeBondOrder(bond);
					if (!changed)
						changed = tryChangeTautomer(mMol.getBondAtom(0, bond), mMol.getBondAtom(1, bond), bond, tautomerWalker)
							   || tryChangeTautomer(mMol.getBondAtom(1, bond), mMol.getBondAtom(0, bond), bond, tautomerWalker);
					changedAny |= changed;
				}
			}
			if (!changedAny)
				break;
		}
	}

	/**
	 * Updates hybridisation for the given atom and return whether previous values were changed
	 * @param atom
	 * @param hybridisation
	 * @return
	 */
	private boolean assignAtom(int atom, int hybridisation) {
		if (mHybridisation[atom] == hybridisation)
			return false;

		mHybridisation[atom] = hybridisation;

		return true;
	}

	/**
	 * Updates hybridisation for the given atom and return whether previous values were changed
	 * @param bond
	 * @param order
	 * @param isFinal
	 * @return
	 */
	private boolean assignBond(int bond, int order, boolean isFinal) {
		if (mBondIsFinal[bond] && order != mMol.getBondOrder(bond))
			System.out.println("$$$ WARNING: Assigning different order ("+mMol.getBondOrder(bond)+"->"+order+") to final bond("+bond+"): "+mMol.getName());

		mBondIsFinal[bond] |= isFinal;

		if (mMol.getBondOrder(bond) == order)
			return false;

		mMol.setBondType(bond, order==1 ? Molecule.cBondTypeSingle : order == 2 ? Molecule.cBondTypeDouble : Molecule.cBondTypeTriple);
		return true;
	}

	/**
	 * Considering current pi electron counts and neighbour counts for both bond atoms,
	 * lookup the typical bond length for single and double bonds. Return the order that
	 * better matches real length of this bond.
	 * @param bond
	 * @return most likely order
	 */
	private int guessBondOrderFromBondLength(int bond) {
		int[] atom = new int[2];
		int[] piCount = new int[2];
		int[] conns = new int[2];
		int[] freeValence = new int[2];
		double metalLigandBonus = 0.0;	// if we have a coordinating metal atom, we increase tendency to higher bond order

		for (int i=0; i<2; i++) {
			atom[i] = mMol.getBondAtom(i, bond);
			if (mHybridisation[atom[i]] == 3 && mMol.getAtomicNo(atom[i]) < 14)
				return 1;

			freeValence[i] = getFreeAtomValence(atom[i], bond, true);
			if (freeValence[i] == 0)
				return 1;

			conns[i] = mMol.getConnAtoms(atom[i]);
			for (int k=0; k<conns[i]; k++) {
				int connBond = mMol.getConnBond(atom[i], k);
				if (connBond != bond)
					piCount[i] += mMol.getBondOrder(connBond) - 1;
			}

			if (mMol.getAllConnAtomsPlusMetalBonds(atom[i]) > mMol.getAllConnAtoms(atom[i]))
				metalLigandBonus = METAL_LIGAND_BONUS;
		}

		int doubleBondIndex = BondLengthSet.getBondIndex(2, false, false,
				mMol.getAtomicNo(atom[0]), mMol.getAtomicNo(atom[1]), piCount[0]+1, piCount[1]+1, conns[0], conns[1], false);
		int singleBondIndex = BondLengthSet.getBondIndex(1, false, false,
				mMol.getAtomicNo(atom[0]), mMol.getAtomicNo(atom[1]), piCount[0], piCount[1], conns[0], conns[1], false);

		if (singleBondIndex == -1 && doubleBondIndex == -1)
			return 1;

		double singleBondLength = BondLengthSet.getBondLength(singleBondIndex);
		double doubleBondLength = BondLengthSet.getBondLength(doubleBondIndex);

		if (singleBondIndex == -1)
			singleBondLength = doubleBondLength / SINGLE_DOUBLE_BOND_LENGTH_FACTOR;
		if (doubleBondIndex == -1)
			doubleBondLength = SINGLE_DOUBLE_BOND_LENGTH_FACTOR * singleBondLength;

		if (freeValence[0] >= 2 && freeValence[1] >= 2
		 && (conns[0] == 1 || mHybridisation[atom[0]] == 1)
		 && (conns[1] == 1 || mHybridisation[atom[1]] == 1)) {
			int tripleBondIndex = BondLengthSet.getBondIndex(3, false, false,
					mMol.getAtomicNo(atom[0]), mMol.getAtomicNo(atom[1]), 2, 2, conns[0], conns[1], false);
			double tripleBondLength = (tripleBondIndex == -1) ?
					DOUBLE_TRIPLE_BOND_LENGTH_FACTOR * doubleBondLength : BondLengthSet.getBondLength(tripleBondIndex);
			if (mBondLength[bond] < (0.5 * (doubleBondLength + tripleBondLength)) + metalLigandBonus)
				return 3;
		}

		return mBondLength[bond] < (0.5 * (doubleBondLength + singleBondLength) + metalLigandBonus) ? 2 : 1;
	}

	private void correctAtomCharges() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			int occupiedValence = mMol.getOccupiedValence(atom);
			if (mMol.getAtomCharge(atom) == 0
			 && ((occupiedValence >= 4 && mMol.getAtomicNo(atom) == 7)
			  || (occupiedValence == 3 && mMol.getAtomicNo(atom) == 16))) {
				// Correct quarternary uncharged nitrogen, tertiary uncharged sulfur, by adding a charge and,
				// if possible, add a negative counter charge on a single bonded neighbour oxygen.
				mMol.setAtomCharge(atom, 1);
				if (mHybridisation[atom] == 1 && mMol.getAtomicNo(atom) == 7) {
					for (int i=0; i<mMol.getConnAtoms(atom); i++) {
						int connAtom = mMol.getConnAtom(atom, i);
						// isonitrile: add a negative charge to the carbon
						if (mMol.getAtomicNo(connAtom) == 6
						 && mMol.getBondOrder(mMol.getConnBond(atom, i)) == 3
						 && mMol.getConnAtoms(connAtom) == 1
						 && mMol.getAtomCharge(connAtom) == 0) {
							mMol.setAtomCharge(connAtom, -1);
							break;
						}
						// azide: add a negative charge to the nitrogen neighbour
						if (mMol.getAtomicNo(connAtom) == 7
						 && mMol.getBondOrder(mMol.getConnBond(atom, i)) == 2
						 && mMol.getConnAtoms(connAtom) == 1
						 && mMol.getAtomCharge(connAtom) == 0) {
							mMol.setAtomCharge(connAtom, -1);
							break;
						}
						// amine oxid: add a negative charge to the oxygen neighbour
						if (mMol.getAtomicNo(connAtom) == 8
						 && mMol.getBondOrder(mMol.getConnBond(atom, i)) == 1
						 && mMol.getConnAtoms(connAtom) == 1
						 && mMol.getAtomCharge(connAtom) == 0) {
							mMol.setAtomCharge(connAtom, -1);
							break;
						}
					}
				}
				else {
					boolean chargedOxygenFound = false;
					int unchargedOxygen = -1;
					for (int i=0; i<mMol.getConnAtoms(atom); i++) {
						int connAtom = mMol.getConnAtom(atom, i);
						if (mMol.getAtomicNo(connAtom) == 8
						 && mMol.getBondOrder(mMol.getConnBond(atom, i)) == 1
						 && mMol.getConnAtoms(connAtom) == 1) {
							if (mMol.getAtomCharge(connAtom) == 0)
								unchargedOxygen = connAtom;
							else
								chargedOxygenFound = true;
						}
					}
					if (!chargedOxygenFound && unchargedOxygen != -1)
						mMol.setAtomCharge(unchargedOxygen, -1);
				}
			}
		else if (mMol.isElectronegative(atom) && mMol.getAtomCharge(atom) == 1 && occupiedValence <= 3)
			mMol.setAtomCharge(atom, 0);
		}
	}

	private void correctExceededValences() {
		int[] freeValence = new int[mMol.getAtoms()];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			freeValence[atom] = getFreeAtomValence(atom, -1, false);

			// we also correct wrong cumulated double bonds here
			if (freeValence[atom] >= 0
			 && getAtomPi(atom) == 2
			 && mMol.getConnAtoms(atom) == 2
			 && mHybridisation[atom] != 1
			 && mMol.getAtomicNo(atom) <= 7)
				freeValence[atom] = -1;
		}

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			while (freeValence[atom] < 0) {
				int bestConnBond = -1;
				double bestDelta = -10.0;	// make sure to even handle much too short bonds
				for (int j=0; j<mMol.getConnAtoms(atom); j++) {
					int connBond = mMol.getConnBond(atom, j);
					if (mMol.getBondOrder(connBond) == 2) {
						int connAtom = mMol.getConnAtom(atom, j);
//						if (mMol.getLowestFreeValence(connAtom) < 0) {
//							bestConnBond = connBond;
//							break;
//						}
						double delta = mBondLength[connBond] - predictedNonAromaticDoubleBondLength(connBond);
						if (freeValence[connAtom] < 0)
							delta += 10;	// make sure to prefer exceeded-valence-neighbours

						if (mMol.getAtomicNo(connAtom) == 8 && mMol.getConnAtoms(connAtom) == 1 && getAtomPi(atom) == 2)
							delta -= ENOL_TAUTOMER_MALUS;

						if (bestDelta < delta) {
							bestDelta = delta;
							bestConnBond = connBond;
						}
					}
				}
				if (bestConnBond != -1) {
					// Override final bonds if necessary
					mMol.setBondOrder(bestConnBond, 1);
					mBondIsFinal[bestConnBond] = true;
					freeValence[mMol.getBondAtom(0, bestConnBond)]++;
					freeValence[mMol.getBondAtom(1, bestConnBond)]++;
				}
				else {
					break;
				}
			}
		}
	}

	private void convertObviousHydroxyIminesToAmides() {
		for (int oxygen=0; oxygen<mMol.getAtoms(); oxygen++) {
			if (mMol.getAtomicNo(oxygen) == 8
			 && mMol.getConnAtoms(oxygen) == 1) {
				int oxygenBond = mMol.getConnBond(oxygen, 0);
				if (!mBondIsFinal[oxygenBond]
				 && mMol.getBondType(oxygenBond) == Molecule.cBondTypeSingle) {
					int carbon = mMol.getConnAtom(oxygen, 0);
					if (!mIsAromaticAtom[carbon]) {
						Neighbour[] neighbour = getNeighbourAtoms(carbon, -1, 0, 0, 0, true);
						if (neighbour.length >= 2) {
							if (neighbour[0].getAtom() == oxygen) {
								for (int i=1; i<neighbour.length; i++) {
									if (!mBondIsFinal[neighbour[i].getBond()]
									 && mMol.getBondType(neighbour[i].getBond()) == Molecule.cBondTypeDouble) {
										assignBond(neighbour[i].getBond(), 1, true);
										assignBond(oxygenBond, 2, true);
										break;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	private void correctNonLinearTripleBonds() {
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondOrder(bond) == 3) {
				for (int i=0; i<2; i++) {
					int atom = mMol.getBondAtom(i, bond);
					if (mMol.getConnAtoms(atom) == 2
					 && mHybridisation[atom] != 1
					 && mMol.getAtomicNo(atom) <= 7) {
						mMol.setBondOrder(bond, 2);
						break;
					}
				}
			}
		}
	}

	private void correctMetalLigandCharges() {
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondType(bond) == Molecule.cBondTypeMetalLigand) {
				for (int i=0; i<2; i++) {
					int metalAtom = mMol.getBondAtom(i, bond);
					if (mMol.isMetalAtom(metalAtom)) {
						if (mMol.getAtomicNo(mMol.getBondAtom(1-i, bond)) == 6
						 && prefersCovalentCarbonBond(metalAtom)) {
							mMol.setBondType(bond, Molecule.cBondTypeSingle);
						}
						else {
							int atom = mMol.getBondAtom(1-i, bond);
							if (qualifiesAsNegativeMetalLigand(atom, metalAtom)) {
								mMol.setAtomCharge(atom, -1);
								oxydateMetalAtom(metalAtom);
							}
						}
					}
				}
			}
		}
	}

	private boolean prefersCovalentCarbonBond(int metalAtom) {
		int atomicNo = mMol.getAtomicNo(metalAtom);
		return atomicNo == 80	// Hg
			|| atomicNo == 82	// Pb
			|| atomicNo == 50	// Sn
			|| atomicNo == 48	// Cd
			|| atomicNo == 30	// Zn
			|| atomicNo == 13;	// Al
	}

	private boolean qualifiesAsNegativeMetalLigand(int atom, int metal) {
		if (mMol.getImplicitHydrogens(atom) == 0
		 || mMol.getAtomCharge(atom) != 0)
			return false;

		if (mMol.isHalogene(atom)
		 || AtomFunctionAnalyzer.isAcidicOxygen(mMol, atom))	// any electronegative atom with implicit hydrogen
			return true;

		if (implicitHydrogenCollidesWithMetalBond(atom, metal))
			return true;

		// cyanide
		if (mMol.getAtomicNo(atom) == 6
		 && mMol.getConnAtoms(atom) == 1
		 && mMol.getConnBondOrder(atom, 0) == 3
		 && mMol.getAtomicNo(mMol.getConnAtom(atom, 0)) == 7
		 && mMol.getConnAtoms(mMol.getConnAtom(atom, 0)) == 1) {
			return true;
		}

		return false;
	}

	private void oxydateMetalAtom(int metalAtom) {
		byte[] oxidationState = Molecule.cCommonOxidationState[mMol.getAtomicNo(metalAtom)];
		if (oxidationState != null) {
			int maxOxidationState = oxidationState[oxidationState.length-1];
			if (mMol.getAtomCharge(metalAtom) < maxOxidationState)
				mMol.setAtomCharge(metalAtom, mMol.getAtomCharge(metalAtom)+1);
		}
	}

	private boolean implicitHydrogenCollidesWithMetalBond(int atom, int metal) {
		if (mMol.getConnAtoms(atom) != 0) {
			Coordinates metalDirection = mMol.getAtomCoordinates(metal).subC(mMol.getAtomCoordinates(atom));
			Coordinates neighbourDirection = new Coordinates();
			for (int i=0; i<mMol.getConnAtoms(atom); i++)
				neighbourDirection.add(mMol.getAtomCoordinates(atom).subC(mMol.getAtomCoordinates(mMol.getConnAtom(atom, i))));

			return metalDirection.getAngle(neighbourDirection) < METAL_HYDROGEN_COLLISION_ANGLE;
		}
		return false;
	}

	private double predictedNonAromaticDoubleBondLength(int bond) {
		int[] atom = new int[2];
		int[] piCount = new int[2];
		int[] conns = new int[2];

		for (int i=0; i<2; i++) {
			atom[i] = mMol.getBondAtom(i, bond);
			conns[i] = mMol.getConnAtoms(atom[i]);
			for (int k=0; k<conns[i]; k++) {
				int connBond = mMol.getConnBond(atom[i], k);
				if (connBond != bond)
					piCount[i] += mMol.getBondOrder(connBond) - 1;
			}
		}

		int index = BondLengthSet.getBondIndex(2, false, false,
				mMol.getAtomicNo(atom[0]), mMol.getAtomicNo(atom[1]), piCount[0], piCount[1], conns[0], conns[1], false);
		return index == -1 ? 0 : BondLengthSet.getBondLength(index);
	}

	/**
	 * Calculates the atom's free valence optionally neglecting neglectBond and optionally considering final bonds only.
	 * In case of an uncharged nitrogen, we add 1 to the returned free valence considering that we can add a charge later.
	 * @param atom
	 * @param neglectBond -1 or bond to count as single bond
	 * @param finalBondsOnly if true, then non-final bonds are counted as single bonds
	 * @return
	 */
	private int getFreeAtomValence(int atom, int neglectBond, boolean finalBondsOnly) {
		int occupiedValence = mMol.getConnAtoms(atom);
		for (int j=0; j<mMol.getConnAtoms(atom); j++) {
			int connBond = mMol.getConnBond(atom, j);
			if (connBond != neglectBond && (!finalBondsOnly || mBondIsFinal[connBond]))
				occupiedValence += mMol.getBondOrder(connBond) - 1;
		}
		byte[] valenceList = Molecule.getAllowedValences(mMol.getAtomicNo(atom));
		int j=0;
		while ((j<valenceList.length-1) && (occupiedValence > valenceList[j]))
			j++;
		int potentialChargeIncrease = (mMol.getAtomicNo(atom) == 7) ? 1 : 0;
		return valenceList[j] - occupiedValence + potentialChargeIncrease;
	}

	/**
	 * Determines the number of neighbour atoms that match the given conditions. Then, sorts matching neighbour atoms
	 * based on ascending bond lengths and returns them.
	 * @param rootAtom Root atom at which to look for neighbours
	 * @param connAtomicNo 0 or -1 (any electronegative) or atomic number of neighbour atom directly connected to the root atom
	 * @param connNeighbours 0 or neighbour count of the neighbour atom that is directly to the root atom
	 * @param farAtomicNo 0 or atomic number of a required atom being connected to the neighbour atom
	 * @param farNeighbours 0 or neighbour count of the required neighbour being connected to the neighbour atom
	 * @param subtractExpectedBondLengths if true, then neighbours are sorted by deviation from the expected bond length
	 * @return matching neighbour atoms sorted by increasing bond length, or deviation from an expected bond length
	 * and -1 otherwise.
	 */
	private Neighbour[] getNeighbourAtoms(int rootAtom, int connAtomicNo, int connNeighbours, int farAtomicNo, int farNeighbours, boolean subtractExpectedBondLengths) {
		mNeighbourList.clear();
		loop:
		for (int i=0; i<mMol.getConnAtoms(rootAtom); i++) {
			int connAtom = mMol.getConnAtom(rootAtom, i);

			// Checking if the atom that is connected to rootAtom
			// is of the requested atomicNo and has the requested neighbour count
			if (connAtomicNo > 0 && mMol.getAtomicNo(connAtom) != connAtomicNo) continue;
			if (connAtomicNo == -1 && !mMol.isElectronegative(connAtom)) continue;
			if (connNeighbours > 0 && mMol.getConnAtoms(connAtom) != connNeighbours) continue;

			// Checking if another far atom that is connected to a direct neighbour atom
			// is of the requested chemical type and neighbour count
			if (farAtomicNo != 0 || farNeighbours != 0) {
				for (int j=0; j<mMol.getAllConnAtoms(connAtom); j++) {
					int otherAtm = mMol.getConnAtom(connAtom, j);
					if (otherAtm == rootAtom) continue loop;
					if (farAtomicNo != 0 && mMol.getAtomicNo(otherAtm) != farAtomicNo) continue loop;
					if (farNeighbours != 0 && mMol.getAllConnAtoms(otherAtm) != farNeighbours) continue loop;
				}
			}

			int connBond = mMol.getConnBond(rootAtom, i);
			mNeighbourList.add(new Neighbour(rootAtom, connAtom, connBond, mBondLength[connBond]));
		}

		if (subtractExpectedBondLengths)
			calculateBondLengthDeviationsForNeighbourList();

		Neighbour[] neighbours = mNeighbourList.toArray(new Neighbour[0]);
		Arrays.sort(neighbours);
		return neighbours;
	}

	/**
	 * If mNeighbourList contains neighbours of a different kind and if we want to sort the list using bond lengths,
	 * then this method should be called once to convert absolute neighbour bond lengths to relative ones that allows
	 * meaningful sorting based on the bond length deviation from the expected bonds lengths.
	 * This method subtracts for every neighbour the expected single-bond-length from the real lengths
	 */
	private void calculateBondLengthDeviationsForNeighbourList() {
		double[] expectedBondLength = new double[mNeighbourList.size()];
		for (int i=0; i<mNeighbourList.size(); i++) {
			int root = mNeighbourList.get(i).getRootAtom();
			int atom = mNeighbourList.get(i).getAtom();
			int index = BondLengthSet.getBondIndex(1, false,
					false, mMol.getAtomicNo(root), mMol.getAtomicNo(atom),
					0, 0, mMol.getConnAtoms(root), mMol.getConnAtoms(atom), false);
			if (index == -1) {
				expectedBondLength = null;
				break;
			}
			expectedBondLength[i] = BondLengthSet.getBondLength(index);
		}
		if (expectedBondLength != null)
			for (int i=0; i<mNeighbourList.size(); i++)
				mNeighbourList.get(i).subtractExpectedBondLength(expectedBondLength[i]);
	}

}

class Neighbour implements Comparable<Neighbour> {
	private double mBondLengthDeviation;
	private final int mRootAtom,mAtom,mBond;

	public Neighbour(int rootAtom, int atom, int bond, double bondLength) {
		mRootAtom = rootAtom;
		mAtom = atom;
		mBond = bond;
		mBondLengthDeviation = bondLength;
	}

	public int getAtom() {
		return mAtom;
	}

	public int getBond() {
		return mBond;
	}

	public int getRootAtom() {
		return mRootAtom;
	}

	public double getBondLengthDeviation() {
		return mBondLengthDeviation;
	}

	// Use this method if neighbours are of different kinds. Then, instead of sorting based on bond lengths,
	// we need to sort based on deviations from expected bond lengths.
	public void subtractExpectedBondLength(double bondLength) {
		mBondLengthDeviation -= bondLength;
	}

	@Override
	public int compareTo(Neighbour o) {
		return Double.compare(mBondLengthDeviation, o.mBondLengthDeviation);
	}
}
