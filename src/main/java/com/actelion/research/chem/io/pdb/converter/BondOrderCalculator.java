package com.actelion.research.chem.io.pdb.converter;

import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.BondLengthSet;
import com.actelion.research.chem.conf.TorsionDB;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * This class calculates bond orders from atom coordinates only.
 */
public class BondOrderCalculator {
	private static final double HALF_PI = 0.5 * Math.PI;
	private static final double SP_BOND_ANGLE_LIMIT = 0.85 * Math.PI;
	private static final double SP2_ANGLE_SUM = 2.000 * Math.PI;
	private static final double SP3_ANGLE_SUM = 1.833 * Math.PI;
	private static final double OUT_OF_PLANE_TORSION_LIMIT = 0.20 * Math.PI;	// relaxed setting, because PDB errors are sometimes substantial
	private static final double SINGLE_DOUBLE_BOND_LENGTH_FACTOR = 0.92;
	private static final double DOUBLE_TRIPLE_BOND_LENGTH_FACTOR = 0.88;
	private static final double AROMATIC_BOND_LENGTH_TOLERANCE = 0.1;
	private static final double PYRIDINONE_CN_CUT_OFF = 1.370;
	private static final double PYRIDINONE_CO_CUT_OFF = 1.285;
	private static final double QUINONE_CO_CUT_OFF = 1.310;
	private static final double QUINONE_CC_CUT_OFF = 1.420;
	private static final int MAX_ASSIGNMENT_ROUNDS = 8;

	private final Molecule3D mMol;
	private double[] mBondLength;
	private final int[] mHybridisation;
	private final boolean[] mIsFlatBond,mBondIsFinal;
	private final Random mRandom;
	private final ArrayList<Neighbour> mNeighbourList;

	/**
	 * This class calculates bond orders from atom coordinates only.
	 * Neither bonds nor hydrogen atoms are required, but may be present.
	 * If the molecule contains bonds, then these are expected to be metal and single bonds only.
	 * For all atoms the molecule must contain 3D-coordinates, e.g. taken from a PDB or MMCIF file.
	 * @param mol
	 */
	public BondOrderCalculator(Molecule3D mol) throws Exception {
		mMol = mol;
		mRandom = new Random();
		mNeighbourList = new ArrayList<>();

		if (mMol.getAllBonds() == 0)
			BondsCalculator.createBonds(mMol, true, null);

		mMol.ensureHelperArrays(Molecule.cHelperRings);
		mHybridisation = new int[mMol.getAtoms()];
		mIsFlatBond = new boolean[mMol.getBonds()];
		mBondIsFinal = new boolean[mMol.getBonds()];
	}

	public void calculateBondOrders() {
		if (mMol.getBonds() == 0)
			return;

		// Don't touch metal bonds!
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mMol.getBondType(bond) == Molecule.cBondTypeMetalLigand)
				mBondIsFinal[bond] = true;

		// Calculate lengths of all bonds
		mBondLength = new double[mMol.getBonds()];
		for (int bond=0; bond<mMol.getBonds(); bond++)
			mBondLength[bond] = mMol.getCoordinates(mMol.getBondAtom(0, bond)).subC(mMol.getCoordinates(mMol.getBondAtom(1, bond))).getLength();

		determineObviousHybridisation();

		// Determine all but terminal SP1 by bond angles
		// Determine out-of-plane torsions for all bonds
		boolean[] isOutOfPlaneTorsionBond = determineOutOfPlaneTorsions();
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (isOutOfPlaneTorsionBond[bond])
				assignBond(bond, 1, true);

		boolean[] isAromaticBond = determineAromaticBonds(isOutOfPlaneTorsionBond);
		boolean success = new AromaticityResolver(mMol, mBondLength).locateDelocalizedDoubleBonds(isAromaticBond.clone(), true, false);
		if (!success)
			System.out.println("WARNING: Assignment of aromatic ring bonds failed.");

		// declare all assigned aromatic bonds as final
		boolean[] isAromaticAtom = new boolean[mMol.getAtoms()];
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (isAromaticBond[bond]) {
				mBondIsFinal[bond] = true;
				isAromaticAtom[mMol.getBondAtom(0, bond)] = true;
				isAromaticAtom[mMol.getBondAtom(1, bond)] = true;
			}
		}

		// Determine final bond orders for known groups outside of aromatic rings
		determineKnownGroups(isAromaticAtom);

		assignNonAromaticBondOrders();

		correctExceededValences();

		assignForgottenDoubleBonds(isOutOfPlaneTorsionBond);

		correctMetalLigandCharges();

//		try { mMol.canonizeCharge(true, true); } catch(Exception e) {}
	}

	private void determineObviousHybridisation() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getConnAtoms(atom) == 2) {
				if (GeometryCalculator.getAngle(mMol, mMol.getConnAtom(atom, 0), atom, mMol.getConnAtom(atom, 1)) > SP_BOND_ANGLE_LIMIT)
					assignAtom(atom, 1);
			}
			else if (mMol.getAllConnAtoms(atom) == 3) {
				double angleSum = GeometryCalculator.getAngle(mMol, mMol.getConnAtom(atom, 0), atom, mMol.getConnAtom(atom, 1))
								+ GeometryCalculator.getAngle(mMol, mMol.getConnAtom(atom, 0), atom, mMol.getConnAtom(atom, 2))
								+ GeometryCalculator.getAngle(mMol, mMol.getConnAtom(atom, 1), atom, mMol.getConnAtom(atom, 2));
				double planarity = (angleSum - SP3_ANGLE_SUM) / (SP2_ANGLE_SUM - SP3_ANGLE_SUM);

/*				// for carbon atoms we also consider relative bond length of a potential double bond as 50% or the criterion
				if (mMol.getAtomicNo(atom) == 6) {
					double doubleBondLikelyhood = -1.0;
					int conns = mMol.getConnAtoms(atom);
					for (int i=0; i<conns; i++) {
						int connAtom = mMol.getConnAtom(atom, i);
						int connConns = mMol.getConnAtoms(connAtom);
						int doubleBondIndex = BondLengthSet.getBondIndex(2, false, false,
								mMol.getAtomicNo(atom), mMol.getAtomicNo(connAtom), 1, 1, conns, connConns, false);
						int singleBondIndex = BondLengthSet.getBondIndex(1, false, false,
								mMol.getAtomicNo(atom), mMol.getAtomicNo(connAtom), 0, 0, conns, connConns, false);
						if (singleBondIndex != -1 && doubleBondIndex != -1) {
							int connBond = mMol.getConnBond(atom, i);
							double singleBondLength = BondLengthSet.getBondLength(singleBondIndex);
							double doubleBondLength = BondLengthSet.getBondLength(doubleBondIndex);
							if (singleBondLength > doubleBondLength) {
								double likelyhood = (singleBondLength - mBondLength[connBond]) / (singleBondLength - doubleBondLength);
								if (doubleBondLikelyhood < likelyhood)
									doubleBondLikelyhood = likelyhood;
							}
						}
					}
					if (doubleBondLikelyhood != -1.0)
						planarity = (planarity + doubleBondLikelyhood) / 2.0;
				}*/

				if (planarity > 0.80)
					assignAtom(atom, 2);
				else if (planarity < 0.5)
					assignAtom(atom, 3);
			}
			else if (mMol.getAllConnAtoms(atom) == 4) {
				assignAtom(atom, 3);
			}
		}
	}

	/**
	 * Determine all but terminal SP1 by bond angles and determine out-of-plane torsions for all bonds.
	 * @return
	 */
	private boolean[] determineOutOfPlaneTorsions() {
		boolean[] isOutOfPlaneTorsionBond = new boolean[mMol.getBonds()];
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

			if (conns1 > 1 && conns2 > 1) {
				for (int i=0; i<conns1 && !isOutOfPlaneTorsionBond[bond]; i++) {
					tAtom[0] = mMol.getConnAtom(tAtom[1], i);
					if (tAtom[0] != tAtom[2]) {
						for (int j=0; j<conns2 && !isOutOfPlaneTorsionBond[bond]; j++) {
							tAtom[3] = mMol.getConnAtom(tAtom[2], j);
							if (tAtom[3] != tAtom[1]) {
								double torsion = TorsionDB.calculateTorsionExtended(mMol, tAtom);
								isOutOfPlaneTorsionBond[bond] = (HALF_PI - Math.abs(Math.abs(torsion) - HALF_PI)) > OUT_OF_PLANE_TORSION_LIMIT;
							}
						}
					}
				}

				mIsFlatBond[bond] = conns1 + conns2 >= 5 & !isOutOfPlaneTorsionBond[bond];
			}
		}

		return isOutOfPlaneTorsionBond;
	}

	private void determineKnownGroups(boolean[] isAromaticAtom ) {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (isAromaticAtom[atom])
				continue;

			// X-C(=X)-X
			if (mMol.getAtomicNo(atom) == 6 && mMol.getConnAtoms(atom) == 3 && mHybridisation[atom] == 2) {
				Neighbour[] electronegative = getNeighbourAtoms(atom, -1, 0, 0, 0, true);
				if (electronegative.length == 3)
					for (int i=0; i<electronegative.length; i++)
						assignBond(electronegative[i].getBond(), i == 0 ? 2 : 1, true);
			}

			// any-C(=O)-O
			if (mMol.getAtomicNo(atom) == 6 && mMol.getConnAtoms(atom) == 3 && mHybridisation[atom] == 2) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 2)
					for (int i=0; i<oxygen.length; i++)
						assignBond(oxygen[i].getBond(), i == 0 ? 2 : 1, true);
			}

				// any-N(=O)-O
			if (mMol.getAtomicNo(atom) == 7 && mMol.getConnAtoms(atom) == 3) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 2) {
					mMol.setAtomCharge(atom, 1);
					for (int i=0; i<oxygen.length; i++) {
						if (i == 1)
							mMol.setAtomCharge(oxygen[i].getAtom(), -1);
						assignBond(oxygen[i].getBond(), i == 0 ? 2 : 1, true);
					}
				}
			}

			// any-P(=O)(-any)-any
			if (mMol.getAtomicNo(atom) == 15 && mMol.getConnAtoms(atom) == 4) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 1)
					for (int i=0; i<oxygen.length; i++)
						assignBond(oxygen[i].getBond(), i == 0 ? 2 : 1, true);
			}

			// any-S(=O)-any
			if (mMol.getAtomicNo(atom) == 16 && mMol.getConnAtoms(atom) == 3) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 1)
					for (int i=0; i<oxygen.length; i++)
						assignBond(oxygen[i].getBond(), i == 0 ? 2 : 1, true);
			}

			// any-S(=O)(=O)-any
			if (mMol.getAtomicNo(atom) == 16 && mMol.getConnAtoms(atom) == 4) {
				Neighbour[] oxygen = getNeighbourAtoms(atom, 8, 1, 0, 0, false);
				if (oxygen.length >= 2)
					for (int i=0; i<oxygen.length; i++)
						assignBond(oxygen[i].getBond(), i < 2 ? 2 : 1, true);
			}
		}
	}

	private boolean[] determineAromaticBonds(boolean[] isOutOfPlaneTorsionBond) {
		RingCollection ringSet = mMol.getRingSet();
		boolean[] isAromaticRing = new boolean[ringSet.getSize()];
		Arrays.fill(isAromaticRing, true);	// default

		boolean[] isDelocalizedBond = new boolean[mMol.getBonds()];

		for (int r=0; r<ringSet.getSize(); r++) {
			for (int bond : ringSet.getRingBonds(r)) {
				if (isOutOfPlaneTorsionBond[bond]) {
					isAromaticRing[r] = false;
					break;
				}
			}
			if (!isAromaticRing[r])
				continue;

			int[] ringAtom = ringSet.getRingAtoms(r);

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
			}

			if (!isAromaticRing[r])
				continue;

			if (ringAtom.length != 3) {    // 3-rings are always flat
				double angleSum = 0.0;
				for (int i = 0; i<ringAtom.length; i++)
					angleSum += GeometryCalculator.getAngle(mMol, ringAtom[(i + 1) % ringAtom.length], ringAtom[i], ringAtom[(i + ringAtom.length - 1) % ringAtom.length]);

				double angleSumDeviation = Math.abs(angleSum - Math.PI * (ringAtom.length - 2));
				if (angleSumDeviation > ((ringAtom.length == 4) ? Math.PI / 20
						: (ringAtom.length == 5) ? Math.PI / 30    // tiny margin: 108 degrees is flat 5-ring; 109.5 degrees is perfect sp3
						: (ringAtom.length == 6) ? Math.PI / 16 : Math.PI / 15)) {
					isAromaticRing[r] = false;
					continue;
				}
			}

			// Change hydroxy-imine to amide tautomer where bond lengths suggest that within aromatic 6-rings
//			if (ringAtom.length == 6) {
				for (int i=0; i<ringAtom.length; i++) {
					int[] neighbourAtom = ringNeighbourAtoms(ringAtom, i);
					int atom = ringAtom[i];
					int oxo = -1;
					int oxoBond = -1;
					for (int j=0; j<mMol.getConnAtoms(atom); j++) {
						int connAtom = mMol.getConnAtom(atom, j);
						if (connAtom != neighbourAtom[0] && connAtom != neighbourAtom[1] && mMol.getAtomicNo(connAtom) == 8 && mMol.getConnAtoms(connAtom) == 1) {
							oxo = connAtom;
							oxoBond = mMol.getConnBond(atom, j);
						}
					}
					if (oxo != -1) {
						int[] neighbourBond = ringNeighbourBonds(ringSet.getRingBonds(r), i);
						for (int j=0; j<2; j++) {
//							if (mMol.getAtomicNo(neighbourAtom[j]) == 6
//							 && (QUINONE_CO_CUT_OFF - mBondLength[oxoBond] + mBondLength[neighbourBond[j]] - QUINONE_CC_CUT_OFF > 0)) {
//								assignBond(oxoBond, 2, true);
//								isAromaticRing[r] = false;
//								break;
//							}
							if (mMol.getAtomicNo(neighbourAtom[j]) == 7
									&& (PYRIDINONE_CO_CUT_OFF - mBondLength[oxoBond] + mBondLength[neighbourBond[j]] - PYRIDINONE_CN_CUT_OFF > 0)) {
								assignBond(oxoBond, 2, true);
								isAromaticRing[r] = false;
								break;
							}
						}
					}
				}
//			}

			if (ringSet.getRingSize(r) == 6)
				for (int bond : ringSet.getRingBonds(r))
					isDelocalizedBond[bond] = true;
		}

		// flat 5-rings may not be aromatic; we check bond lengths and set ring as non-aromatic in case of too long bond
		for (int r=0; r<ringSet.getSize(); r++) {
			int[] ringBond = ringSet.getRingBonds(r);
			if (ringBond.length == 5) {
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
							if (!mIsFlatBond[bond] && mBondLength[bond] - bondLength > AROMATIC_BOND_LENGTH_TOLERANCE) {
								isAromaticRing[r] = false;
								break;
							}
						}
					}
				}
			}
		}

		boolean[] isAromaticBond = new boolean[mMol.getBonds()];
		for (int r=0; r<mMol.getRingSet().getSize(); r++) {
			if (isAromaticRing[r]) {
				for (int atom : ringSet.getRingAtoms(r))
					if (!isSPSeTe(atom) || mMol.getConnAtoms(atom) == 2)
						assignAtom(atom, 2);

				for (int bond : mMol.getRingSet().getRingBonds(r))
					isAromaticBond[bond] = true;
			}
		}

		return isAromaticBond;
	}

	private int[] ringNeighbourAtoms(int[] ringAtom, int i) {
		int[] atom = new int[2];
		atom[0] = ringAtom[i==0 ? ringAtom.length-1 : i-1];
		atom[1] = ringAtom[i==ringAtom.length-1 ? 0 : i+1];
		return atom;
	}

	private int[] ringNeighbourBonds(int[] ringBond, int i) {
		int[] bond = new int[2];
		bond[0] = ringBond[i==0 ? ringBond.length-1 : i-1];
		bond[1] = ringBond[i];
		return bond;
	}

	private boolean isSPSeTe(int atom) {
		return mMol.getAtomicNo(atom) == 15
			|| mMol.getAtomicNo(atom) == 16
			|| mMol.getAtomicNo(atom) == 34
			|| mMol.getAtomicNo(atom) == 52;
	}

	/**
	 * If we have sp2 carbons without pi-electrons after assigning all double bonds to match bond lengths,
	 * then we probably have an issue with bonds lengths. We check for potential neighbours
	 * @param isOutOfPlaneTorsionBond
	 */
	private void assignForgottenDoubleBonds(boolean[] isOutOfPlaneTorsionBond) {
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);	// we need pi
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mHybridisation[atom] == 2 && mMol.getAtomicNo(atom) == 6 && mMol.getAtomPi(atom) == 0) {
				mNeighbourList.clear();
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					int connBond = mMol.getConnBond(atom, i);
					if (!mBondIsFinal[connBond]
					 && !isOutOfPlaneTorsionBond[connBond]) {
						int connAtom = mMol.getConnAtom(atom, i);
						if (mHybridisation[connAtom] != 3 && mMol.getAtomPi(connAtom) == 0)
							mNeighbourList.add(new Neighbour(connAtom, connBond, mBondLength[connBond]));
					}
				}

				if (!mNeighbourList.isEmpty()) {
					calculateBondLengthsForNeighbourList(atom);
					Neighbour[] neighbour = mNeighbourList.toArray(new Neighbour[0]);
					Arrays.sort(neighbour);
					assignBond(neighbour[0].getBond(), 2, true);
				}
			}
		}
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
		int[] randomAndBondIndex = new int[mMol.getBonds()];
		for (int i=0; i<MAX_ASSIGNMENT_ROUNDS; i++) {
			for (int bond=0; bond<mMol.getBonds(); bond++)
				randomAndBondIndex[bond] = bond + (mRandom.nextInt(0x00007FFF) << 16);
			Arrays.sort(randomAndBondIndex);

			boolean changed = false;
			for (int rabi : randomAndBondIndex) {
				int bond = rabi & 0x0000FFFF;
				if (!mBondIsFinal[bond])
					changed |= assignBond(bond, guessBondOrderFromBondLength(bond), false);
			}
			if (!changed)
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
		if (mBondIsFinal[bond])
			System.out.println("WARNING: bond assignment change after final: "+bond + ": from "+mMol.getBondOrder(bond) + " to "+order);

		mBondIsFinal[bond] |= isFinal;

		if (mMol.getBondOrder(bond) == order)
			return false;

		mMol.setBondType(bond, order==1 ? Molecule.cBondTypeSingle : order == 2 ? Molecule.cBondTypeDouble : Molecule.cBondTypeTriple);
		return true;
	}

	/**
	 * Considering current pi electron counts and neighbour counts for both bond atoms,
	 * lookup the typical bond length for single and double bonds. Return the order that
	 * better matchtes real length of this bond.
	 * @param bond
	 * @return most likely order
	 */
	private int guessBondOrderFromBondLength(int bond) {
		int[] atom = new int[2];
		int[] piCount = new int[2];
		int[] conns = new int[2];
		int[] freeValence = new int[2];

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
			if (mBondLength[bond] < (0.5 * (doubleBondLength + tripleBondLength)))
				return 3;
		}

		return mBondLength[bond] < (0.5 * (doubleBondLength + singleBondLength)) ? 2 : 1;
	}

	private void correctExceededValences() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			int freeValence = mMol.getLowestFreeValence(atom);
			if (freeValence < 0) {
				// Correct quarternary uncharged nitrogen by adding a charge and,
				// if possible, add a negative counter charge on a single bonded neighbour oxygen.
				if (mMol.getAtomicNo(atom) == 7 && mMol.getAtomCharge(atom) == 0) {
					mMol.setAtomCharge(atom, 1);
					freeValence++;
					for (int i=0; i<mMol.getConnAtoms(atom); i++) {
						int connAtom = mMol.getConnAtom(atom, i);
						if (mMol.getAtomicNo(connAtom) == 8
						 && mMol.getBondOrder(mMol.getConnBond(atom, i)) == 1
						 && mMol.getConnAtoms(connAtom) == 1
						 && mMol.getAtomCharge(connAtom) == 0) {
							mMol.setAtomCharge(connAtom, -1);
							break;
						}
					}
				}
			}
			for (int i=0; i<-freeValence; i++) {
				int bestConnBond = -1;
				double bestDelta = -10.0;	// make sure to even handle much too short bonds
				for (int j=0; j<mMol.getConnAtoms(atom); j++) {
					int connBond = mMol.getConnBond(atom, j);
					if (mMol.getBondOrder(connBond) == 2) {
						int connAtom = mMol.getConnAtom(atom, j);
						if (mMol.getLowestFreeValence(connAtom) < 0) {
							bestConnBond = connBond;
							break;
						}
						double delta = mBondLength[connBond] - predictedNonAromaticBondLength(connBond);
						if (bestDelta < delta) {
							bestDelta = delta;
							bestConnBond = connBond;
						}
					}
				}
				if (bestConnBond != -1) {
					assignBond(bestConnBond, 1, true);
				}
			}
		}

	}

	private void correctMetalLigandCharges() {
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondType(bond) == Molecule.cBondTypeMetalLigand) {
				for (int i=0; i<2; i++) {
					int metal = mMol.getBondAtom(i, bond);
					if (mMol.isMetalAtom(metal)) {
						int atom = mMol.getBondAtom(1-i, bond);
						if (qualifiesAsNegativeMetalLigand(atom, metal))
							mMol.setAtomCharge(atom,-1);
					}
				}
			}
		}
	}

	private boolean qualifiesAsNegativeMetalLigand(int atom, int metal) {
		if (mMol.getImplicitHydrogens(atom) == 0
		 || mMol.getAtomCharge(atom) != 0)
			return false;

		if (AtomFunctionAnalyzer.isAcidicOxygen(mMol, atom))
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

	private boolean implicitHydrogenCollidesWithMetalBond(int atom, int metal) {
		if (mMol.getConnAtoms(atom) != 0) {
			Coordinates metalDirection = mMol.getCoordinates(metal).subC(mMol.getCoordinates(atom));
			Coordinates neighbourDirection = new Coordinates();
			for (int i=0; i<mMol.getConnAtoms(atom); i++)
				neighbourDirection.add(mMol.getCoordinates(atom).subC(mMol.getCoordinates(mMol.getConnAtom(atom, i))));

			return metalDirection.getAngle(neighbourDirection) < 0.25 * Math.PI;
		}
		return false;
	}

	private double predictedNonAromaticBondLength(int bond) {
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
		return valenceList[j] - occupiedValence;
	}

	/**
	 * Determines the number of neighbour atoms that match the given conditions. Then, sorts matching neighbour atoms
	 * based on ascending bond lengths and returns them.
	 * @param rootAtom Root atom at which to look for neighbours
	 * @param connAtomicNo 0 or -1 (any electronegative) or atomic number of neighbour atom directly connected to the root atom
	 * @param connNeighbours 0 or neighbour count of the neighbour atom that is directly to the root atom
	 * @param farAtomicNo 0 or atomic number of a required atom being connected to the neighbour atom
	 * @param farNeighbours 0 or neighbour count of the required neighbour being connected to the neighbour atom
	 * @param useBondLengthStatistics if true, then neighbours are sorted by deviation from the expected bond length
	 * @return matching neighbour atoms sorted by increasing bond length, or deviation from an expected bond length
	 * and -1 otherwise.
	 */
	private Neighbour[] getNeighbourAtoms(int rootAtom, int connAtomicNo, int connNeighbours, int farAtomicNo, int farNeighbours, boolean useBondLengthStatistics) {
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
			mNeighbourList.add(new Neighbour(connAtom, connBond, mBondLength[connBond]));
		}

		if (useBondLengthStatistics)
			calculateBondLengthsForNeighbourList(rootAtom);

		Neighbour[] neighbours = mNeighbourList.toArray(new Neighbour[0]);
		Arrays.sort(neighbours);
		return neighbours;
	}

	/**
	 * If mNeighbourList contains neighbours of a different kind and if we want to sort the list using bond lengths,
	 * then this method should be called once to convert absolute neighbour bond lengths to relative ones that allows
	 * meaningful sorting based on the bond length deviation from the expected bonds lengths.
	 * @param rootAtom
	 */
	private void calculateBondLengthsForNeighbourList(int rootAtom) {
		double[] expectedBondLength = new double[mNeighbourList.size()];
		for (int i=0; i<mNeighbourList.size(); i++) {
			int index = BondLengthSet.getBondIndex(1, false,
					false, mMol.getAtomicNo(rootAtom), mMol.getAtomicNo(mNeighbourList.get(i).getAtom()),
					0, 0, mMol.getConnAtoms(rootAtom), mMol.getConnAtoms(mNeighbourList.get(i).getAtom()), false);
			if (index == -1) {
				expectedBondLength = null;
				break;
			}
			expectedBondLength[i] = BondLengthSet.getBondLength(index);
		}
		if (expectedBondLength != null)
			for (int i=0; i<mNeighbourList.size(); i++)
				mNeighbourList.get(i).setExpectedBondLength(expectedBondLength[i]);
	}

}

class Neighbour implements Comparable<Neighbour> {
	private double mBondLength;
	private final int mAtom,mBond;

	public Neighbour(int atom, int bond, double bondLength) {
		mAtom = atom;
		mBond = bond;
		mBondLength = bondLength;
	}

	public int getAtom() {
		return mAtom;
	}

	public int getBond() {
		return mBond;
	}

	// Use this method if neighbours are of different kinds. Then, instead of sorting based on bond lengths,
	// we need to sort based on deviations from expected bond lengths.
	public void setExpectedBondLength(double bondLength) {
		mBondLength -= bondLength;
	}

	@Override
	public int compareTo(Neighbour o) {
		return Double.compare(mBondLength, o.mBondLength);
	}
}
