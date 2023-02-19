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

package org.openmolecules.chem.conf.so;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.BondAngleSet;
import com.actelion.research.chem.conf.BondLengthSet;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.VDWRadii;
import com.actelion.research.util.DoubleFormat;

import java.util.ArrayList;

public class DistanceRule extends ConformationRule {
	private static final double VDW_RADIUS_CORRECTION = 1.02;   // We increase reported VDW radii to increase Lennard-Jones energies

	private static final int PRIORITY_ONE_BOND = 10;    // We also use priorities to distinguish cases
	private static final int PRIORITY_TWO_BONDS = 5;
	private static final int PRIORITY_THREE_BONDS = 3;
	private static final int PRIORITY_FOUR_AND_MORE_BONDS = 1;
	private static final int PRIORITY_DISCONNECTED = 0;

	private double[] mDistance;
	private int[] mNotList;
	private int mPriority;

	/**
	 * Creates a dummy rule to be neglected
	 */
	public DistanceRule() {
		super(null);
		}

	/**
	 * Constructor for 2 bonds in between.
	 * @param atom
	 * @param notList atoms (direct neighbors of atom[0] and atom[1]) to be excluded from movement
	 * @param distance
	 */
	public DistanceRule(int[] atom, int[] notList, double distance, int priority) {
		super(atom);
		mDistance = new double[1];
		mDistance[0] = distance;
		mNotList = notList;
		mPriority = priority;
		}

	public DistanceRule(int[] atom, double minDistance, double maxDistance, int priority) {
		super(atom);
		mDistance = new double[2];
		mDistance[0] = minDistance;
		mDistance[1] = maxDistance;
		mPriority = priority;
		}

	public DistanceRule(int[] atom, int[] notList, double minDistance, double maxDistance, int priority) {
		super(atom);
		mDistance = new double[2];
		mDistance[0] = minDistance;
		mDistance[1] = maxDistance;
		mNotList = notList;
		mPriority = priority;
		}

	@Override
	public int getRuleType() {
		return RULE_TYPE_DISTANCE;
		}

	public boolean isFixedDistance() {
		return mDistance.length == 1;
		}

	public static void calculateRules(ArrayList<ConformationRule> ruleList, StereoMolecule mol) {
		BondLengthSet bondLengthSet = new BondLengthSet(mol);
		BondAngleSet bondAngleSet = new BondAngleSet(mol, bondLengthSet);

		DistanceRule[][] rule = new DistanceRule[mol.getAllAtoms()][];
		for (int i=1; i<mol.getAllAtoms(); i++)
			rule[i] = new DistanceRule[i];

		// distances with 1 bond between both atoms
		for (int bond=0; bond<mol.getAllBonds(); bond++) {
			int[] atom = combineAtoms(mol.getBondAtom(0, bond), mol.getBondAtom(1, bond));
			setFixedDistance(rule, atom, atom, bondLengthSet.getLength(bond), PRIORITY_ONE_BOND);
			}

		for (int atom=0; atom<mol.getAtoms(); atom++) {
			for (int i=1; i<mol.getAllConnAtoms(atom); i++) {
				int connAtom1 = mol.getConnAtom(atom, i);
				int connBond1 = mol.getConnBond(atom, i);
				double bondLength1 = bondLengthSet.getLength(connBond1);

					// distances with 2 bonds between both atoms
				for (int j=0; j<i; j++) {
					int connAtom2 = mol.getConnAtom(atom, j);
					int connBond2 = mol.getConnBond(atom, j);

					double angle = bondAngleSet.getConnAngle(atom, i, j);

					double bondLength2 = bondLengthSet.getLength(connBond2);
					double distance = Math.sqrt(bondLength1*bondLength1+bondLength2*bondLength2
												-2*bondLength1*bondLength2*Math.cos(angle));
					int[] notAtom = new int[1];
					notAtom[0] = atom;
					setFixedDistance(rule, combineAtoms(connAtom1, connAtom2), notAtom, distance, PRIORITY_TWO_BONDS);
					}
				}
			}

		// distances with 3 bonds between both atoms (special cases only)
//		int[] bondRingSize = calculateBondRingSizes(mol);
		for (int bond=0; bond<mol.getAllBonds(); bond++) {
			if (mol.isAromaticBond(bond)
			 || (mol.isRingBond(bond) && mol.getBondRingSize(bond) <= 5))
				continue;

			int[] atom = new int[2];
			for (int i=0; i<2; i++)
				atom[i] = mol.getBondAtom(i, bond);

			if (mol.getAllConnAtoms(atom[0]) > 1
			 && mol.getAllConnAtoms(atom[1]) > 1) {

					// triple bonds
				if (mol.getBondOrder(bond) == 3) {
					double distance = bondLengthSet.getLength(bond);
					int[] outerAtom = new int[2];
					for (int i = 0; i < 2; i++) {
						for (int j = 0; j < mol.getAllConnAtoms(atom[i]); j++) {
							int connBond = mol.getConnBond(atom[i], j);
							if (connBond != bond) {
								distance += bondLengthSet.getLength(connBond);
								outerAtom[i] = mol.getConnAtom(atom[i], j);
								break;
								}
							}
						}
					int[] notAtom = new int[2];
					notAtom[0] = atom[0];
					notAtom[1] = atom[1];
					setFixedDistance(rule, combineAtoms(outerAtom[0], outerAtom[1]), notAtom, distance, PRIORITY_THREE_BONDS);
					}
				else if (mol.getBondOrder(bond) == 2
					  && mol.getAtomPi(atom[0]) == 1
					  && mol.getAtomPi(atom[1]) == 1
					  && mol.getBondParity(bond) != Molecule.cBondParityUnknown) {
					// strainless double bond with stereo information
					// (including symmetrical ones with parityNone)
					int[][] connAtom = new int[2][];
					int[][] connBond = new int[2][];
					double[][] connAngle = new double[2][];
					for (int i=0; i<2; i++) {
						connAtom[i] = new int[mol.getAllConnAtoms(atom[i])-1];
						connBond[i] = new int[mol.getAllConnAtoms(atom[i])-1];
						connAngle[i] = new double[mol.getAllConnAtoms(atom[i])-1];

						int doubleBondOpponentIndex = -1;
						for (int j=0; j<mol.getAllConnAtoms(atom[i]); j++) {
							if (mol.getConnAtom(atom[i], j) == atom[1-i]) {
								doubleBondOpponentIndex = j;
								break;
								}
							}

						int connIndex = 0;
						for (int j=0; j<mol.getAllConnAtoms(atom[i]); j++) {
							if (j != doubleBondOpponentIndex) {
								connAtom[i][connIndex] = mol.getConnAtom(atom[i], j);
								connBond[i][connIndex] = mol.getConnBond(atom[i], j);
								connAngle[i][connIndex] = bondAngleSet.getConnAngle(atom[i], doubleBondOpponentIndex, j);
								connIndex++;
								}
							}
						}

					for (int i=0; i<connAtom[0].length; i++) {
						for (int j=0; j<connAtom[1].length; j++) {
							boolean isE = (mol.getBondParity(bond) == Molecule.cBondParityEor1);
							if (connAtom[0].length == 2 && connAtom[0][i] > connAtom[0][1-i])
								isE = !isE;
							if (connAtom[1].length == 2 && connAtom[1][j] > connAtom[1][1-j])
								isE = !isE;
							setDoubleBondDistance(connAtom[0][i], connAtom[1][j],
												  connBond[0][i], connBond[1][j],
												  connAngle[0][i], connAngle[1][j],
												  bond, isE, bondLengthSet, rule, mol);
							}
						}
					}
				else if (mol.getBondOrder(bond) != 0) {
					int[] opponentIndex = new int[2];
					for (int i = 0; i < 2; i++) {
						for (int j = 0; j < mol.getAllConnAtoms(atom[i]); j++) {
							if (mol.getConnAtom(atom[i], j) == atom[1 - i]) {
								opponentIndex[i] = j;
								break;
								}
							}
						}
					for (int i = 0; i < mol.getAllConnAtoms(atom[0]); i++) {
						if (i != opponentIndex[0]) {
							for (int j = 0; j < mol.getAllConnAtoms(atom[1]); j++) {
								if (j != opponentIndex[1]) {
									if (mol.getAtomPi(atom[0]) == 0
											&& mol.getAtomPi(atom[1]) == 0)
										setSingleBondConnAtomDistance(
												mol.getConnAtom(atom[0], i), mol.getConnAtom(atom[1], j),
												mol.getConnBond(atom[0], i), mol.getConnBond(atom[1], j),
												bondAngleSet.getConnAngle(atom[0], opponentIndex[0], i),
												bondAngleSet.getConnAngle(atom[1], opponentIndex[1], j),
												bond, bondLengthSet, rule, mol);
									else
										setAnyBondConnAtomDistance(
												mol.getConnAtom(atom[0], i), mol.getConnAtom(atom[1], j),
												mol.getConnBond(atom[0], i), mol.getConnBond(atom[1], j),
												bondAngleSet.getConnAngle(atom[0], opponentIndex[0], i),
												bondAngleSet.getConnAngle(atom[1], opponentIndex[1], j),
												bond, bondLengthSet, rule, mol);
									}
								}
							}
						}
					}
				}
			}

			// distances over 4 bonds in allenes
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.getAtomPi(atom) == 2
			 && mol.getConnAtoms(atom) == 2
			 && mol.getConnBondOrder(atom, 0) == 2
			 && mol.getConnBondOrder(atom, 1) == 2) {
				int atom1 = mol.getConnAtom(atom, 0);
				int atom2= mol.getConnAtom(atom, 1);
				for (int i=0; i<mol.getAllConnAtoms(atom1); i++) {
					int conn1 = mol.getConnAtom(atom1, i);
					if (conn1 != atom) {
						for (int j=0; j<mol.getAllConnAtoms(atom2); j++) {
							int conn2 = mol.getConnAtom(atom2, j);
							if (conn2 != atom) {
								double angle1 = bondAngleSet.getAngle(atom1, atom, conn1);
								double angle2 = bondAngleSet.getAngle(atom2, atom, conn2);
								double bondLength1 = bondLengthSet.getLength(mol.getConnBond(atom1, i));
								double bondLength2 = bondLengthSet.getLength(mol.getConnBond(atom2, j));
								double dx = bondLengthSet.getLength(mol.getConnBond(atom, 0))
										  + bondLengthSet.getLength(mol.getConnBond(atom, 1))
										  - bondLength1*Math.cos(angle1)
										  - bondLength2*Math.cos(angle2);
								double dy = bondLength1*Math.sin(angle1);
								double dz = bondLength2*Math.sin(angle2);
								int[] notAtom = new int[2];
								notAtom[0] = atom1;
								notAtom[1] = atom2;
								setFixedDistance(rule, combineAtoms(conn1, conn2), notAtom, Math.sqrt(dx*dx+dy*dy+dz*dz), PRIORITY_THREE_BONDS);
								}
							}
						}
					}
				}
			}

		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			calculateLongDistanceRules(rule, atom, mol, bondLengthSet);

		calculateDisconnectedDistanceRules(rule, mol);

		for (int i=1; i<mol.getAllAtoms(); i++)
			for (int j=0; j<i; j++)
				if (rule[i][j].mDistance != null)	// skip dummy rules
					ruleList.add(rule[i][j]);
		}

	private static void calculateLongDistanceRules(DistanceRule[][] rule, int rootAtom, StereoMolecule mol, BondLengthSet bondLengthSet) {
		int[] bondCount = new int[mol.getAllAtoms()];
		int[] graphAtom = new int[mol.getAllAtoms()];
		double[] distanceToRoot = new double[mol.getAllAtoms()];

		graphAtom[0] = rootAtom;
		int current = 0;
		int highest = 0;
	 	while (current <= highest) {
			int parent = graphAtom[current];
			for (int i=0; i<mol.getAllConnAtoms(parent); i++) {
				int candidate = mol.getConnAtom(parent, i);
				if (bondCount[candidate] == 0 && candidate != rootAtom) {
					graphAtom[++highest] = candidate;
					bondCount[candidate] = bondCount[parent] + 1;

					if (bondCount[candidate] == 2)
						distanceToRoot[candidate] = (candidate < rootAtom) ?
													rule[rootAtom][candidate].mDistance[0]
												  : rule[candidate][rootAtom].mDistance[0];
						// distances with 3 or more bonds in between
					else if (bondCount[candidate] > 2) {
						distanceToRoot[candidate] = distanceToRoot[parent]
												  + bondLengthSet.getLength(mol.getConnBond(parent, i));

						if (candidate < rootAtom && rule[rootAtom][candidate] == null) {
							if (bondCount[candidate] == 3) {
								rule[rootAtom][candidate] = new DistanceRule();	// add dummy rule to block this atom combination
								}
							else {
								int[] notList = new int[2];
								notList[0] = graphAtom[1];
								notList[1] = parent;
								rule[rootAtom][candidate] = new DistanceRule(combineAtoms(rootAtom, candidate), notList,
										getVDWRadius(rootAtom, mol) + getVDWRadius(candidate, mol), distanceToRoot[candidate], 0);
								}
							}
						}
					}
				}
			current++;
			}
		}

	private static void calculateDisconnectedDistanceRules(DistanceRule[][] rule, StereoMolecule mol) {
		for (int atom1=1; atom1<mol.getAllAtoms(); atom1++) {
			for (int atom2=0; atom2<atom1; atom2++) {
				if (rule[atom1][atom2] == null) {
					rule[atom1][atom2] = new DistanceRule(combineAtoms(atom1, atom2),
							getVDWRadius(atom1, mol) + getVDWRadius(atom2, mol), Double.MAX_VALUE, -1);
					}
				}
			}
		}

	private static double getVDWRadius(int atom, StereoMolecule mol) {
		int atomicNo = mol.getAtomicNo(atom);
		return (atomicNo < VDWRadii.VDW_RADIUS.length) ? VDWRadii.VDW_RADIUS[atomicNo] : 2.0f;
		}

	private static void setDoubleBondDistance(int atom1, int atom2, int bond1, int bond2,
	                                          double angle1, double angle2, int ezBond, boolean isE,
									   		  BondLengthSet bondLengthSet, DistanceRule[][] rule,
									   		  StereoMolecule mol) {
		double s1 = bondLengthSet.getLength(ezBond)
				  - bondLengthSet.getLength(bond1) * Math.cos(angle1)
				  - bondLengthSet.getLength(bond2) * Math.cos(angle2);
		double s2 = bondLengthSet.getLength(bond1) * Math.sin(angle1);
		if (isE)
			s2 += bondLengthSet.getLength(bond2) * Math.sin(angle2);
		else
			s2 -= bondLengthSet.getLength(bond2) * Math.sin(angle2);
		int[] notAtom = new int[2];
		notAtom[0] = mol.getBondAtom(0, ezBond);
		notAtom[1] = mol.getBondAtom(1, ezBond);
		setFixedDistance(rule, combineAtoms(atom1, atom2), notAtom, Math.sqrt(s1*s1+s2*s2), PRIORITY_THREE_BONDS);
		}

	/**
	 * Adds a rules for two atoms with 3 bonds in between, where the central
	 * bond is rotatable. We assume the geometry with the central bond torsion
	 * in 60 degrees as minimum distance and with 180 degrees as maximum distance.
	 * @param atom1
	 * @param atom2
	 * @param bond1
	 * @param bond2
	 * @param angle1
	 * @param angle2
	 * @param centralBond
	 * @param bondLengthSet
	 * @param rule
	 * @param mol
	 */
	private static void setSingleBondConnAtomDistance(int atom1, int atom2, int bond1, int bond2,
	                                                  double angle1, double angle2, int centralBond,
	                                                  BondLengthSet bondLengthSet, DistanceRule[][] rule,
	                                                  StereoMolecule mol) {
		if (atom1 == atom2)
			return;	// possible in 3-membered ring

		int[] atom = combineAtoms(atom1, atom2);
		if (rule[atom[0]][atom[1]] != null && rule[atom[0]][atom[1]].isFixedDistance())
			return;

		double sinTorsion = 0.866;	// default: 60 degrees, i.e. we assume atoms come not closer than in gauche-position
		double cosTorsion = 0.5;	// default: 60 degrees
		if (mol.isRingBond(centralBond)) {
			int ringSize = mol.getBondRingSize(centralBond);
			if (ringSize < 6) {
				sinTorsion = 0.0f;	// for strained rings we allow syn-position without strain
				cosTorsion = 1.0f;
			}
		}

		// distance along central bond, which is independent of central bond torsion
		double dx = bondLengthSet.getLength(centralBond)
				- bondLengthSet.getLength(bond1) * Math.cos(angle1)
				- bondLengthSet.getLength(bond2) * Math.cos(angle2);
		double s1 = bondLengthSet.getLength(bond1) * Math.sin(angle1);
		double s2 = bondLengthSet.getLength(bond2) * Math.sin(angle2);
		double dyMax = s1 + s2;
		double dyMin = s1 - s2 * cosTorsion;
		double dz = s2 * sinTorsion;
		double min = Math.sqrt(dx*dx+dyMin*dyMin+dz*dz);
		double max = Math.sqrt(dx*dx+dyMax*dyMax);
		int[] notAtom = new int[2];
		notAtom[0] = mol.getBondAtom(0, centralBond);
		notAtom[1] = mol.getBondAtom(1, centralBond);

		DistanceRule currentRule = rule[atom[0]][atom[1]];
		if (currentRule == null) {
			rule[atom[0]][atom[1]] = new DistanceRule(atom, notAtom, min, max, PRIORITY_THREE_BONDS);
		}
		else {
			currentRule.mDistance[0] = Math.min(currentRule.mDistance[0], min);
			currentRule.mDistance[1] = Math.min(currentRule.mDistance[1], max);
		}
	}

	/**
	 * Adds a rules for two atoms with 3 bonds in between, where the central
	 * bond is rotatable. We assume the geometry with the central bond torsion
	 * in 60 degrees as minimum distance and with 180 degrees as maximum distance.
	 * @param atom1
	 * @param atom2
	 * @param bond1
	 * @param bond2
	 * @param angle1
	 * @param angle2
	 * @param centralBond
	 * @param bondLengthSet
	 * @param rule
	 * @param mol
	 */
	private static void setAnyBondConnAtomDistance(int atom1, int atom2, int bond1, int bond2,
												   double angle1, double angle2, int centralBond,
												   BondLengthSet bondLengthSet, DistanceRule[][] rule,
												   StereoMolecule mol) {
		if (atom1 == atom2)
			return;	// possible in 3-membered ring

		int[] atom = combineAtoms(atom1, atom2);
		if (rule[atom[0]][atom[1]] != null && rule[atom[0]][atom[1]].isFixedDistance())
			return;

		// distance along central bond, which is independent of central bond torsion
		double dx = bondLengthSet.getLength(centralBond)
				  - bondLengthSet.getLength(bond1) * Math.cos(angle1)
				  - bondLengthSet.getLength(bond2) * Math.cos(angle2);
		double s1 = bondLengthSet.getLength(bond1) * Math.sin(angle1);
		double s2 = bondLengthSet.getLength(bond2) * Math.sin(angle2);
		double dyMax = s1 + s2;
		double dyMin = s1 - s2;
		double min = Math.sqrt(dx*dx+dyMin*dyMin);
		double max = Math.sqrt(dx*dx+dyMax*dyMax);
		int[] notAtom = new int[2];
		notAtom[0] = mol.getBondAtom(0, centralBond);
		notAtom[1] = mol.getBondAtom(1, centralBond);

		DistanceRule currentRule = rule[atom[0]][atom[1]];
		if (currentRule == null) {
			rule[atom[0]][atom[1]] = new DistanceRule(atom, notAtom, min, max, PRIORITY_THREE_BONDS);
			}
		else {
			currentRule.mDistance[0] = Math.min(currentRule.mDistance[0], min);
			currentRule.mDistance[1] = Math.min(currentRule.mDistance[1], max);
			}
		}

	/**
	 * Defines a fixed distance rule between the given atoms.
	 * If there is already a distance rule defined for these atoms,
	 * then the priority decides, which one gets precedence. If priorities are
	 * equal, then the distance value and not-lists are merged.
	 * The distance value will be a mean one and the notLists are combined.
	 * If there is a distance range already defined, then it is replaced by
	 * a new fixed distance rule.
	 * @param rule distance rule matrix between all atoms
	 * @param atom both atoms in array with first atom being the larger one
	 * @param notList
	 * @param distance
	 * @param priority
	 */
	private static void setFixedDistance(DistanceRule[][] rule, int[] atom, int[] notList, double distance, int priority) {
		if (rule[atom[0]][atom[1]] == null) {
			rule[atom[0]][atom[1]] = new DistanceRule(atom, notList, distance, priority);
			}
		else if (rule[atom[0]][atom[1]].mDistance.length == 2
			  || rule[atom[0]][atom[1]].mPriority < priority) {
			rule[atom[0]][atom[1]] = new DistanceRule(atom, notList, distance, priority);
			}
		else if (rule[atom[0]][atom[1]].mPriority == priority) {
			rule[atom[0]][atom[1]].mDistance[0] = (rule[atom[0]][atom[1]].mDistance[0] + distance) / 2f;
			rule[atom[0]][atom[1]].mNotList = mergeNotLists(rule[atom[0]][atom[1]].mNotList, notList);
			}
		}

	private static final int[] mergeNotLists(int[] nl1, int[] nl2) {
		if (nl1 == null)
			return nl2;
		if (nl2 == null)
			return nl1;
		int[] nl = new int[nl1.length+nl2.length];
		int index = 0;
		for (int atom:nl1)
			nl[index++] = atom;
		for (int atom:nl2)
			nl[index++] = atom;
		return nl;
		}

	/**
	 * puts atom1 and atom2 into an array, such that the first atom is the larger one
	 * @param atom1
	 * @param atom2
	 * @return
	 */
	private static int[] combineAtoms(int atom1, int atom2) {
		int[] atom = new int[2];
		if (atom1 > atom2) {
			atom[0] = atom1;
			atom[1] = atom2;
			}
		else {
			atom[0] = atom2;
			atom[1] = atom1;
			}
		return atom;
		}

	private static int[] calculateBondRingSizes(StereoMolecule mol) {
		int[] bondRingSize = new int[mol.getBonds()];
		RingCollection ringSet = mol.getRingSet();
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			int[] ringBond = ringSet.getRingBonds(ring);
			for (int i=0; i<ringBond.length; i++) {
				if (bondRingSize[ringBond[i]] == 0
				 || bondRingSize[ringBond[i]] > ringBond.length) {
					bondRingSize[ringBond[i]] = ringBond.length;
					}
				}
			}
		return bondRingSize;
		}

	@Override
	public boolean apply(Conformer conformer, double cycleFactor) {
		double dx = conformer.getX(mAtom[1]) - conformer.getX(mAtom[0]);
		double dy = conformer.getY(mAtom[1]) - conformer.getY(mAtom[0]);
		double dz = conformer.getZ(mAtom[1]) - conformer.getZ(mAtom[0]);
		double distance = Math.sqrt(dx*dx+dy*dy+dz*dz);

		double distanceFactor = 0.0f;
		if (mDistance.length == 2) {	// is min and max
			if (distance < mDistance[0]) {
				distanceFactor = (distance-mDistance[0]) / distance;
				}
			else if (distance > mDistance[1]) {
				distanceFactor = (distance-mDistance[1]) / distance;
				}
			}
		else {	// exact distance
			if (distance < mDistance[0]) {
				distanceFactor = (distance-mDistance[0]) / distance;
				}
			else if (distance > mDistance[0]) {
				distanceFactor = (distance-mDistance[0]) / distance;
				}
			}

		if (Math.abs(distanceFactor) < 0.001)
			return false;

		double factor = cycleFactor * distanceFactor;

		StereoMolecule mol = conformer.getMolecule();

		if (mPriority == PRIORITY_ONE_BOND) {
			if (mol.getAllConnAtoms(mAtom[0]) == 1
			 && mol.getAllConnAtoms(mAtom[1]) != 1) {
				conformer.getCoordinates(mAtom[0]).add(dx*factor, dy*factor, dz*factor);
				return true;
				}
			if (mol.getAllConnAtoms(mAtom[0]) != 1
			 && mol.getAllConnAtoms(mAtom[1]) == 1) {
				conformer.getCoordinates(mAtom[1]).add(-dx*factor, -dy*factor, -dz*factor);
				return true;
				}
			}

		factor /= 2f;
		moveGroup(conformer, mAtom[0], mNotList, dx*factor, dy*factor, dz*factor);
		moveGroup(conformer, mAtom[1], mNotList, -dx*factor, -dy*factor, -dz*factor);

		return true;
		}

	@Override
	public double addStrain(Conformer conformer, double[] atomStrain) {
		double strain = getStrain(conformer);
		if (atomStrain != null && strain > 0.0) {
			atomStrain[mAtom[0]] += strain / 2;
			atomStrain[mAtom[1]] += strain / 2;
			}
		return strain;
		}

	private double getStrain(Conformer conformer) {
		double distance = conformer.getCoordinates(mAtom[1]).distance(conformer.getCoordinates(mAtom[0]));
		if (mDistance.length == 2) {
			if (distance < VDW_RADIUS_CORRECTION * mDistance[0]) {
				return calculateVDWStrain(conformer.getMolecule(), distance);
				}
			else if (distance > mDistance[1]) {
				return calculateRelaxedStrain((distance - mDistance[1]) / distance);
				}
			}
		else {
			double strain = (mPriority == PRIORITY_ONE_BOND) ?
					calculateDirectConnectionStrain(Math.abs(mDistance[0] - distance) / Math.max(mDistance[0], distance))
						  : (mPriority == PRIORITY_TWO_BONDS) ?
					calculateTwoBondStrain(Math.abs(mDistance[0] - distance) / Math.max(mDistance[0], distance))
					: calculateRelaxedStrain(Math.abs(mDistance[0] - distance) / Math.max(mDistance[0], distance));
			if (Math.abs(strain) > 0.01f)
				return strain;
			}
		return 0.0;
		}

	public void printStrain(Conformer conformer) {
		double distance = conformer.getCoordinates(mAtom[1]).distance(conformer.getCoordinates(mAtom[0]));
		double strain = 0.0;
		double must = 0.0;
		if (mDistance.length == 2) {
			if (distance < VDW_RADIUS_CORRECTION * mDistance[0]) {
				must = mDistance[0];
				strain = calculateVDWStrain(conformer.getMolecule(), distance);
				}
			else if (distance > mDistance[1]) {
				must = mDistance[1];
				strain = calculateRelaxedStrain((distance - mDistance[1]) / distance);
				}
			}
		else {
			must = mDistance[0];
			strain = (mPriority == PRIORITY_ONE_BOND) ?
					calculateDirectConnectionStrain(Math.abs(mDistance[0] - distance) / Math.max(mDistance[0], distance))
					: (mPriority == PRIORITY_TWO_BONDS) ?
					calculateTwoBondStrain(Math.abs(mDistance[0] - distance) / Math.max(mDistance[0], distance))
					: calculateRelaxedStrain(Math.abs(mDistance[0] - distance) / Math.max(mDistance[0], distance));
			}
		if (strain > 0.001) {
			System.out.print("Atoms("+Molecule.cAtomLabel[conformer.getMolecule().getAtomicNo(mAtom[0])]+mAtom[0]
					+(mPriority == PRIORITY_ONE_BOND? "-" : mPriority == PRIORITY_TWO_BONDS? "-?-" : mPriority == PRIORITY_THREE_BONDS? "-?-?-" : "-...-")
					+Molecule.cAtomLabel[conformer.getMolecule().getAtomicNo(mAtom[1])]+mAtom[1]+") distance:"+DoubleFormat.toString(distance, 4));
			System.out.println(" must:"+DoubleFormat.toString(must,4)+" strain:"+DoubleFormat.toString(strain, 4));
			}
		}

	private double calculateDirectConnectionStrain(double deltaDistance) {
		// We use a quadratic potential with an energy penalty of 100.0 kcal/mol when 10% off
		return 10000.0 * deltaDistance * deltaDistance;
		}

	private double calculateTwoBondStrain(double deltaDistance) {
		// We use a quadratic potential with an energy penalty of 50 kcal/mol when 10% off
		return 8000.0 * deltaDistance * deltaDistance;
		}

	private double calculateRelaxedStrain(double deltaDistance) {
		// We use a quadratic potential with an energy penalty of 50 kcal/mol when 10% off
		return 4000.0 * deltaDistance * deltaDistance;
		}

	private double calculateVDWStrain(StereoMolecule mol, double distance) {
		// We use the repulsion part of the Lennard-Jones potential
		double vdwradii = VDW_RADIUS_CORRECTION * getVDWRadius(mAtom[1], mol) + getVDWRadius(mAtom[0], mol);
		double reldist = distance / vdwradii;
		double reldist6 = Math.pow(reldist, -6);
		double constant = 2.0f; // constant=1.0 causes a return value of 1.66 kcal/mol with reldist=0.9
		return (reldist >= 1.0) ? 0.0 : constant * (reldist6 * reldist6 - reldist6);
		}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("distance rule:");
		super.addAtomList(sb);
		if (mDistance.length == 1)
			sb.append(" distance:"+DoubleFormat.toString(mDistance[0]));
		else
			sb.append(" min:"+DoubleFormat.toString(mDistance[0])+" max:"+DoubleFormat.toString(mDistance[1]));
		if (mNotList != null) {
			sb.append(" not:"+mNotList[0]);
			for (int i=1; i<mNotList.length; i++)
				sb.append(","+mNotList[i]);
			}
		return sb.toString();
		}
	}
