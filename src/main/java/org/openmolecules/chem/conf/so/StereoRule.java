/*
 * @(#)StereoRule.java
 *
 * Copyright 2013 openmolecules.org, Inc. All Rights Reserved.
 *
 * NOTICE: All information contained herein is, and remains the property
 * of openmolecules.org.  The intellectual and technical concepts contained
 * herein are proprietary to openmolecules.org.
 * Actelion Pharmaceuticals Ltd. is granted a non-exclusive, non-transferable
 * and timely unlimited usage license.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.so;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

import java.util.ArrayList;

public class StereoRule extends ConformationRule {
	public StereoRule(int[] atom) {
		super(atom);
		}

    public static void calculateRules(ArrayList<ConformationRule> ruleList, StereoMolecule mol) {
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			int parity = mol.getAtomParity(atom);
			if (parity == Molecule.cAtomParity1 || parity == Molecule.cAtomParity2) {
				int[] atomList = new int[5];
				for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
					int connAtom = mol.getConnAtom(atom, i);
					int index = 0;
					while (index < i && connAtom > atomList[index])
						index++;
					for (int j=i-1; j>=index; j--)
						atomList[j+1] = atomList[j];
					atomList[index] = connAtom;
					}

				if (mol.getAllConnAtoms(atom) == 3)
					atomList[3] = -1;

				atomList[4] = atom;

				if (parity == Molecule.cAtomParity1) {
					int temp = atomList[2];
					atomList[2] = atomList[1];
					atomList[1] = temp;
					}

				ruleList.add(new StereoRule(atomList));
				}
			}
		}

	@Override
	public boolean apply(Conformer conformer, double cycleFactor) {
		double[] n = getNormalVector(conformer, mAtom);
		boolean invertAtom3 = (mAtom[3] != -1
							&& getStereoAngleCosine(conformer, n, mAtom[0], mAtom[3]) > 0.0);
		boolean invertAtom4 = (getStereoAngleCosine(conformer, n, mAtom[0], mAtom[4]) > 0.0);

		if (!invertAtom3 && !invertAtom4)
			return false;

		// invert stereocenter by moving atom[3] and/or atom[4] through plane of other atoms
		double d = n[0]*conformer.getX(mAtom[0])
				 + n[1]*conformer.getY(mAtom[0])
				 + n[2]*conformer.getZ(mAtom[0]);	// distance of Hesse equation of plane

		double distance3 = (mAtom[3] == -1) ? 0.0
						 : n[0]*conformer.getX(mAtom[3])
						 + n[1]*conformer.getY(mAtom[3])
						 + n[2]*conformer.getZ(mAtom[3]) - d;

		double distance4 = n[0]*conformer.getX(mAtom[4])
						 + n[1]*conformer.getY(mAtom[4])
						 + n[2]*conformer.getZ(mAtom[4]) - d;

		if (mAtom[3] == -1 || conformer.getMolecule().getConnAtoms(mAtom[3]) == 1 || !invertAtom4) {
			// keep atoms 0,1,2 in plane and move atom3 and/or atom4 through plane
			if (invertAtom3)
				conformer.getCoordinates(mAtom[3]).add(-2.0*distance3*n[0],
													  -2.0*distance3*n[1],
													  -2.0*distance3*n[2]);
			if (invertAtom4)
				conformer.getCoordinates(mAtom[4]).add(-2.0*distance4*n[0],
													  -2.0*distance4*n[1],
													  -2.0*distance4*n[2]);
			}
		else {	// for neighbors and central atom needs to be inverted: keep center of gravity
			for (int i=0; i<3; i++)
				conformer.getCoordinates(mAtom[i]).add(0.5*distance4*n[0],
													  0.5*distance4*n[1],
													  0.5*distance4*n[2]);
			if (invertAtom3)
				conformer.getCoordinates(mAtom[3]).add(0.5*distance4*n[0] - 2.0*distance3*n[0],
													  0.5*distance4*n[1] - 2.0*distance3*n[1],
													  0.5*distance4*n[2] - 2.0*distance3*n[2]);
			conformer.getCoordinates(mAtom[4]).add(-1.5*distance4*n[0],
												  -1.5*distance4*n[1],
												  -1.5*distance4*n[2]);
			}

		return true;
		}

	@Override
	public double addStrain(Conformer conformer, double[] atomStrain) {
		double totalStrain = 0;
		double[] n = getNormalVector(conformer, mAtom);
		if (getStereoAngleCosine(conformer, n, mAtom[0], mAtom[4]) > 0.0
		 || (mAtom[3] != -1 && getStereoAngleCosine(conformer, n, mAtom[0], mAtom[3]) > 0.0)) {
			for (int i=0; i<mAtom.length; i++) {
				if (mAtom[i] != -1) {
					double panalty = 0.25; // arbitrary value 0.5 * 0.5;
					atomStrain[mAtom[i]] += panalty;
					totalStrain += panalty;
					}
				}
			}
		return totalStrain;
		}

	/**
	 * Returns the angle between vector n and vector (atom0->atom)
	 * @param conformer
	 * @param n vector n (typically the normal vector of the plane a0,a1,a2)
	 * @param atom0	root of second vector
	 * @param atom end of second vector
	 * @return
	 */
	private double getStereoAngleCosine(Conformer conformer, double[] n, int atom0, int atom) {
		// calculate the three vectors leading from atom[0] to the other three atoms
		double[] v = new double[3];
		v[0] = conformer.getX(atom) - conformer.getX(atom0);
		v[1] = conformer.getY(atom) - conformer.getY(atom0);
		v[2] = conformer.getZ(atom) - conformer.getZ(atom0);

		// calculate cos(angle) of coords[2] to normal vector
		return (v[0]*n[0]+v[1]*n[1]+v[2]*n[2]) / Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
		}

	/**
	 * Calculates the normal vector of plane (a0->a1, a0->a2)
	 * @param conformer
	 * @param atom
	 * @return normal vector
	 */
	private double[] getNormalVector(Conformer conformer, int[] atom) {
		double[] n = new double[3];
		getNormalVector(conformer, atom, n);
		return n;
		}

	/**
	 * Calculates the normal vector of plane (a0->a1, a0->a2)
	 * @param conformer
	 * @param atom
	 * @param n receives the normal vector of the plane a0,a1,a2
	 * @return
	 */
	private void getNormalVector(Conformer conformer, int[] atom, double[] n) {
		// calculate the three vectors leading from atom[0] to the other three atoms
		double[][] coords = new double[2][3];
		for (int i=0; i<2; i++) {
			coords[i][0] = conformer.getX(atom[i+1]) - conformer.getX(atom[0]);
			coords[i][1] = conformer.getY(atom[i+1]) - conformer.getY(atom[0]);
			coords[i][2] = conformer.getZ(atom[i+1]) - conformer.getZ(atom[0]);
			}

		// calculate the normal vector (vector product of coords[0] and coords[1])
		n[0] = coords[0][1]*coords[1][2]-coords[0][2]*coords[1][1];
		n[1] = coords[0][2]*coords[1][0]-coords[0][0]*coords[1][2];
		n[2] = coords[0][0]*coords[1][1]-coords[0][1]*coords[1][0];

		double l = Math.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
		n[0] /= l;
		n[1] /= l;
		n[2] /= l;
		}

	@Override
	public int getRuleType() {
		return RULE_TYPE_STEREO;
		}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("stereo rule:");
		super.addAtomList(sb);
		return sb.toString();
		}
	}
