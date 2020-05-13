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

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

import java.util.ArrayList;

public class TetrahedralStereoRule extends ConformationRule {
	private int[] mRotatableAtom;
	private Coordinates mAxisOfRotation;

	public TetrahedralStereoRule(StereoMolecule mol, int[] atom) {
		super(atom);
		mRotatableAtom = getRotatableAtoms(mol, atom[4]);
		}

    public static void calculateRules(ArrayList<ConformationRule> ruleList, StereoMolecule mol) {
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.getAllConnAtoms(atom) >= 3) {
				int parity = mol.getAtomParity(atom);
				if ((parity == Molecule.cAtomParity1 || parity == Molecule.cAtomParity2)
				 && !mol.isCentralAlleneAtom(atom)) {
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

					ruleList.add(new TetrahedralStereoRule(mol, atomList));
					}
				}
			}
	    }

	@Override
	public boolean apply(Conformer conformer, double cycleFactor) {
		double[] n = getPlaneVector(conformer);
		boolean invert = (mAtom[3] == -1 && isOnSameSide(conformer, n, mAtom[0], mAtom[4]))
					  || (mAtom[3] != -1 && isOnSameSide(conformer, n, mAtom[0], mAtom[3]));

		if (!invert)
			return false;

		for (int atom:mRotatableAtom)
			rotateAtom(conformer, atom, mAtom[1], mAxisOfRotation, Math.PI);

		return true;




		/* invert stereocenter by moving atom[3] and/or atom[4] through plane of other atoms
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

		return true;*/
		}

	@Override
	public double addStrain(Conformer conformer, double[] atomStrain) {
		double totalStrain = 0;
		double[] n = getPlaneVector(conformer);
		if ((mAtom[3] == -1 && isOnSameSide(conformer, n, mAtom[0], mAtom[4]))
		 || (mAtom[3] != -1 && isOnSameSide(conformer, n, mAtom[0], mAtom[3]))) {
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
	 * Returns whether atom1 is on the same side of the plane defined by atom0 and n
	 * @param conformer
	 * @param n vector of the plane a0,a1,a2
	 * @param atom0	atom on plane
	 * @param atom1 atom not on plane
	 * @return
	 */
	private boolean isOnSameSide(Conformer conformer, double[] n, int atom0, int atom1) {
		// calculate the three vectors leading from atom[0] to the other three atoms
		double[] v = new double[3];
		v[0] = conformer.getX(atom1) - conformer.getX(atom0);
		v[1] = conformer.getY(atom1) - conformer.getY(atom0);
		v[2] = conformer.getZ(atom1) - conformer.getZ(atom0);

		// calculate cos(angle) of a0->a1 to normal vector
		return v[0]*n[0]+v[1]*n[1]+v[2]*n[2] > 0;
		}

	/**
	 * Calculates the vector of plane (a0->a1, a0->a2)
	 * @param conformer
	 * @return
	 */
	private double[] getPlaneVector(Conformer conformer) {
		double[][] coords = new double[2][3];
		for (int i=0; i<2; i++) {
			coords[i][0] = conformer.getX(mAtom[i+1]) - conformer.getX(mAtom[0]);
			coords[i][1] = conformer.getY(mAtom[i+1]) - conformer.getY(mAtom[0]);
			coords[i][2] = conformer.getZ(mAtom[i+1]) - conformer.getZ(mAtom[0]);
			}

		// calculate the vector of the plane (vector product of coords[0] and coords[1])
		double[] v = new double[3];
		v[0] = coords[0][1]*coords[1][2]-coords[0][2]*coords[1][1];
		v[1] = coords[0][2]*coords[1][0]-coords[0][0]*coords[1][2];
		v[2] = coords[0][0]*coords[1][1]-coords[0][1]*coords[1][0];
		return v;
		}

	private int[] getRotatableAtoms(StereoMolecule mol, int atom) {
		int neighbours = mol.getAllConnAtoms(atom);
		boolean[][] isMemberAtom = new boolean[neighbours][mol.getAllAtoms()];
		int[] substituentSize = new int[neighbours];

		// find sizes of the 3 or 4 substituents (if not in a ring)
		for (int i=0; i<neighbours; i++)
			if (!mol.isRingBond(mol.getConnBond(atom, i)))
				substituentSize[i] = mol.getSubstituent(atom, mol.getConnAtom(atom, i), isMemberAtom[i], null, null);

		int rotatableAtomCount = 0;
		boolean[] isRotatableAtom = null;
		for (int s=2; s<neighbours; s++) {
			int smallestIndex = -1;
			int smallestSize = Integer.MAX_VALUE;
			for (int i=0; i<neighbours; i++) {
				if (smallestSize > substituentSize[i] && substituentSize[i] != 0) {
					smallestSize = substituentSize[i];
					smallestIndex = i;
					}
				}
			if (smallestIndex != -1) {
				rotatableAtomCount += substituentSize[smallestIndex];
				substituentSize[smallestIndex] = 0; // mark to not find it again
				if (isRotatableAtom == null)
					isRotatableAtom = isMemberAtom[smallestIndex];
				else
					for (int a=0; a<mol.getAllAtoms(); a++)
						isRotatableAtom[a] |= isMemberAtom[smallestIndex][a];
				}
			else {
				if (isRotatableAtom == null)
					isRotatableAtom = new boolean[mol.getAllAtoms()];
				for (int i=mol.getAllConnAtoms(atom)-1; i>=0; i--) {
					int connAtom = mol.getConnAtom(atom, i);
					if (!isRotatableAtom[connAtom]) {
						isRotatableAtom[connAtom] = true;
						rotatableAtomCount++;
						break;
						}
					}
				}
			}

		int index = 0;
		int[] rotatableAtom = new int[rotatableAtomCount];

		for (int i=0; i<isRotatableAtom.length; i++)
			if (isRotatableAtom[i])
				rotatableAtom[index++] = i;

		mAxisOfRotation = new Coordinates();
		Coordinates b = new Coordinates();
		for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
			int connAtom = mol.getConnAtom(atom, i);
			if (!isRotatableAtom[connAtom])
				mAxisOfRotation.add(mol.getCoordinates(connAtom));
			else if (neighbours == 4)
				b.add(mol.getCoordinates(connAtom));
		}
		mAxisOfRotation.scale(0.5);
		if (neighbours == 4)
			b.scale(0.5);
		else
			b.add(mol.getCoordinates(atom));
		mAxisOfRotation.sub(b).unit();

		return rotatableAtom;
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
