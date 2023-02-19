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

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.Angle;
import com.actelion.research.util.DoubleFormat;

public class TorsionDescriptor implements Comparable<TorsionDescriptor> {
	private static final float TORSION_EQUIVALENCE_TOLERANCE = (float)Math.PI / 12f;

	private float[] mTorsion;
	private float[] mMaxTorsion;

	/**
	 *
	 * @param torsion torsion angles of rotatable bonds from 0 to maxTorsion
	 * @param maxTorsion maximum torsions of rotatable bonds defining individual symmetry dependent range
	 */
	public TorsionDescriptor(float[] torsion, float[] maxTorsion) {
		mTorsion = torsion;
		mMaxTorsion = maxTorsion;
		}

	/**
	 * Creates a TorsionDescriptor from conformer's coordinates
	 * drawing atom connectivity from the conformer's molecule.
	 * The torsion descriptor is not canonical, unless the passed molecule is canonical.
	 * Rotatable bonds need to carry at least one external non-hydrogen neighbor on each side.
	 * @param conformer
	 *
	public TorsionDescriptor(Conformer conformer, int[] rotatableBond) {
		mTorsion = new float[rotatableBond.length];
		StereoMolecule mol = conformer.getMolecule();
		mol.ensureHelperArrays(Molecule.cHelperSymmetrySimple);
		int[] atom = new int[4];
		for (int i=0; i<rotatableBond.length; i++) {
			for (int j=0; j<2; j++) {
				atom[j+1] = mol.getBondAtom(j, rotatableBond[i]);
				for (int k=0; k<mol.getConnAtoms(atom[j+1]); k++) {
					int refAtom = mol.getConnAtom(atom[j+1], k);
					if (refAtom != mol.getBondAtom(1-j, rotatableBond[i])) {
						atom[3*j] = refAtom;
						break;
						}
					}
				}
			mTorsion[i] = (float)getNormalizedTorsion(conformer.calculateTorsion(atom), getSymmetryClass());
			}
		}*/

	/**
	 * Returns true, if none of the torsion angles are more different
	 * than TORSION_EQUIVALENCE_TOLERANCE;
	 * @param td
	 * @return
	 */
	public boolean equals(TorsionDescriptor td) {
		for (int i=0; i<mTorsion.length; i++) {
			float dif = Math.abs(mTorsion[i] - td.mTorsion[i]);
			if (dif > TORSION_EQUIVALENCE_TOLERANCE
			 && dif < mMaxTorsion[i] - TORSION_EQUIVALENCE_TOLERANCE)
				return false;
			}
		return true;
		}

	/**
	 * Returns 0, if none of the torsion angles are more different than TORSION_EQUIVALENCE_TOLERANCE;
	 * Returns -1 if the first non-equivalent torsion angle is smaller for this than for the passed TorsionDescriptor td.
	 * @param td
	 * @return
	 */
	@Override
	public int compareTo(TorsionDescriptor td) {
		for (int i=0; i<mTorsion.length; i++) {
			float dif = Math.abs(mTorsion[i] - td.mTorsion[i]);
			if (dif > TORSION_EQUIVALENCE_TOLERANCE
			 && dif < mMaxTorsion[i] - TORSION_EQUIVALENCE_TOLERANCE)
				return (dif < mMaxTorsion[i]/2) ^ (mTorsion[i] < td.mTorsion[i]) ? 1 : -1;
			}
		return 0;
		}

	/**
	 * The relevance of a rotatable bond and its torsion angle for creating substantially different conformers
	 * depends on how close the bond is to the center of the molecule. Bond relevance values range from
	 * 1.0/atomCount (e.g. bond to methyl group) to 1.0 (bond dividing molecule into two equally large parts).
	 * Ring bonds are assigned a relevance value of 0.33 independent of their location.
	 * @param mol
	 * @param rotatableBond
	 * @return
	 */
	public static float[] getRotatableBondWeights(StereoMolecule mol, int[] rotatableBond) {
		return TorsionRelevanceHelper.getRelevance(mol, rotatableBond);
		}

	/**
	 * Calculates a similarity value between td and this considering
	 * individual torsion values, the importance of the rotatable bond,
	 * and the ratio of rotatable/non-rotatable bonds.
	 * @param td
	 * @return
	 */
	public float getDissimilarity(TorsionDescriptor td, float[] bondWeight) {
		assert(mTorsion.length == td.mTorsion.length);

		if (mTorsion.length == 0)
			return 0.0f;

		float meanAngleDiff = 0f;

		float weightSum = 0f;
		for (int i=0; i<mTorsion.length; i++) {
			meanAngleDiff += bondWeight[i] * Math.abs(Angle.difference(mTorsion[i], td.mTorsion[i]));
			weightSum += bondWeight[i];
			}
		return meanAngleDiff / ((float)Math.PI * weightSum);
		}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<mTorsion.length; i++) {
			sb.append(i == 0? "Torsions: " : ", ");
			sb.append(DoubleFormat.toString(mTorsion[i], 3)+"("+Math.round(mMaxTorsion[i]+0.5f)+")");
			}
		return sb.toString();
		}
	}
