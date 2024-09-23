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
import com.actelion.research.chem.StereoMolecule;

public class MolecularFlexibilityCalculator {
	public MolecularFlexibilityCalculator() {
		TorsionDB.initialize(TorsionDB.MODE_ANGLES);
		}

	/**
	 * Calculates a molecular flexibility as a value from 0.0 to 1.0 considering torsion statistics
	 * derived from the CSD database (torsion maxima, frequencies, and 50% intervals) of rotatable bonds.
	 * It also considers the effect that a specific torsion has on the overall geometry of the
	 * molecule, i.e. central bonds are weighted higher than those in the periphery.
	 * Currently it is assumed that the mol <b>consists of one(!!!) fragment only</b>.
	 * <br><br>
	 * In more detail the following steps are done:<br>
	 * -- relevant rotatable bonds are determined, i.e. single bonds, which are not in a ring with less than 6 members,
	 *   where both atoms are sp2 or sp3 and carry at least one more non-hydrogen neighbor,
	 *   which are not symmetrically redundant, and those where a torsion change modifies the relative location of atoms.
	 *   For chains of conjugated triple bonds the following applies:
	 *   If at least one terminal sp2/sp3 atom has no external neighbor, then no single bond is considered rotatable.
	 *   Otherwise that terminal single bond connecting the smaller substituent is considered the only rotatable bond
	 *   of the linear atom strand.<br>
	 * -- the local environment of any rotatable bonds is characterized by its first and second shell of neighbours plus
	 *   various properties as ring membership, aromaticity, stereo configuration, if applicable, etc.<br>
	 * -- a bond specific flexibility value is calculated from the torsion angle histogram of equivalent bonds
	 *   within any high-resolution x-ray structures of the CSD database. (in the rare cases with no CSD precedents
	 *   the histogram is roughly predicted). Frequency distributions with wide and multiple distribution maxima of
	 *   similar heights receive the local flexibility values around close to 1.0 while histograms with one narrow single
	 *   peak are close to 0.0.<br>
	 * -- a weighting factor is assigned to every rotatable bond as follows:<br>
	 *   - for ring bonds: factor=0.33, since ring bonds cannot be changed without affecting typically two other rings bonds<br>
	 *   - other bonds: factor=sqrt(2.0 * smallerSideNonHydrogenAtomCount / moleculeNonHydrogenAtomCount)<br>
	 * -- from the number of all bonds, rotatabale bonds, their specific flexibility values and weighting factor an
	 *   overall flexibility value is calculated with a non-linear incremental approach.<br>
	 * @param mol one molecular fragment(!)
	 * @return
	 */
	public float calculateMolecularFlexibility(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		if (mol.getAllAtoms() == 0)
			return 0f;

		boolean[] isRotatableBond = new boolean[mol.getBonds()];
		int rotatableBondCount = TorsionDB.findRotatableBonds(mol, false, isRotatableBond);
		if (rotatableBondCount == 0)
			return 0f;

		float[] bondFlexibility = calculateBondFlexibilities(mol, isRotatableBond);
		float[] bondImportance = TorsionRelevanceHelper.getRelevance(mol, (boolean[])null);

/*
System.out.println();
System.out.println("Molecular flexibility:"+calculateMoleculeFlexibility(bondFlexibility, bondImportance));
for (int bond=0; bond<mol.getBonds(); bond++) {
 System.out.print("bond:"+bond+" flexibility"+DoubleFormat.toString(bondFlexibility[bond])+" importance"+DoubleFormat.toString(bondImportance[bond]));
 if (isRotatableBond[bond]) {
  int[] torsionAtom = new int[4];
  String torsionID = TorsionDB.getTorsionID(mol, bond, torsionAtom, null);
  short[] torsion = TorsionDB.getTorsions(torsionID);
  short[] frequency = null;
  short[][] range = null;
  if (torsion != null) {
   frequency = TorsionDB.getTorsionFrequencies(torsionID);
   range = TorsionDB.getTorsionRanges(torsionID);
  }
  else {
   TorsionPrediction prediction = new TorsionPrediction(mol, torsionAtom);
   torsion = prediction.getTorsions();
   frequency = prediction.getTorsionFrequencies();
   range = prediction.getTorsionRanges();
  }
  System.out.print(" torsions:");
  for (int i=0; i<torsion.length; i++)
   System.out.print(""+torsion[i]+"(f:"+frequency[i]+", r:"+range[i][0]+"-"+range[i][1]+"); ");
 }
 System.out.println();
}
*/
		return calculateMoleculeFlexibility(bondFlexibility, bondImportance);
		}

	public float[] calculateBondFlexibilities(StereoMolecule mol, boolean[] isRotatableBond) {
		float[] bondFlexibility = new float[mol.getBonds()];
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (isRotatableBond[bond]) {
				int[] torsionAtom = new int[4];
				String torsionID = TorsionDB.getTorsionID(mol, bond, torsionAtom, null);
				short[] torsion = TorsionDB.getTorsions(torsionID);
				short[] frequency = null;
				short[][] range = null;
				if (torsion != null) {
					frequency = TorsionDB.getTorsionFrequencies(torsionID);
					range = TorsionDB.getTorsionRanges(torsionID);
				}
				else {
					TorsionPrediction prediction = new TorsionPrediction(mol, torsionAtom);
					torsion = prediction.getTorsions();
					frequency = prediction.getTorsionFrequencies();
					range = prediction.getTorsionRanges();
				}

				if (!mol.isAromaticBond(bond)) {
					bondFlexibility[bond] = (torsion == null) ? 1f : calculateBondFlexibility(torsion, frequency, range);

					// we need to reduce in rings increasingly with ring size, because flexibility is hampered
					if (mol.isRingBond(bond))
						bondFlexibility[bond] *= 0.5f * (1f - 3f / mol.getBondRingSize(bond));
				}
			}
		}
		return bondFlexibility;
	}

	private float calculateBondFlexibility(short[] torsion, short[] frequency, short[][] range) {
//		for (int i=0; i<torsion.length; i++)
//			System.out.println("torsion:"+torsion[i]+" frequency:"+frequency[i]+" range:"+range[i][0]+"-"+range[i][1]);

		int count = torsion.length;

		float maxFrequency = 0f;
		int maxFrequencyIndex = -1;
		float maxEmptySpace = -1f;
		int maxEmptySpaceIndex = -1;
		for (int i=0; i<count; i++) {
			if (maxFrequency < frequency[i]) {
				maxFrequency = frequency[i];
				maxFrequencyIndex = i;
				}
			float r = rightEmptySpace(i, count, range);
			if (maxEmptySpace < r) {
				maxEmptySpace = r;
				maxEmptySpaceIndex = i;
				}
			}

		float emptyRangeSum = maxEmptySpace;

		float factor = 1f;
		int ti = maxEmptySpaceIndex;
		while (right(ti, count) != maxFrequencyIndex) {
			ti = right(ti, count);
			factor *= (float)(1f - Math.pow(frequency[ti] / maxFrequency, 0.3));
			emptyRangeSum += factor * (rightEmptySpace(ti, count, range) + range[ti][1] - range[ti][0]);
			}
		factor = 1f;
		ti = maxEmptySpaceIndex;
		while (ti != maxFrequencyIndex) {
			int oi = ti;
			ti = left(ti, count);
			factor *= (1f - (float)Math.sqrt(frequency[oi] / maxFrequency));
			emptyRangeSum += factor * (rightEmptySpace(ti, count, range) + range[oi][1] - range[oi][0]);
			}

//System.out.println("emptyRangeSum:"+emptyRangeSum+" flexibility:"+(0.5f+0.5f*(float)Math.cos(emptyRangeSum*Math.PI/360)));
		return 0.5f+0.5f*(float)Math.cos(emptyRangeSum*Math.PI/360);
		}

	private int left(int i, int count) {
		return (i == 0) ? count-1 : i-1;
		}

	private int right(int i, int count) {
		return (i == count-1) ? 0 : i+1;
		}

	private int rightEmptySpace(int i, int count, short[][] range) {
		return (i==count-1) ? 360+range[0][0]-range[i][1] : range[i+1][0]-range[i][1];

		}

	private float calculateMoleculeFlexibility(float[] bondFlexibility, float[] bondWeight) {
		final float CORRECTION_FACTOR = 0.7f;

		float meanFlexibility = 0f;
		float weightSum = 0f;
		for (int i=0; i<bondFlexibility.length; i++) {
			meanFlexibility += bondWeight[i] * bondFlexibility[i];
			weightSum += bondWeight[i];
			}
		float rawFlexibility = (weightSum == 0f) ? 0f : meanFlexibility / weightSum;
		return (float)(1.0-Math.pow(1-Math.pow(rawFlexibility, CORRECTION_FACTOR) ,1.0/CORRECTION_FACTOR));
		}
	}
