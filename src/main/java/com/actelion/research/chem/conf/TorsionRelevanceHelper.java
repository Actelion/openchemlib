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

public class TorsionRelevanceHelper {
	/**
	 * The relevance of a rotatable bond and its torsion angle for creating substantially different conformers
	 * depends on how close the bond is to the center of the molecule. Bond relevance values range from
	 * 1.0/atomCount (e.g. bond to methyl group) to 1.0 (bond dividing molecule into two equally large parts).
	 * Ring bonds are assigned a relevance value of 0.33 independent of their location.
	 * @param mol
	 * @param rotatableBond array containing bond indexes for which to calculate relevance values
	 * @return
	 */
	public static final float[] getRelevance(StereoMolecule mol, int[] rotatableBond) {
		mol.ensureHelperArrays(Molecule.cHelperRings);

		float[] bondWeight = new float[rotatableBond.length];
		for (int i=0; i<bondWeight.length; i++) {
			if (mol.isRingBond(rotatableBond[i])) {
				bondWeight[i] = 0.33f;
				}
			else {
				int atom1 = mol.getBondAtom(0, rotatableBond[i]);
				int atom2 = mol.getBondAtom(1, rotatableBond[i]);
				if (mol.getConnAtoms(atom1) == 1 || mol.getConnAtoms(atom2) == 1) {
					bondWeight[i] = 1f / mol.getAtoms();	// rotates hydrogens only
					}
				else {
					int atomCount = mol.getSubstituentSize(atom1, atom2);
					bondWeight[i] = 2f * Math.min(atomCount, mol.getAtoms() - atomCount) / mol.getAtoms();
					}
				}
			}
		return bondWeight;
	}

	/**
	 * The relevance of a rotatable bond and its torsion angle for creating substantially different conformers
	 * depends on how close the bond is to the center of the molecule. Bond relevance values range from
	 * 1.0/atomCount (e.g. bond to methyl group) to 1.0 (bond dividing molecule into two equally large parts).
	 * @param mol
	 * @param isRotatableBond if null, then the relevance is calculated for every non-H-bond
	 * @return array with bond relevance values for all rotatable or all non-H-bonds
	 */
	public static final float[] getRelevance(StereoMolecule mol, boolean[] isRotatableBond) {
		mol.ensureHelperArrays(Molecule.cHelperRings);

		float[] bondWeight = new float[mol.getBonds()];
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (isRotatableBond == null || isRotatableBond[bond]) {
				if (mol.isRingBond(bond)) {
					bondWeight[bond] = 0.33f;
					}
				else {
					int atom1 = mol.getBondAtom(0, bond);
					int atom2 = mol.getBondAtom(1, bond);
					if (mol.getConnAtoms(atom1) == 1 || mol.getConnAtoms(atom2) == 1) {
						bondWeight[bond] = 1f / mol.getAtoms();	// rotates hydrogens only
						}
					else {
						int atomCount = mol.getSubstituentSize(atom1, atom2);
						bondWeight[bond] = 2f * Math.min(atomCount, mol.getAtoms() - atomCount) / mol.getAtoms();
						}
					}
				}
			}
		return bondWeight;
		}
	}
