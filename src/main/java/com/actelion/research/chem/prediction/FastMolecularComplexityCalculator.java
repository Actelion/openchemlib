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
 * @author Thomas Sander, Modest v. Korff
 */

package com.actelion.research.chem.prediction;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.TreeSet;

/**
 * This is a fast version of the MolecularComplexityCalculator in the com.actelion.chem.properties.complexity
 * package. It determines the number of distinct fragments within the given molecule up to a bod count of 7.
 */
public class FastMolecularComplexityCalculator {
	protected static final int MAX_BOND_COUNT = 7;

	/**
	 * Ambiguous bonds are normalized.
	 * @param mol
	 * @return
	 */
	public static float assessComplexity(StereoMolecule mol) {
		int bondCount = Math.min(mol.getBonds()/2, MAX_BOND_COUNT);

		if (bondCount < 2)
			return 0;

		mol.ensureHelperArrays(Molecule.cHelperSymmetrySimple);
		StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());
		TreeSet<String> fragmentSet = new TreeSet<>();
		int[] atomMap = new int[mol.getAllAtoms()];

		boolean[][] bondsTouch = new boolean[mol.getBonds()][mol.getBonds()];
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			for (int i=1; i<mol.getConnAtoms(atom); i++) {
				for (int j=0; j<i; j++) {
					int bond1 = mol.getConnBond(atom, i);
					int bond2 = mol.getConnBond(atom, j);
					bondsTouch[bond1][bond2] = true;
					bondsTouch[bond2][bond1] = true;
				}
			}
		}

		boolean[] bondIsMember = new boolean[mol.getBonds()];
		int maxLevel = bondCount - 2;
		int[] levelBond = new int[maxLevel+1];
		for (int rootBond=0; rootBond<mol.getBonds(); rootBond++) {
			bondIsMember[rootBond] = true;
			int level = 0;
			levelBond[0] = rootBond;
			while (true) {
				boolean levelBondFound = false;
				while (!levelBondFound && levelBond[level] < mol.getBonds()-1) {
					levelBond[level]++;
					if (!bondIsMember[levelBond[level]]) {
						for (int bond=rootBond; bond<mol.getBonds(); bond++) {
							if (bondIsMember[bond] && bondsTouch[bond][levelBond[level]]) {
								levelBondFound = true;
								break;
							}
						}
					}
				}

				if (levelBondFound) {
					bondIsMember[levelBond[level]] = true;
					if (level == maxLevel) {
						mol.copyMoleculeByBonds(fragment, bondIsMember, true, atomMap);
						fragmentSet.add(new Canonizer(fragment).getIDCode());
						bondIsMember[levelBond[level]] = false;
					}
					else {
						level++;
						levelBond[level] = rootBond;
					}
				}
				else {
					if (--level < 0)
						break;
					bondIsMember[levelBond[level]] = false;
				}
			}
			bondIsMember[rootBond] = false;
		}

		return (float)Math.log(fragmentSet.size()) / bondCount;
	}
}
