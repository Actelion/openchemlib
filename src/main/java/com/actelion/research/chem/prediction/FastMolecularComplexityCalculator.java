/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 *
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander,Modest von Korff
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
	/**
	 * Ambiguous bonds are normalized.
	 * @param mol
	 * @return
	 */
	public static float assessComplexity(StereoMolecule mol) {
		final int MAX_BOND_COUNT = 7;
		int bondCount = Math.min(mol.getBonds()/2, MAX_BOND_COUNT);

		if (bondCount < 2)
			return 0;

		mol.ensureHelperArrays(Molecule.cHelperRings);
		StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());
		TreeSet<String> fragmentSet = new TreeSet<String>();
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
		}

		return (float)Math.log(fragmentSet.size()) / bondCount;
		}
	}
