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
 * @author Thomas Sander
 */

package com.actelion.research.chem.prediction;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

public class MolecularShapeCalculator {
	/**
	 * Returns the number of bonds of the shortest path between those two atoms
	 * with the largest topological distance divided by the number of non-hydrogen atoms.
	 * @param mol
	 * @return
	 */
	public static float assessShape(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		if (mol.getAtoms() == 0)
			return -1;
		if (mol.getBonds() == 0)
			return 0;

		int maxLength = 0;
		for (int atom=0; atom<mol.getAtoms(); atom++)
			if (mol.getConnAtoms(atom) == 1 || mol.isRingAtom(atom))
				maxLength = Math.max(maxLength, findHighestAtomDistance(mol, atom));

		return (float)(maxLength+1) / (float)mol.getAtoms();
	}

	/**
	 * Calculates the topological distance to the topologically most remote atom.
	 * @param mol
	 * @param startAtom
	 * @return number of bonds from startAtom to remote atom
	 */
	private static int findHighestAtomDistance(StereoMolecule mol, int startAtom) {
		int[] graphLevel = new int[mol.getAtoms()];
		int[] graphAtom = new int[mol.getAtoms()];

		graphAtom[0] = startAtom;
		graphLevel[startAtom] = 1;

		int current = 0;
		int highest = 0;
		while (current <= highest /* && graphLevel[current] <= maxLength */) {
			int parent = graphAtom[current];
			for (int i=0; i<mol.getConnAtoms(parent); i++) {
				int candidate = mol.getConnAtom(parent, i);
				if (graphLevel[candidate] == 0) {
					graphAtom[++highest] = candidate;
					graphLevel[candidate] = graphLevel[parent]+1;
				}
			}
			current++;
		}
		return graphLevel[graphAtom[highest]] - 1;
	}
}
