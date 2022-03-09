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
