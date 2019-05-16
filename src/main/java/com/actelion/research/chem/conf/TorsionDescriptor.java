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
*/

package com.actelion.research.chem.conf;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.Angle;
import com.actelion.research.util.DoubleFormat;

public class TorsionDescriptor implements Comparable<TorsionDescriptor> {
	private static final float TORSION_EQUIVALENCE_TOLERANCE = (float)Math.PI / 12f;

	private float[] mTorsion;

	/**
	 * Calculates an array of all rotatable bonds that can be used
	 * multiple times as parameter to calculateDescriptor().
	 * If the molecule contains marked atoms, these are not considered
	 * part of the molecule, when detecting rotatable bonds. Any non-aromatic
	 * single bond with at least one non-H, non-marked neighbor to either side
	 * is considered a rotatable bond, if none of the bond atoms is marked.
	 * @param mol the molecule behind multiple conformers
	 */
	public static int[] getRotatableBonds(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);

		int count = 0;
		for (int bond=0; bond<mol.getBonds(); bond++)
			if (qualifiesAsDescriptorBond(mol, bond))
				count++;
		int[] rotatableBond = new int[count];
		count = 0;
		for (int bond=0; bond<mol.getBonds(); bond++)
			if (qualifiesAsDescriptorBond(mol, bond))
				rotatableBond[count++] = bond;
		return rotatableBond;
		}

	private static boolean qualifiesAsDescriptorBond(StereoMolecule mol, int bond) {
		if (mol.getBondOrder(bond) != 1 || mol.isAromaticBond(bond))
			return false;

		for (int i=0; i<2; i++) {
			int bondAtom = mol.getBondAtom(i, bond);
			if (mol.isMarkedAtom(bondAtom))
				return false;

			if (mol.getAtomPi(bondAtom) == 2) {
				int otherBondAtom = mol.getBondAtom(1-i, bond);
				if (mol.getAtomPi(otherBondAtom) == 2)
					return false;
				bondAtom = getFirstNonSPAtom(mol, otherBondAtom, bondAtom);
				if (bondAtom == -1 || bondAtom < otherBondAtom) // we only take one bond (that with the lower atom index)
					return false;
				}

			int connAtoms = mol.getConnAtoms(bondAtom);
			if (connAtoms == 1)
				return false;
			int qualifiedConnAtoms = 0;
			for (int j=0; j<connAtoms; j++) {
				int connAtom = mol.getConnAtom(bondAtom, j);
				if (!mol.isMarkedAtom(connAtom))
					qualifiedConnAtoms++;
				}
			if (qualifiedConnAtoms < 2)
				return false;
			}

		return true;
		}

	private static int getFirstNonSPAtom(StereoMolecule mol, int rearAtom, int atom) {
		if (mol.getConnAtoms(atom) == 2) {
			for (int i=0; i<2; i++) {
				int nextAtom = mol.getConnAtom(atom, i);
				if (nextAtom != rearAtom) {
					if (mol.getAtomPi(nextAtom) == 2)
						return getFirstNonSPAtom(mol, atom, nextAtom);
					else
						return nextAtom;
					}
				}
			}
		return -1;
		}

	/**
	 * Creates a TorsionDescriptor the coordinates of the passed molecule using the default method to detect rotatable bonds.
	 * The torsion descriptor is not canonical, unless the passed molecule is canonical.
	 * Rotatable bonds need to carry at least one external non-hydrogen neighbor on each side.
	 * @param mol
	 */
	public TorsionDescriptor(StereoMolecule mol) {
		this(mol, getRotatableBonds(mol));
		}

	/**
	 * Creates a TorsionDescriptor from conformer's coordinates drawing atom connectivity from the conformer's molecule
	 * and using the default method to detect rotatable bonds.
	 * The torsion descriptor is not canonical, unless the passed molecule is canonical.
	 * Rotatable bonds need to carry at least one external non-hydrogen neighbor on each side.
	 * @param conformer
	 */
	public TorsionDescriptor(Conformer conformer) {
		this(conformer, getRotatableBonds(conformer.getMolecule()));
		}

	/**
	 * Creates a TorsionDescriptor from the coordinates of the passed molecule.
	 * The torsion descriptor is not canonical, unless the passed molecule is canonical.
	 * Rotatable bonds need to carry at least one external non-hydrogen neighbor on each side.
	 * @param mol
	 * @param rotatableBond those bonds considered to be rotatable
	 */
	public TorsionDescriptor(StereoMolecule mol, int[] rotatableBond) {
		mTorsion = new float[rotatableBond.length];
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		int[] atom = new int[4];
		for (int i=0; i<rotatableBond.length; i++) {
			atom[1] = mol.getBondAtom(0, rotatableBond[i]);
			if (mol.getAtomPi(atom[1]) == 2) {
				atom[1] = getFirstNonSPAtom(mol, atom[2], atom[1]);
				atom[0] = mol.getConnAtom(atom[1], mol.getAtomPi(mol.getConnAtom(atom[1], 0)) == 2 ? 1 : 0);
				}
			else {
				atom[0] = mol.getConnAtom(atom[1], mol.getConnAtom(atom[1], 0) == atom[2] ? 1 : 0);
				}
			atom[2] = mol.getBondAtom(1, rotatableBond[i]);
			if (mol.getAtomPi(atom[1]) == 2) {
				atom[2] = getFirstNonSPAtom(mol, atom[1], atom[2]);
				atom[3] = mol.getConnAtom(atom[2], mol.getAtomPi(mol.getConnAtom(atom[2], 0)) == 2 ? 1 : 0);
				}
			else {
				atom[3] = mol.getConnAtom(atom[2], mol.getConnAtom(atom[2], 0) == atom[1] ? 1 : 0);
				}
			mTorsion[i] = (float)mol.calculateTorsion(atom);
			}
		}

	/**
	 * Creates a TorsionDescriptor from conformer's coordinates
	 * drawing atom connectivity from the conformer's molecule.
	 * The torsion descriptor is not canonical, unless the passed molecule is canonical.
	 * Rotatable bonds need to carry at least one external non-hydrogen neighbor on each side.
	 * @param conformer
	 */
	public TorsionDescriptor(Conformer conformer, int[] rotatableBond) {
		mTorsion = new float[rotatableBond.length];
		StereoMolecule mol = conformer.getMolecule();
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
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
			mTorsion[i] = (float)conformer.calculateTorsion(atom);
			}
		}

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
			 && dif < 2*(float)Math.PI - TORSION_EQUIVALENCE_TOLERANCE)
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
			 && dif < 2*(float)Math.PI - TORSION_EQUIVALENCE_TOLERANCE)
				return (dif < Math.PI) ^ (mTorsion[i] < td.mTorsion[i]) ? 1 : -1;
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
			sb.append(DoubleFormat.toString(mTorsion[i], 3));
			}
		return sb.toString();
		}
	}
