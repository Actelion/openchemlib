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

public class TorsionDescriptorHelper {
	private static final int SYMMETRY_360 = 0;  // 0 -> 359 degrees
	private static final int SYMMETRY_180 = 1;  // 0 -> 179 (equal to 180 -> 359)
	private static final int SYMMETRY_120 = 2;  // 0 -> 119 (equals 120 -> 239 and 240 -> 359)
	private static final int SYMMETRY_60  = 3;  // 0 -> 59 (60 -> 119, 120 -> 179, 180 -> 239, 240 -> 299, 300 -> 359)

	private static double[] SYMMETRY_REDUNDANCY_ANGLE = { 2*Math.PI, Math.PI, 2*Math.PI/3, Math.PI/3};

	private static final int HALF_SYMMETRY_C1 = 0;  // three distinct terminal neighbors
	private static final int HALF_SYMMETRY_D1 = 1;  // e.g. single terminal neighbor or two equal sp3 neighbors
	private static final int HALF_SYMMETRY_D2 = 2;  // two equal sp2 neighbors
	private static final int HALF_SYMMETRY_D3 = 3;  // for simplicity reasons this is covered by D1

	private static final int[][] SYMMETRY =
		  { { SYMMETRY_360, SYMMETRY_360, SYMMETRY_360, SYMMETRY_120 },
			{ SYMMETRY_360, SYMMETRY_360, SYMMETRY_180, SYMMETRY_120 },
			{ SYMMETRY_360, SYMMETRY_180, SYMMETRY_180, SYMMETRY_60  },
			{ SYMMETRY_120, SYMMETRY_120, SYMMETRY_60 , SYMMETRY_120 } };

	private StereoMolecule mMol;
	private int[] mRotatableBond;
	private int[][] mAtomSequence;
	private int[][] mRearAtom;
	private int[] mSymmetryClass;

	public TorsionDescriptorHelper(StereoMolecule mol) {
		mMol = mol;
		mMol.ensureHelperArrays(Molecule.cHelperSymmetrySimple);
		mRotatableBond = findRotatableBonds(mol);
		findAtomSequences();
		}

	public TorsionDescriptorHelper(StereoMolecule mol, int[] rotatableBond) {
		mMol = mol;
		mMol.ensureHelperArrays(Molecule.cHelperSymmetrySimple);
		mRotatableBond = rotatableBond;
		findAtomSequences();
		}

	/**
	 * Creates a TorsionDescriptor from the coordinates of the molecule passed to the constructor
	 * using the default method to detect rotatable bonds.
	 * The torsion descriptor is not canonical, unless the passed molecule is canonical.
	 * Rotatable bonds need to carry at least one external non-hydrogen neighbor on each side.
	 */
	public TorsionDescriptor getTorsionDescriptor() {
		float[] torsion = new float[mRotatableBond.length];
		float[] maxTorsion = new float[mRotatableBond.length];

		for (int i=0; i<mRotatableBond.length; i++) {
			torsion[i] = (float) getNormalizedTorsion(TorsionDB.calculateTorsionExtended(mMol, mAtomSequence[i]), mSymmetryClass[i]);
			maxTorsion[i] = (float)SYMMETRY_REDUNDANCY_ANGLE[mSymmetryClass[i]];
			}

		return new TorsionDescriptor(torsion, maxTorsion);
		}

	/**
	 * Creates a TorsionDescriptor from the coordinates of the passed conformer assuming that its
	 * underlying molecule is the same that was passed to this TorsionDescriptorHelper's constructor.
	 * This TorsionDescriptorHelper uses the default method to detect rotatable bonds.
	 * The torsion descriptor is not canonical, unless the passed molecule is canonical.
	 * Rotatable bonds need to carry at least one external non-hydrogen neighbor on each side.
	 */
	public TorsionDescriptor getTorsionDescriptor(Conformer conformer) {
		float[] torsion = new float[mRotatableBond.length];
		float[] maxTorsion = new float[mRotatableBond.length];

		for (int i=0; i<mRotatableBond.length; i++) {
			torsion[i] = (float) getNormalizedTorsion(TorsionDB.calculateTorsionExtended(conformer, mAtomSequence[i]), mSymmetryClass[i]);
			maxTorsion[i] = (float)SYMMETRY_REDUNDANCY_ANGLE[mSymmetryClass[i]];
			}

		return new TorsionDescriptor(torsion, maxTorsion);
		}

	/**
	 * Calculates an array of all rotatable bonds that can be used
	 * multiple times as parameter to calculateDescriptor().
	 * If the molecule contains marked atoms, these are not considered
	 * part of the molecule, when detecting rotatable bonds. Any non-aromatic,
	 * non-3-ring single bond with at least one non-H, non-marked neighbor to either side
	 * is considered a rotatable bond, if none of the bond atoms is marked.
	 * An exception are linear chains of sp-hybridized atoms, of which not more than one
	 * bond is considered rotatable.
	 */
	public static int[] findRotatableBonds(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);

		boolean[] isRotatableBond = new boolean[mol.getBonds()];
		int count = 0;
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (qualifiesAsDescriptorBond(mol, bond)) {
				isRotatableBond[bond] = true;
				count++;
				}
			}

		int[] rotatableBond = new int[count];

		count = 0;
		for (int bond=0; bond<mol.getBonds(); bond++)
			if (isRotatableBond[bond])
				rotatableBond[count++] = bond;

		return rotatableBond;
		}

	private static boolean qualifiesAsDescriptorBond(StereoMolecule mol, int bond) {
		if (mol.getBondOrder(bond) != 1
		 || mol.isAromaticBond(bond)
		 || mol.getBondRingSize(bond) == 3)
			return false;

		int[] bondAtom = new int[2];
		for (int i=0; i<2; i++) {
			bondAtom[i] = mol.getBondAtom(i, bond);
			if (mol.isMarkedAtom(bondAtom[i])
			 || mol.getNonHydrogenNeighbourCount(bondAtom[i]) <= 1)
				return false;
			}

		if (isLinearAtom(mol, bondAtom[0])
		 || isLinearAtom(mol, bondAtom[1])) {
			int[] maxBond = new int[1];
			maxBond[0] = bond;

			int[][] chainEndAtom = new int[2][];
			for (int i=0; i<2; i++) {
				if (isLinearAtom(mol, bondAtom[i])) {
					chainEndAtom[i] = new int[2];
					chainEndAtom[i][0] = bondAtom[1-i];
					chainEndAtom[i][1] = bondAtom[i];
					if (!getFirstNonSPAtom(mol, bondAtom[1-i], chainEndAtom[i], maxBond))
						return false;
					}
				}

			if (maxBond != null && bond != maxBond[0])	// in case of linear chains, we only consider the bond with the highest index rotatable
				return false;

			// now we replace the bondAtom by the end-of-chain atom
			for (int i=0; i<2; i++)
				if (chainEndAtom[i] != null)
					bondAtom[i] = chainEndAtom[i][1];
			}

		for (int i=0; i<2; i++) {
			int connAtoms = mol.getConnAtoms(bondAtom[i]);
			if (connAtoms == 1)
				return false;
			int qualifiedConnAtoms = 0;
			for (int j=0; j<connAtoms; j++) {
				int connAtom = mol.getConnAtom(bondAtom[i], j);
				if (!mol.isMarkedAtom(connAtom))
					qualifiedConnAtoms++;
				}
			if (qualifiedConnAtoms < 2)
				return false;
			}

		return true;
		}

	/**
	 * @param mol
	 * @param atom
	 * @return whether the atom is a sp-hybridized atom with exactly two non-hydrogen neighbours
	 */
	private static boolean isLinearAtom(StereoMolecule mol, int atom) {
		return mol.getConnAtoms(atom) == 2 && mol.getAtomPi(atom) == 2 && mol.getAtomicNo(atom) <= 7;
	}

	/**
	 * Stepwise walks along the sp-atom chain starting from the connected atoms in the atom array
	 * away from atom[0] and updates the atom[] array with the atom indexes of the neighbor bond
	 * until atom[1] is a non-sp-atom and atom[0] is the last sp-atom seen of the chain,
	 * (which, of course, is a direct neighbor of atom[1]).
	 * @param mol
	 * @param atom contains two connected atoms, of which atom[1] is sp-hybridized
	 * @param maxBond null or contains the highest non-marked rotatable bond in the chain
	 * @return false if the end of the chain is a sp-hybridized atom or in case of a cyclic chain
	 */
	private static boolean getFirstNonSPAtom(StereoMolecule mol, int startAtom, int[] atom, int[] maxBond) {
		for (int i=0; i<2; i++) {
			int nextAtom = mol.getConnAtom(atom[1], i);
			if (nextAtom != atom[0]) {
				if (nextAtom == startAtom)  // we found a cycle
					return false;

				int nextBond = mol.getConnBond(atom[1], i);

				atom[0] = atom[1];
				atom[1] = nextAtom;

				if (mol.getConnAtoms(nextAtom) == 1)
					return false;

				if (maxBond != null && !mol.isMarkedAtom(atom[0]) && !mol.isMarkedAtom(atom[1]))
					maxBond[0] = Math.max(maxBond[0], nextBond);

				if (!isLinearAtom(mol, nextAtom))
					return true;

				return getFirstNonSPAtom(mol, startAtom, atom, maxBond);
				}
			}
		return false;	// should never happen
		}

	private void findAtomSequences() {
		mAtomSequence = new int[mRotatableBond.length][4];
		mRearAtom = new int[mRotatableBond.length][2];
		mSymmetryClass = new int[mRotatableBond.length];

		int[] atom = new int[2];	// rear bond atom, front bond atom, first conn atom
		for (int i=0; i<mRotatableBond.length; i++) {
			for (int j=0; j<2; j++) {
				atom[0] = mMol.getBondAtom(1-j, mRotatableBond[i]);
				atom[1] = mMol.getBondAtom(j, mRotatableBond[i]);

				if (isLinearAtom(mMol, atom[1]))
					getFirstNonSPAtom(mMol, atom[0], atom, null);

				mAtomSequence[i][1+j] = atom[1];
				mRearAtom[i][j] = atom[0];
				}

			int halfSymmetry1 = getHalfSymmetry(mAtomSequence[i][1], mRearAtom[i][0]);
			int halfSymmetry2 = getHalfSymmetry(mAtomSequence[i][2], mRearAtom[i][1]);
			mSymmetryClass[i] = SYMMETRY[halfSymmetry1][halfSymmetry2];

			// update atom in atom sequence, where we should consider virtual torsions or need to select a particular neighbour atom
			mAtomSequence[i][0] = chooseSequenceEndAtom(mAtomSequence[i][1], mRearAtom[i][0], halfSymmetry1);
			mAtomSequence[i][3] = chooseSequenceEndAtom(mAtomSequence[i][2], mRearAtom[i][1], halfSymmetry2);
			}
		}

	/**
	 * Determines the symmetry of one end of the 4-atom sequence,
	 * which may be one of:
	 * HALF_SYMMETRY_C1: not symmetric due to stereo center or tetrahedral nitrogen.
	 * HALF_SYMMETRY_D1: mirror plane only due to
	 * - one terminal atom only
	 * - two distinct atoms at sp2 center
	 * - exactly 2 symmetrical atoms at sp3 center.
	 * HALF_SYMMETRY_D2: two symmetrical atoms at sp2 center.
	 * HALF_SYMMETRY_D3: three symmetrical atoms at sp3 center.
	 * The fragment's helper array level should be cHelperSymmetrySimple.
	 * @param atom one of the bond atoms of the rotatable bond
	 * @param rearAtom the rear atom of the rotatable bond
	 * @return
	 */
	private int getHalfSymmetry(int atom, int rearAtom) {
		if (mMol.getConnAtoms(atom) == 2)
			return HALF_SYMMETRY_D1;

		int[] connAtom = getTerminalAtoms(atom, rearAtom);

		if (mMol.getConnAtoms(atom) == 3) {
			if (mMol.getSymmetryRank(connAtom[0]) == mMol.getSymmetryRank(connAtom[1]))
				return isFlatAtom(atom) ? HALF_SYMMETRY_D2 : HALF_SYMMETRY_D1;
			else
				return isFlatAtom(atom) ? HALF_SYMMETRY_D1 : HALF_SYMMETRY_C1;
		}

		if (mMol.getConnAtoms(atom) == 4) {
			// two equal ranks with additional neighbor that will serve as reference atom
			if (mMol.getSymmetryRank(connAtom[0]) == mMol.getSymmetryRank(connAtom[1])
			 && mMol.getSymmetryRank(connAtom[0]) == mMol.getSymmetryRank(connAtom[2]))
				return HALF_SYMMETRY_D3;

			if (mMol.getSymmetryRank(connAtom[0]) == mMol.getSymmetryRank(connAtom[1])
			 || mMol.getSymmetryRank(connAtom[0]) == mMol.getSymmetryRank(connAtom[2])
			 || mMol.getSymmetryRank(connAtom[1]) == mMol.getSymmetryRank(connAtom[2]))
				return HALF_SYMMETRY_D1;
		}

		return HALF_SYMMETRY_C1;
	}

	/**
	 * Checks whether atom is sp2 hybridized.
	 * Amide nitrogens are also considered to be sp2.
	 * @param atom
	 * @return
	 */
	private boolean isFlatAtom(int atom) {
		if ((mMol.getAtomPi(atom) == 1 && mMol.getAtomicNo(atom) < 10)
		  || mMol.isAromaticAtom(atom)
		  || mMol.isFlatNitrogen(atom))
			return true;

		return false;
		}

	private int[] getTerminalAtoms(int atom, int rearAtom) {
		int index = 0;
		int[] connAtom = new int[mMol.getConnAtoms(atom)-1];
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (mMol.getConnAtom(atom, i) != rearAtom)
				connAtom[index++] = mMol.getConnAtom(atom, i);
		return connAtom;
		}

	private int chooseSequenceEndAtom(int rootAtom, int rearAtom, int halfSymmetry) {
		if (halfSymmetry == HALF_SYMMETRY_D1 && !isFlatAtom(rootAtom)) {
			if (mMol.getConnAtoms(rootAtom) == 3)	// two symmetrical neighbours at sp3
				return -1;

			if (mMol.getConnAtoms(rootAtom) == 3) {	// three neighbours of which two are symmetrical
				int[] connAtom = getTerminalAtoms(rootAtom, rearAtom);
				if (mMol.getSymmetryRank(connAtom[0]) == mMol.getSymmetryRank(connAtom[1]))
					return connAtom[2];
				if (mMol.getSymmetryRank(connAtom[0]) == mMol.getSymmetryRank(connAtom[2]))
					return connAtom[1];
				else
					return connAtom[0];
				}
			}

		// return highest ranking neighbor or rootAtom that is not rearAtom
		int maxRank = -1;
		int maxRankAtom = -1;
		for (int i=0; i<mMol.getConnAtoms(rootAtom); i++) {
			int connAtom = mMol.getConnAtom(rootAtom, i);
			if (connAtom != rearAtom && maxRank < mMol.getSymmetryRank(connAtom)) {
				maxRank = mMol.getSymmetryRank(connAtom);
				maxRankAtom = connAtom;
				}
			}

		return maxRankAtom;
		}

	/**
	 * Normalizes a torsion angle considering the rotatable bonds symmetry type
	 * by returning the lowest symmetrically equivalent torsion that is >= 0.
	 * @param angle
	 * @param symmetryClass
	 * @return angle within native range of symmetry type
	 */
	private static double getNormalizedTorsion(double angle, int symmetryClass) {
		double limit = SYMMETRY_REDUNDANCY_ANGLE[symmetryClass] / 2;

		while (angle < -limit)
			angle += 2*Math.PI;

		while (angle >= limit)
			angle -= SYMMETRY_REDUNDANCY_ANGLE[symmetryClass];

		return angle;
		}
	}
