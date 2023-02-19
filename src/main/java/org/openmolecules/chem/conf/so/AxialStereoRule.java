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

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Handles atropisomers as well as chiral allenes
 */
public class AxialStereoRule extends ConformationRule {
	private int[] mRotatableAtom;
	private boolean mInverse, mPositiveTorsionWanted;

	public AxialStereoRule(StereoMolecule mol, int[] atom, int[] rearArom, boolean isInRing, boolean positiveTorsionWanted) {
		super(atom);
		mRotatableAtom = getRotatableAtoms(mol, isInRing, atom, rearArom);
		mPositiveTorsionWanted = positiveTorsionWanted;
		}

    public static void calculateRules(ArrayList<ConformationRule> ruleList, StereoMolecule mol) {
		for (int centralAtom=0; centralAtom<mol.getAtoms(); centralAtom++) {
			if (mol.isCentralAlleneAtom(centralAtom)) {
			    int parity = mol.getAtomParity(centralAtom);
			    if (parity == Molecule.cAtomParity1 || parity == Molecule.cAtomParity2) {
				    int[] atomList = new int[4];
				    int[] rearAtom = new int[2];
				    for (int i=0; i<mol.getConnAtoms(centralAtom); i++) {
					    int atom = mol.getConnAtom(centralAtom, i);
					    int lowConn = Integer.MAX_VALUE;
					    for (int j=0; j<mol.getConnAtoms(atom); j++) {
						    int connAtom = mol.getConnAtom(atom, j);
						    if (connAtom != centralAtom && lowConn > connAtom)
							    lowConn = connAtom;
					    }
					    atomList[3*i] = lowConn;
					    atomList[1+i] = atom;
					    rearAtom[i] = centralAtom;
				    }
if (atomList[0] == Integer.MAX_VALUE || atomList[3] == Integer.MAX_VALUE) { // TODO remove this
 System.out.println("Unexpected MAX_VALUE idcode:"+new Canonizer(mol).getIDCode());
 break;
}
				    ruleList.add(new AxialStereoRule(mol, atomList, rearAtom, mol.isRingAtom(centralAtom), parity == Molecule.cAtomParity1));
			    }
		    }
	    }

		for (int bond=0; bond<mol.getBonds(); bond++) {
			int parity = mol.getBondParity(bond);
			if ((parity == Molecule.cBondParityEor1 || parity == Molecule.cBondParityZor2) && mol.getBondOrder(bond) != 2) {
				int[] atomList = new int[4];
				int[] rearAtom = new int[2];
				for (int i=0; i<2; i++) {
					int atom = mol.getBondAtom(i, bond);
					int rear = mol.getBondAtom(1-i, bond);
					int lowConn = Integer.MAX_VALUE;
					for (int j=0; j<mol.getConnAtoms(atom); j++) {
						int connAtom = mol.getConnAtom(atom, j);
						if (connAtom != rear && lowConn > connAtom)
							lowConn = connAtom;
						}
					atomList[3*i] = lowConn;
					atomList[1+i] = atom;
					rearAtom[i] = mol.getBondAtom(1-i, bond);
					}

				ruleList.add(new AxialStereoRule(mol, atomList, rearAtom, mol.isRingBond(bond), parity == Molecule.cBondParityZor2));
				}
			}
		}

	/**
	 * Detects all conne
	 * @param mol
	 * @param firstAtom
	 * @param rearAtom
	 * @return
	 */
	private int[] getAromaticFragmentAtoms(StereoMolecule mol, int firstAtom, int rearAtom) {
		int[] aromAtom = new int[mol.getAtoms()];
		boolean[] isMember = new boolean[mol.getAtoms()];
		isMember[rearAtom] = true;  // to prevent including the rearAtom
		isMember[firstAtom] = true;
		aromAtom[0] = firstAtom;
		int current = 0;
		int highest = 1;
		while (current <= highest) {
			for (int i=0; i<mol.getConnAtoms(aromAtom[current]); i++) {
				if (mol.isAromaticBond(mol.getConnBond(aromAtom[current], i))) {
					int candidate = mol.getConnAtom(aromAtom[current], i);
					if (!isMember[candidate]) {
						aromAtom[++highest] = candidate;
						isMember[candidate] = true;
						}
					}
				}
			current++;
			}
		for (int i=0; i<=highest; i++) {
			for (int j=0; j<mol.getConnAtoms(aromAtom[i]); j++) {
				int candidate = mol.getConnAtom(aromAtom[i], j);
				if (!isMember[candidate] && !mol.isRingBond(mol.getConnBond(aromAtom[i], j)))
					mol.getSubstituent(aromAtom[i], candidate, isMember, null, null);
				}
			}
		isMember[rearAtom] = false;
		for (int i=0; i<=highest; i++)
			isMember[aromAtom[i]] = false;
		for (int atom=0; atom<isMember.length; atom++)
			if (isMember[atom])
				aromAtom[++highest] = atom;

		return Arrays.copyOf(aromAtom, highest+1);
		}

	private int[] getRotatableAtoms(StereoMolecule mol, boolean isInRing, int[] atomList, int[] rearAtom) {
		// If both aromatic systems are connected by another chain of atoms, then we rotate just one aromatic system
		// with all connected substituents that don't connect to the other side.
		if (isInRing)
			return getAromaticFragmentAtoms(mol, atomList[1], atomList[2]);

		// Otherwise we rotate the smaller side of the molecule seen from the binap bond perspective
		boolean[][] isMemberAtom = new boolean[2][mol.getAllAtoms()];
		int[] count = new int[2];
		for (int i=0; i<2; i++)
			count[i] = mol.getSubstituent(atomList[1+i], rearAtom[i], isMemberAtom[i], null, null);

		int smallerSubstituentIndex = (count[0] < count[1]) ? 0 : 1;
		int[] rotatableAtom = new int[count[smallerSubstituentIndex]];
		int index = 0;
		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			if (isMemberAtom[smallerSubstituentIndex][atom])
				rotatableAtom[index++] = atom;

		mInverse = (smallerSubstituentIndex == 1);

		// Torsion atomList[1] is bondAtom 0 and atomList[2] is bondAtom 1.
		// Default is to rotate atomList[3] and its connected atom environment.
		// If atomList[0] and environment is smaller, then mInverse is true and
		// we need to rotate with inverted torsion angles.
		return rotatableAtom;
		}

	@Override
	public boolean apply(Conformer conformer, double cycleFactor) {
		double torsion = conformer.calculateTorsion(mAtom);
		if (torsion == 0 || (mPositiveTorsionWanted == (torsion > 0)))
			return false;

		Coordinates unit = conformer.getCoordinates(mAtom[2]).subC(conformer.getCoordinates(mAtom[1]));
		unit.unit();

		double rotation = Math.abs(torsion) < Math.PI/2 ? -2.0*torsion : 2.0*(Math.PI - torsion);
		if (mInverse)
			rotation = -rotation;

		for (int atom:mRotatableAtom)
			rotateAtom(conformer, atom, conformer.getCoordinates(mAtom[1]), unit, rotation);

		return true;
		}

	@Override
	public double addStrain(Conformer conformer, double[] atomStrain) {
		double strain = conformer.calculateTorsion(mAtom) < 0 ? 1.0 : 0.0;
		if (atomStrain != null)
			for (int atom:mAtom)
				atomStrain[atom] += strain / 4;
		return strain;
		}

	@Override
	public int getRuleType() {
		return RULE_TYPE_BINAP;
		}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("binap rule ("+(mPositiveTorsionWanted?"positive":"negative")+" torsion):");
		super.addAtomList(sb);
		return sb.toString();
		}
	}
