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

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

import java.util.ArrayList;

public class PlaneRule extends ConformationRule {
	private int[] mPlaneAtom;	// these are the atoms that define the plane

	public PlaneRule(int[] atom, StereoMolecule mol) {
		super(atom);
		int count = 0;
		for (int i=0; i<atom.length; i++)
			if (mol.getConnAtoms(atom[i]) != 1)
				count++;
		if (count > 2) {
			mPlaneAtom = new int[count];
			count = 0;
			for (int i=0; i<atom.length; i++)
				if (mol.getConnAtoms(atom[i]) != 1)
					mPlaneAtom[count++] = atom[i];
			}
		else {
			mPlaneAtom = atom;
			}
		}

	@Override
	public int getRuleType() {
		return RULE_TYPE_PLANE;
		}

    public static void calculateRules(ArrayList<ConformationRule> ruleList, StereoMolecule mol) {
		boolean[] isFlatBond = new boolean[mol.getBonds()];

		RingCollection ringSet = mol.getRingSet();
		boolean[] isAromaticRing = new boolean[ringSet.getSize()];
		ringSet.determineAromaticity(isAromaticRing, new boolean[ringSet.getSize()], new int[ringSet.getSize()], true);
		for (int ring=0; ring<ringSet.getSize(); ring++)
			if (isAromaticRing[ring])
				for (int i=0; i<ringSet.getRingSize(ring); i++)
					isFlatBond[ringSet.getRingBonds(ring)[i]] = true;

		int[] atomicNo = new int[2];
		for (int bond=0; bond<mol.getBonds(); bond++) {
			int bondAtom1 = mol.getBondAtom(0, bond);
			int bondAtom2 = mol.getBondAtom(1, bond);
		    atomicNo[0] = mol.getAtomicNo(bondAtom1);
		    atomicNo[1] = mol.getAtomicNo(bondAtom2);
			isFlatBond[bond] |= (mol.isAromaticBond(bond)
							|| (mol.getBondOrder(bond) == 2
							 && atomicNo[0] <= 8 && atomicNo[1] <= 8
							 && mol.getAtomPi(bondAtom1) == 1
							 && mol.getAtomPi(bondAtom2) == 1
							 && mol.getAllConnAtoms(bondAtom1) > 1
							 && mol.getAllConnAtoms(bondAtom2) > 1));
			if (!isFlatBond[bond]) {
					// check if bond is an amide or ester bond
				if (mol.getBondOrder(bond) == 1) {
					for (int i=0; i<2; i++) {
						if ((atomicNo[i] == 7 || atomicNo[i] == 8) && atomicNo[1-i] == 6) {
							int carbon = mol.getBondAtom(1-i, bond);
							for (int j=0; j<mol.getConnAtoms(carbon); j++) {
								if (mol.getConnBondOrder(carbon, j) == 2) {
									int connAtom = mol.getConnAtom(carbon, j);
									if (mol.getAtomicNo(connAtom) == 7
									 || mol.getAtomicNo(connAtom) == 8
									 || mol.getAtomicNo(connAtom) == 16) {
										isFlatBond[bond] = true;
										break;
										}
									}
								}
							}
						}
					}
				}
			}

		boolean[] isFlatBondAtom = new boolean[mol.getAtoms()];
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (isFlatBond[bond]) {
				isFlatBondAtom[mol.getBondAtom(0, bond)] = true;
				isFlatBondAtom[mol.getBondAtom(1, bond)] = true;
				}
			}

		int[] fragmentAtom = new int[mol.getAllAtoms()];
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (isFlatBond[bond]) {
				fragmentAtom[0] = mol.getBondAtom(0, bond);
				int count = getFlatFragmentAtoms(fragmentAtom, isFlatBond, mol);
				int[] atomList = new int[count];
				for (int i=0; i<count; i++)
					atomList[i] = fragmentAtom[i];
				ruleList.add(new PlaneRule(atomList, mol));
				}
			}

		// this covers flat atoms that are not part of flat bonds, e.g. C in C=S, B if not B-, etc.
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (!isFlatBondAtom[atom]) {
				if ((mol.getAtomicNo(atom) == 5 && mol.getAtomCharge(atom) == 0 && mol.getAllConnAtoms(atom) <= 3)
				 || (mol.getAtomicNo(atom) <= 8 && mol.getAtomPi(atom) == 1 && mol.getAllConnAtoms(atom) > 1)
				 || (mol.isFlatNitrogen(atom) && mol.getAtomPi(atom) != 2 && mol.getAllConnAtoms(atom) > 1)) {
					int[] atomList = new int[1+mol.getAllConnAtoms(atom)];
					for (int i=0; i<mol.getAllConnAtoms(atom); i++)
						atomList[i] = mol.getConnAtom(atom, i);
					atomList[mol.getAllConnAtoms(atom)] = atom;
					ruleList.add(new PlaneRule(atomList, mol));
					}
				}
			}
    	}

	private static int getFlatFragmentAtoms(int[] fragmentAtom, boolean[] isFlatBond, StereoMolecule mol) {
		// locate all atoms connected directly via flat bonds
		boolean[] isFragmentMember = new boolean[mol.getAllAtoms()];
		isFragmentMember[fragmentAtom[0]] = true;
		int current = 0;
		int highest = 0;
	 	while (current <= highest && mol.getAtomPi(fragmentAtom[current]) < 2) {
			for (int i=0; i<mol.getConnAtoms(fragmentAtom[current]); i++) {
				int candidateAtom = mol.getConnAtom(fragmentAtom[current], i);
				int candidateBond = mol.getConnBond(fragmentAtom[current], i);
				if (isFlatBond[candidateBond]) {
					if (!isFragmentMember[candidateAtom]) {
						fragmentAtom[++highest] = candidateAtom;
						isFragmentMember[candidateAtom] = true;
						}
					isFlatBond[candidateBond] = false;
					}
				}
			current++;
			}
	
			// attach first sphere of atoms connected via non-flat bonds
		for (int i=highest; i>=0; i--) {
			if (mol.getAtomicNo(fragmentAtom[i]) <= 8) {
				for (int j=0; j<mol.getAllConnAtoms(fragmentAtom[i]); j++) {
					int connAtom = mol.getConnAtom(fragmentAtom[i], j);
					if (!isFragmentMember[connAtom] && mol.getConnBondOrder(fragmentAtom[i], j) != 0) {
						fragmentAtom[++highest] = connAtom;
						isFragmentMember[connAtom] = true;
						}
					}
				}
			}
	
		return highest+1;
		}

	@Override
	public boolean apply(Conformer conformer, double cycleFactor) {
		Coordinates cog = new Coordinates();
		Coordinates n = new Coordinates();
		double[][] coords = new double[mPlaneAtom.length][3];
		ConformationRule.calculateNearestPlane(conformer, mPlaneAtom, cog, n, coords);

		for (int i=0; i<mAtom.length; i++) {
			double distance = -(n.x*(conformer.getX(mAtom[i]) - cog.x)
							  + n.y*(conformer.getY(mAtom[i]) - cog.y)
							  + n.z*(conformer.getZ(mAtom[i]) - cog.z));
			moveGroup(conformer, mAtom[i], mAtom, 0.5 * distance * cycleFactor * n.x,
												  0.5 * distance * cycleFactor * n.y,
												  0.5 * distance * cycleFactor * n.z);
			}

		return true;
		}

	@Override
	public double addStrain(Conformer conformer, double[] atomStrain) {
		Coordinates cog = new Coordinates();
		Coordinates n = new Coordinates();
		double[][] coords = new double[mAtom.length][3];
		ConformationRule.calculateNearestPlane(conformer, mAtom, cog, n, coords);

		double totalStrain = 0;
		for (int i=0; i<mAtom.length; i++) {
			double distance = -(n.x*coords[i][0] + n.y*coords[i][1] + n.z*coords[i][2]);
			double panalty = 20.0 * distance*distance;
			if (atomStrain != null)
				atomStrain[mAtom[i]] += panalty;
			totalStrain += panalty;
			}

		return totalStrain;
		}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("plane rule:");
		super.addAtomList(sb);
		return sb.toString();
		}
	}
