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

package com.actelion.research.chem;

public class MoleculeNeutralizer {
	private static final int cInitialMultiChargeFragmentMembers = 16;

	/**
	 * Tries to neutralize all charged atoms of the molecule unless a charged atom
	 * has a neighbour atom with opposite charge. Quarternary ammonium or B(-) won't
	 * be touched. If positive charges remain, this method tries to deprotonate
	 * halogenes and acidic oxygen atoms to neutralize the molecule.
	 * @param mol
	 * @return remaining overall charge, which may be different from 0 if complete neutralization cannot be achieved
	 */
	public static int neutralizeChargedMolecule(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);

		int overallCharge = 0;
		boolean[] hasOppositeChargedNeighbour = new boolean[mol.getAtoms()];
		boolean[] hydrogensForDeletion = null;
		int[] valence = new int[mol.getAtoms()];

		for (int bond = 0; bond < mol.getBonds(); bond++) {
			int atom1 = mol.getBondAtom(0, bond);
			int atom2 = mol.getBondAtom(1, bond);
			int charge1 = mol.getAtomCharge(atom1);
			int charge2 = mol.getAtomCharge(atom2);
			if (charge1 != 0 && charge2 != 0 && ((charge1 < 0) ^ (charge2 < 0))) {
				hasOppositeChargedNeighbour[atom1] = true;
				hasOppositeChargedNeighbour[atom2] = true;
			}
			valence[atom1] += mol.getBondOrder(bond);
			valence[atom2] += mol.getBondOrder(bond);
		}

		for (int chargedAtom = 0; chargedAtom < mol.getAtoms(); chargedAtom++) {
			overallCharge += mol.getAtomCharge(chargedAtom);

			if (mol.getAtomCharge(chargedAtom) == 1) {
				if (mol.getAtomicNo(chargedAtom) == 7) {
					if (!hasOppositeChargedNeighbour[chargedAtom]) {
						if (valence[chargedAtom] <= 3) {
							overallCharge -= 1;
							mol.setAtomCharge(chargedAtom, 0);
							if (mol.getConnAtoms(chargedAtom) != mol.getAllConnAtoms(chargedAtom)) {
								mol.deleteAtom(mol.getConnAtom(chargedAtom, mol.getAllConnAtoms(chargedAtom) - 1));
								mol.ensureHelperArrays(Molecule.cHelperRings);
							}
						}
						else if (mol.isAromaticAtom(chargedAtom)) {
							// For quarternary aromatic charged nitrogen we check the complete annelated
							// aromatic system, whether we can remove a hydrogen from another non-basic
							// aromatic nitrogen. If this is the case, then we need to change single and double
							// bonds to get a valid resonance structure again.
							boolean[] isAromaticAtom = new boolean[mol.getAtoms()];
							boolean[] isAromaticBond = new boolean[mol.getBonds()];
							mol.findRingSystem(chargedAtom, true, isAromaticAtom, isAromaticBond);

							// now we try to locate a pyrrole nitrogen that may serve as an electron pair donor by loosing a hydrogen
							for (int donorAtom = 0; donorAtom < mol.getAtoms(); donorAtom++) {
								if (isAromaticAtom[donorAtom]
										&& mol.getAtomicNo(donorAtom) == 7
										&& mol.getAtomCharge(donorAtom) == 0
										&& valence[donorAtom] == 2) {    // a potential donor for the electron pair is found
									if (removeHydrogenAndDelocalize(mol, isAromaticBond, chargedAtom, donorAtom)) {
										overallCharge -= 1;
										break;
									}
								}
							}
						}
					}
				}
			}
			else if (mol.getAtomCharge(chargedAtom) < 0) {
				if (mol.getAtomicNo(chargedAtom) == 6
						|| mol.getAtomicNo(chargedAtom) == 7
						|| mol.getAtomicNo(chargedAtom) == 8
						|| mol.getAtomicNo(chargedAtom) == 16) {
					if (!hasOppositeChargedNeighbour[chargedAtom]) {
						overallCharge -= mol.getAtomCharge(chargedAtom);
						mol.setAtomCharge(chargedAtom, 0);
					}
					else {
						int[] member = new int[cInitialMultiChargeFragmentMembers];
						member[0] = chargedAtom;
						int noOfMembers = 1;
						int memberToProcess = 0;

						while (memberToProcess < noOfMembers) {
							for (int bond = 0; bond < mol.getAllBonds(); bond++) {
								int atom = -1;
								if (mol.getBondAtom(0, bond) == member[memberToProcess])
									atom = mol.getBondAtom(1, bond);
								else if (mol.getBondAtom(1, bond) == member[memberToProcess])
									atom = mol.getBondAtom(0, bond);

								if (atom == -1)
									continue;

								if (mol.getAtomCharge(atom) != 0) {
									boolean found = false;
									for (int i = 0; i < noOfMembers; i++) {
										if (atom == member[i]) {
											found = true;
											break;
										}
									}

									if (!found) {
										if (noOfMembers == member.length) {
											// for cartridge compatibility use old Java 5 method rather than Arrays.copy():
											int[] copy = new int[2 * member.length];
											System.arraycopy(member, 0, copy, 0, member.length);
											member = copy;
										}

										member[noOfMembers++] = atom;
									}
								}
							}

							memberToProcess++;
						}

						int fragmentsCharge = 0;
						for (int i = 0; i < noOfMembers; i++) {
							fragmentsCharge += mol.getAtomCharge(member[i]);
						}

						if (fragmentsCharge < 0) {
							int leastElectronegative = -1;
							int leastElectronegativity = 99;
							for (int i = 0; i < noOfMembers; i++) {
								if (mol.getAtomCharge(member[i]) < 0) {
									if (leastElectronegativity >
											electronegativity(mol.getAtomicNo(member[i]))) {
										leastElectronegativity =
												electronegativity(mol.
														getAtomicNo(member[i]));
										leastElectronegative = member[i];
									}
								}
							}

							if (leastElectronegative != -1) {
								overallCharge -= mol.getAtomCharge(leastElectronegative);
								mol.setAtomCharge(leastElectronegative, 0);
							}
						}
					}
				}
			}
		}

		if (hydrogensForDeletion != null)
			mol.deleteAtoms(hydrogensForDeletion);

		// if we cannot fully neutralize by proton abstraction (R4N+, pyrylium etc)
		if (overallCharge > 0)
			overallCharge = removeAcidicHydrogens(mol, overallCharge);

		return overallCharge;
	}

	/**
	 * Tries to remove a hydrogen from donorAtom and freshly localize the double bonds
	 * of the flagged delocalized molecule part in order to get rid of the charge on
	 * chargedAtom. If successful the molecule is changed accordingly and true is returned.
	 * Otherwise mol is not changed.
	 *
	 * @param mol
	 * @param isAromaticBond
	 * @param chargedAtom
	 * @param donorAtom
	 * @return true if mol was successfully changed
	 */
	private static boolean removeHydrogenAndDelocalize(StereoMolecule mol, boolean[] isAromaticBond,
												int chargedAtom, int donorAtom)
	{

		// memorize original bond types and copy bond delocalization mask
		int[] oldBondType = new int[mol.getBonds()];
		boolean[] bondMask = new boolean[mol.getBonds()];
		for (int bond = 0; bond < mol.getBonds(); bond++) {
			oldBondType[bond] = mol.getBondType(bond);
			if (isAromaticBond[bond]) {
				bondMask[bond] = true;
				mol.setBondType(bond, Molecule.cBondTypeSingle);
			}
		}

		// The new donor is now part of the delocalized system, but all
		// remaining aromatic atoms without pi-bonds need to be un-flagged
		// to prevent their connected bonds to be touched by the AromaticityResolver.
		for (int atom = 0; atom < mol.getAtoms(); atom++) {
			if (atom != donorAtom
					&& mol.getAtomPi(atom) == 0
					&& bondMask[mol.getConnBond(donorAtom, 0)]) {
				for (int i = 0; i < mol.getConnAtoms(atom); i++) {
					int bond = mol.getConnBond(atom, i);
					mol.setBondType(bond, Molecule.cBondTypeSingle);
					bondMask[bond] = false;
				}
			}
		}

		// The formerly charged nitrogen looses a double bond and
		// is not part of the delocalized system, either.
		for (int i = 0; i < mol.getConnAtoms(chargedAtom); i++) {
			int bond = mol.getConnBond(chargedAtom, i);
			mol.setBondType(bond, Molecule.cBondTypeSingle);
			bondMask[bond] = false;
		}

		if (new AromaticityResolver(mol).locateDelocalizedDoubleBonds(bondMask)) {
			mol.setAtomCharge(chargedAtom, 0);

			// remove explicit hydrogen from donor, if present
			if (mol.getConnAtoms(donorAtom) != mol.getAllConnAtoms(donorAtom)) {
				mol.deleteAtom(mol.getConnAtom(donorAtom, mol.getAllConnAtoms(donorAtom) - 1));
			}

			mol.ensureHelperArrays(Molecule.cHelperRings);
			return true;
		}

		for (int bond = 0; bond < mol.getBonds(); bond++) {
			mol.setBondType(bond, oldBondType[bond]);
		}

		return false;
	}

	private static int removeAcidicHydrogens(StereoMolecule mol, int maxHydrogens) {
		if (maxHydrogens > 0)	// HF
			maxHydrogens = removeHydrogensFromHalogene(mol, maxHydrogens, 9);
		if (maxHydrogens > 0)	// HCl
			maxHydrogens = removeHydrogensFromHalogene(mol, maxHydrogens, 17);
		if (maxHydrogens > 0)	// HBr
			maxHydrogens = removeHydrogensFromHalogene(mol, maxHydrogens, 35);
		if (maxHydrogens > 0)	// HI
			maxHydrogens = removeHydrogensFromHalogene(mol, maxHydrogens, 53);

		if (maxHydrogens > 0) {	// first handle all acidic oxygens vicinal to positively charged atoms
			for (int atom = 0; atom < mol.getAtoms(); atom++) {
				if (mol.getAtomCharge(atom) > 0) {
					boolean found = false;
					for (int i=0; i<mol.getConnAtoms(atom); i++) {
						int connAtom = mol.getConnAtom(atom, i);
						if (mol.getAtomCharge(connAtom) == 0
						 && mol.isElectronegative(connAtom)
						 && mol.getImplicitHydrogens(connAtom) > 0) {
							mol.setAtomCharge(connAtom, -1);
							maxHydrogens--;
							if (maxHydrogens == 0)
								return 0;
							found = true;
							break;
						}
					}
					if (found)
						continue;
				}
			}
		}

		if (maxHydrogens > 0)	// sulfonic & sulfinic acid
			maxHydrogens = removeAcidicHydrogensFromAcid(mol, maxHydrogens, 8, 16);
		if (maxHydrogens > 0)	// sulfonic & sulfinic acid
			maxHydrogens = removeAcidicHydrogensFromAcid(mol, maxHydrogens, 8, 15);
		if (maxHydrogens > 0)	// carboxylic acid
			maxHydrogens = removeAcidicHydrogensFromAcid(mol, maxHydrogens, 8, 6);
		if (maxHydrogens > 0)	// sulfonamids
			maxHydrogens = removeAcidicHydrogensFromAcid(mol, maxHydrogens, 7, 16);

		return maxHydrogens;
	}

	private static int removeHydrogensFromHalogene(StereoMolecule mol, int maxHydrogens, int atomicNo) {
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.getAtomicNo(atom) == atomicNo
			 && mol.getAtomCharge(atom) == 0
			 && mol.getConnAtoms(atom) == 0) {
				mol.setAtomCharge(atom, -1);
				maxHydrogens--;
				if (maxHydrogens == 0)
					return 0;
			}
		}
		return maxHydrogens;
	}

	private static int removeAcidicHydrogensFromAcid(StereoMolecule mol, int maxHydrogens, int atomicNo, int centralAtomicNo) {
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.getAtomicNo(atom) == atomicNo
			 && mol.getAtomCharge(atom) == 0
			 && mol.getImplicitHydrogens(atom) > 0) {
				boolean deprotonated = false;
				for (int i=0; i<mol.getConnAtoms(atom); i++) {
					if (mol.getConnBondOrder(atom, i) == 1) {
						int centralAtom = mol.getConnAtom(atom, i);
						if (mol.getAtomicNo(centralAtom) == centralAtomicNo) {
							boolean oxoFound = false;
							boolean deprotonatedFound = false;
							for (int j=0; j<mol.getConnAtoms(centralAtom); j++) {
								int connAtom = mol.getConnAtom(centralAtom, j);
								if (mol.getAtomCharge(connAtom) < 0) {
									deprotonatedFound = true;
									break;
								}
								if (connAtom != atom
								 && mol.getAtomicNo(connAtom) == 8
								 && mol.getConnBondOrder(centralAtom, j) == 2) {
									oxoFound = true;
								}
							}
							if (!deprotonatedFound && oxoFound) {
								mol.setAtomCharge(atom, -1);
								maxHydrogens--;
								if (maxHydrogens == 0)
									return 0;
								deprotonated = true;
							}
						}
					}
				if (deprotonated)
					break;
				}
			}
		}
		return maxHydrogens;
	}

	private static int electronegativity(int atomicNo) {
		switch (atomicNo) {
			case 6:
				return 1; // C
			case 53:
				return 2; // I
			case 33:
				return 3; // As
			case 34:
				return 4; // Se
			case 35:
				return 5; // Br
			case 15:
				return 6; // P
			case 16:
				return 7; // S
			case 17:
				return 8; // Cl
			case 7:
				return 9; // N
			case 8:
				return 10; // O
			case 9:
				return 11; // F
		}
		return 0;
	}
}
