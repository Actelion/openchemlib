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

package com.actelion.research.chem;

import java.util.Arrays;

public class ScaffoldHelper {

	/**
	 * Modifies mol in-place to represent the Murcko scaffold from the original mol.
	 * If skeletonOnly is true, then all bond types are reduced to non-stereo-single bonds
	 * and all atoms are converted to carbon atoms.
	 * @param mol may be an empty molecule if there are no ring systems
	 * @param skeletonOnly reduce to simple skeleton
	 */
	public static void createMurckoScaffold(StereoMolecule mol, boolean skeletonOnly) {
		boolean[] isScaffoldMember = findMurckoScaffold(mol);
		if (isScaffoldMember == null) {
			mol.deleteMolecule();
			return;
			}
		for (int atom=0; atom<mol.getAtoms(); atom++)
			if (!isScaffoldMember[atom])
				mol.markAtomForDeletion(atom);
		mol.deleteMarkedAtomsAndBonds();
		if (skeletonOnly)
			makeSkeleton(mol);
		}

	/**
	 * Creates Murcko scaffold as new StereoMolecule from mol without touching mol.
	 * If skeletonOnly is true, then all bond types are reduced to non-stereo-single bonds
	 * and all atoms are converted to carbon atoms.
	 * @param mol
	 * @param skeletonOnly reduce to simple skeleton
	 * @return Murcko scaffold or null, if mol does not contain any rings
	 */
	public static StereoMolecule getMurckoScaffold(StereoMolecule mol, boolean skeletonOnly) {
		boolean[] isScaffoldMember = findMurckoScaffold(mol);
		if (isScaffoldMember == null)
			return null;
		StereoMolecule scaffold = new StereoMolecule(mol.getAtoms(), mol.getBonds());
		mol.copyMoleculeByAtoms(scaffold, isScaffoldMember, false, null);
		if (skeletonOnly)
			makeSkeleton(scaffold);
		return scaffold;
		}

	private static void makeSkeleton(StereoMolecule mol) {
		for (int bond=0; bond<mol.getAllBonds(); bond++)
			mol.setBondType(bond, Molecule.cBondTypeSingle);
		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			mol.setAtomicNo(atom, 6);
		}

	/**
	 * This method flags all atoms that are a member of the Murcko fragment of mol.
	 * The Murcko fragment is detected in three steps:
	 * - flag all ring atoms
	 * - recursively flag all atoms that are on a direct path connecting any two flagged atoms
	 * - recursively flag all atoms that are connected to any flagged atom via double bond
	 * Hydrogen atoms are neglected.
	 * @param mol
	 * @return null if mol doesn't contain any rings
	 */
	public static boolean[] findMurckoScaffold(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		if (mol.getRingSet().getSize() == 0)
			return null;

		boolean[] isMember = new boolean[mol.getAtoms()];

		// add all ring members
		for (int atom=0; atom<mol.getAtoms(); atom++)
			if (mol.isRingAtom(atom))
				isMember[atom] = true;

		// recursively add all exocyclic chains that lead to a member atom
		boolean found;
		do {
			found = false;
			for (int atom=0; atom<mol.getAtoms(); atom++) {
				if (isMember[atom] && mol.getConnAtoms(atom) > 2) {
					for (int i=0; i<mol.getConnAtoms(atom); i++) {
						int connAtom = mol.getConnAtom(atom, i);
						if (!isMember[connAtom]) {
							addShortestPathToMember(mol, atom, connAtom, isMember);
							}
						}
					}
				}
			} while (found);

		// recursively add all exocyclic and exochain double bonds
		do {
			found = false;
			for (int atom=0; atom<mol.getAtoms(); atom++) {
				if (isMember[atom] && mol.getAtomPi(atom) != 0) {
					for (int i=0; i<mol.getConnAtoms(atom); i++) {
						int connAtom = mol.getConnAtom(atom, i);
						if (!isMember[connAtom]
						 && mol.getConnBondOrder(atom, i) > 1) {
							isMember[connAtom] = true;
							found |= (connAtom < atom);
							}
						}
					}
				}
			} while (found);

		return isMember;
		}

	private static void addShortestPathToMember(StereoMolecule mol, int atom1, int atom2, boolean[] isMember) {
		int[] parentAtom = new int[mol.getAtoms()];	// serves as indication of being a path member and allows backtracking the path
		int[] graphAtom = new int[mol.getAtoms()];
		Arrays.fill(parentAtom, -1);

		graphAtom[0] = atom2;
		parentAtom[atom1] = -2;	// != -1 to mark it being a part of the graph
		parentAtom[atom2] = atom1;
		int current = 0;
		int highest = 0;
		while (current <= highest) {
			for (int i=0; i<mol.getConnAtoms(graphAtom[current]); i++) {
				int candidate = mol.getConnAtom(graphAtom[current], i);
				if (parentAtom[candidate] == -1) {
					if (isMember[candidate]) {	// add path and return
						int atom = graphAtom[current];
						while (!isMember[atom]) {
							isMember[atom] = true;
							atom = parentAtom[atom];
							}
						return;
						}
					graphAtom[++highest] = candidate;
					parentAtom[candidate] = graphAtom[current];
					}
				}
			current++;
			}
		}

	/**
	 * Modifies mol in-place to represent the most central ring system of the original mol.
	 * @param mol may be an empty molecule if there are no ring systems
	 */
	public static void createMostCentralRingSystem(StereoMolecule mol) {
		boolean[] isScaffoldMember = findMostCentralRingSystem(mol);
		if (isScaffoldMember == null) {
			mol.deleteMolecule();
			return;
			}
		for (int atom=0; atom<mol.getAtoms(); atom++)
			if (!isScaffoldMember[atom])
				mol.markAtomForDeletion(atom);
		mol.deleteMarkedAtomsAndBonds();
		}

	/**
	 * Creates most central ring system as new StereoMolecule from mol without touching mol.
	 * @param mol
	 * @return most central ring system or null if there is none
	 */
	public static StereoMolecule getMostCentralRingSystem(StereoMolecule mol) {
		boolean[] isScaffoldMember = findMostCentralRingSystem(mol);
		if (isScaffoldMember == null)
			return null;
		StereoMolecule scaffold = new StereoMolecule(mol.getAtoms(), mol.getBonds());
		mol.copyMoleculeByAtoms(scaffold, isScaffoldMember, false, null);
		return scaffold;
		}

	/**
	 * This method flags all atoms that are a member of that ring system in mol,
	 * which is most central within this molecule graph.
	 * @param mol
	 * @return null if mol does not contain any rings
	 */
	public static boolean[] findMostCentralRingSystem(StereoMolecule mol) {
		boolean[] isCuttableBond = findCuttableBonds(mol);
		if (isCuttableBond == null)
			return null;

		int[] fragmentNo = new int[mol.getAllAtoms()];
		int fragmentCount = mol.getFragmentNumbers(fragmentNo, isCuttableBond, false);
		float[] fragmentDistanceSum = new float[fragmentCount];
		int[] fragmentAtomCount = new int[fragmentCount];
		boolean[] fragmentContainsRing = new boolean[fragmentCount];
		float[] meanAtomDistance = mol.getAverageTopologicalAtomDistance();
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			fragmentDistanceSum[fragmentNo[atom]] += meanAtomDistance[atom];
			fragmentAtomCount[fragmentNo[atom]]++;
			if (mol.isRingAtom(atom))
				fragmentContainsRing[fragmentNo[atom]] = true;
			}
		double minSum = Double.MAX_VALUE;
		int minFragment = -1;
		for (int i=0; i<fragmentCount; i++) {
			if (fragmentContainsRing[i]
			 && minSum > (double)fragmentDistanceSum[i]/fragmentAtomCount[i]) {
				minSum = (double)fragmentDistanceSum[i]/fragmentAtomCount[i];
				minFragment = i;
				}
			}
		boolean[] isFragmentMember = new boolean[mol.getAtoms()];
		for (int atom=0; atom<mol.getAtoms(); atom++)
			isFragmentMember[atom] = (fragmentNo[atom] == minFragment);
		return isFragmentMember;
		}

	private static boolean[] findCuttableBonds(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		if (mol.getRingSet().getSize() == 0)
			return null;
		boolean[] isCuttableBond = new boolean[mol.getBonds()];
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (!mol.isRingBond(bond)
			 && (mol.isRingAtom(mol.getBondAtom(0, bond)) || mol.isRingAtom(mol.getBondAtom(1, bond))))
				isCuttableBond[bond] = true;
			}
		return isCuttableBond;
		}
	}
