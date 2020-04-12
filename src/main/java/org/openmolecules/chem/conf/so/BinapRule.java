/*
 * @(#)StereoRule.java
 *
 * Copyright 2013 openmolecules.org, Inc. All Rights Reserved.
 *
 * NOTICE: All information contained herein is, and remains the property
 * of openmolecules.org.  The intellectual and technical concepts contained
 * herein are proprietary to openmolecules.org.
 * Actelion Pharmaceuticals Ltd. is granted a non-exclusive, non-transferable
 * and timely unlimited usage license.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.so;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

import java.util.ArrayList;
import java.util.Arrays;

public class BinapRule extends ConformationRule {
	private int[] mRotatableAtom;
	private boolean mInverse;

	public BinapRule(StereoMolecule mol, int[] atom, int bond) {
		super(atom);
		mRotatableAtom = getRotatableAtoms(mol, bond, atom);
		}

    public static void calculateRules(ArrayList<ConformationRule> ruleList, StereoMolecule mol) {
		for (int bond=0; bond<mol.getBonds(); bond++) {
			int parity = mol.getBondParity(bond);
			if ((parity == Molecule.cBondParityEor1 || parity == Molecule.cBondParityZor2) && mol.getBondOrder(bond) != 2) {
				int[] atomList = new int[4];
				for (int i=0; i<2; i++) {
					int atom = mol.getBondAtom(i, bond);
					int rearAtom = mol.getBondAtom(1-i, bond);
					int lowConn = Integer.MAX_VALUE;
					for (int j=0; j<mol.getConnAtoms(atom); j++) {
						int connAtom = mol.getConnAtom(atom, j);
						if (connAtom != rearAtom && lowConn > connAtom)
							lowConn = connAtom;
						}
					atomList[3*i] = lowConn;
					atomList[1+i] = atom;
					}
				if (parity == Molecule.cBondParityEor1) {
					int atom = mol.getBondAtom(0, bond);
					int rearAtom = mol.getBondAtom(1, bond);
					for (int i=0; i<mol.getConnAtoms(atom); i++) {
						int connAtom = mol.getConnAtom(atom, i);
						if (connAtom != rearAtom && connAtom != atomList[0]) {
							atomList[0] = connAtom;
							break;
							}
						}
					}

				ruleList.add(new BinapRule(mol, atomList, bond));
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

	private int[] getRotatableAtoms(StereoMolecule mol, int bond, int[] atomList) {
		// If both aromatic systems are connected by another chain of atoms, then we rotate just one aromatic system
		// with all connected substituents that don't connect to the other side.
		if (mol.isRingBond(bond))
			return getAromaticFragmentAtoms(mol, atomList[1], atomList[2]);

		// Otherwise we rotate the smaller side of the molecule seen from the binap bond perspective
		boolean[][] isMemberAtom = new boolean[2][mol.getAllAtoms()];
		int[] count = new int[2];
		for (int i=0; i<2; i++)
			count[i] = mol.getSubstituent(atomList[1+i], atomList[2-i], isMemberAtom[i], null, null);

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
		if (torsion >= 0)  // positive torsions are desired
			return false;

		Coordinates unit = conformer.getCoordinates(mAtom[2]).subC(conformer.getCoordinates(mAtom[1]));
		unit.unit();

		double rotation = torsion > -Math.PI/2 ? -2.0*torsion : 2.0*(Math.PI - torsion);
		if (mInverse)
			rotation = -rotation;

		for (int atom:mRotatableAtom)
			rotateAtom(conformer, atom, mAtom[1], unit, rotation);

		return true;
		}

	@Override
	public double addStrain(Conformer conformer, double[] atomStrain) {
		return conformer.calculateTorsion(mAtom) < 0 ? 1.0 : 0.0;
		}

	private void rotateAtom(Conformer conformer, int atom, int refAtom, Coordinates unit, double theta) {
		double x = unit.x;
		double y = unit.y;
		double z = unit.z;
		double c = Math.cos(theta);
		double s = Math.sin(theta);
		double t = 1-c;
		double mx = conformer.getX(atom) - conformer.getX(refAtom);
		double my = conformer.getY(atom) - conformer.getY(refAtom);
		double mz = conformer.getZ(atom) - conformer.getZ(refAtom);
		conformer.setX(atom, conformer.getX(refAtom) + (t*x*x+c)*mx + (t*x*y+s*z)*my + (t*x*z-s*y)*mz);
		conformer.setY(atom, conformer.getY(refAtom) + (t*x*y-s*z)*mx + (t*y*y+c)*my + (t*y*z+s*x)*mz);
		conformer.setZ(atom, conformer.getZ(refAtom) + (t*x*z+s*y)*mx + (t*z*y-s*x)*my + (t*z*z+c)*mz);
		}

	@Override
	public int getRuleType() {
		return RULE_TYPE_BINAP;
		}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("binap rule:");
		super.addAtomList(sb);
		return sb.toString();
		}
	}
