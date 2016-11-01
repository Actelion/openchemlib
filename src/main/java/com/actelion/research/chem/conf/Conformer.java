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

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;

public class Conformer {
	private Coordinates[] mCoordinates;
	private StereoMolecule mMol;
	private short[] mBondTorsion;

	/**
	 * Creates a Conformer, i.e. basically a new set of atom coordinates for the given Molecule.
	 * The conformer is initialized with an exact copy of the molecules internal coordinates.
	 * @param mol
	 */
	public Conformer(StereoMolecule mol) {
		mMol = mol;
		mCoordinates = new Coordinates[mol.getAllAtoms()];
		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			mCoordinates[atom] = new Coordinates(mol.getCoordinates(atom));
		}

	/**
	 * Creates a new conformer as an exact copy of the given one.
	 * @param c
	 */
	public Conformer(Conformer c) {
		mMol = c.mMol;

		mCoordinates = new Coordinates[c.mCoordinates.length];
		for (int i = 0; i< mCoordinates.length; i++)
			mCoordinates[i] = new Coordinates(c.mCoordinates[i]);

		if (c.mBondTorsion != null) {
			mBondTorsion = new short[c.mBondTorsion.length];
			for (int i=0; i<c.mBondTorsion.length; i++)
				mBondTorsion[i] = c.mBondTorsion[i];
			}
		}

	public int getSize() {
		return mCoordinates.length;
		}

	public double getX(int atom) {
		return mCoordinates[atom].x;
		}
	public double getY(int atom) {
		return mCoordinates[atom].y;
		}
	public double getZ(int atom) {
		return mCoordinates[atom].z;
		}

	public Coordinates getCoordinates(int atom) {
		return mCoordinates[atom];
		}

	/**
	 * Copies x,y,z from coords into this atom's coordinates
	 * @param atom
	 * @param coords
	 */
	public void setCoordinates(int atom, Coordinates coords) {
		mCoordinates[atom].set(coords);
		}

	/**
	 * Replaces the atom's Coordinates object by the given one
	 * @param atom
	 * @param coords
	 */
	public void setCoordinatesReplace(int atom, Coordinates coords) {
		mCoordinates[atom] = coords;
	}

	public void setX(int atom, double x) {
		mCoordinates[atom].x = x;
		}
	public void setY(int atom, double y) {
		mCoordinates[atom].y = y;
		}
	public void setZ(int atom, double z) {
		mCoordinates[atom].z = z;
		}

	/**
	 * Returns the current bond torsion angle in degrees, it is was set before.
	 * @param bond
	 * @return -1 or previously set torsion angle in the range 0 ... 359
	 */
	public int getBondTorsion(int bond) {
		return mBondTorsion == null ? -1 : mBondTorsion[bond];
		}

	/**
	 * Sets the current bond torsion to be retrieved later.
	 * @param bond
	 * @param torsion in degrees
	 */
	public void setBondTorsion(int bond, short torsion) {
		if (mBondTorsion == null)
			mBondTorsion = new short[mMol.getAllBonds()];
		mBondTorsion[bond] = torsion;
		}

	/**
	 * @return reference to the original molecule
	 */
	public StereoMolecule getMolecule() {
		return mMol;
		}

	/**
	 * Copies the molecule's atom coordinates to this Conformer.
	 * @param mol molecule identical to the original molecule passed in the Constructor
	 */
	public void copyFrom(StereoMolecule mol) {
		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			mCoordinates[atom].set(mol.getAtomX(atom), mol.getAtomY(atom), mol.getAtomZ(atom));
		}

	/**
	 * Copies the conformer's atom coordinates to this Conformer.
	 * @param conformer other conformer of the molecule passed in the Constructor
	 */
	public void copyFrom(Conformer conformer) {
		for (int atom=0; atom<conformer.getSize(); atom++)
			mCoordinates[atom].set(conformer.getX(atom), conformer.getY(atom), conformer.getZ(atom));
		}

	/**
	 * Copies this Conformer's atom coordinates to the given molecule.
	 * If no molecule is given, then a compact copy of this conformer's
	 * molecule is created and returned with this conformer's coordinates.
	 * @param mol null or original molecule passed in the Constructor or identical copy
	 * @return this conformer as StereoMolecule
	 */
	public StereoMolecule toMolecule(StereoMolecule mol) {
		if (mol == null)
			mol = mMol.getCompactCopy();
		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			mol.getCoordinates(atom).set(mCoordinates[atom]);
		return mol;
		}
	}
