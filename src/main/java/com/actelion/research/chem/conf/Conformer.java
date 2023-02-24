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

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;

import java.util.Arrays;

public class Conformer implements Comparable<Conformer> {
	private Coordinates[] mCoordinates;
	private StereoMolecule mMol;
	private String mName;
	private short[] mBondTorsion;
	private double mEnergy;
	private double mLikelihood;
	private TorsionDescriptor mTorsionDescriptor;

	/**
	 * Creates a Conformer, i.e. basically a new set of atom coordinates for the given Molecule.
	 * The conformer is initialized with an exact copy of the molecules internal coordinates.
	 *
	 * @param mol
	 */
	public Conformer(StereoMolecule mol) {
		mMol = mol;
		mCoordinates = new Coordinates[mol.getAllAtoms()];
		for (int atom = 0; atom < mol.getAllAtoms(); atom++)
			mCoordinates[atom] = new Coordinates(mol.getCoordinates(atom));

		mEnergy = Double.NaN;
	}

	/**
	 * Creates a new conformer as an exact copy of the given one.
	 *
	 * @param c
	 */
	public Conformer(Conformer c) {
		this(c, c.mMol);
	}

	/**
	 * Use at own risk!!
	 *
	 * Initializes the conformer as a copy of c but replaces the
	 * StereoMolecule by the supplied molecule.
	 *
	 * @param c
	 * @param mol
	 */
	public Conformer(Conformer c, StereoMolecule mol) {
		mMol = mol;

		mCoordinates = new Coordinates[c.mCoordinates.length];
		for (int i = 0; i < mCoordinates.length; i++)
			mCoordinates[i] = new Coordinates(c.mCoordinates[i]);

		if (c.mBondTorsion != null)
			mBondTorsion = Arrays.copyOf(c.mBondTorsion, c.mBondTorsion.length);

		mName = (c.mName == null || c.mName.endsWith(" (copy)")) ? c.mName : c.mName.concat(" (copy)");

		mEnergy = Double.NaN;
	}


	/**
	 * Translate this conformer's coordinates such that its center of gravity
	 * is moved to P(0,0,0) assuming all atoms have the same mass.
	 *
	 * @return this conformer with centered coordinates
	 */
	public Conformer center() {
		Coordinates cog = new Coordinates();
		for (int atom = 0; atom < mMol.getAllAtoms(); atom++)
			cog.add(mCoordinates[atom]);
		cog.scale(1.0 / mMol.getAllAtoms());

		for (int atom = 0; atom < mMol.getAllAtoms(); atom++)
			mCoordinates[atom].sub(cog);

		return this;
	}

	/**
	 * Translate this conformer's coordinates by adding the dx,dy,dz shifts
	 * to all atom coordinates.
	 *
	 * @return this conformer with translated coordinates
	 */
	public Conformer translate(double dx, double dy, double dz) {
		for (int atom = 0; atom < mMol.getAllAtoms(); atom++)
			mCoordinates[atom].add(dx, dy, dz);

		return this;
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

	public Coordinates[] getCoordinates() {
		return mCoordinates;
	}

	public Coordinates getCoordinates(int atom) {
		return mCoordinates[atom];
	}

	/**
	 * Copies x,y,z from coords into this atom's coordinates
	 *
	 * @param atom
	 * @param coords
	 */
	public void setCoordinates(int atom, Coordinates coords) {
		mCoordinates[atom].set(coords);
	}

	/**
	 * Replaces the atom's Coordinates object by the given one
	 *
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
	 * Removes atoms Coordinates objects from those atoms that are marked in the given array.
	 * Make sure to also remove those atoms from the underlying Molecule.
	 * When this method is used, bond torsion angles are lost.
	 * @param isToBeDeleted
	 * @return
	 */
	public int deleteAtoms(boolean[] isToBeDeleted) {
		int count = 0;
		for (int i = 0; i < mCoordinates.length; i++)
			if (isToBeDeleted[i])
				count++;

		if (count != 0) {
			Coordinates[] newCoords = new Coordinates[mCoordinates.length - count];
			int newIndex = 0;
			for (int i = 0; i < mCoordinates.length; i++)
				if (!isToBeDeleted[i])
					newCoords[newIndex++] = mCoordinates[i];

			mCoordinates = newCoords;
			mBondTorsion = null;
		}

		return count;
	}

	/**
	 * Calculates a signed torsion as an exterior spherical angle
	 * from a valid 4-atom strand.
	 * Looking along the central bond, the torsion angle is 0.0, if the
	 * projection of front and rear bonds point in the same direction.
	 * If the front bond is rotated in the clockwise direction, the angle
	 * increases, i.e. has a positive value.
	 * http://en.wikipedia.org/wiki/Dihedral_angle
	 *
	 * @param atom 4 valid atom indices defining a connected atom sequence
	 * @return torsion in the range: -pi <= torsion <= pi
	 */
	public double calculateTorsion(int[] atom) {
		Coordinates c1 = getCoordinates(atom[0]);
		Coordinates c2 = getCoordinates(atom[1]);
		Coordinates c3 = getCoordinates(atom[2]);
		Coordinates c4 = getCoordinates(atom[3]);

		Coordinates v1 = c2.subC(c1);
		Coordinates v2 = c3.subC(c2);
		Coordinates v3 = c4.subC(c3);

		Coordinates n1 = v1.cross(v2);
		Coordinates n2 = v2.cross(v3);

		return -Math.atan2(v2.getLength() * v1.dot(n2), n1.dot(n2));
	}

	/**
	 * Returns the current bond torsion angle in degrees, as it was set before.
	 * @param bond
	 * @return -1 or previously set torsion angle in the range 0 ... 359
	 */
	public short getBondTorsion(int bond) {
		return mBondTorsion == null ? -1 : mBondTorsion[bond];
	}

	/**
	 * Sets the current bond torsion to be retrieved later.
	 * @param bond
	 * @param torsion in degrees
	 */
	public void setBondTorsion(int bond, short torsion) {
		if (mBondTorsion == null) {
			mBondTorsion = new short[mMol.getAllBonds()];
			Arrays.fill(mBondTorsion, (short)-1);
			}

		while (torsion < 0)
			torsion += 360;
		while (torsion >= 360)
			torsion -= 360;

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
	 *
	 * @param mol molecule identical to the original molecule passed in the Constructor
	 */
	public void copyFrom(StereoMolecule mol) {
		for (int atom = 0; atom < mol.getAllAtoms(); atom++)
			mCoordinates[atom].set(mol.getCoordinates(atom));
	}

	/**
	 * Copies this conformer's atom coordinates to mol.
	 *
	 * @param mol molecule identical to the original molecule passed in the Constructor
	 */
	public void copyTo(StereoMolecule mol) {
		for (int atom = 0; atom < mol.getAllAtoms(); atom++)
			mol.getCoordinates(atom).set(mCoordinates[atom]);
	}

	/**
	 * Copies the conformer's atom coordinates and torsion values to this Conformer.
	 *
	 * @param conformer other conformer of the molecule passed in the Constructor
	 */
	public void copyFrom(Conformer conformer) {
		for (int atom = 0; atom < conformer.getSize(); atom++)
			mCoordinates[atom].set(conformer.mCoordinates[atom]);

		if (conformer.mBondTorsion == null)
			mBondTorsion = null;
		else
			mBondTorsion = Arrays.copyOf(conformer.mBondTorsion, conformer.mBondTorsion.length);
	}

	/**
	 * Copies this Conformer's atom coordinates to the associated molecule.
	 *
	 * @return this conformer's associated StereoMolecule with updated coordinates
	 */
	public StereoMolecule toMolecule() {
		return toMolecule(mMol);
	}

	/**
	 * Copies this Conformer's atom coordinates to the given molecule.
	 * If no molecule is given, then a compact copy of this conformer's
	 * molecule is created and returned with this conformer's coordinates.
	 *
	 * @param mol null or original molecule passed in the Constructor or identical copy
	 * @return this conformer as StereoMolecule
	 */
	public StereoMolecule toMolecule(StereoMolecule mol) {
		if (mol == null)
			mol = mMol.getCompactCopy();
		for (int atom = 0; atom < mol.getAllAtoms(); atom++)
			mol.getCoordinates(atom).set(mCoordinates[atom]);
		if (mName != null)
			mol.setName(mName);
		return mol;
	}

	public double getEnergy() {
		return mEnergy;
	}

	public void setEnergy(double energy) {
		mEnergy = energy;
	}

	public double getLikelihood() {
		return mLikelihood;
	}

	public void setLikelihood(double likelihood) {
		mLikelihood = likelihood;
	}


	public String getName() {
		return mName == null ? mMol.getName() : mName;
	}

	public void setName(String name) {
		mName = name;
	}

	/**
	 * Returns true, if none of the torsion angles between both conformers
	 * are more different than TorsionDescriptor.TORSION_EQUIVALENCE_TOLERANCE;
	 * Calling this method requires that calculateDescriptor() has been called earlier.
	 * @param conformer
	 * @return true if all torsions are similar
	 */
	public boolean equals(Conformer conformer) {
		ensureDescriptors(conformer);
		return mTorsionDescriptor.equals(conformer.mTorsionDescriptor);
	}

	@Override
	public int compareTo(Conformer conformer) {
		ensureDescriptors(conformer);
		return mTorsionDescriptor.compareTo(conformer.mTorsionDescriptor);
	}

	/**
	 * Calculates the torsion descriptor for the current coordinates using the passed set of rotatable bonds.
	 * You may use TorsionDescriptorHelper.findRotatableBonds() once and pass it for every new conformer to this method.
	 * @param rotatableBond set of rotatable bonds to be considered
	 */
	public void calculateDescriptor(int[] rotatableBond) {
		mTorsionDescriptor = new TorsionDescriptorHelper(getMolecule(), rotatableBond).getTorsionDescriptor(this);
	}

	private void ensureDescriptors(Conformer c) {
		if (mTorsionDescriptor == null || c.mTorsionDescriptor == null) {
			TorsionDescriptorHelper tdh = new TorsionDescriptorHelper(mMol);
			if (mTorsionDescriptor == null)
				mTorsionDescriptor = tdh.getTorsionDescriptor(this);
			if (c.mTorsionDescriptor == null)
				c.mTorsionDescriptor = tdh.getTorsionDescriptor(c);
		}
	}


	/**
	 * Serializes the given conformer to a string.
	 * <p>
	 * Use the constructor Conformer(String serialized) to deserialize
	 * the conformer.
	 *
	 * @return
	 *
	public String serialize() {
		StringBuilder sb = new StringBuilder();

		Canonizer c = new Canonizer(mMol);
		String s_idcode = c.getIDCode();

		List<String> s_coordinates = new ArrayList<>();
		for (Coordinates coords:mCoordinates)
			s_coordinates.add(encodeCoordinates(coords));

		String s_bondTorsion = null;
		if (mBondTorsion == null) {
			s_bondTorsion = "N";
		} else {
			s_bondTorsion = "X" + SimpleEncoder.encodeShortArray(mBondTorsion);
		}

		sb.append(s_idcode);
		sb.append(" ");
		sb.append(String.join(";", s_coordinates));
		sb.append(" ");
		sb.append(String.join(";", s_bondTorsion));
		sb.append(" ");
		sb.append("" + mEnergy);
		sb.append(" ");
		sb.append("" + mName);

		return sb.toString();
	}*/

	/**
	 * Constructor for deserialization.
	 *
	 * @param serialized
	 *
	public Conformer(String serialized) {
		mTorsionDescriptor = null;

		String[] splits = serialized.split(" ");

		String idcode   = splits[0];
		String coords   = splits[1];
		String torsions = splits[2];
		String energy   = splits[3];
		String name     = splits[4];

		IDCodeParser p = new IDCodeParser();
		mMol = new StereoMolecule();
		p.parse(mMol,idcode);
		String[] c_split = coords.split(";");
		mCoordinates = new Coordinates[c_split.length];
		for(int zi=0;zi<c_split.length;zi++) { mCoordinates[zi] = decodeCoordinates(c_split[zi]); }

		if(torsions.startsWith("N")) {
			mBondTorsion = null;
		}
		else {
			mBondTorsion = SimpleEncoder.decodeShortArray(torsions.substring(1));
		}

		mEnergy = Double.parseDouble(energy);
		mName = name;
	}

	public static String encodeCoordinates(Coordinates c) {
		StringBuilder sb = new StringBuilder();

		byte[] data = ByteBuffer.allocate(4*3).putFloat((float) c.x)
				.putFloat((float) c.y)
				.putFloat((float) c.z).array();
		return SimpleEncoder.convertByteArrayToS3IDCode(data);
	}

	public static Coordinates decodeCoordinates(String s) {
		byte[] data = SimpleEncoder.convertS3IDCodeToByteArray(s);
		byte[] fx = new byte[]{data[0],data[1],data[2],data[3]};
		byte[] fy = new byte[]{data[4],data[5],data[6],data[7]};
		byte[] fz = new byte[]{data[8],data[9],data[10],data[11]};
		float cx = ByteBuffer.wrap(fx).getFloat();
		float cy = ByteBuffer.wrap(fy).getFloat();
		float cz = ByteBuffer.wrap(fz).getFloat();
		return new Coordinates(cx,cy,cz);
	}
*/
}
