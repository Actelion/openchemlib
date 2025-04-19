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
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

public class AtomAssembler {
	private final StereoMolecule mMol;

	public AtomAssembler(StereoMolecule mol) {
		mMol = mol;
		}

	public int addImplicitHydrogens() {
		mMol.ensureHelperArrays(Molecule.cHelperRings);
		int total = 0;

		for (int atom=0; atom<mMol.getAtoms(); atom++)
			total += addImplicitHydrogens(atom, true);

		for (int atom=0; atom<mMol.getAtoms(); atom++)
			total += addImplicitHydrogens(atom, false);

		return total;
		}


	public int addImplicitHydrogens(int atom) {
		return addImplicitHydrogens(atom, false);
	}

	/**
	 * @param atom
	 * @param skipRotatableHydrogens if true and if atom is rotatable such that new hydrogen positions depend on atom rear bond rotation, then skip adding these atoms for now
	 * @return
	 */
	private int addImplicitHydrogens(int atom, boolean skipRotatableHydrogens) {
		int atomicNo = mMol.getAtomicNo(atom);
		if (atomicNo == 0)
			return 0;

		int count = mMol.getImplicitHydrogens(atom);
		if (count == 0)
			return 0;

		if (mMol.getConnAtoms(atom) == 0)
			return addHydrogensToSingleAtom(atom, count);

		int pi = mMol.getAtomPi(atom);
		int sp = (atomicNo >= 15) ? 3 : (pi == 2) ? 1 : (pi == 1 || atomicNo == 5 || mMol.isFlatNitrogen(atom)) ? 2 : 3;

		Coordinates croot = mMol.getCoordinates(atom);

		double bondLength = getHydrogenBondLength(atom);

		if (sp == 1) {	// simple case, where we need to extend linearly
			Coordinates cconn = mMol.getCoordinates(mMol.getConnAtom(atom, 0));
			Coordinates cnew = croot.addC(croot.subC(cconn).unit().scale(bondLength));
			int newAtom = mMol.addAtom(1);
			mMol.setAtomX(newAtom, cnew.x);
			mMol.setAtomY(newAtom, cnew.y);
			mMol.setAtomZ(newAtom, cnew.z);
			mMol.addBond(atom, newAtom, Molecule.cBondTypeSingle);
			return 1;
			}

		// find an 3-atom sequence that will be the basis for dihedral angle positioning
		int[] atomSequence = new int[4];
		atomSequence[2] = atom;
		for (int i1=0; i1<mMol.getConnAtoms(atom); i1++) {
			atomSequence[1] = mMol.getConnAtom(atom, i1);
			for (int i0=0; i0<mMol.getConnAtoms(atomSequence[1]); i0++) {
				atomSequence[0] = mMol.getConnAtom(atomSequence[1], i0);
				if (atomSequence[0] != atom) {
					// sequence found. Now we check whether certain dihedrals are already blocked

					if (mMol.getConnAtoms(atom) == 3) {
						// must be sp3 with already two of three positions occupied; we calculate the dihedral of the missing H
						atomSequence[3] = -1;
						double dihedral = TorsionDB.calculateTorsionExtended(mMol, atomSequence);
						addHydrogen(mMol.getCoordinates(atomSequence[0]), mMol.getCoordinates(atomSequence[1]),
								mMol.getCoordinates(atom), atom, Math.PI*109/180, dihedral, bondLength);
						return 1;
						}

					double angle = (sp == 2) ? Math.PI*2/3 : Math.PI*109/180;

					if (mMol.getConnAtoms(atom) == 2) {
						// we have one position already occupied

						if (count == 1 && skipRotatableHydrogens)
							return 0;

						for (int i3=0; i3<mMol.getConnAtoms(atom); i3++) {
							atomSequence[3] = mMol.getConnAtom(atom, i3);
							if (atomSequence[3] != atomSequence[1]) {
								double dihedral = TorsionDB.calculateTorsionExtended(mMol, atomSequence);
								dihedral += (sp == 2) ? Math.PI : Math.PI*2/3;
								if (dihedral > Math.PI)
									dihedral -= 2*Math.PI;
								// TODO in case of count==1 we might avoid a potentially crowded side
								addHydrogen(mMol.getCoordinates(atomSequence[0]), mMol.getCoordinates(atomSequence[1]),
										mMol.getCoordinates(atom), atom, angle, dihedral, bondLength);

								if (count != 1) {
									dihedral += Math.PI*2/3;
									if (dihedral > Math.PI)
										dihedral -= 2*Math.PI;
									addHydrogen(mMol.getCoordinates(atomSequence[0]), mMol.getCoordinates(atomSequence[1]),
											mMol.getCoordinates(atom), atom, angle, dihedral, bondLength);
									}
								return count;
								}
							}
						}

					if (skipRotatableHydrogens)
						return 0;

					// no competing atoms
					// if we have a single bonded option for atomSequence[0], then take that
					for (int i=i0+1; i<mMol.getConnAtoms(atomSequence[1]); i++) {
						int alternative = mMol.getConnAtom(atomSequence[1], i);
						if (alternative != atom && mMol.getConnBondOrder(atomSequence[1], i) == 1) {
							atomSequence[0] = alternative;
							break;
							}
						}

					// Here we have one of more hydrogen atoms at a terminal atom, i.e. an atom with one non-H neighbour!

					// If we have a double bond or flat nitrogen, we cannot really rotate. Preferred is then anti (PI) to single bonded atomSequence[0],
					// e.g. we prefer -OH syn to =O in carboxylic acid, which is anti to alpha carbon
					double dihedral = Math.PI;

					// If we have a rotatable bond, we look for the closest colliding atom and choose a rotation state to avoid that
					if (sp == 3) {
						dihedral = getMostBlockedDihedral(atomSequence, bondLength);
						if (count == 3)
							dihedral += Math.PI/3;
						else if (count == 2)
							dihedral += Math.PI*2/3;
						else if (count == 1)
							dihedral += Math.PI;
					}
					for (int i=0; i<count; i++) {
						addHydrogen(mMol.getCoordinates(atomSequence[0]), mMol.getCoordinates(atomSequence[1]),
								mMol.getCoordinates(atom), atom, angle, dihedral, bondLength);

						dihedral += (sp == 2) ? Math.PI : Math.PI*2/3;
						}
					return count;
					}
				}
			}

		// from here are atoms with only one shell of neighbours where we cannot apply the dihedral procedure

		if (count == 1 && sp == mMol.getConnAtoms(atom)) {
			// simple case, where we can just invert the sum of all neighbor bond vectors
			Coordinates v = new Coordinates();
			for (int i=0; i<mMol.getConnAtoms(atom); i++)
				v.add(croot.subC(mMol.getCoordinates(mMol.getConnAtom(atom, i))).unit());
			v.unit();
			int hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, croot.x + bondLength * v.x);
			mMol.setAtomY(hydrogen, croot.y + bondLength * v.y);
			mMol.setAtomZ(hydrogen, croot.z + bondLength * v.z);
			return 1;
			}

		double angle = (sp == 2) ? Math.PI*2/3 : Math.PI*109.5/180;
		double rotationDif = (sp == 2) ? Math.PI : Math.PI*2/3;

		int rearAtom = mMol.getConnAtom(atom, 0);
		Coordinates v = croot.subC(mMol.getCoordinates(rearAtom));
		boolean vIsParallelToZ = (v.z != 0 && (Math.abs(v.x)+Math.abs(v.y))/Math.abs(v.z) < 0.01);

		// to ensure staggered hydrogens on opposite side of bond we use an upsideDown value to distinguish bond directions
		boolean upsideDown = false;
		if (vIsParallelToZ) {
			if (v.z < 0)
				upsideDown = true;
			}
		else {
			if (v.x == 0) {
				if (v.y<0)
					upsideDown = true;
				}
			else {
				if (v.x<0)
					upsideDown = true;
				}
			}

		double rotation = 0.0;

		if (mMol.getConnAtoms(atom) == 2) {
			// we have another neighbour that determines the rotation of the hydrogen atom(s); atom must be sp3

			// if v is reasonably parallel to the z-axis, we swap v and v2
			int otherAtom = mMol.getConnAtom(atom, 1);
			if (vIsParallelToZ) {
				rearAtom = otherAtom;
				otherAtom = mMol.getConnAtom(atom, 0);
				v = croot.subC(mMol.getCoordinates(rearAtom));
				vIsParallelToZ = false;
				}

			Coordinates dummyCoords = new Coordinates(mMol.getCoordinates(rearAtom));
			dummyCoords.add(0, 0, 1);

			rotation = calculateDihedral(dummyCoords, mMol.getCoordinates(rearAtom), mMol.getCoordinates(atom), mMol.getCoordinates(otherAtom));
			rotation += rotationDif;
			}

		if (vIsParallelToZ) {
			double dz = -Math.cos(angle) * bondLength;
			if (upsideDown) {
				rotation += Math.PI;
				dz = -dz;
				}
			double r = bondLength * Math.sin(Math.PI - angle);
			for (int i=0; i<count; i++) {
				int hydrogen = mMol.addAtom(1);
				mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
				mMol.setAtomX(hydrogen, croot.x + r * Math.sin(rotation));
				mMol.setAtomY(hydrogen, croot.y + r * Math.cos(rotation));
				mMol.setAtomZ(hydrogen, croot.z + dz);

				rotation += rotationDif;
				}
			return count;
			}
		else {
			for (int i=0; i<count; i++) {
				// now we need to transform coordinates
				Coordinates c1 = new Coordinates(mMol.getCoordinates(rearAtom));
				c1.add(0, 0, upsideDown ? -1 : 1);
				Coordinates c2 = mMol.getCoordinates(rearAtom);
				Coordinates c3 = mMol.getCoordinates(atom);

				addHydrogen(c1, c2, c3, atom, angle, rotation, bondLength);

				rotation += rotationDif;
				}

			return count;
			}
		}

	private int addHydrogensToSingleAtom(int atom, int count) {
		Coordinates p = mMol.getCoordinates(atom);
		double length = getHydrogenBondLength(atom);
		switch (count) {
		case 1:
			int hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, p.x+length);
			mMol.setAtomY(hydrogen, p.y);
			mMol.setAtomZ(hydrogen, p.z);
			return 1;
		case 2:
			double angle = Math.PI/180*(mMol.getAtomicNo(atom) == 8 ? 104.45 : mMol.getAtomicNo(atom) == 16 ? 92.1 : 109.5);
			hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, p.x+length);
			mMol.setAtomY(hydrogen, p.y);
			mMol.setAtomZ(hydrogen, p.z);
			hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, p.x+length*Math.cos(angle));
			mMol.setAtomY(hydrogen, p.y+length*Math.sin(angle));
			mMol.setAtomZ(hydrogen, p.z);
			return 2;
		case 3:
			angle = Math.PI/180*21;
			double dx = length*Math.cos(angle);
			double dz = length*Math.sin(angle);
			hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, p.x+dx);
			mMol.setAtomY(hydrogen, p.y);
			mMol.setAtomZ(hydrogen, p.z-dz);
			hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, p.x-0.5*dx);
			mMol.setAtomY(hydrogen, p.y+0.866*dx);
			mMol.setAtomZ(hydrogen, p.z-dz);
			hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, p.x-0.5*dx);
			mMol.setAtomY(hydrogen, p.y-0.866*dx);
			mMol.setAtomZ(hydrogen, p.z-dz);
			return 3;
		case 4:
			angle = Math.PI/180*19.5;
			dx = length*Math.cos(angle);
			dz = length*Math.sin(angle);
			hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, p.x-dx);
			mMol.setAtomY(hydrogen, p.y);
			mMol.setAtomZ(hydrogen, p.z+dz);
			hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, p.x+0.5*dx);
			mMol.setAtomY(hydrogen, p.y+0.866*dx);
			mMol.setAtomZ(hydrogen, p.z+dz);
			hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, p.x+0.5*dx);
			mMol.setAtomY(hydrogen, p.y-0.866*dx);
			mMol.setAtomZ(hydrogen, p.z+dz);
			hydrogen = mMol.addAtom(1);
			mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
			mMol.setAtomX(hydrogen, p.x);
			mMol.setAtomY(hydrogen, p.y);
			mMol.setAtomZ(hydrogen, p.z-length);
			return 4;
			}

		// water 104.45 degrees, 95.84 pm
		// ammonia: 107.8 degrees,  101.7 pm
		// methane: 109.5 degrees, 108.7 pm
		// borane: 120 degrees, 119 pm
		// HF 91.7 pm
		// HCl 127.4 pm
		// HBr 141.4 pm
		// HI 160.9 pm
		// H2S 92.1 degrees, 133.6 pm

		return 0;
		}

	private double getHydrogenBondLength(int atom) {
		int index = BondLengthSet.getBondIndex(1, false, false, mMol.getAtomicNo(atom), 1, mMol.getAtomPi(atom), 0);
		return (index == -1) ? 1.09 : BondLengthSet.getBondLength(index);
	}

	private double calculateDihedral(Coordinates c1, Coordinates c2, Coordinates c3, Coordinates c4) {
		Coordinates v1 = c2.subC(c1);
		Coordinates v2 = c3.subC(c2);
		Coordinates v3 = c4.subC(c3);

		Coordinates n1 = v1.cross(v2);
		Coordinates n2 = v2.cross(v3);

		return -Math.atan2(v2.getLength() * v1.dot(n2), n1.dot(n2));
		}

	/**
	 * Adds a new single bonded atom to rootAtom such that bond length c3->newAtom, angle c2->c3->newAtom and dihedral c1->c2->c3->newAtom are met.
	 * @param atomSequence
	 * @param bondLength
	 */
	private double getMostBlockedDihedral(int[] atomSequence, double bondLength) {
		Coordinates croot = mMol.getCoordinates(atomSequence[2]);
		Coordinates cconn = mMol.getCoordinates(atomSequence[2]);
		Coordinates ctest = croot.addC(croot.subC(cconn).unit().scale(0.5*bondLength));
		final double angle = 109*Math.PI/180;
		double orbitRadius = bondLength * Math.sin(angle);

		final double FACTOR = 0.85;
		double maxCollision = 0.0;
		double maxTorsion = Double.NaN;
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			if (atom != atomSequence[1] && atom != atomSequence[2]) {
				Coordinates catom = mMol.getCoordinates(atom);
				// Simplified approach: we keep catom probe in center of the orbit of its potential positions
				// and add the orbit radius to the accepted VDW-radii distance.
				double minDist = orbitRadius + FACTOR * VDWRadii.getVDWRadius(1)+VDWRadii.getVDWRadius(mMol.getAtomicNo(atom));
				double dx = Math.abs(ctest.x - catom.x);
				if (dx < minDist) {
					double dy = Math.abs(ctest.y - catom.y);
					if (dy < minDist) {
						double dz = Math.abs(ctest.z - catom.z);
						if (dz < minDist) {
							if (Math.sqrt(dx*dx + dy*dy + dz*dz) < minDist) {	// we have a potential collision and need to calculate precisely
								atomSequence[3] = atom;
								double dihedral = (atom == atomSequence[0]) ? 0.0 : TorsionDB.calculateTorsionExtended(mMol, atomSequence);
								Coordinates p = getCoordinatesWithConstraints(mMol.getCoordinates(atomSequence[0]), mMol.getCoordinates(atomSequence[1]),
										mMol.getCoordinates(atomSequence[2]), angle, dihedral, bondLength);

								dx = Math.abs(p.x - catom.x);
								dy = Math.abs(p.y - catom.y);
								dz = Math.abs(p.z - catom.z);
								double distance = Math.sqrt(dx*dx + dy*dy + dz*dz);
								double collision = minDist - orbitRadius - distance;
								if (maxCollision < collision) {
									maxCollision = collision;
									maxTorsion = dihedral;
								}
							}
						}
					}
				}
			}
		}

		return maxTorsion;
	}

	/**
	 * Adds a new single bonded atom to rootAtom such that bond length c3->newAtom, angle c2->c3->newAtom and dihedral c1->c2->c3->newAtom are met.
	 * @param c1
	 * @param c2
	 * @param c3
	 * @param rootAtom
	 * @param angle
	 * @param dihedral
	 * @param bondLength
	 */
	private void addHydrogen(Coordinates c1, Coordinates c2, Coordinates c3, int rootAtom, double angle, double dihedral, double bondLength) {
		Coordinates p = getCoordinatesWithConstraints(c1, c2, c3, angle, dihedral, bondLength);
		int hydrogen = mMol.addAtom(1);
		mMol.addBond(rootAtom, hydrogen, Molecule.cBondTypeSingle);
		mMol.setAtomX(hydrogen, p.x);
		mMol.setAtomY(hydrogen, p.y);
		mMol.setAtomZ(hydrogen, p.z);
	}

	/**
	 * Adds a new single bonded atom to rootAtom such that bond length c3->newAtom, angle c2->c3->newAtom and dihedral c1->c2->c3->newAtom are met.
	 * @param c1
	 * @param c2
	 * @param c3
	 * @param angle
	 * @param dihedral
	 * @param bondLength
	 */
	private Coordinates getCoordinatesWithConstraints(Coordinates c1, Coordinates c2, Coordinates c3, double angle, double dihedral, double bondLength) {
		double r = bondLength * Math.sin(Math.PI - angle);
		double x = -r * Math.sin(dihedral);
		double y = r * Math.cos(dihedral);
		double z = bondLength * Math.cos(Math.PI - angle);

		Coordinates axisZ = c3.subC(c2).unit();
		Coordinates axisX = c1.subC(c2).cross(axisZ).unit();    // needs to be mapped to x-axis
		Coordinates axisY = axisX.cross(axisZ).unit();			// needs to be mapped to y-axis

		double[][] m = new double[3][3];
		m[0][0] = -axisX.x;
		m[0][1] = -axisX.y;
		m[0][2] = -axisX.z;
		m[1][0] = -axisY.x;
		m[1][1] = -axisY.y;
		m[1][2] = -axisY.z;
		m[2][0] = axisZ.x;
		m[2][1] = axisZ.y;
		m[2][2] = axisZ.z;
		return new Coordinates(x, y, z).rotate(m).add(c3);
		}
	}
