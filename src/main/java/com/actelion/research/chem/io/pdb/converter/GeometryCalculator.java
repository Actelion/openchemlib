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
 * @author Joel Freyss
 */

package com.actelion.research.chem.io.pdb.converter;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;


/**
 * Utility class to perform 3D geometry calculations on molecules
 */
public class GeometryCalculator {
	
	public final static Coordinates getCoordinates(StereoMolecule mol, int atm) {
		return new Coordinates(mol.getAtomX(atm), mol.getAtomY(atm), mol.getAtomZ(atm));
	}

	
	/**
	 * Gets the Angle between 3 atoms
	 * @param mol
	 * @param a1
	 * @param a2
	 * @param a3
	 * @return the angle
	 */
	public final static double getAngle(StereoMolecule mol, int a1, int a2, int a3) {
		Coordinates c1 = mol.getCoordinates(a1);
		Coordinates c2 = mol.getCoordinates(a2);
		Coordinates c3 = mol.getCoordinates(a3);

		return c1.subC(c2).getAngle(c3.subC(c2));
	}
	
	public final static double getAngle(Coordinates c1, Coordinates c2, Coordinates c3) {
		return c1.subC(c2).getAngle(c3.subC(c2));
	}
	
  	/**
	 * Gets the Dihedral Angle between 4 atoms
	 * @param mol
     * @param a1
     * @param a2
     * @param a3
	 * @param a4
	 * @return the angle
	 */
	public final static double getDihedral(StereoMolecule mol, int a1, int a2, int a3, int a4) {
		Coordinates c1 = mol.getCoordinates(a1);
		Coordinates c2 = mol.getCoordinates(a2);
		Coordinates c3 = mol.getCoordinates(a3);
		Coordinates c4 = mol.getCoordinates(a4);
		return c1.getDihedral(c2, c3, c4);
	}
		
	/**
	 * Gets the center of Gravity of a molecule
	 * @param mol
	 * @return
	 */
	public final static Coordinates getCenterGravity(StereoMolecule mol) {
		Coordinates c = new Coordinates();
		for(int i=0; i<mol.getAllAtoms(); i++) {
			c.x += mol.getAtomX(i); 
			c.y += mol.getAtomY(i); 
			c.z += mol.getAtomZ(i); 
		}
		c.x /= mol.getAllAtoms();
		c.y /= mol.getAllAtoms();
		c.z /= mol.getAllAtoms();
		
		return c;
	}
	
	/**
	 * Gets the Bounds of a molecule
	 * @param molecule
	 * @return an Array of Coordinares [lowerBounds, upperbounds]
	 */
	public final static Coordinates[] getBounds(StereoMolecule molecule) {
		if(molecule.getAllAtoms()==0) return new Coordinates[]{new Coordinates(0, 0, 0), new Coordinates(0, 0, 0)};
		Coordinates[] coords = new Coordinates[]{new Coordinates(Float.MAX_VALUE, Float.MAX_VALUE, Float.MAX_VALUE), new Coordinates(-Float.MAX_VALUE, -Float.MAX_VALUE, -Float.MAX_VALUE)};
		for(int i=0; i<molecule.getAllAtoms(); i++) {
			coords[0].x = Math.min(coords[0].x, molecule.getAtomX(i)); 
			coords[0].y = Math.min(coords[0].y, molecule.getAtomY(i)); 
			coords[0].z = Math.min(coords[0].z, molecule.getAtomZ(i)); 

			coords[1].x = Math.max(coords[1].x, molecule.getAtomX(i)); 
			coords[1].y = Math.max(coords[1].y, molecule.getAtomY(i)); 
			coords[1].z = Math.max(coords[1].z, molecule.getAtomZ(i)); 
		}
		return coords;
	}		

	/**
	 * Translate a Molecule
	 * @param molecule
	 * @param c
	 */
	public final static void translate(StereoMolecule molecule, Coordinates c) {
		for(int i=0; i<molecule.getAllAtoms(); i++) {
			molecule.setAtomX(i, molecule.getAtomX(i)+c.x);
			molecule.setAtomY(i, molecule.getAtomY(i)+c.y);
			molecule.setAtomZ(i, molecule.getAtomZ(i)+c.z);
		}
	}
	
}
