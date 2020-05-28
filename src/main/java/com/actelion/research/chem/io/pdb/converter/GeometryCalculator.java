/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
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
