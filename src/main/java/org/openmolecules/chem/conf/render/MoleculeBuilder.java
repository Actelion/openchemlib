package org.openmolecules.chem.conf.render;

import com.actelion.research.chem.Coordinates;

/**
 * Method set to be implemented by any class that uses MoleculeArchitect
 * to construct a molecule for a 3D environment for rendering.
 */
public interface MoleculeBuilder {
	public void init();
	public void addSphere(int atom, int bond, Coordinates c, double radius, int argb);
	public void addCylinder(int bond, double radius, double length, Coordinates c, double rotationY, double rotationZ, int argb);
	public void done();
	}
