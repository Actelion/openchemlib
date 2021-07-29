package com.actelion.research.chem.contrib;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

public class HydrogenHandler {
	final static double cos30 = Math.cos(30.0 / 360.0 * Math.PI * 2);
	final static double sin30 = Math.sin(30.0 / 360.0 * Math.PI * 2);
	final static double cos330 = Math.cos(-30.0 / 360.0 * Math.PI * 2);
	final static double sin330 = Math.sin(-30.0 / 360.0 * Math.PI * 2);
	final static double cos60 = Math.cos(60.0 / 360.0 * Math.PI * 2);
	final static double sin60 = Math.sin(60.0 / 360.0 * Math.PI * 2);
	final static double cos300 = Math.cos(-60.0 / 360.0 * Math.PI * 2);
	final static double sin300 = Math.sin(-60.0 / 360.0 * Math.PI * 2);
	final static double cos45 = Math.cos(45.0 / 360.0 * Math.PI * 2);
	final static double sin45 = Math.sin(45.0 / 360.0 * Math.PI * 2);
	final static double cos315 = Math.cos(-45.0 / 360.0 * Math.PI * 2);
	final static double sin315 = Math.sin(-45.0 / 360.0 * Math.PI * 2);
	final static double H_BOND_RATIO=0.7;

	// changed from method name from expandHydrogensAtAtom(); TLS 9.Nov.2015
	public static void addImplicitHydrogens(StereoMolecule molecule, int iAtom) {
		final int atomicNo = 1;
		// Following line let crash the "cleanup"
		// molecule.setAtomMarker(iAtom, true);
		int nHydrogens = molecule.getImplicitHydrogens(iAtom);
		double x = molecule.getAtomX(iAtom);
		double y = molecule.getAtomY(iAtom);
		switch (nHydrogens) {
			case 1: {
				int nNeighbours = molecule.getConnAtoms(iAtom);
				// this is for salts like HCL
				double dx, dy;
				int iNew;
				if (nNeighbours == 0) {
					double distMin = Double.MAX_VALUE;
					int iMin = -1;
					for (int atm = 0; atm < molecule.getAllAtoms(); atm++) {
						if (atm == iAtom)
							continue;
						double deltaX = x - molecule.getAtomX(atm);
						double deltaY = y - molecule.getAtomY(atm);
						double currDist = Math.sqrt((deltaX * deltaX) + (deltaY * deltaY));
						if (distMin > currDist) {
							distMin = currDist;
							iMin = atm;
						}
					}
					dx = x - molecule.getAtomX(iMin);
					dy = y - molecule.getAtomY(iMin);
				} else {
					dx = x - molecule.getAtomX(molecule.getConnAtom(iAtom, 0));
					dy = y - molecule.getAtomY(molecule.getConnAtom(iAtom, 0));
				}
				if (nNeighbours == 1) {
					// rotate by 45
					iNew = molecule.addAtom((float)(x + cos45 * dx + sin45 * dy),(float)(y - sin45 * dx + cos45 * dy));
				} else if (nNeighbours == 2) {
					dx = x - 0.5 * (molecule.getAtomX(molecule.getConnAtom(iAtom, 0)) + molecule.getAtomX(molecule.getConnAtom(iAtom, 1)));
					dy = y - 0.5 * (molecule.getAtomY(molecule.getConnAtom(iAtom, 0)) + molecule.getAtomY(molecule.getConnAtom(iAtom, 1)));
					iNew = molecule.addAtom((float)(x + dx),(float)(y + dy));
				} else if (nNeighbours == 3) {  // case of existing 3 atoms. We need to put the hydrogen next to the chiral bond !
					// looking for a stereobond ...
					int stereoAtom1=molecule.getConnAtom(iAtom, 0);
					for (int i=1; i<3; i++) {
						int bond=molecule.getConnBond(iAtom, i);
						if ((molecule.getBondType(bond)== Molecule.cBondTypeDown) ||
								(molecule.getBondType(bond)==Molecule.cBondTypeUp)) {
							stereoAtom1=molecule.getConnAtom(iAtom, i);
						}
					}

					// should we put the hydrogen next to the chiral bond ? It depends ot the sum of the 3 vectors
					double angle1=Math.abs(Molecule.getAngleDif(
							Molecule.getAngle(molecule.getAtomX(iAtom),
									molecule.getAtomY(iAtom),
									molecule.getAtomX(molecule.getConnAtom(iAtom, 0)),
									molecule.getAtomY(molecule.getConnAtom(iAtom, 0))
							),
							Molecule.getAngle(molecule.getAtomX(iAtom),
									molecule.getAtomY(iAtom),
									molecule.getAtomX(molecule.getConnAtom(iAtom, 1)),
									molecule.getAtomY(molecule.getConnAtom(iAtom, 1))
							)
					));
					double angle2=Math.abs(Molecule.getAngleDif(
							Molecule.getAngle(molecule.getAtomX(iAtom),
									molecule.getAtomY(iAtom),
									molecule.getAtomX(molecule.getConnAtom(iAtom, 0)),
									molecule.getAtomY(molecule.getConnAtom(iAtom, 0))
							),
							Molecule.getAngle(molecule.getAtomX(iAtom),
									molecule.getAtomY(iAtom),
									molecule.getAtomX(molecule.getConnAtom(iAtom, 2)),
									molecule.getAtomY(molecule.getConnAtom(iAtom, 2))
							)
					));

					double angle3=Math.abs(Molecule.getAngleDif(
							Molecule.getAngle(molecule.getAtomX(iAtom),
									molecule.getAtomY(iAtom),
									molecule.getAtomX(molecule.getConnAtom(iAtom, 1)),
									molecule.getAtomY(molecule.getConnAtom(iAtom, 1))
							),
							Molecule.getAngle(molecule.getAtomX(iAtom),
									molecule.getAtomY(iAtom),
									molecule.getAtomX(molecule.getConnAtom(iAtom, 2)),
									molecule.getAtomY(molecule.getConnAtom(iAtom, 2))
							)
					));

					// System.out.println(angle1+" - "+angle2+" - "+angle3);
					boolean normal=true;
					if ((angle1>angle2) && (angle1>angle3)) {
						if ((angle2+angle3) < Math.PI) {
							normal=false;
							dx = x - 0.5 * (molecule.getAtomX(molecule.getConnAtom(iAtom, 0)) + molecule.getAtomX(molecule.getConnAtom(iAtom, 1)));
							dy = y - 0.5 * (molecule.getAtomY(molecule.getConnAtom(iAtom, 0)) + molecule.getAtomY(molecule.getConnAtom(iAtom, 1)));
						}
					} else if ((angle2>angle1) && (angle2>angle3)) {
						if ((angle1+angle3) < Math.PI) {
							normal=false;
							dx = x - 0.5 * (molecule.getAtomX(molecule.getConnAtom(iAtom, 0)) + molecule.getAtomX(molecule.getConnAtom(iAtom, 2)));
							dy = y - 0.5 * (molecule.getAtomY(molecule.getConnAtom(iAtom, 0)) + molecule.getAtomY(molecule.getConnAtom(iAtom, 2)));
						}
					} else {
						if ((angle1+angle2) < Math.PI) {
							normal=false;
							dx = x - 0.5 * (molecule.getAtomX(molecule.getConnAtom(iAtom, 1)) + molecule.getAtomX(molecule.getConnAtom(iAtom, 2)));
							dy = y - 0.5 * (molecule.getAtomY(molecule.getConnAtom(iAtom, 1)) + molecule.getAtomY(molecule.getConnAtom(iAtom, 2)));
						}
					}

					if (normal) {
						// what is the closest atom to the atom 1 ?
						int stereoAtom2=molecule.getConnAtom(iAtom, 0);
						double distance=Double.MAX_VALUE;
						for (int i=0; i<3; i++) {
							int atom=molecule.getConnAtom(iAtom, i);
							if (atom!=stereoAtom1) {
								double currentDistance=Math.pow(molecule.getAtomX(iAtom)-molecule.getAtomX(atom),2)+Math.pow(molecule.getAtomY(iAtom)-molecule.getAtomY(atom),2);
								if (currentDistance<distance) {
									stereoAtom2=atom;
									distance=currentDistance;
									System.out.println("Minimal distance: "+stereoAtom1+" - "+stereoAtom2+" - "+distance+" - "+molecule.getAtomicNo(stereoAtom2));
								}
							}
						}
						iNew = molecule.addAtom((molecule.getAtomX(stereoAtom1) + molecule.getAtomX(stereoAtom2))/2, (molecule.getAtomY(stereoAtom1) + molecule.getAtomY(stereoAtom2))/2);
					} else { // 3 bonds on the same side ...
						// oposite of middle bond ...
						iNew = molecule.addAtom((float)(x + dx),(float)(y + dy));
					}
				} else {
					iNew = molecule.addAtom((float)(x + dx),(float)(y + dy));
				}
				molecule.setAtomicNo(iNew, atomicNo);
				molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);

			}
			break;
			case 2:
				// these are the most complex ones
				// diastereotopic protons, amines, amides
				int nNeighbours = molecule.getConnAtoms(iAtom);

				if (nNeighbours == 1) {
					// amines ,mides, etc
					double dx = x - molecule.getAtomX(molecule.getConnAtom(iAtom, 0));
					double dy = y - molecule.getAtomY(molecule.getConnAtom(iAtom, 0));
					// rotate by 60 degrees forwards and backwards in the x,y-plane
					int iNew;
					// backward
					iNew = molecule.addAtom((float)(x + (cos60 * dx - sin60 * dy)*H_BOND_RATIO),(float)(y + (sin60 * dx + cos60 * dy)*H_BOND_RATIO));
					molecule.setAtomicNo(iNew, atomicNo);
					molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
					// forward
					iNew = molecule.addAtom((float)(x + (cos300 * dx - sin300 * dy)*H_BOND_RATIO),(float)(y + (sin300 * dx + cos300 * dy)*H_BOND_RATIO));
					molecule.setAtomicNo(iNew, atomicNo);
					molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
				} else if (nNeighbours == 2) {
					int iNew;
					double dx1 = x - molecule.getAtomX(molecule.getConnAtom(iAtom, 0));
					double dy1 = y - molecule.getAtomY(molecule.getConnAtom(iAtom, 0));
					double dx2 = x - molecule.getAtomX(molecule.getConnAtom(iAtom, 1));
					double dy2 = y - molecule.getAtomY(molecule.getConnAtom(iAtom, 1));

					double length1=Math.sqrt(dx1*dx1+dy1*dy1)*H_BOND_RATIO;
					double length2=Math.sqrt(dx2*dx2+dy2*dy2)*H_BOND_RATIO;
					double dx=dx1+dx2;
					double dy=dy1+dy2;
					double length=Math.sqrt(dx*dx+dy*dy);
					double averageLength=(length1+length2)/2;
					dx=dx/length*averageLength;
					dy=dy/length*averageLength;

					int stereoBond=molecule.getStereoBond(iAtom);
					iNew = molecule.addAtom((float)(x + cos30 * dx - sin30 * dy),(float)(y + sin30 * dx + cos30 * dy));
					molecule.setAtomicNo(iNew, atomicNo);
					if (stereoBond>-1) {
						molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
					} else {
						molecule.addBond(iAtom, iNew, Molecule.cBondTypeUp);
					}

					iNew = molecule.addAtom((float)(x + cos330 * dx - sin330 * dy),(float)(y + sin330 * dx + cos330 * dy));
					molecule.setAtomicNo(iNew, atomicNo);
					molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);

				/*
				iNew = molecule.addAtom(x + dx1, y + dy1);
				molecule.setAtomicNo(iNew, atomicNo);
				molecule.addBond(iAtom, iNew, Molecule.cBondTypeUp);

				iNew = molecule.addAtom(x + dx2, y + dy2);
				molecule.setAtomicNo(iNew, atomicNo);
				molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
				*/
				} else {
					for (int iHydrogen = 0; iHydrogen < nHydrogens; iHydrogen++) {
						int iNew = molecule.addAtom((float)x, (float)y);
						molecule.setAtomicNo(iNew, atomicNo);
						molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
					}
				}
				break;
			case 3:
				// CH3, NH3 and the like
			{
				int iNew;
				double dx;
				double dy;
				if (molecule.getConnAtoms(iAtom)>0) {
					dx = (x - molecule.getAtomX(molecule.getConnAtom(iAtom, 0)))*H_BOND_RATIO;
					dy = (y - molecule.getAtomY(molecule.getConnAtom(iAtom, 0)))*H_BOND_RATIO;
					// backwards
					iNew = molecule.addAtom((float)(x + dx),(float)(y + dy));
					molecule.setAtomicNo(iNew, atomicNo);
					molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
					// to the left
					iNew = molecule.addAtom((float)(x - dy),(float)(y + dx));
					molecule.setAtomicNo(iNew, atomicNo);
					molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
					// to the right
					iNew = molecule.addAtom((float)(x + dy),(float)(y - dx));
					molecule.setAtomicNo(iNew, atomicNo);
					molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
				} else {
					dx = molecule.getAverageBondLength(true);
					dy = molecule.getAverageBondLength(true);
					// forward
					iNew = molecule.addAtom((float)(x + dx),(float)(y + dy));
					molecule.setAtomicNo(iNew, atomicNo);
					molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
					// backwards top
					iNew = molecule.addAtom((float)(x - dy*cos60),(float)(y + dx*sin60));
					molecule.setAtomicNo(iNew, atomicNo);
					molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
					// backwards bottom
					iNew = molecule.addAtom((float)(x - dy*cos60),(float)(y - dx*sin60));
					molecule.setAtomicNo(iNew, atomicNo);
					molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
				}
			}
			break;
			default:
				//
			{
				for (int iHydrogen = 0; iHydrogen < nHydrogens; iHydrogen++) {
					int iNew = molecule.addAtom((float)x,(float)y);
					molecule.setAtomicNo(iNew, atomicNo);
					molecule.addBond(iAtom, iNew, Molecule.cBondTypeSingle);
				}
				break;
			}
		}
	}

	// changed from method name from expandHydrogens(); TLS 9.Nov.2015
	public static void addImplicitHydrogens(StereoMolecule molecule) {
		molecule.ensureHelperArrays(Molecule.cHelperNeighbours);
		int nAtomsBefore = molecule.getAtoms();
		for (int iAtom = 0; iAtom < nAtomsBefore; iAtom++) {
			addImplicitHydrogens(molecule, iAtom);
		}
	}

	/* This seems strange: explicit hydrogens are counted twice.
	 * I suggest using one of the following two functions:	TLS 9.Nov.2015
	 *
	public static int getNumberOfExplicitHydrogens(StereoMolecule molecule) {
		molecule.ensureHelperArrays(Molecule.cHelperNeighbours);
		return molecule.getAllAtoms() - molecule.getAtoms();
	}
	public static int getNumberOfAllHydrogens(StereoMolecule molecule) {
		molecule.ensureHelperArrays(Molecule.cHelperNeighbours);
		int nbHydrogens = 0;
		for (int iAtom = 0; iAtom < molecule.getAtoms(); iAtom++)
			nbHydrogens += molecule.getAllHydrogens(iAtom);
		return nbHydrogens;
	}
	*/

	public static int getNumberOfHydrogens(StereoMolecule molecule) {
		molecule.ensureHelperArrays(Molecule.cHelperNeighbours);
		int nbHydrogens = 0;
		for (int iAtom = 0; iAtom < molecule.getAllAtoms(); iAtom++) {
			if (molecule.getAtomicNo(iAtom) == 1)
				nbHydrogens++;
			else
				nbHydrogens += molecule.getPlainHydrogens(iAtom);
		}
		return nbHydrogens;
	}
}
