package org.openmolecules.chem.interaction.rf;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;

public class RFInteraction {
	private int mPAtom, mLAtom, mPType, mLType;
	private double mDistance;
	private InteractionGeometry mL2PGeometry,mP2LGeometry;

	RFInteraction(StereoMolecule protein, StereoMolecule ligand, int pAtom, int lAtom, int pType, int lType, double distance) {
		mPAtom = pAtom;
		mLAtom = lAtom;
		mPType = pType;
		mLType = lType;
		mDistance = distance;
		mP2LGeometry = new InteractionGeometry(protein, ligand, pAtom, lAtom);
		mL2PGeometry = new InteractionGeometry(ligand, protein, lAtom, pAtom);
	}

	public int getPAtom() {
		return mPAtom;
	}

	public int getLAtom() {
		return mLAtom;
	}

	public int getPType() {
		return mPType;
	}

	public int getLType() {
		return mLType;
	}

	public double getDistance() {
		return mDistance;
	}

	public double getP2LAngle() {
		return mP2LGeometry.mAngle;
	}

	public double getL2PAngle() {
		return mL2PGeometry.mAngle;
	}

	public double getP2LTorsion() {
		return mP2LGeometry.mTorsion;
	}

	public double getL2PTorsion() {
		return mL2PGeometry.mTorsion;
	}

	public String getP2LGeometryType() {
		return mP2LGeometry.mGeometryType;
	}

	public String getL2PGeometryType() {
		return mL2PGeometry.mGeometryType;
	}

	static class InteractionGeometry {
		private double mTorsion,mAngle;
		private String mGeometryType = null;

		public InteractionGeometry(StereoMolecule mol, StereoMolecule remoteMol, int atom, int remoteAtom) {
			Coordinates p = mol.getAtomCoordinates(atom);
			Coordinates remoteP = remoteMol.getAtomCoordinates(remoteAtom);
			Coordinates rearP = null;
			Coordinates farP = null;	// defines plane with p and rearP

			if (mol.getConnAtoms(atom) == 1) {
				int rearAtom = mol.getConnAtom(atom, 0);
				rearP = mol.getAtomCoordinates(rearAtom);
				if (mol.getAtomPi(rearAtom) == 2) {
					farP = new Coordinates(1, 1, 1).add(rearP);	// create an arbitrary plane atom
					mGeometryType = "n1-sp";
				}
				else if (mol.getAtomPi(rearAtom) == 1 || mol.isAromaticAtom(rearAtom)) {
					int farAtom = mol.getConnAtom(rearAtom, mol.getConnAtom(rearAtom, 0) == atom ? 1 : 0);
					farP = mol.getAtomCoordinates(farAtom);
					mGeometryType = "n1-sp2";
				}
				else if (mol.getAtomPi(rearAtom) == 0) {
					if (mol.getConnAtoms(rearAtom) == 1) {
						farP = new Coordinates(1, 1, 1).add(rearP);	// create an arbitrary plane atom
						mGeometryType = "n1-sp3n1";
					}
					else if (mol.getConnAtoms(rearAtom) == 2) {
						int farAtom = mol.getConnAtom(rearAtom, mol.getConnAtom(rearAtom, 0) == atom ? 1 : 0);
						farP = mol.getAtomCoordinates(farAtom);
						mGeometryType = "n1-sp3n2";
					}
					else if (mol.getConnAtoms(rearAtom) == 3) {
						farP = new Coordinates(0, 0, 0);
						for (int i=0; i<mol.getConnAtoms(rearAtom); i++)
							if (mol.getConnAtom(rearAtom, i) != atom)
								farP.add(mol.getAtomCoordinates(mol.getConnAtom(rearAtom, i)));
						farP.scale(0.5);
						mGeometryType = "n1-sp3n3";
					}
					else if (mol.getConnAtoms(rearAtom) == 4) {
						int maxAtomicNo = -1;
						int farAtom = -1;
						for (int i=0; i<mol.getConnAtoms(rearAtom); i++) {
							int connAtom = mol.getConnAtom(rearAtom, i);
							if (connAtom != atom && maxAtomicNo < mol.getAtomicNo(connAtom)) {
								farAtom = connAtom;
								maxAtomicNo = mol.getAtomicNo(connAtom);
							}
						}
						farP = mol.getAtomCoordinates(farAtom);
						mGeometryType = "n1-sp3n4";
					}
				}
			}
			else if (mol.getConnAtoms(atom) == 2) {
				if (mol.getAtomPi(atom) == 2) {
					rearP = mol.getAtomCoordinates(mol.getConnAtom(atom, 0));
					farP = new Coordinates(1, 1, 1).add(rearP);	// create an arbitrary plane atom
					mGeometryType = "n2sp";
				}
				else if (mol.getAtomPi(atom) == 1 || mol.isAromaticAtom(atom)) {
					rearP = mol.getAtomCoordinates(mol.getConnAtom(atom, 0)).centerC(mol.getAtomCoordinates(mol.getConnAtom(atom, 1)));
					farP = mol.getAtomCoordinates(mol.getConnAtom(atom, 0));
					mGeometryType = "n2sp2";
				}
				else if (mol.getAtomPi(atom) == 0) {
					rearP = mol.getAtomCoordinates(mol.getConnAtom(atom, 0)).centerC(mol.getAtomCoordinates(mol.getConnAtom(atom, 1)));
					farP = mol.getAtomCoordinates(mol.getConnAtom(atom, 0));
					mGeometryType = "n2sp3";
				}
			}
			else if (mol.getConnAtoms(atom) == 3) {
				if (mol.getAtomPi(atom) == 1 || mol.isAromaticAtom(atom)) {
					int minScore = Integer.MAX_VALUE;	// higher score of aromatic, ring, atomicNo, neighbour count
					int minAtom = -1;
					for (int i=0; i<mol.getConnAtoms(atom); i++) {
						int connAtom = mol.getConnAtom(atom, i);
						int connBond = mol.getConnBond(atom, i);
						int score = mol.getAtomicNo(atom)
								+ 256 * mol.getConnAtoms(atom)
								+ (mol.isRingBond(connBond) ? 1024 : 0)
								+ (mol.isAromaticBond(connBond) ? 2048 : 0);
						if (minScore > score) {
							minScore = score;
							minAtom = connAtom;
						}
					}
					rearP = new Coordinates(0, 0, 0);
					for (int i=0; i<mol.getConnAtoms(atom); i++)
						if (mol.getConnAtom(atom, i) != minAtom)
							rearP.add(mol.getAtomCoordinates(mol.getConnAtom(atom, i)));
					rearP.scale(0.5);
					int farAtom = mol.getConnAtom(atom, mol.getConnAtom(atom, 0) == minAtom ? 1 : 0);
					farP = mol.getAtomCoordinates(farAtom);
					mGeometryType = "n3sp2";
				}
				else if (mol.getAtomPi(atom) == 0) {
					rearP = mol.getAtomCoordinates(mol.getConnAtom(atom, 0)).addC(
							mol.getAtomCoordinates(mol.getConnAtom(atom, 1))).add(
							mol.getAtomCoordinates(mol.getConnAtom(atom, 2))).scale(1.0/3.0);
					int maxAtomicNo = -1;
					int farAtom = -1;
					for (int i=0; i<mol.getConnAtoms(atom); i++) {
						int connAtom = mol.getConnAtom(atom, i);
						if (maxAtomicNo < mol.getAtomicNo(connAtom)) {
							farAtom = connAtom;
							maxAtomicNo = mol.getAtomicNo(connAtom);
						}
					}
					farP = mol.getAtomCoordinates(farAtom);
					mGeometryType = "n3sp3";
				}
			}
			else if (mol.getConnAtoms(atom) == 4) {
				int minScore = Integer.MAX_VALUE;	// higher score of atomicNo, neighbour count
				int minAtom = -1;
				for (int i=0; i<mol.getConnAtoms(atom); i++) {
					int connAtom = mol.getConnAtom(atom, i);
					int score = mol.getAtomicNo(atom)
							+ 256 * mol.getConnAtoms(atom);
					if (minScore > score) {
						minScore = score;
						minAtom = connAtom;
					}
				}
				rearP = new Coordinates(0, 0, 0);
				for (int i=0; i<mol.getConnAtoms(atom); i++)
					if (mol.getConnAtom(atom, i) != minAtom)
						rearP.add(mol.getAtomCoordinates(mol.getConnAtom(atom, i)));
				rearP.scale(1.0/3.0);
				int farAtom = mol.getConnAtom(atom, mol.getConnAtom(atom, 0) == minAtom ? 1 : 0);
				farP = mol.getAtomCoordinates(farAtom);
				mGeometryType = "n4";
			}

			if (rearP == null || farP == null)
				return;

			mAngle = p.subC(rearP).getAngle(remoteP.subC(p));

			Coordinates v1 = rearP.subC(farP);
			Coordinates v2 = p.subC(rearP);
			Coordinates v3 = remoteP.subC(p);
			Coordinates n1 = v1.cross(v2);
			Coordinates n2 = v2.cross(v3);
			mTorsion = Math.abs(Math.atan2(v2.getLength() * v1.dot(n2), n1.dot(n2)));
		}
	}
}
