package org.openmolecules.chem.interaction.rf;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

public class InteractionGeometry {
	public static final String[] TYPE_NAME = {"n<1|n>4", "n1-sp/n1", "n1-sp2", "n1-sp3n2", "n1-sp3n3", "n1-sp3n4",
			"n2sp", "n2sp2", "n2sp3", "n3sp2", "n3sp3", "n4"};
	private static final int TYPE_UNSUPPORTED = 0;
	private static final int TYPE_N1SP_N1N1 = 1;
	private static final int TYPE_N1SP2 = 2;
	private static final int TYPE_N1_SP3N2 = 3;
	private static final int TYPE_N1_SP3N3 = 4;
	private static final int TYPE_N1_SP3N4 = 5;
	private static final int TYPE_N2SP = 6;
	private static final int TYPE_N2SP2 = 7;
	private static final int TYPE_N2SP3 = 8;
	private static final int TYPE_N3SP2 = 9;
	private static final int TYPE_N3SP3 = 10;
	private static final int TYPE_N4 = 11;

	private double mTorsion,mAngle;
	private int mType;
	private final Coordinates mConnectionP;
	private Coordinates mRearP1,mRearP2,mFarP;

	public InteractionGeometry(StereoMolecule mol, int atom) {
		mol.ensureHelperArrays(Molecule.cHelperRings);

		mConnectionP = mol.getAtomCoordinates(atom);
		mRearP2 = null;	// default

		if (mol.getConnAtoms(atom) == 0 || mol.getConnAtoms(atom) > 4) {
			mType = TYPE_UNSUPPORTED;
		}
		else if (mol.getConnAtoms(atom) == 1) {
			int rearAtom = mol.getConnAtom(atom, 0);
			mRearP1 = mol.getAtomCoordinates(rearAtom);
			if (mol.getConnAtoms(rearAtom) == 1 || mol.getAtomPi(rearAtom) == 2) {
				mFarP = new Coordinates(1, 1, 1).add(mRearP1);	// create an arbitrary plane atom
				mType = TYPE_N1SP_N1N1;
			}
			else if (mol.getAtomPi(rearAtom) == 1 || mol.isAromaticAtom(rearAtom) || mol.isFlatNitrogen(rearAtom)) {
				int farIndex = getMostBulkyNeighbourIndex(mol, rearAtom, atom);
				mFarP = mol.getAtomCoordinates(mol.getConnAtom(rearAtom, farIndex));
				mType = TYPE_N1SP2;
			}
			else {
				if (mol.getConnAtoms(rearAtom) == 2) {
					int farAtom = mol.getConnAtom(rearAtom, mol.getConnAtom(rearAtom, 0) == atom ? 1 : 0);
					mFarP = mol.getAtomCoordinates(farAtom);
					mType = TYPE_N1_SP3N2;
				}
				else if (mol.getConnAtoms(rearAtom) == 3) {
					mFarP = new Coordinates(0, 0, 0);
					for (int i=0; i<mol.getConnAtoms(rearAtom); i++)
						if (mol.getConnAtom(rearAtom, i) != atom)
							mFarP.add(mol.getAtomCoordinates(mol.getConnAtom(rearAtom, i)));
					mFarP.scale(0.5);
					mType = TYPE_N1_SP3N3;
				}
				else {
					int farAtom = mol.getConnAtom(rearAtom, getMostBulkyNeighbourIndex(mol, rearAtom, atom));
					mFarP = mol.getAtomCoordinates(farAtom);
					mType = TYPE_N1_SP3N4;
				}
			}
		}
		else if (mol.getConnAtoms(atom) == 2) {
			if (mol.getAtomPi(atom) == 2) {
				int rearIndex = getMostBulkyNeighbourIndex(mol, atom, -1);
				mRearP1 = mol.getAtomCoordinates(mol.getConnAtom(atom, rearIndex));
				mFarP = new Coordinates(1, 1, 1).add(mRearP1);	// create an arbitrary plane atom
				mType = TYPE_N2SP;
			}
			else {
				mRearP1 = new Coordinates(0, 0, 0);
				for (int i=0; i<mol.getConnAtoms(atom); i++)
					mRearP1.add(mol.getAtomCoordinates(mol.getConnAtom(atom, i)));
				mRearP1.scale(0.5);

				int farIndex = getMostBulkyNeighbourIndex(mol, atom, -1);
				mFarP = mol.getAtomCoordinates(mol.getConnAtom(atom, farIndex));
				mType = mol.getAtomPi(atom) == 1 || mol.isAromaticAtom(atom) || mol.isFlatNitrogen(atom) ? TYPE_N2SP2 : TYPE_N2SP3;
			}
		}
		else if (mol.getConnAtoms(atom) == 3) {
			if (mol.getAtomPi(atom) == 1 || mol.isAromaticAtom(atom) || mol.isFlatNitrogen(atom)) {
				// Rear point is a point perpendicular to the plane and behind the interaction atom in the plane.
				// Since at this time the other interaction atom is not known, we have two potential solutions,
				// one on each side of the plane. We calculate both and decide later which one to use.
				Coordinates p1 = mol.getAtomCoordinates(mol.getConnAtom(atom, 0));
				Coordinates p2 = mol.getAtomCoordinates(mol.getConnAtom(atom, 1));
				Coordinates v = mConnectionP.subC(p1).cross(mConnectionP.subC(p2));
				mRearP1 = mConnectionP.addC(v);
				mRearP2 = mConnectionP.subC(v);
				int farIndex = getMostBulkyNeighbourIndex(mol, atom, -1);
				mFarP = mol.getAtomCoordinates(mol.getConnAtom(atom, farIndex));
				mType = TYPE_N3SP2;
			}
			else if (mol.getAtomPi(atom) == 0) {
				mRearP1 = new Coordinates(0, 0, 0);
				for (int i=0; i<mol.getConnAtoms(atom); i++)
					mRearP1.add(mol.getAtomCoordinates(mol.getConnAtom(atom, i)));
				mRearP1.scale(0.333333333333);
				int farIndex = getMostBulkyNeighbourIndex(mol, atom, -1);
				mFarP = mol.getAtomCoordinates(mol.getConnAtom(atom, farIndex));
				mType = TYPE_N3SP3;
			}
		}
		else {
			int rearIndex = getMostBulkyNeighbourIndex(mol, atom, -1);
			mRearP1 = mol.getAtomCoordinates(mol.getConnAtom(atom, rearIndex));
			int farIndex = getMostBulkyNeighbourIndex(mol, atom, mol.getConnAtom(atom, rearIndex));
			mFarP = mol.getAtomCoordinates(mol.getConnAtom(atom, farIndex));
			mType = TYPE_N4;
		}
	}

	public InteractionGeometry(StereoMolecule mol, StereoMolecule remoteMol, int atom, int remoteAtom) {
		this(mol, atom);

		Coordinates remoteP = remoteMol.getAtomCoordinates(remoteAtom);

		if (mType == TYPE_UNSUPPORTED) {
			mAngle = 0;
			mTorsion = 0;
			return;
		}

		Coordinates rearP = (mRearP2 == null
				|| mRearP2.distanceSquared(remoteP) < mRearP1.distanceSquared(remoteP)) ? mRearP1 : mRearP2;

		mAngle = mConnectionP.subC(rearP).getAngle(remoteP.subC(mConnectionP));

		Coordinates v1 = rearP.subC(mFarP);
		Coordinates v2 = mConnectionP.subC(rearP);
		Coordinates v3 = remoteP.subC(mConnectionP);
		Coordinates n1 = v1.cross(v2);
		Coordinates n2 = v2.cross(v3);
		mTorsion = Math.abs(Math.atan2(v2.getLength() * v1.dot(n2), n1.dot(n2)));
	}

	public double getAngle() {
		return mAngle;
	}

	public double getTorsion() {
		return mTorsion;
	}

	public int getType() {
		return mType;
	}

	public String getTypeName() {
		return TYPE_NAME[mType];
	}

	private int getMostBulkyNeighbourIndex(StereoMolecule mol, int atom, int vetoAtom) {
		int maxPenalty = 0;
		int index = -1;
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int connAtom = mol.getConnAtom(atom, i);
			if (connAtom != vetoAtom) {
				int penalty = 16384 * mol.getConnAtoms(connAtom)
						+ (mol.isAromaticBond(mol.getConnBond(atom, i)) ? 8192 : 0)
						+ 1024 * mol.getAllHydrogens(connAtom)
						+ (mol.getAtomicNo(connAtom) > 10 ? 512 : 0)
						+ 255 - mol.getAtomicNo(connAtom);
				if (maxPenalty < penalty) {
					maxPenalty = penalty;
					index = i;
				}
			}
		}
		return index;
	}

	private Coordinates getPointAbovePlane(Coordinates p, Coordinates p1, Coordinates p2, Coordinates remoteP) {
		Coordinates v = p.subC(p1).cross(p.subC(p2));
		Coordinates r1 = p.addC(v);
		if (remoteP == null)
			return r1;

		// if an interacting remote atom is given, then return that point, which is on the other side of the plane
		Coordinates r2 = p.subC(v);
		return r1.distanceSquared(remoteP) < r2.distanceSquared(remoteP) ? r2 : r1;
	}
}
