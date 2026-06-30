package org.openmolecules.chem.interaction.rf;

import com.actelion.research.chem.StereoMolecule;

public class RFInteraction {
	private final int mPAtom,mLAtom,mPType,mLType;
	private final double mDistance,mRelDistance;
	private final InteractionGeometry mL2PGeometry,mP2LGeometry;

	RFInteraction(StereoMolecule protein, StereoMolecule ligand, int pAtom, int lAtom, int pType, int lType, double distance, double relDistance) {
		mPAtom = pAtom;
		mLAtom = lAtom;
		mPType = pType;
		mLType = lType;
		mDistance = distance;
		mRelDistance = relDistance;
		mP2LGeometry = new InteractionGeometry(protein, ligand, pAtom, lAtom);
		mL2PGeometry = new InteractionGeometry(ligand, protein, lAtom, pAtom);
	}

	/**
	 * Calculates the geometry dependent RF-value for the given interaction.
	 * @param uncertaintyHolder null or double[1] to receive the uncertainty value
	 * @return RF-value or NaN in case of unknown atom types
	 */
	public double getRFValue(double[] uncertaintyHolder) {
		return RFKnowledgeBase.getInstance().getRFValue(this, uncertaintyHolder);
	}

	/**
	 * Determines the geometry independent (raw) RF-value for the given interaction.
	 * @return RF-value or NaN in case of unknown atom types
	 */
	public double getRawRFValue() {
		return RFKnowledgeBase.getInstance().getRawRFValue(mLType, mPType);
	}

	/**
	 * Determines the uncertainty of the geometry independent (raw) RF-value for the given interaction.
	 * @return RF-value or NaN in case of unknown atom types
	 */
	public double getRawUncertainty() {
		return RFKnowledgeBase.getInstance().getRawUncertainty(mLType, mPType);
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

	public double getRelDistance() {
		return mRelDistance;
	}

	public double getP2LAngle() {
		return mP2LGeometry.getAngle();
	}

	public double getL2PAngle() {
		return mL2PGeometry.getAngle();
	}

	public double getP2LTorsion() {
		return mP2LGeometry.getTorsion();
	}

	public double getL2PTorsion() {
		return mL2PGeometry.getTorsion();
	}

	public int getP2LGeometryType() {
		return mP2LGeometry.getType();
	}

	public int getL2PGeometryType() {
		return mL2PGeometry.getType();
	}

	public String getP2LGeometryName() {
		return mP2LGeometry.getTypeName();
	}

	public String getL2PGeometryName() {
		return mL2PGeometry.getTypeName();
	}
}
