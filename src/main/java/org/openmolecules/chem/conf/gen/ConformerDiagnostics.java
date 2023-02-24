package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.conf.Conformer;

import java.util.ArrayList;

public class ConformerDiagnostics {
	private String mIDCode,mCoords,mName;
	private BaseConformer mBaseConformer;
	private double mCollisionIntensity,mLikelihood;
	private int[] mTorsionIndex,mConformerIndex,mCollisionAtoms,mFixedTorsion;
	private ArrayList<String> mEliminationRules;
	private StringBuilder mCollisionLog;
	private boolean mSuccess;

	protected ConformerDiagnostics(TorsionSet ts) {
		mTorsionIndex = ts.getTorsionIndexes();
		mConformerIndex = ts.getConformerIndexes();
		mCollisionLog = new StringBuilder();
		mEliminationRules = new ArrayList<>();
	}

	protected void setConformer(BaseConformer bc, Conformer c) {
		Canonizer can = new Canonizer(c.toMolecule(null));
		mIDCode = can.getIDCode();
		mCoords = can.getEncodedCoordinates();
		mBaseConformer = bc;
		mName = c.getName();
		mFixedTorsion = new int[bc.getRotatableBonds().length];
		mLikelihood = c.getLikelihood();
		for (int i=0; i<mFixedTorsion.length; i++)
			mFixedTorsion[i] = c.getBondTorsion(bc.getRotatableBonds()[i].getBond());
	}

	protected BaseConformer getBaseConformer() {
		return mBaseConformer;
	}

	protected void setCollisionStrain(double v) {
		mCollisionIntensity = v;
	}

	protected void setCollisionAtoms(int[] atoms) {
		mCollisionAtoms = atoms;
	}

	protected void setSuccess(boolean b) {
		mSuccess = b;
	}

	protected void addEliminationRule(String rule) {
		mEliminationRules.add(rule);
	}

	protected void writeCollisionLog(String log) {
		if (mCollisionLog.length() != 0)
			mCollisionLog.append("<NL>");
		mCollisionLog.append(log);
	}

	public String getIDCode() {
		return mIDCode;
	}

	public String getCoords() {
		return mCoords;
	}

	public String getName() {
		return mName;
	}

	public boolean isSuccess() {
		return mSuccess;
	}

	public double getCollisionIntensity() {
		return mCollisionIntensity;
	}

	public double getLikelihood() {
		return mLikelihood;
	}

	public int[] getCollisionAtoms() {
		return mCollisionAtoms;
	}

	public int[] getRigidFragmentIndexes() {
		return mConformerIndex;
	}

	public int[] getTorsionIndexes() {
		return mTorsionIndex;
	}

	public int[] getFixedTorsions() {
		return mFixedTorsion;
	}

	public String getCollisionLog() {
		return mCollisionLog.toString();
	}

	public ArrayList<String> getEliminationRules() {
		return mEliminationRules;
	}
}