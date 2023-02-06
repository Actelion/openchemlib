package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.conf.Conformer;

import java.util.ArrayList;

public class ConformerDiagnostics {
	private String mIDCode,mCoords;
	private double mCollisionIntensity;
	private int[] mTorsion,mTorsionIndex,mConformerIndex,mCollisionAtoms;
	private ArrayList<String> mEliminationRules;
	private StringBuilder mCollisionLog;
	private boolean mSuccess;

	protected ConformerDiagnostics(TorsionSet ts) {
		mTorsionIndex = ts.getTorsionIndexes();
		mConformerIndex = ts.getConformerIndexes();
		mCollisionLog = new StringBuilder();
		mEliminationRules = new ArrayList<>();
	}

	protected void setConformer(Conformer c) {
		Canonizer can = new Canonizer(c.toMolecule(null));
		mIDCode = can.getIDCode();
		mCoords = can.getEncodedCoordinates();
	}

	protected void setCollisionIntensity(double v) {
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
			mCollisionLog.append("\n");
		mCollisionLog.append(log);
	}

	public String getIDCode() {
		return mIDCode;
	}

	public String getCoords() {
		return mCoords;
	}

	public boolean isSuccess() {
		return mSuccess;
	}

	public double getCollisionIntensity() {
		return mCollisionIntensity;
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

	public String getCollisionLog() {
		return mCollisionLog.toString();
	}

	public ArrayList<String> getEliminationRules() {
		return mEliminationRules;
	}
}