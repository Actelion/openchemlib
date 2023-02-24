package org.openmolecules.chem.conf.gen;

public class TorsionSetWithEliminationRuleStrain implements Comparable<TorsionSetWithEliminationRuleStrain> {
	private TorsionSet mTorsionSet;
	private double mRuleStrain;

	public TorsionSetWithEliminationRuleStrain(TorsionSet torsionSet, double ruleStrain) {
		mTorsionSet = torsionSet;
		mRuleStrain = ruleStrain;
	}

	public TorsionSet getTorsionSet() {
		return mTorsionSet;
	}

	public double getRuleStrain() {
		return mRuleStrain;
	}

	@Override
	public int compareTo(TorsionSetWithEliminationRuleStrain ts) {
		return mRuleStrain < ts.mRuleStrain ? -1 : mRuleStrain > ts.mRuleStrain ? 1 : 0;
	}

	public boolean equals(TorsionSetWithEliminationRuleStrain ts) {
		return mRuleStrain == ts.mRuleStrain;
	}
}
