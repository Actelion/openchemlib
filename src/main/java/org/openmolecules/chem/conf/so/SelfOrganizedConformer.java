/*
 * @(#)SelfOrganizedConformer.java
 *
 * Copyright 2014 openmolecules.org, Inc. All Rights Reserved.
 *
 * NOTICE: All information contained herein is, and remains the property
 * of openmolecules.org.  The intellectual and technical concepts contained
 * herein are proprietary to openmolecules.org.
 * Actelion Pharmaceuticals Ltd. is granted a non-exclusive, non-transferable
 * and timely unlimited usage license.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.so;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.TorsionDescriptor;
import com.actelion.research.chem.conf.TorsionDescriptorHelper;

import java.util.ArrayList;

public class SelfOrganizedConformer extends Conformer {
	private static final double	MAX_ATOM_STRAIN = 0.01;

	// adjust such that a molecule strain of atom count times MAX_AVERARAGE_ATOM_STRAIN reduces frequency by factor 100
	private static final double	MAX_AVERARAGE_ATOM_STRAIN = 0.001;

	private double	    mMaxAtomStrain,mTotalStrain;
	private double[]	mAtomStrain,mRuleStrain;
	private boolean 	mIsUsed;
	private TorsionDescriptor mTorsionDescriptor;

	public SelfOrganizedConformer(StereoMolecule mol) {
		super(mol);
		}

	/**
	 * Checks whether the total strain of this Conformer is larger than that of conformer,
	 * assuming that the calculated strain values are up-to-date.
	 * @param conformer
	 * @return
	 */
	public boolean isWorseThan(SelfOrganizedConformer conformer) {
		return mTotalStrain > conformer.mTotalStrain;
		}

	public void calculateStrain(ArrayList<ConformationRule> ruleList) {
		if (mAtomStrain != null)
			return;

		mAtomStrain = new double[getMolecule().getAllAtoms()];
		mRuleStrain = new double[ConformationRule.RULE_NAME.length];

		for (ConformationRule rule:ruleList)
			if (rule.isEnabled())
				mRuleStrain[rule.getRuleType()] += rule.addStrain(this, mAtomStrain);

		mMaxAtomStrain = 0f;
		mTotalStrain = 0f;
		for (int atom=0; atom<getMolecule().getAllAtoms(); atom++) {
			mTotalStrain += mAtomStrain[atom];
			if (mMaxAtomStrain < mAtomStrain[atom])
				mMaxAtomStrain = mAtomStrain[atom];
			}
		}

	public double getAtomStrain(int atom) {
		return mAtomStrain[atom];
		}

	public double getRuleStrain(int rule) {
		return mRuleStrain[rule];
		}

	public double getHighestAtomStrain() {
		return mMaxAtomStrain;
		}

	public double getTotalStrain() {
		return mTotalStrain;
		}

	/**
	 * Tries to estimate the relative likelyhood of this conformer from atom strains
	 * considering an unstrained conformer to have a likelyhood of 1.0.
	 * @return conformer likelyhood
	 */
	public double getLikelyhood() {
		return Math.pow(100, -mTotalStrain / (SelfOrganizedConformer.MAX_AVERARAGE_ATOM_STRAIN*getMolecule().getAllAtoms()));
		}

	/**
	 * @param ruleList may be null, if isAcceptable() was called earlier and neither ruleList not conformer were changes since
	 * @return
	 */
	protected boolean isAcceptable(ArrayList<ConformationRule> ruleList) {
		if (ruleList != null)
			calculateStrain(ruleList);
		return (mMaxAtomStrain < MAX_ATOM_STRAIN
			 && mTotalStrain < MAX_AVERARAGE_ATOM_STRAIN * getMolecule().getAllAtoms());
		}

	public void invalidateStrain() {
		mAtomStrain = null;
		mRuleStrain = null;
		}

	/**
	 * Calculates the torsion descriptor for the current coordinates.
	 * Use calculateRotatableBondsForDescriptor() once and pass it
	 * for every new conformer to this method.
	 * @param rotatableBond set of rotatable bonds to be considered
	 */
	public void calculateDescriptor(int[] rotatableBond) {
		mTorsionDescriptor = new TorsionDescriptorHelper(getMolecule(), rotatableBond).getTorsionDescriptor(this);
		}

	/**
	 * Returns true, if none of the torsion angles between both conformers
	 * are more different than TorsionDescriptor.TORSION_EQUIVALENCE_TOLERANCE;
	 * Calling this method requires that calculateDescriptor() has been called earlier.
	 * @param conformer
	 * @return true if all torsions are similar
	 */
	public boolean equals(SelfOrganizedConformer conformer) {
		return mTorsionDescriptor.equals(conformer.mTorsionDescriptor);
		}

	public boolean isUsed() {
		return mIsUsed;
		}

	public void setUsed(boolean isUsed) {
		mIsUsed = isUsed;
		}
	}
