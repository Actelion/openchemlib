/*
 * Copyright 2013-2020 Thomas Sander, openmolecules.org
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.so;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

import java.util.ArrayList;

public class SelfOrganizedConformer extends Conformer {
	private static final double	STRAIN_FOR_LIKELIHOOD_FACTOR_10 = 1.36;   // we try to match kcal/mol with strain

	private double	    mTotalStrain,mLikelihood;
	private double[]	mAtomStrain,mRuleStrain;
	private boolean 	mIsUsed;

	public SelfOrganizedConformer(StereoMolecule mol) {
		super(mol);
		mTotalStrain = Double.NaN;
		}

	public SelfOrganizedConformer(SelfOrganizedConformer conformer) {
		super(conformer);
		copyStrainFrom(conformer);
	}

	/**
	 * Copies the conformer's atom coordinates, torsion and strains values to this Conformer.
	 * @param conformer other conformer of the molecule passed in the Constructor
	 */
	public void copyFrom(SelfOrganizedConformer conformer) {
		super.copyFrom(conformer);
		copyStrainFrom(conformer);
		}

	private void copyStrainFrom(SelfOrganizedConformer conformer) {
		mTotalStrain = conformer.mTotalStrain;
		mLikelihood = conformer.mLikelihood;
		mIsUsed = conformer.mIsUsed;
		mAtomStrain = conformer.mAtomStrain == null ? null : conformer.mAtomStrain.clone();
		mRuleStrain = conformer.mRuleStrain == null ? null : conformer.mRuleStrain.clone();
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

		mTotalStrain = 0f;
		for (int atom=0; atom<getMolecule().getAllAtoms(); atom++)
			mTotalStrain += mAtomStrain[atom];

		mLikelihood = -1.0;
		}

	public double getAtomStrain(int atom) {
		return mAtomStrain == null ? Double.NaN : mAtomStrain[atom];
		}

	public double getTotalStrain() {
		return mTotalStrain;
		}

	public double getRuleStrain(int rule) {
		return mRuleStrain == null ? Double.NaN : mRuleStrain[rule];
		}

	/**
	 * Tries to estimate the relative likelihood of this conformer from atom strains
	 * considering an unstrained conformer to have a likelihood of 1.0.
	 * @return conformer likelihood
	 */
	public double getLikelihood() {
		if (mLikelihood == -1)
			mLikelihood = Math.pow(10, -mTotalStrain / STRAIN_FOR_LIKELIHOOD_FACTOR_10);

		return mLikelihood;
		}

	public void invalidateStrain() {
		mTotalStrain = Double.MAX_VALUE;
		mAtomStrain = null;
		mRuleStrain = null;
		}

	public boolean isUsed() {
		return mIsUsed;
		}

	public void setUsed(boolean isUsed) {
		mIsUsed = isUsed;
		}
	}
