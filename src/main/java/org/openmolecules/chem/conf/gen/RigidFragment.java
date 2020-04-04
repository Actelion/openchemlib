/*
 * @(#)Rigid3DFragment.java
 *
 * Copyright 2013 openmolecules.org, Inc. All Rights Reserved.
 *
 * NOTICE: All information contained herein is, and remains the property
 * of openmolecules.org.  The intellectual and technical concepts contained
 * herein are proprietary to openmolecules.org.
 * Actelion Pharmaceuticals Ltd. is granted a non-exclusive, non-transferable
 * and timely unlimited usage license.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;

public class RigidFragment {
	int   mCoreAtomCount;
	int[] mCoreToFragmentAtom;
	int[] mFragmentToOriginalAtom;
	int[] mExtendedToFragmentAtom;
	int[] mOriginalToExtendedAtom;

	Conformer[] mConformerList;
	double[] mConformerLikelihood;

	public RigidFragment(int coreAtomCount,
	                     int[] coreToFragmentAtom,
	                     int[] fragmentToOriginalAtom,
	                     int[] extendedToFragmentAtom,
	                     int[] originalToExtendedAtom,
	                     Conformer[] conformerList,
	                     double[] conformerLikelyhood) {
		this.mCoreAtomCount = coreAtomCount;
		this.mCoreToFragmentAtom = coreToFragmentAtom;
		this.mFragmentToOriginalAtom = fragmentToOriginalAtom;
		this.mExtendedToFragmentAtom = extendedToFragmentAtom;
		this.mOriginalToExtendedAtom = originalToExtendedAtom;
		this.mConformerList          = conformerList;
		this.mConformerLikelihood    = conformerLikelyhood;
	}

	public Coordinates getCoreCoordinates(int conformer, int atom) {
		return mConformerList[conformer].getCoordinates(mCoreToFragmentAtom[atom]);
	}

	public Coordinates getExtendedCoordinates(int conformer, int atom) {
		return mConformerList[conformer].getCoordinates(mExtendedToFragmentAtom[atom]);
	}

	public int getConformerCount() {
		return mConformerList.length;
	}

	public double getConformerLikelihood(int i) {
		return mConformerLikelihood[i];
	}

	/**
	 * Calculates a random conformer index giving conformers with lower strains
	 * a higher chance to be selected. With a progress value of 0.0 selection
	 * likelyhoods are proportional to conformer likelyhoods due to lower strains.
	 * With increasing progress value the higher frequent conformers get less
	 * and less preferred until 1.0 without any preference.
	 *
	 * @param random
	 * @param progress 0...1
	 */
	public int getLikelyRandomConformerIndex(double random, double progress) {
		double sum = 0;
		for (int t = 0; t < mConformerLikelihood.length; t++) {
			double contribution = (1f - progress) * mConformerLikelihood[t] + progress / mConformerLikelihood.length;
			sum += contribution;
			if (random <= sum)
				return t;
		}
		return mConformerLikelihood.length - 1;  // should never reach this
	}

	public int originalToExtendedAtom(int originalAtom) {
		return mOriginalToExtendedAtom[originalAtom];
	}

	public int coreToOriginalAtom(int atom) {
		return mFragmentToOriginalAtom[mCoreToFragmentAtom[atom]];
	}

	public int extendedToOriginalAtom(int atom) {
		return mFragmentToOriginalAtom[mExtendedToFragmentAtom[atom]];
	}

	public int getConnectionPointCount() {
		return mExtendedToFragmentAtom.length - mCoreToFragmentAtom.length;
	}

	/**
	 * @return count of core atoms, i.e. atoms inside of rotatable bonds
	 */
	public int getCoreSize() {
		return mCoreToFragmentAtom.length;
	}

	/**
	 * @return count of core and 1st shell atoms, i.e. core and rotatable bond atoms
	 */
	public int getExtendedSize() {
		return mExtendedToFragmentAtom.length;
	}
}
