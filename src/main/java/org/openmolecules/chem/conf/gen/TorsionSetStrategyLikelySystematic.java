/*
 * @(#)TorsionSetStrategyLikelySystematic.java
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

import java.util.Arrays;
import java.util.Comparator;

/**
 * TorsionSetStrategy that systematically creates all possible TorsionSets in
 * batches while keeping a focus on the likelyhood of individual torsions.
 * @author Thomas Sander
 */
public class TorsionSetStrategyLikelySystematic extends TorsionSetStrategy {
	private int[]			mCurrentMaxTorsionIndex,mCurrentMaxConformerIndex;
	private TorsionSet[]	mAvailableTorsionSet;
	private int				mAvailableTorsionSetIndex;

	public TorsionSetStrategyLikelySystematic(RotatableBond[] rotatableBond, Rigid3DFragment[] fragment) {
		super(rotatableBond, fragment);
		mCurrentMaxTorsionIndex = new int[rotatableBond.length];
		mCurrentMaxConformerIndex = new int[fragment.length];
		mAvailableTorsionSet = new TorsionSet[1];
		mAvailableTorsionSet[0] = createTorsionSet(new int[rotatableBond.length], new int[fragment.length]);
		mAvailableTorsionSetIndex = -1;
		}

	/**
	 * Determine the next likely set of torsion angles and return the TorsionSet
	 * containing an array of indexes referring to the torsion value in the
	 * RotatableBond and containing the contribution of this conformer to all
	 * non colliding ones.
	 * @return
	 */
	@Override
	public TorsionSet createTorsionSet(TorsionSet previousTorsionSet) {
		if (mAvailableTorsionSet == null)
			return null;

		if (mAvailableTorsionSetIndex + 1 == mAvailableTorsionSet.length)
			createNextAvailableSet();

		if (mAvailableTorsionSet == null)
			return null;

		mAvailableTorsionSetIndex++;
		return mAvailableTorsionSet[mAvailableTorsionSetIndex];
		}

	private void createNextAvailableSet() {
		int bestIndex = -1;
		boolean isConformerIndex = false;
		double minLoss = Double.MAX_VALUE;
		for (int i=0; i<mCurrentMaxTorsionIndex.length; i++) {
			if (mCurrentMaxTorsionIndex[i] < mRotatableBond[i].getTorsionCount()-1) {
				double loss = mRotatableBond[i].getTorsionLikelyhood(mCurrentMaxTorsionIndex[i])
							/ mRotatableBond[i].getTorsionLikelyhood(mCurrentMaxTorsionIndex[i]+1);
				if (minLoss > loss) {
					minLoss = loss;
					bestIndex = i;
					}
				}
			}
		for (int i=0; i<mCurrentMaxConformerIndex.length; i++) {
			if (mCurrentMaxConformerIndex[i] < mRigidFragment[i].getConformerCount()-1) {
				double loss = mRigidFragment[i].getConformerLikelyhood(mCurrentMaxConformerIndex[i])
						    / mRigidFragment[i].getConformerLikelyhood(mCurrentMaxConformerIndex[i]+1);
				if (minLoss > loss) {
					minLoss = loss;
					bestIndex = i;
					isConformerIndex = true;
					}
				}
			}

		if (bestIndex == -1) {	// no unused torsion set remaining
			mAvailableTorsionSet = null;
			return;
			}

		int count = 1;
		for (int i=0; i<mCurrentMaxTorsionIndex.length; i++)
			if (isConformerIndex || i != bestIndex)
				count *= (mCurrentMaxTorsionIndex[i] + 1);
		for (int i=0; i<mCurrentMaxConformerIndex.length; i++)
			if (!isConformerIndex || i != bestIndex)
				count *= (mCurrentMaxConformerIndex[i] + 1);

		mAvailableTorsionSet = new TorsionSet[count];

		if (isConformerIndex)
			mCurrentMaxConformerIndex[bestIndex]++;
		else
			mCurrentMaxTorsionIndex[bestIndex]++;

		for (int i=0; i<count; i++) {
			int n = i;
			int[] torsionIndex = new int[mCurrentMaxTorsionIndex.length];
			for (int j=0; j<mCurrentMaxTorsionIndex.length; j++) {
				if (isConformerIndex || j != bestIndex) {
					int torsionCount = mCurrentMaxTorsionIndex[j] + 1;
					if (torsionCount != 1) {
						torsionIndex[j] = n % torsionCount;
						n /= torsionCount;
						}
					}
				}
			int[] conformerIndex = new int[mCurrentMaxConformerIndex.length];
			for (int j=0; j<mCurrentMaxConformerIndex.length; j++) {
				if (!isConformerIndex || j != bestIndex) {
					int conformerCount = mCurrentMaxConformerIndex[j] + 1;
					if (conformerCount != 1) {
						conformerIndex[j] = n % conformerCount;
						n /= conformerCount;
						}
					}
				}

			if (isConformerIndex)
				conformerIndex[bestIndex] = mCurrentMaxConformerIndex[bestIndex];
			else
				torsionIndex[bestIndex] = mCurrentMaxTorsionIndex[bestIndex];

			mAvailableTorsionSet[i] = createTorsionSet(torsionIndex, conformerIndex);
			}

		Arrays.sort(mAvailableTorsionSet, new Comparator<TorsionSet>() {
			@Override
			public int compare(TorsionSet ts1, TorsionSet ts2) {
				return ts1.getLikelyhood() == ts2.getLikelyhood() ? 0
					 : ts1.getLikelyhood() > ts2.getLikelyhood() ? -1 : 1;	// we want the highest likelyhood first
				}
			} );
		mAvailableTorsionSetIndex = -1;
		}
	}
