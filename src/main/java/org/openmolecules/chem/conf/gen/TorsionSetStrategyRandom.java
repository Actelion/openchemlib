/*
 * @(#)TorsionSetStrategyRandom.java
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

import java.util.Random;

public class TorsionSetStrategyRandom extends TorsionSetStrategy {
	private static final int MAX_TRIES_FOR_NEW = 64;
	private Random mRandom;
	private boolean mPreferLikelyTorsions;

	/**
	 * This simple TorsionSetStrategy produces random sets of torsion indexes in a loop.
	 * If it finds a new torsion set within a reasonable number of cycles, the set is returned.
	 * If this is not successful or if all permutations of torsion indexes was already delivered,
	 * then it stops creating new permutations and returns null.<br>
	 * Torsion indices are chosen either pure randomly or optionally with a bias
	 * towards those torsion angles, that are more frequently populated in the CSD
	 * and cause less collision strain.
	 * @param rotatableBond
	 * @param preferLikelyTorsions if set then more frequent angles are picked with higher probability
	 * @param seed
	 */
	public TorsionSetStrategyRandom(RotatableBond[] rotatableBond, Rigid3DFragment[] fragment, boolean preferLikelyTorsions, long seed) {
		super(rotatableBond, fragment);
		mPreferLikelyTorsions = preferLikelyTorsions;
		mRandom = (seed == 0) ? new Random() : new Random(seed);
		}

	public Random getRandom() {
		return mRandom;
		}

	/**
	 * Build a new set of torsion angles by pure or biased random picking.
	 * @param previousTorsionSet is not used
	 * @return a random and new torsion index set.
	 */
	@Override
	public TorsionSet createTorsionSet(TorsionSet previousTorsionSet) {
		if (getTorsionSetCount() == getPermutationCount())
			return null;	// if we have found all permutations

		int[] torsionIndex = new int[mRotatableBond.length];
		int[] conformerIndex = new int[mRigidFragment.length];
		TorsionSet ts = null;
		int count = 0;

		do {
			if (count++ == MAX_TRIES_FOR_NEW)
				return null;

			if (mPreferLikelyTorsions) {
				double progress = (double)count/(double)MAX_TRIES_FOR_NEW;
				for (int j=0; j<mRotatableBond.length; j++)
					torsionIndex[j] = mRotatableBond[j].getLikelyRandomTorsionIndex(mRandom.nextDouble(), progress);
				for (int j=0; j<mRigidFragment.length; j++)
					conformerIndex[j] = mRigidFragment[j].getLikelyRandomConformerIndex(mRandom.nextDouble(), progress);
				}
			else {
				for (int j=0; j<mRotatableBond.length; j++)
					torsionIndex[j] = mRandom.nextInt(mRotatableBond[j].getTorsionCount());
				for (int j=0; j<mRigidFragment.length; j++)
					conformerIndex[j] = mRandom.nextInt(mRigidFragment[j].getConformerCount());
				}

			ts = createTorsionSet(torsionIndex, conformerIndex);
			} while (!isNewTorsionSet(ts));

		return ts;
		}
	}
