/*
 * @(#)TorsionSetStrategyAdaptiveRandom.java
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

public class TorsionSetStrategyAdaptiveRandom extends TorsionSetStrategyRandom {
	private static final int MAX_TRIES_FOR_NEW = 64;
	private boolean mStartWithMostProbable;

	/**
	 * This torsion set strategy produces random sets of torsion indices until a torsion set
	 * collides. Then it updates individual torsion indices of those rotatable bonds
	 * that connect colliding fragments.
	 * Torsion indices are picked either by pure random or with a twisted likelyhood
	 * towards towards those angles, that show higher frequencies in the CSD.
	 * @param rotatableBond
	 * @param preferLikelyTorsions if set then more frequent torsions are picked with higher probability
	 * @param startWithMostProbable if true then the first torsion set returned contains for every bond the most frequent torsion
	 * @param seed
	 * @return
	 */
	public TorsionSetStrategyAdaptiveRandom(RotatableBond[] rotatableBond, Rigid3DFragment[] fragment, boolean preferLikelyTorsions, boolean startWithMostProbable, long seed) {
		super(rotatableBond, fragment, preferLikelyTorsions, seed);
		mStartWithMostProbable = startWithMostProbable;
		}

	@Override
	public TorsionSet createTorsionSet(TorsionSet previousTorsionSet) {
		if (previousTorsionSet == null) {
			if (mStartWithMostProbable)
				return createTorsionSet(new int[mRotatableBond.length], new int[mRigidFragment.length]);
			else
				return super.createTorsionSet(null);
			}

		if (previousTorsionSet.getCollisionIntensitySum() == 0.0)
			return super.createTorsionSet(previousTorsionSet);

		double[] collisionIntensitySum = getBondAndFragmentCollisionIntensities(previousTorsionSet);
		int[] collidingTorsionIndex = previousTorsionSet.getTorsionIndexes();
		int[] collidingConformerIndex = previousTorsionSet.getConformerIndexes();
		boolean[] indexTried = new boolean[collisionIntensitySum.length];
		for (int count=0; count<MAX_TRIES_FOR_NEW; count++) {

			// build a new weighted random scale not including already tried indexes
			if (indexTried[0] || mRotatableBond[0].getTorsionCount() == 1)
				collisionIntensitySum[0] = 0.0;
			int index = 1;
			for (int i=1; i<mRotatableBond.length; i++) {
				if (indexTried[index] || mRotatableBond[i].getTorsionCount() == 1)
					collisionIntensitySum[index] = collisionIntensitySum[index-1];
				else
					collisionIntensitySum[index] += collisionIntensitySum[index-1];
				index++;
				}
			for (int i=0; i<mRigidFragment.length; i++) {
				if (indexTried[index] || mRigidFragment[i].getConformerCount() == 1)
					collisionIntensitySum[index] = collisionIntensitySum[index-1];
				else
					collisionIntensitySum[index] += collisionIntensitySum[index-1];
				index++;
				}

			// if we cannot rotate any of the colliding bonds nor have alternative fragment conformers, then create a random new torsion set
			if (collisionIntensitySum[collisionIntensitySum.length-1] == 0.0)
				return super.createTorsionSet(previousTorsionSet);

			// find a bond or fragment
			double random = getRandom().nextDouble() * collisionIntensitySum[collisionIntensitySum.length-1];
			for (int i=0; i<collisionIntensitySum.length; i++) {
				if (random < collisionIntensitySum[i]) {
					int[] torsionIndex = Arrays.copyOf(collidingTorsionIndex, collidingTorsionIndex.length);
					int[] conformerIndex = Arrays.copyOf(collidingConformerIndex, collidingConformerIndex.length);
					if (i < mRotatableBond.length) {
						for (int j=1; j<mRotatableBond[i].getTorsionCount(); j++) {
							torsionIndex[i]++;
							if (torsionIndex[i] == mRotatableBond[i].getTorsionCount())
								torsionIndex[i] = 0;
							TorsionSet ts = createTorsionSet(torsionIndex, conformerIndex);
							if (isNewTorsionSet(ts))
								return ts;
							}
						}
					else {
						i -= mRotatableBond.length;
						for (int j=1; j<mRigidFragment[i].getConformerCount(); j++) {
							conformerIndex[i]++;
							if (conformerIndex[i] == mRigidFragment[i].getConformerCount())
								conformerIndex[i] = 0;
							TorsionSet ts = createTorsionSet(torsionIndex, conformerIndex);
							if (isNewTorsionSet(ts))
								return ts;
							}
						}

					indexTried[i] = true;
					break;
					}
				}
			}

		return super.createTorsionSet(previousTorsionSet);
		}
	}
