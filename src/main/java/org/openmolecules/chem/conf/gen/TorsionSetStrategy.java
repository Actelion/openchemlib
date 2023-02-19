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

package org.openmolecules.chem.conf.gen;

import com.actelion.research.util.SortedList;
import com.actelion.research.util.UniqueList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

/**
 * Knowing all rotatable bonds of an underlying molecule and knowing those rigid fragments
 * that are connected by them, the TorsionSetStrategy provides a mechanism to deliver
 * valid and unique torsion sets, each effectively defining an individual conformer.
 * While the strategies of the implementing classes differ (random, biased random, systematic, ...)
 * this parent class provides the following common functionality:<br>
 * <li>maintains a history of already delivered torsion sets.
 * <li>provides a novelty check for torsion sets suggested by derived classes.
 * <li>check new torsion sets for atom collisions.
 * <li>maintains a list of torsion set subsets that cause atom collisions.
 * <li>provides a quick check, whether a new torsion set contains a subset causing collisions.
 * <br><br>A TorsionSetStrategy is employed by a ConformerGenerator, which repeatedly calls
 * getNextTorsionIndexes(), creates the conformer, checks for atom collisions and reports back,
 * whether and how serious atom collisions occurred.
 */
public abstract class TorsionSetStrategy {
	private static final int MAX_TOTAL_COUNT = 10000;   // maximum of distinct torsion sets to be checked

	// Slightly colliding torsion sets, which are initially not considered acceptable, are collected in a cache.
	// After SECOND_CHOICE_START_INDEX torsion sets have been generated we start tolerating such torsion sets.
	// From that index to MAX_TOTAL_COUNT we linearly increase a collision tolerance value below which a
	// torsion set is considered acceptable. Thus, from SECOND_CHOICE_START_INDEX we check the cache, whether
	// it contains a tolerable torsion set before we generate a new one.

	// divisor to calculate count index from which to start collecting second choice conformers
	private static final int SECOND_CHOICE_START_FRACTION = 5;

	// We determine the lowest collision intensity of all conformers before starting to include
	// second choice conformers. Unavoidable inherent molecule strain may cause this lowest value
	// to be above 0.0, e.g. in a,a'-substituted naphtalines. This low value may not be higher than
	// MAX_COLLISION_INTENSITY_BASE.
	private static final double MAX_COLLISION_INTENSITY_BASE = 0.5; // sum of squares of distances below VDW radii

	// tolerated maximum collision intensity difference on top of the lowest collision intensity found
	private static final double SECOND_CHOICE_MAX_TOLERANCE = 0.2; // sum of squares of distances below VDW radii

	public static final double MAX_ALLOWED_COLLISION_INTENSITY = MAX_COLLISION_INTENSITY_BASE + SECOND_CHOICE_MAX_TOLERANCE;

	protected ConformerGenerator mConformerGenerator;
	protected RotatableBond[] mRotatableBond;
	protected RigidFragment[] mRigidFragment;
	private TorsionSetEncoder mTorsionSetEncoder;
	private int mFragmentCount,mCollisionCount,mTotalCount,mMaxTotalCount,mPermutationCount;
	private boolean mUsingSecondChoices;
	private int[][][] mBondsBetweenFragments;
	private int[][] mConnFragmentNo;
	private int[][] mConnRotatableBondNo;
	private int[] mGraphFragment;
	private int[] mGraphBond;
	private int[] mGraphParent;
	private boolean[] mGraphFragmentHandled;
	private double mLowestCollisionStrain;
	private UniqueList<TorsionSet> mTorsionSetList;
	private SortedList<TorsionSet> mSecondChoiceList;

	public TorsionSetStrategy(ConformerGenerator conformerGenerator) {
		mConformerGenerator = conformerGenerator;
		mRotatableBond = conformerGenerator.getRotatableBonds();
		mRigidFragment = conformerGenerator.getRigidFragments();
		mTorsionSetEncoder = new TorsionSetEncoder(mRigidFragment, mRotatableBond);

		// create arrays of neighbor fragment no's
		mFragmentCount = 0;
		for (RotatableBond rb:mRotatableBond)
			mFragmentCount = Math.max(mFragmentCount, Math.max(1+rb.getFragmentNo(0), 1+rb.getFragmentNo(1)));
		int[] count = new int[mFragmentCount];
		for (RotatableBond rb:mRotatableBond) {
			count[rb.getFragmentNo(0)]++;
			count[rb.getFragmentNo(1)]++;
			}
		mConnFragmentNo = new int[count.length][];
		mConnRotatableBondNo = new int[count.length][];
		for (int i=0; i<count.length; i++) {
			mConnFragmentNo[i] = new int[count[i]];
			mConnRotatableBondNo[i] = new int[count[i]];
			}
		Arrays.fill(count, 0);
		for (int i=0; i<mRotatableBond.length; i++) {
			int f1 = mRotatableBond[i].getFragmentNo(0);
			int f2 = mRotatableBond[i].getFragmentNo(1);
			mConnFragmentNo[f1][count[f1]] = f2;
			mConnFragmentNo[f2][count[f2]] = f1;
			mConnRotatableBondNo[f1][count[f1]] = i;
			mConnRotatableBondNo[f2][count[f2]] = i;
			count[f1]++;
			count[f2]++;
			}

		mGraphFragment = new int[mFragmentCount];
		mGraphBond = new int[mFragmentCount];
		mGraphParent = new int[mFragmentCount];
		mGraphFragmentHandled = new boolean[mFragmentCount];

		// initialize array for rotatable bond sequences from fragment pairs
		mBondsBetweenFragments = new int[mFragmentCount][][];
		for (int f1=1; f1<mFragmentCount; f1++) {
			mBondsBetweenFragments[f1] = new int[f1][];
			for (int f2=0; f2<f1; f2++)
				mBondsBetweenFragments[f1][f2] = getRotatableBondsBetween(f1, f2);
			}

		mPermutationCount = 1;

		for (RotatableBond rb:mRotatableBond)
			mPermutationCount *= rb.getTorsionCount();
		for (RigidFragment rf:mRigidFragment)
			mPermutationCount *= rf.getConformerCount();
		if (mPermutationCount <= 0)
			mPermutationCount = Integer.MAX_VALUE;

		mTorsionSetList = new UniqueList<>();
		mSecondChoiceList = new SortedList<>(Comparator.comparingDouble(ts -> ((TorsionSet)ts).getCollisionIntensitySum()));
		mCollisionCount = 0;
		mTotalCount = 0;
		mLowestCollisionStrain = MAX_COLLISION_INTENSITY_BASE;
		mUsingSecondChoices = false;
		mMaxTotalCount = Math.min(MAX_TOTAL_COUNT, mPermutationCount);
		}

	public int[] getBondsBetweenFragments(int f1, int f2) {
		return mBondsBetweenFragments[f1][f2];
		}

	public void setMaxTotalCount(int maxTotalCount) {
		mMaxTotalCount = Math.min(maxTotalCount, mPermutationCount);
		}

	/*	public UniqueList<TorsionSet> getTorsionSetList() {
		return mTorsionSetList;
		}
*/

	/**
	 * @return pre-calculated number of all possible torsion index permutations
	 */
	public int getPermutationCount() {
		return mPermutationCount;
		}

	/**
	 * @return number of generated torsion sets that resulted in collisions
	 */
	public int getFailureCount() {
		return mCollisionCount;
		}

	/**
	 * @return number of generated torsion sets till now
	 */
	public int getTorsionSetCount() {
		return mTotalCount;
		}

	public TorsionSetEncoder getTorsionSetEncoder() {
		return mTorsionSetEncoder;
	}

	/**
	 * Creates a new TorsionSet object from a torsion index array.
	 * An overall likelihood of the new torsion set is calculated
	 * by multiplying individual likelihoods of the specific torsion angles
	 * of all underlying RotatableBonds.
	 * @param torsionIndex
	 * @param conformerIndex null or conformer index for every rigid fragment
	 * @return
	 */
	protected TorsionSet createTorsionSet(int[] torsionIndex, int[] conformerIndex) {
		double likelihood = 1.0;
		if (conformerIndex != null)
			for (int rf=0; rf<mRigidFragment.length; rf++)
				likelihood *= mRigidFragment[rf].getConformerLikelihood(conformerIndex[rf]);
		BaseConformer baseConformer = mConformerGenerator.getBaseConformer(conformerIndex == null ?
				new int[mRigidFragment.length] : conformerIndex);
		for (int rb=0; rb<mRotatableBond.length; rb++)
			likelihood *= baseConformer.getTorsionLikelyhood(rb, torsionIndex[rb]);
		return new TorsionSet(torsionIndex, conformerIndex, mTorsionSetEncoder, likelihood);
		}

	protected boolean isNewTorsionSet(TorsionSet ts) {
		return !mTorsionSetList.contains(ts);
		}

	/**
	 * Creates the next set of torsion indexes to be tried by the ConformerGenerator.
	 * If the previous set obtained by this method resulted in an atom collision,
	 * then the collision intensity matrix among all rigid fragments and the overall
	 * collision intensity sum must have been reported to the torsion set.
	 * Torsion sets returned by this method are guaranteed to avoid torsion sequences that
	 * had previously caused collisions.
	 * To achieve this goal, TorsionSetStrategy repeatedly requests new TorsionSets
	 * from the strategy implementation and returns the first set that is not in conflict
	 * with the collision rules already collected or until no further set exists.
	 * @param previousTorsionSet delivered by this method (or null if it is the first call)
	 * @return torsion index set that adheres to already known collision rules
	 */
	public final TorsionSet getNextTorsionSet(TorsionSet previousTorsionSet, ConformerSetDiagnostics diagnostics) {
		// Some molecules have unavoidable internal strains,
		// which we try to determine until we start returning second choices.
		if (previousTorsionSet != null
		 && mLowestCollisionStrain > previousTorsionSet.getCollisionIntensitySum())
			mLowestCollisionStrain = previousTorsionSet.getCollisionIntensitySum();

		if (mUsingSecondChoices)
			return getBestSecondChoice();

		if (mTotalCount == mMaxTotalCount) {
			if (diagnostics != null)
				diagnostics.setExitReason("maxTotal(" + mMaxTotalCount + ") reached A; collisions:" + mCollisionCount);

			return null;
			}

		if (previousTorsionSet != null
		 && !previousTorsionSet.isUsed()	// if it was already used, then the collision was considered tolerable by the ConformerGenerator
		 && previousTorsionSet.getCollisionIntensitySum() != 0) {
			processCollisions(previousTorsionSet);

			if (previousTorsionSet.getCollisionIntensitySum() < mLowestCollisionStrain + SECOND_CHOICE_MAX_TOLERANCE)
				mSecondChoiceList.add(previousTorsionSet);

			mCollisionCount++;
			}

		double tolerance = calculateCollisionTolerance();

		TorsionSet ts = getSecondChoiceTorsionSet(tolerance);
		if (ts != null)
			return ts;

		ts = createTorsionSet(previousTorsionSet);

		while (ts != null
				&& matchesEliminationRule(ts, tolerance, mConformerGenerator.getBaseConformer(
						ts.getConformerIndexes()).getEliminationRules())) {
			mCollisionCount++;
			mTotalCount++;
			mTorsionSetList.add(ts);

			if (mTotalCount == mMaxTotalCount) {
				if (diagnostics != null)
					diagnostics.setExitReason("maxTotal(\"+mMaxTotalCount+\") reached B; collisions:"+mCollisionCount);
				return null;
				}

			tolerance = calculateCollisionTolerance();
			ts = getSecondChoiceTorsionSet(tolerance);
			if (ts != null)
				return ts;

			ts = createTorsionSet(ts);
			}

		if (ts == null) {
			if (mSecondChoiceList.size() != 0) {
				mUsingSecondChoices = true;

				// remove all second choices beyond smallest found collision intensity plus tolerance
				while (mSecondChoiceList.size() != 0
					&& mSecondChoiceList.get(mSecondChoiceList.size()-1).getCollisionIntensitySum() > mLowestCollisionStrain + SECOND_CHOICE_MAX_TOLERANCE)
					mSecondChoiceList.remove(mSecondChoiceList.size()-1);

				ts = getBestSecondChoice();
				}

			if (ts == null && diagnostics != null)
				diagnostics.setExitReason("Inner strategy stop criterion hit");

			return ts;
			}

		mTotalCount++;
		mTorsionSetList.add(ts);

		return ts;
		}

	private TorsionSet getSecondChoiceTorsionSet(double tolerance) {
		if (tolerance == 0.0)
			return null;

		if (mSecondChoiceList.size() == 0
		 || mSecondChoiceList.get(0).getCollisionIntensitySum() > tolerance)
			return null;

		TorsionSet ts = mSecondChoiceList.get(0);
		mSecondChoiceList.remove(0);
		return ts;
		}

	private TorsionSet getBestSecondChoice() {
		if (mSecondChoiceList.size() == 0)
			return null;

		int index = -1;

		if (this instanceof TorsionSetStrategyRandom) {
			index = ((TorsionSetStrategyRandom)this).getRandom().nextInt(mSecondChoiceList.size());
			}
		else {
			double minCollisionIntensitySum = Double.MAX_VALUE;
			for (int i=0; i<mSecondChoiceList.size(); i++) {
				TorsionSet ts = mSecondChoiceList.get(i);
				if (minCollisionIntensitySum > ts.getCollisionIntensitySum()) {
					minCollisionIntensitySum = ts.getCollisionIntensitySum();
					index = i;
					}
				}
			}

		TorsionSet ts = mSecondChoiceList.get(index);
		mSecondChoiceList.remove(index);
		return ts;
		}

	/**
	 * If no collision free torsion set can be constructed, this method is called
	 * to get the torsion set with the least atom collision strain.
	 * @return
	 */
	public TorsionSet getBestCollidingTorsionIndexes() {
		double bestCollisionIntensity = Double.MAX_VALUE;
		TorsionSet bestTorsionSet = null;
		for (int i=0; i<mTorsionSetList.size(); i++) {
			TorsionSet ts = mTorsionSetList.get(i);
			if (bestCollisionIntensity > ts.getCollisionIntensitySum()) {
				bestCollisionIntensity = ts.getCollisionIntensitySum();
				bestTorsionSet = ts;
				}
			}
		return bestTorsionSet;
		}

	/**
	 * Provide the next set of torsion angles using a specific strategy and
	 * considering, which angle combinations were already tried, which had failed, and
	 * (depending on the strategy) considering the likelyhoods of particular torsions.
	 * The implementation must not return the same torsion set multiple times, either by using
	 * a systematic strategy or by calling isNewTorsionSet() to sort out redundant ones.
	 * @param previousTorsionSet previous torsion set which may or may not be used by strategy
	 * @return
	 */
	public abstract TorsionSet createTorsionSet(TorsionSet previousTorsionSet);

	/**
	 * Must be called if the torsion indexes delivered with getNextTorsionIndexes()
	 * caused an atom collision.
	 * Completes a dictionary of bond sequences with specific torsions that causes
	 * atom collisions. Depending on the strategy implementation, this disctionary
	 * is taken into account, when creating new torsion sets.
	 */
	private void processCollisions(TorsionSet torsionSet) {
		double[][] collisionIntensityMatrix = torsionSet.getCollisionIntensityMatrix();
		int[] torsionIndex = torsionSet.getTorsionIndexes();
		for (int f1=1; f1<collisionIntensityMatrix.length; f1++) {
			if (collisionIntensityMatrix[f1] != null) {
				for (int f2=0; f2<f1; f2++) {
					if (collisionIntensityMatrix[f1][f2] != 0f) {
						int[] rotatableBondIndex = mBondsBetweenFragments[f1][f2];
						TorsionSetEliminationRule rule = new TorsionSetEliminationRule(torsionIndex,
								rotatableBondIndex, collisionIntensityMatrix[f1][f2], mTorsionSetEncoder);
						boolean isCovered = false;
						ArrayList<TorsionSetEliminationRule> obsoleteList = null;
						ArrayList<TorsionSetEliminationRule> currentList = mConformerGenerator.getBaseConformer(
								torsionSet.getConformerIndexes()).getEliminationRules();
						for (TorsionSetEliminationRule er:currentList) {
							if (er.isCovered(rule)) {
								isCovered = true;
								break;
								}
							if (er.isMoreGeneral(rule)) {
								if (obsoleteList == null)
									obsoleteList = new ArrayList<>();
								obsoleteList.add(er);
								}
							}
						if (!isCovered) {
							if (obsoleteList != null)
								currentList.removeAll(obsoleteList);
							currentList.add(rule);
//							eliminateTorsionSets(mask, data);
// currently validity checking is done in getNextTorsionIndexes() against the list of elimination rules
							}
						}
					}
				}
			}
		}

	/**
	 * Calculates for every rotatable bond and for every rigid fragment a collision intensity sum
	 * for the given torsion/conformer state from the collision rules already known.
	 * If the given torsion or conformer index of the collidingTorsionSet is covered by a rule
	 * then the rule's collision intensity is added to all involved bond's intensity sums.
	 * Strategies may focus on those bond torsions or conformers first that have the highest
	 * collision intensity sums.
	 * @param collidingTorsionSet
	 * @return collision intensity sum for every rotatable bond and rigid fragment
	 */
	protected double[] getBondAndFragmentCollisionIntensities(TorsionSet collidingTorsionSet) {
		double[] collisionIntensity = new double[mRotatableBond.length+mRigidFragment.length];
		ArrayList<TorsionSetEliminationRule> elimRules = mConformerGenerator.getBaseConformer(
				collidingTorsionSet.getConformerIndexes()).getEliminationRules();

		for (TorsionSetEliminationRule er:elimRules)
			if (collidingTorsionSet.matches(er, 0.0))
				for (int i=0; i<collisionIntensity.length; i++)
					if (mTorsionSetEncoder.isMaskSet(er, i))
						collisionIntensity[i] += er.getCollisionIntensity();

		return collisionIntensity;
		}

	/**
	 * If the implementation of the TorsionStrategy is caching some kind of a pre-calculated
	 * list of TorsionSets, then those sets should be removed that match the elimination
	 * condition defined by mask and data, i.e. TorsionSets that return true on matches(mask, data).
	 * @param mask
	 * @param data
	 *
	public abstract void eliminateTorsionSets(long[] mask, long[] data);	*/


	/**
	 * With best current knowledge about colliding torsion combinations
	 * and based on the individual frequencies of currently active torsions
	 * this method returns the conformers's overall contribution to the
	 * total set of non colliding conformers.
	 * @return this conformer's contribution to all conformers
	 */
	public double getContribution(TorsionSet torsionSet) {
		double likelyhood = 1.0;
		for (int rf=0; rf<mRigidFragment.length; rf++)
			likelyhood *= mRigidFragment[rf].getConformerLikelihood(torsionSet.getConformerIndexes()[rf]);
		BaseConformer baseConformer = mConformerGenerator.getBaseConformer(torsionSet.getConformerIndexes());
		for (int rb=0; rb<mRotatableBond.length; rb++)
			likelyhood *= baseConformer.getTorsionLikelyhood(rb, torsionSet.getTorsionIndexes()[rb]);
		return likelyhood;
		}

	private boolean matchesEliminationRule(TorsionSet ts, double tolerance, ArrayList<TorsionSetEliminationRule> elimRules) {
		for (TorsionSetEliminationRule er:elimRules)
			if (ts.matches(er, tolerance))
				return true;
		return false;
		}

	public double calculateCollisionTolerance() {
		int secondChoiceStartIndex = mMaxTotalCount / SECOND_CHOICE_START_FRACTION;
		return (mTotalCount <= secondChoiceStartIndex) ? 0.0
			: mLowestCollisionStrain
			  + SECOND_CHOICE_MAX_TOLERANCE * (mTotalCount - secondChoiceStartIndex)
											/ (mMaxTotalCount - secondChoiceStartIndex);
		}

	private int[] getRotatableBondsBetween(int f1, int f2) {
		Arrays.fill(mGraphFragmentHandled, false);
		mGraphFragment[0] = f1;
		mGraphFragmentHandled[f1] = true;
		int current = 0;
		int highest = 0;
		while (current <= highest) {
			for (int i=0; i<mConnFragmentNo[mGraphFragment[current]].length; i++) {
				int candidate = mConnFragmentNo[mGraphFragment[current]][i];
				if (candidate == f2) {
					int bondCount = 1;
					int index = current;
					while (index != 0) {
						bondCount++;
						index = mGraphParent[index];
						}
					int[] rotatableBondIndex = new int[bondCount];
					rotatableBondIndex[0] = mConnRotatableBondNo[mGraphFragment[current]][i];
					bondCount = 1;
					while (current != 0) {
						rotatableBondIndex[bondCount++] = mGraphBond[current];
						current = mGraphParent[current];
						}
					return rotatableBondIndex;
					}
				if (!mGraphFragmentHandled[candidate]) {
					mGraphFragmentHandled[candidate] = true;
					highest++;
					mGraphFragment[highest] = candidate;
					mGraphBond[highest] = mConnRotatableBondNo[mGraphFragment[current]][i];
					mGraphParent[highest] = current;
					}
				}
			current++;
			}
		return null;
		}
	}
