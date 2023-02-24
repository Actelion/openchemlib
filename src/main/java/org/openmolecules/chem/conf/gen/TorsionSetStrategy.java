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

	// The ConformerGenerator repeatedly requests new TorsionSets from the TorsionSetStrategy. It checks any
	// returned TorsionSet for atom collisions: If a TorsionSet's collision strain is smaller than the lowest
	// previously seen strain plus SECOND_CHOICE_TOLERANCE, then a conformer is constructed from the TorsionSet.
	// Independent of whether a conformer was built, the ConformerGenerator returns the previous torsion set
	// when requesting a new one.
	// If a returned TorsionSet was not used to build a conformer, but has a strain lower than MAX_COLLISION_INTENSITY,
	// then it is cached by the TorsionSetStrategy as a second choice for later.
	// TorsionSets with a strain higher than SECOND_CHOICE_TOLERANCE are analyzed by the TorsionSetStrategy for
	// collisions and elimination rules are derived for the future generation of new torsion sets.
	// When no more TorsionSet can be created, the second choice TorsionSets are returned with increasing
	// collision strain. Of those TorsionSets with strains higher than MAX_COLLISION_INTENSITY, the best is kept.
	// If no acceptable TorsionSets could be created then the ConformerGenerator may choose to use the best unacceptable one.

	public static final double COLLISION_STRAIN_FOR_FACTOR_10 = 1.36;   // Collision strains are supposed to be on a kcal/mol scale
	public static final double MAX_COLLISION_STRAIN = 6;
	private static final double MAX_LOWEST_COLLISION_STRAIN = 3;
	private static final double SECOND_CHOICE_TOLERANCE = 3;

	protected ConformerGenerator mConformerGenerator;
	protected RotatableBond[] mRotatableBond;
	protected RigidFragment[] mRigidFragment;
	private TorsionSetEncoder mTorsionSetEncoder;
	private int mFragmentCount,mCollisionCount,mMaxTotalCount,mPermutationCount,mSuccessCount;
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
	private SortedList<TorsionSetWithEliminationRuleStrain> mStrainedTorsionSetCache;
	private TorsionSet mBestUnacceptableTorsionSet;

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
		mSecondChoiceList = new SortedList<>(Comparator.comparingDouble(ts -> ((TorsionSet)ts).getCollisionStrainSum()));
		mStrainedTorsionSetCache = new SortedList<>();
		mCollisionCount = 0;
		mSuccessCount = 0;
		mLowestCollisionStrain = MAX_LOWEST_COLLISION_STRAIN;
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
		return mTorsionSetList.size();
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
		TorsionSet ts = new TorsionSet(torsionIndex, conformerIndex, mTorsionSetEncoder);
		ts.setContribution(getContribution(ts));
		return ts;
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
		if (previousTorsionSet != null && previousTorsionSet.isUsed())
			mSuccessCount++;

		// Some molecules have unavoidable internal strains,
		// which we try to determine until we start returning second choices.
		if (previousTorsionSet != null
		 && mLowestCollisionStrain > previousTorsionSet.getCollisionStrainSum())
			mLowestCollisionStrain = previousTorsionSet.getCollisionStrainSum();

		// If we switched to second choices then return second choices until no one is left.
		if (mUsingSecondChoices)
			return getBestSecondChoice();

		if (mTorsionSetList.size() == mMaxTotalCount) {
			if (diagnostics != null)
				diagnostics.setExitReason("maxTotal(" + mMaxTotalCount + ") reached A; collisions:" + mCollisionCount);

			return null;
			}

		// - Learn from the collisions of the previous torsion set, if the strain is high enough to potentially be a problem later.
		// - If collision is tolerable, keep the torsion set among second choices for later
		if (previousTorsionSet != null
		 && previousTorsionSet.getCollisionStrainSum() > SECOND_CHOICE_TOLERANCE) {

			BaseConformer baseConformer = (diagnostics != null) ? mConformerGenerator.getBaseConformer(previousTorsionSet.getConformerIndexes()) : null;
			ArrayList<TorsionSetEliminationRule> eliminationRules = (diagnostics != null) ? baseConformer.getEliminationRules() : null;
			int elimRuleCount = (diagnostics == null) ? 0 : eliminationRules.size();

			processCollisions(previousTorsionSet);

			if (diagnostics != null)
				for (int i=elimRuleCount; i<eliminationRules.size(); i++)
					diagnostics.get(previousTorsionSet).addEliminationRule(mTorsionSetEncoder.createRuleString(eliminationRules.get(i), baseConformer));

			if (!previousTorsionSet.isUsed()) {
				if (previousTorsionSet.getCollisionStrainSum() <MAX_COLLISION_STRAIN)
					mSecondChoiceList.add(previousTorsionSet);
				else if (mBestUnacceptableTorsionSet == null
					  || mBestUnacceptableTorsionSet.getCollisionStrainSum() > previousTorsionSet.getCollisionStrainSum())
					mBestUnacceptableTorsionSet = previousTorsionSet;
				}

			mCollisionCount++;
			}

		double tolerance = calculateCollisionTolerance();

		// Calculate a tolerable collision threshold that increases with the number of returned torsion sets.
		// And return a tolerable torsion set, if we have one in the list.
		// Note: Second choice torsion sets have valid collision sums
/* Current calculate collision tolerance just get smaller over time. Thus, it is impossible that previous second choice
   sets become now first choice. If we introduce an increase of the tolerance over time, we should uncomment this.
		TorsionSet oldTS = getBestSecondChoiceIfBelow(tolerance);
		if (oldTS != null)
			return oldTS;
*/

		TorsionSet ts = getTorsionSet(tolerance, previousTorsionSet);

		while (ts != null) {
			// With the current knowledge of elimination rule, check if we can rule out the torsion set because of unacceptable strain
			double collisionStrain = calculateKnownCollisionStrainSum(ts, mConformerGenerator.getBaseConformer(ts.getConformerIndexes()).getEliminationRules());
			if (collisionStrain < tolerance)
				break;

			if (collisionStrain <MAX_COLLISION_STRAIN)
				mStrainedTorsionSetCache.add(new TorsionSetWithEliminationRuleStrain(ts, collisionStrain));

			mCollisionCount++;
			mTorsionSetList.add(ts);

			if (mTorsionSetList.size() == mMaxTotalCount) {
				if (diagnostics != null)
					diagnostics.setExitReason("maxTotal(\"+mMaxTotalCount+\") reached B; collisions:"+mCollisionCount);
				return null;
				}

			tolerance = calculateCollisionTolerance();

/* Current calculate collision tolerance just get smaller over time. Thus, it is impossible that previous second choice
   sets become now first choice. If we introduce an increase of the tolerance over time, we should uncomment this.
			TorsionSet oldTS = getBestSecondChoiceIfBelow(tolerance);
			if (oldTS != null)
				return oldTS;
*/

			ts = getTorsionSet(tolerance, ts);
			}

		if (ts == null) {
			if (mSecondChoiceList.size() != 0) {
				mUsingSecondChoices = true;

				// remove all second choices beyond smallest found collision intensity plus tolerance
				while (mSecondChoiceList.size() != 0
					&& mSecondChoiceList.get(mSecondChoiceList.size()-1).getCollisionStrainSum() > mLowestCollisionStrain + SECOND_CHOICE_TOLERANCE)
					mSecondChoiceList.remove(mSecondChoiceList.size()-1);

				ts = getBestSecondChoice();
				}

			if (ts == null && diagnostics != null)
				diagnostics.setExitReason("Inner strategy stop criterion hit");

			return ts;
			}

		mTorsionSetList.add(ts);

		return ts;
		}

	/**
	 * If we have an earlier created unused torsion set in the cache, return that.
	 * Otherwise creates a new torsion set using the strategy implementation.
	 * @param tolerance
	 * @param previousTorsionSet
	 * @return cached or new torsion set
	 */
	private TorsionSet getTorsionSet(double tolerance, TorsionSet previousTorsionSet) {
		// If we have an earlier generated torsion set elimination rules at time of creation
		// didn't indicate a higher strain than tolerance, the return that...
		if (mStrainedTorsionSetCache.size() != 0) {
			TorsionSetWithEliminationRuleStrain tswers = mStrainedTorsionSetCache.get(0);
			if (tswers.getRuleStrain() < tolerance) {
				mStrainedTorsionSetCache.remove(0);
				return tswers.getTorsionSet();
				}
			}

		return createTorsionSet(previousTorsionSet);
		}

	private TorsionSet getBestSecondChoiceIfBelow(double tolerance) {
		if (mSecondChoiceList.size() == 0 || tolerance == 0.0)
			return null;

		if (mSecondChoiceList.get(0).getCollisionStrainSum() > tolerance)
			return null;

		TorsionSet ts = mSecondChoiceList.get(0);
		mSecondChoiceList.remove(0);
		return ts;
		}

	private TorsionSet getBestSecondChoice() {
		if (mSecondChoiceList.size() == 0)
			return null;

		int index = (this instanceof TorsionSetStrategyRandom) ?
					((TorsionSetStrategyRandom)this).getRandom().nextInt(mSecondChoiceList.size()) : 0;

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
		return mBestUnacceptableTorsionSet;
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
		double[][] collisionIntensityMatrix = torsionSet.getCollisionStrainMatrix();
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
			if (collidingTorsionSet.getCollisionStrainIfMatches(er) != 0)
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

	private double calculateKnownCollisionStrainSum(TorsionSet ts, ArrayList<TorsionSetEliminationRule> elimRules) {
		double sum = 0.0;
		for (TorsionSetEliminationRule er:elimRules)
			sum += ts.getCollisionStrainIfMatches(er);

		return sum;
		}

	public double calculateCollisionTolerance() {
		return mLowestCollisionStrain + SECOND_CHOICE_TOLERANCE;
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
