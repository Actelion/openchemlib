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

import java.util.Arrays;

/**
 * TorsionSetStrategy that systematically creates all possible TorsionSets in
 * batches while keeping a focus on the likelyhood of individual torsions.
 * @author Thomas Sander
 */
public class TorsionSetStrategyLikelySystematic extends TorsionSetStrategy {
	private int[]			mCurrentMaxTorsionIndex,mCurrentMaxConformerIndex;
	private TorsionSet[]	mAvailableTorsionSet;
	private int				mAvailableTorsionSetIndex;

	public TorsionSetStrategyLikelySystematic(ConformerGenerator conformerGenerator) {
		super(conformerGenerator);
		mCurrentMaxTorsionIndex = new int[mRotatableBond.length];
		mCurrentMaxConformerIndex = new int[mRigidFragment.length];
		mAvailableTorsionSet = new TorsionSet[1];
		mAvailableTorsionSet[0] = createTorsionSet(new int[mRotatableBond.length], new int[mRigidFragment.length]);
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
		BaseConformer baseConformer = mConformerGenerator.getBaseConformer(mCurrentMaxConformerIndex);
		for (int rb=0; rb<mCurrentMaxTorsionIndex.length; rb++) {
			if (mCurrentMaxTorsionIndex[rb] < mRotatableBond[rb].getTorsionCount()-1) {
				double loss = baseConformer.getTorsionLikelyhood(rb, mCurrentMaxTorsionIndex[rb])
							/ baseConformer.getTorsionLikelyhood(rb, mCurrentMaxTorsionIndex[rb]+1);
				if (minLoss > loss) {
					minLoss = loss;
					bestIndex = rb;
					}
				}
			}
		for (int rf=0; rf<mCurrentMaxConformerIndex.length; rf++) {
			if (mCurrentMaxConformerIndex[rf] < mRigidFragment[rf].getConformerCount()-1) {
				double loss = mRigidFragment[rf].getConformerLikelihood(mCurrentMaxConformerIndex[rf])
						    / mRigidFragment[rf].getConformerLikelihood(mCurrentMaxConformerIndex[rf]+1);
				if (minLoss > loss) {
					minLoss = loss;
					bestIndex = rf;
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

		Arrays.sort(mAvailableTorsionSet, (ts1, ts2) -> {
			return Double.compare(ts2.getContribution(), ts1.getContribution());    // we want the highest likelyhood first
			} );
		mAvailableTorsionSetIndex = -1;
		}
	}
