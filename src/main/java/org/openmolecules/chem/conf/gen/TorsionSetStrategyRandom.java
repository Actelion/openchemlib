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
	 * @param conformerGenerator
	 * @param preferLikelyTorsions if set then more frequent angles are picked with higher probability
	 * @param seed
	 */
	public TorsionSetStrategyRandom(ConformerGenerator conformerGenerator, boolean preferLikelyTorsions, long seed) {
		super(conformerGenerator);
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
		TorsionSet ts;
		int count = 0;

		do {
			if (count++ == MAX_TRIES_FOR_NEW)
				return null;

			if (mPreferLikelyTorsions) {
				double progress = (double)count/(double)MAX_TRIES_FOR_NEW;
				for (int rf=0; rf<mRigidFragment.length; rf++)
					conformerIndex[rf] = mRigidFragment[rf].getLikelyRandomConformerIndex(mRandom.nextDouble(), progress);
				BaseConformer baseConformer = mConformerGenerator.getBaseConformer(conformerIndex);
				for (int rb=0; rb<mRotatableBond.length; rb++)
					torsionIndex[rb] = baseConformer.getLikelyRandomTorsionIndex(rb, mRandom.nextDouble(), progress);
				}
			else {
				for (int j=0; j<mRigidFragment.length; j++)
					conformerIndex[j] = mRandom.nextInt(mRigidFragment[j].getConformerCount());
				for (int j=0; j<mRotatableBond.length; j++)
					torsionIndex[j] = mRandom.nextInt(mRotatableBond[j].getTorsionCount());
				}

			ts = createTorsionSet(torsionIndex, conformerIndex);
			} while (!isNewTorsionSet(ts));

		return ts;
		}
	}
