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

import com.actelion.research.chem.Coordinates;
import org.openmolecules.chem.conf.so.SelfOrganizedConformer;

public class RigidFragment {
	int   mCoreAtomCount;
	int[] mCoreToFragmentAtom;
	int[] mFragmentToOriginalAtom;
	int[] mExtendedToFragmentAtom;
	int[] mOriginalToExtendedAtom;

	SelfOrganizedConformer[] mConformerList;
	double[] mConformerLikelihood;

	public RigidFragment(int coreAtomCount,
	                     int[] coreToFragmentAtom,
	                     int[] fragmentToOriginalAtom,
	                     int[] extendedToFragmentAtom,
	                     int[] originalToExtendedAtom,
	                     SelfOrganizedConformer[] conformerList,
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

	public SelfOrganizedConformer getConformer(int i) {
		return mConformerList[i];
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
