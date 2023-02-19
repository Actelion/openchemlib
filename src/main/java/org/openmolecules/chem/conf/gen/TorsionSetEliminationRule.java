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

public class TorsionSetEliminationRule {
	private long[] mMask;
	private long[] mData;
	private double mCollisionIntensity;

	public TorsionSetEliminationRule(int[] torsionIndex, int[] rotatableBondIndex, double collisionIntensity, TorsionSetEncoder encoder) {
		mMask = new long[encoder.getEncodingLongCount()];
		mData = new long[encoder.getEncodingLongCount()];
		encoder.encodeRule(torsionIndex, rotatableBondIndex, mMask, mData);
		mCollisionIntensity = collisionIntensity;
		}

	public long[] getMask() {
		return mMask;
		}

	public long[] getData() {
		return mData;
		}

	public double getCollisionIntensity() {
		return mCollisionIntensity;
		}

	/**
	 * Checks whether mask and data as elimination rule are covered by this rule.
	 * @param rule
	 * @return true if mask & data are already covered and don't need to be considered
	 */
	public boolean isCovered(TorsionSetEliminationRule rule) {
		boolean isCovered = true;
		for (int i=0; i<mMask.length; i++)
			if ((~rule.mMask[i] & mMask[i]) != 0l || (rule.mData[i] & mMask[i]) != mData[i])
				isCovered = false;

		return isCovered;
		}

	/**
	 * Checks whether mask and data constitute a more general rule than this
	 * and therefore include this.
	 * @param rule
	 * @return true if mask & data are more general than this
	 */
	public boolean isMoreGeneral(TorsionSetEliminationRule rule) {
		boolean isMoreGeneral = true;
		for (int i=0; i<mMask.length; i++)
			if ((~mMask[i] & rule.mMask[i]) != 0l || (mData[i] & rule.mMask[i]) != rule.mData[i])
				 isMoreGeneral = false;

		return isMoreGeneral;
		}
	}
