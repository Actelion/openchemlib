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

	public TorsionSetEliminationRule(long[] mask, long[] data, double collisionIntensity) {
		mMask = mask;
		mData = data;
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
	 * @param mask
	 * @param data
	 * @return true if mask & data are already covered and don't need to be considered
	 */
	public boolean isCovered(long[] mask, long[] data) {
		boolean isCovered = true;
		for (int i=0; i<mMask.length; i++)
			if ((~mask[i] & mMask[i]) != 0l || (data[i] & mMask[i]) != mData[i])
				isCovered = false;

		return isCovered;
		}

	/**
	 * Checks whether mask and data constitute a more general rule than this
	 * and therefore include this.
	 * @param mask
	 * @param data
	 * @return true if mask & data are more general than this
	 */
	public boolean isMoreGeneral(long[] mask, long[] data) {
		boolean isMoreGeneral = true;
		for (int i=0; i<mMask.length; i++)
			if ((~mMask[i] & mask[i]) != 0l || (mData[i] & mask[i]) != data[i])
				 isMoreGeneral = false;

		return isMoreGeneral;
		}
	}
