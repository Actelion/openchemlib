/*
 * Copyright (c) 1997 - 2016
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the the copyright holder nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package com.actelion.research.chem.coords;

/**
 * Created by thomas on 9/23/16.
 */
public class FragmentAssociation {
	private InventorFragment[] mFragment;
	private double[] mX,mY;
	private int[] mCount;

	public FragmentAssociation(InventorFragment f1, InventorFragment f2, int atomIndex1, int atomIndex2) {
		mFragment = new InventorFragment[2];
		mFragment[0] = f1;
		mFragment[1] = f2;
		mX = new double[2];
		mY = new double[2];
		mX[0] = mFragment[0].getAtomX(atomIndex1);
		mY[0] = mFragment[0].getAtomY(atomIndex1);
		mX[1] = mFragment[1].getAtomX(atomIndex2);
		mY[1] = mFragment[1].getAtomY(atomIndex2);
		mCount = new int[2];
		mCount[0] = 1;
		mCount[1] = 1;
		}

	/**
	 * Uses center of gravity of the fragments as anchor points
	 * @param f1
	 * @param f2
	 */
	public FragmentAssociation(InventorFragment f1, InventorFragment f2) {
		mFragment = new InventorFragment[2];
		mFragment[0] = f1;
		mFragment[1] = f2;
		mX = new double[2];
		mY = new double[2];
		mCount = new int[2];
		for (int i=0; i<2; i++) {
			for (int j=0; j<mFragment[i].size(); j++) {
				mX[i] += mFragment[i].getAtomX(j);
				mY[i] += mFragment[i].getAtomY(j);
				}
			mCount[i] = mFragment[i].size();
			}
		}

	public void add(int atomIndex1, int atomIndex2) {
		mX[0] += mFragment[0].getAtomX(atomIndex1);
		mY[0] += mFragment[0].getAtomY(atomIndex1);
		mX[1] += mFragment[1].getAtomX(atomIndex2);
		mY[1] += mFragment[1].getAtomY(atomIndex2);
		mCount[0]++;
		mCount[1]++;
		}

	public InventorFragment getFragment(int i) {
		return mFragment[i];
		}

	public int getPriority() {
		return mFragment[0].size() * mFragment[1].size();
		}

	public void arrange(double minDistance, boolean keepFirstFragment) {
		double[] angle = new double[2];
		for (int i=0; i<2; i++) {
			mX[i] /= mCount[i];
			mY[i] /= mCount[i];
			angle[i] = mFragment[i].calculatePreferredAttachmentAngle(mX[i], mY[i], mFragment[1-i].size(), minDistance);
			}

		mFragment[0].rotate(mX[0], mY[0], Math.PI/2.0 - angle[0]);
		mFragment[1].rotate(mX[1], mY[1], Math.PI*3.0/2.0 - angle[1]);

		double yMin = Double.MAX_VALUE;
		double yMax = -Double.MAX_VALUE;
		double dy = mY[0] - mY[1];
		for (int i=0; i<mFragment[1].mAtomY.length; i++) {
			mFragment[1].mAtomY[i] += dy;
			if (yMin > mFragment[1].mAtomY[i])
				yMin = mFragment[1].mAtomY[i];
			if (yMax < mFragment[1].mAtomY[i])
				yMax = mFragment[1].mAtomY[i];
			}

		double range = yMax - yMin + 2*minDistance;
		int binCount = (int)Math.ceil(range);
		yMin += (range - binCount) / 2 - minDistance;
		double[] leftX = new double[binCount];
		for (int i=0; i<binCount; i++)
			leftX[i] = mX[1] + minDistance;
		for (int i=0; i<mFragment[1].mAtomY.length; i++) {
			double relY = mFragment[1].mAtomY[i] - yMin;
			int low = (int)(relY - minDistance);
			int high = Math.min((int)(relY + minDistance), binCount-1);
			for (int j=low; j<=high; j++) {
				if (leftX[j] > mFragment[1].mAtomX[i])
					leftX[j] = mFragment[1].mAtomX[i];
				}
			}

		for (int i=0; i<binCount; i++)
			leftX[i] -= minDistance;

		double dx = mX[0] - mX[1];
		for (int i=0; i<mFragment[0].mAtomX.length; i++) {
			int index = (int)(mFragment[0].mAtomY[i] - yMin);
			if (index >= 0 && index < leftX.length
			 && dx < mFragment[0].mAtomX[i] - leftX[index])
				dx = mFragment[0].mAtomX[i] - leftX[index];
			}

		for (int i=0; i<mFragment[1].mAtomX.length; i++)
			mFragment[1].mAtomX[i] += dx;

		if (keepFirstFragment) {
			mFragment[0].rotate(mX[0], mY[0], angle[0] - Math.PI / 2.0);
			mFragment[1].rotate(mX[0], mY[0], angle[0] - Math.PI / 2.0);
			}
		}
	}
