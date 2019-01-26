/*
 * @(#)TorsionSetEliminationRule.java
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
