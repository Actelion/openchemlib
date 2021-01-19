/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.chem.descriptor;

/**
 * This is a descriptor handler, where the input object is an integer array
 * that typically represents counts of some sort. This class may be used
 * if input objects are more complex and the descriptors derived from them
 * can be represented as integer vectors.<br>
 * This class provides similarity calculation and descriptor en- and decoding.
 */

public class DescriptorHandlerIntVector<U extends Object> implements DescriptorHandler<int[], U> {
    protected static final int[] FAILED_OBJECT = new int[0];
    private static final double CORRECTION_FACTOR = 0.42;

    private double correctionFactor;

	public DescriptorHandlerIntVector() {
		correctionFactor = CORRECTION_FACTOR;
	}

	public void setCorrectionFactor(double correctionFactor) {
		this.correctionFactor = correctionFactor;
	}

	@Override
	public float getSimilarity(int[] d1, int[] d2) {
        if (d1 == null || d2 == null)
            return Float.NaN;

        int total = 0;
        int matching = 0;
        for (int i=0; i<d1.length; i++) {
            total += Math.max(d1[i], d2[i]);
            matching += Math.min(d1[i], d2[i]);
            }

        return normalizeValue((float)matching/(float)total);
		}

	private float normalizeValue(double value) {
		return value <= 0.0f ? 0.0f
			 : value >= 1.0f ? 1.0f
			 : (float)(1.0-Math.pow(1-Math.pow(value, correctionFactor) ,1.0/correctionFactor));
		}

	@Override
	public DescriptorInfo getInfo() {
        return DescriptorConstants.DESCRIPTOR_IntegerVector;
		}

	@Override
	public String getVersion() {
		return DescriptorConstants.DESCRIPTOR_IntegerVector.version;
		}

	@Override
	public String encode(int[] d) {
        return calculationFailed(d) ? FAILED_STRING : new String(new DescriptorEncoder().encodeIntArray(d));
		}

	@Override
	public int[] decode(String s) {
		return s == null ? null : new DescriptorEncoder().decodeIntArray(s);
		}

	@Override
	public int[] decode(byte[] bytes) {
		return bytes == null ? null : new DescriptorEncoder().decodeIntArray(bytes);
		}

	@Override
	public int[] createDescriptor(U o) {
		return (int[])o;	// the object is the descriptor
		}

	@Override
	public boolean calculationFailed(int[] d) {
        return d==null || d.length==0;
		}

	@Override
	public DescriptorHandler<int[], U> getThreadSafeCopy() {
		return this;
		}
	}
