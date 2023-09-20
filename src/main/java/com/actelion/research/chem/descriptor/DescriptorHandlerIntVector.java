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

package com.actelion.research.chem.descriptor;

import java.nio.charset.StandardCharsets;

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
        return calculationFailed(d) ? FAILED_STRING : new String(new DescriptorEncoder().encodeIntArray(d), StandardCharsets.UTF_8);
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
