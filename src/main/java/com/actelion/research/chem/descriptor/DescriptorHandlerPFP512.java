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
*/

package com.actelion.research.chem.descriptor;

import com.actelion.research.chem.StereoMolecule;

public class DescriptorHandlerPFP512 extends AbstractDescriptorHandlerFP<StereoMolecule> implements DescriptorConstants {
    private static final double CORRECTION_FACTOR = 0.85;

    private static DescriptorHandlerPFP512 sDefaultInstance;

	public static DescriptorHandlerPFP512 getDefaultInstance() {
		synchronized (DescriptorHandlerPFP512.class) {
			if (sDefaultInstance == null) {
				sDefaultInstance = new DescriptorHandlerPFP512();
			}
		}
		return sDefaultInstance;
	}

	public DescriptorInfo getInfo() {
		return DESCRIPTOR_PFP512;
	}

	public String getVersion() {
		return DescriptorConstants.DESCRIPTOR_PFP512.version;
	}

	public int[] createDescriptor(StereoMolecule mol) {
		if (mol ==null)
			return null;

		java.util.BitSet bitset = new FingerPrintGenerator().getFingerprint(mol);

		if (bitset == null)
			return FAILED_OBJECT;

		int[] fp = new int[16];
		int mask = 1;
		for (int i = 0; i < 32; i++) {
			for (int j = 0; j < 16; j++)
				if (bitset.get(16 * i + j))
					fp[j] += mask;
			mask <<= 1;
		}

		return fp;
	}

	public DescriptorHandler<int[], StereoMolecule> getThreadSafeCopy() {
		return this;
	}

	@Override
    public float getSimilarity(int[] o1, int[] o2) {
		return normalizeValue(super.getSimilarity(o1, o2));
    }

	private float normalizeValue(double value) {
		return value <= 0.0f ? 0.0f
			 : value >= 1.0f ? 1.0f
			 : (float)(1.0-Math.pow(1-Math.pow(value, CORRECTION_FACTOR) ,1.0/CORRECTION_FACTOR));
	}
}
