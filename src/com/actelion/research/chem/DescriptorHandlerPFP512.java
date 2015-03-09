/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
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

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.FingerPrintGenerator;

public class DescriptorHandlerPFP512 extends AbstractDescriptorHandlerFP<StereoMolecule>
		implements DescriptorConstants {
    private static final double CORRECTION_FACTOR = 0.8;

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
		return "1.0";
	}

	public int[] createDescriptor(StereoMolecule mol) {
		java.util.BitSet bitset = new FingerPrintGenerator().getFingerprint(mol);

		if (bitset == null)
			return FAILED_OBJECT;

		int[] fp = new int[16];
		int mask = 1;
		for (int i = 0; i < 32; i++) {
			for (int j = 0; j < 16; j++)
				if (bitset.get(32 * i + j))
					fp[j] += mask;
			mask <<= 1;
		}

		return fp;
	}

	public DescriptorHandler<int[], StereoMolecule> getDeepCopy() {
		return new DescriptorHandlerPFP512();
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
