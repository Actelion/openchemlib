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

import com.actelion.research.chem.*;

public class DescriptorHandlerFFP512 extends AbstractDescriptorHandlerFP<StereoMolecule> {
	public static final String VERSION = SSSearcherWithIndex.cIndexVersion;
	private static DescriptorHandlerFFP512 sDefaultInstance;
	private static final int sIntCount = (SSSearcherWithIndex.getNoOfKeys() + 31) / 32;

	public static DescriptorHandlerFFP512 getDefaultInstance() {
		synchronized(DescriptorHandlerFFP512.class) {
			if (sDefaultInstance == null)
				sDefaultInstance = new DescriptorHandlerFFP512();
			}
		return sDefaultInstance;
		}

	public DescriptorInfo getInfo() {
		return DescriptorConstants.DESCRIPTOR_FFP512;
		}

	public String getVersion() {
		return VERSION;
		}

	public int[] decode(String s) {
		int[] descriptor = (s != null && s.length() == 128) ?
			SSSearcherWithIndex.getIndexFromHexString(s) : super.decode(s);
		return (descriptor != null && descriptor.length == sIntCount) ? descriptor : null;
		}

	public int[] decode(byte[] bytes) {
		int[] descriptor = (bytes != null && bytes.length == 128) ?
			SSSearcherWithIndex.getIndexFromHexString(bytes) : super.decode(bytes);
		return (descriptor != null && descriptor.length == sIntCount) ? descriptor : null;
		}

	public int[] createDescriptor(StereoMolecule mol) {
		int[] descriptor = new SSSearcherWithIndex().createIndex(mol);
		return (descriptor == null) ? FAILED_OBJECT : descriptor;
		}
	
	
	public DescriptorHandler<int[], StereoMolecule> getDeepCopy() {
		return new DescriptorHandlerFFP512();
		}
	}
