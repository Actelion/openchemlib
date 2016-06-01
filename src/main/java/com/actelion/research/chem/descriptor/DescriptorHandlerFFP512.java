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
	
	public DescriptorHandler<int[], StereoMolecule> getThreadSafeCopy() {
		return this;
		}
	}
