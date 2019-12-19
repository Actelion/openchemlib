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

public class DescriptorHelper implements DescriptorConstants {
	
	public static final String TAG_SIMILARITY = "Similarity";
	
    public static int getDescriptorType(String shortName) {
        DescriptorInfo descriptorInfo = getDescriptorInfo(unifyShortName(shortName));
        return (descriptorInfo == null) ? DESCRIPTOR_TYPE_UNKNOWN
                                        : descriptorInfo.type;
    }

    public static DescriptorInfo getDescriptorInfo(String shortName) {
        for (int i=0; i<DESCRIPTOR_EXTENDED_LIST.length; i++)
            if (DESCRIPTOR_EXTENDED_LIST[i].shortName.equals(unifyShortName(shortName)))
                return DESCRIPTOR_EXTENDED_LIST[i];
        return null;
    }

    public static boolean isBinaryFingerprint(String shortName) {
        DescriptorInfo descriptorInfo = getDescriptorInfo(unifyShortName(shortName));
        return (descriptorInfo == null) ? false
                                        : descriptorInfo.isBinary;
    }

    public static boolean isDescriptorShortName(String shortName) {
        for (int i=0; i<DESCRIPTOR_EXTENDED_LIST.length; i++)
            if (DESCRIPTOR_EXTENDED_LIST[i].shortName.equals(unifyShortName(shortName)))
                return true;
        return false;
    }

    public static String shortNameToName(String shortName) {
		for (int i=0; i<DESCRIPTOR_EXTENDED_LIST.length; i++)
			if (DESCRIPTOR_EXTENDED_LIST[i].shortName.equals(shortName))
				return DESCRIPTOR_EXTENDED_LIST[i].name;
		return null;
	}

	public static String nameToShortName(String name) {
		for (int i=0; i<DESCRIPTOR_EXTENDED_LIST.length; i++)
			if (DESCRIPTOR_EXTENDED_LIST[i].name.equals(name))
				return DESCRIPTOR_EXTENDED_LIST[i].shortName;
		return null;
	}

	private static String unifyShortName(String shortname) {
        return "PP3DMM2".equals(shortname) ? "Flexophore" : shortname;
    }

	/**
	 * Creates a header tag name from the descriptor short name.
	 * The tag is used to store virtual screening scores. 
	 * @param dh
	 * @return
	 */
	public static String getTagDescriptorSimilarity(ISimilarityCalculator<?> dh) {
		return TAG_SIMILARITY + dh.getInfo().shortName;
		}

	public static String getTagDescriptorSimilarity(SimilarityCalculatorInfo info) {
		return TAG_SIMILARITY + info.shortName;
		}

	public static String getTagDescriptorSimilarity(DescriptorInfo dh){
		return TAG_SIMILARITY + dh.shortName;
		}
	}
