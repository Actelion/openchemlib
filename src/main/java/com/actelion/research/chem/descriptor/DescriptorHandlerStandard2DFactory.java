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

public class DescriptorHandlerStandard2DFactory implements DescriptorConstants,DescriptorHandlerFactory {
	protected static DescriptorHandlerStandard2DFactory sFactory;

	public static DescriptorHandlerFactory getFactory() {
		if (sFactory == null) {
			synchronized(DescriptorHandlerStandard2DFactory.class) {
				if (sFactory == null)
					sFactory = new DescriptorHandlerStandard2DFactory();
				}
			}
		return sFactory;
		}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public DescriptorHandler getDefaultDescriptorHandler(String shortName) {

		if (DESCRIPTOR_FFP512.shortName.equals(shortName))
			return DescriptorHandlerLongFFP512.getDefaultInstance();
		if (DESCRIPTOR_PFP512.shortName.equals(shortName))
			return DescriptorHandlerLongPFP512.getDefaultInstance();
//		if (DESCRIPTOR_SSSFP.shortName.equals(shortName))
//			return DescriptorHandlerSSSPathFp.getDefaultInstance();
		if (DESCRIPTOR_ALLFRAG.shortName.equals(shortName))
			return DescriptorHandlerAllFragmentsFP.getDefaultInstance();
		if (DESCRIPTOR_HashedCFp.shortName.equals(shortName))
			return DescriptorHandlerLongCFP.getDefaultInstance();
		if (DESCRIPTOR_SkeletonSpheres.shortName.equals(shortName))
			return DescriptorHandlerSkeletonSpheres.getDefaultInstance();
		if (DESCRIPTOR_OrganicFunctionalGroups.shortName.equals(shortName))
			return DescriptorHandlerFunctionalGroups.getDefaultInstance();

        if (DESCRIPTOR_ReactionFP.shortName.equals(shortName))
            return DescriptorHandlerReactionFP.getDefaultInstance();

		return null;
		}

    @SuppressWarnings({ "unchecked", "rawtypes" })
	public DescriptorHandler create(String shortName) {

		if (DESCRIPTOR_FFP512.shortName.equals(shortName))
			return new DescriptorHandlerLongFFP512();
		if (DESCRIPTOR_PFP512.shortName.equals(shortName))
			return new DescriptorHandlerLongPFP512();
//	    if (DESCRIPTOR_SSSFP.shortName.equals(shortName))
//		    return new DescriptorHandlerSSSPathFp();
	    if (DESCRIPTOR_ALLFRAG.shortName.equals(shortName))
		    return new DescriptorHandlerAllFragmentsFP();
		if (DESCRIPTOR_HashedCFp.shortName.equals(shortName))
			return new DescriptorHandlerLongCFP();
		if (DESCRIPTOR_SkeletonSpheres.shortName.equals(shortName))
			return new DescriptorHandlerSkeletonSpheres();
		if (DESCRIPTOR_OrganicFunctionalGroups.shortName.equals(shortName))
			return new DescriptorHandlerFunctionalGroups();

        if (DESCRIPTOR_ReactionFP.shortName.equals(shortName))
            return new DescriptorHandlerReactionFP();

		return null;
		}
	}
