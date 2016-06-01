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

import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionSearcher;

public class DescriptorHandlerReactionIndex extends AbstractDescriptorHandlerFP<Reaction> {
    private static DescriptorHandlerReactionIndex sDefaultInstance;

    public static DescriptorHandlerReactionIndex getDefaultInstance() {
        if (sDefaultInstance == null) {
        	synchronized(DescriptorHandlerReactionIndex.class) {
        		sDefaultInstance = new DescriptorHandlerReactionIndex();
        	}
        }
        return sDefaultInstance;
    }

    public DescriptorInfo getInfo() {
        return DescriptorConstants.DESCRIPTOR_ReactionIndex;
    }

    public String getVersion() {
        return ReactionSearcher.cIndexVersion;
    }

    public int[] createDescriptor(Reaction rxn) {
	    if (rxn ==null)
		    return null;

        int[] descriptor = new ReactionSearcher().createIndex(rxn);
        return (descriptor == null) ? FAILED_OBJECT : descriptor;
    }
    
	public DescriptorHandler<int[], Reaction> getThreadSafeCopy() {
		return this;
	}

}
