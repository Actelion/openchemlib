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
        int[] descriptor = new ReactionSearcher().createIndex(rxn);
        return (descriptor == null) ? FAILED_OBJECT : descriptor;
    }
    
	public DescriptorHandler<int[], Reaction> getDeepCopy() {
		return new DescriptorHandlerReactionIndex();
	}

}
