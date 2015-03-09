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
 * @author Christian Rufener
 */

package com.actelion.research.chem;

import com.actelion.research.chem.reaction.ReactionArrow;

public class DrawingObjectFactory
{
	private DrawingObjectFactory() {}

	public static AbstractDrawingObject createObject(String descriptor) {
		final String START = AbstractDrawingObject.DESCRIPTOR_START+AbstractDrawingObject.DESCRIPTOR_TYPE;
	    if (!descriptor.startsWith(START)
	     || !descriptor.endsWith(AbstractDrawingObject.DESCRIPTOR_END))
	        return null;

	    int index = descriptor.indexOf('\"', START.length());
	    if (index == -1)
	    	return null;

	    final String type = descriptor.substring(START.length(), index);
	    final String detail = descriptor.substring(START.length()+type.length()+1,
				 								   descriptor.length()-AbstractDrawingObject.DESCRIPTOR_END.length());
	    
	    if (type.equals(ReactionArrow.TYPE_STRING))
	        return new ReactionArrow(detail);
	    if (type.equals(TextDrawingObject.TYPE_STRING))
	        return new TextDrawingObject(detail);
	
	    return null;
	    }
	}
