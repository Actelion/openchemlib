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

package com.actelion.research.chem;

import java.util.ArrayList;

public class DrawingObjectList extends ArrayList<AbstractDrawingObject> {
    static final long serialVersionUID = 0x20060724;

    public DrawingObjectList() {
		super();
		}

	public DrawingObjectList(DrawingObjectList l) {
		super();
		try {
			if (l != null) {
				for (int i = 0; i < l.size(); i++)
					add((AbstractDrawingObject) l.get(i).clone());
			}
		} catch (Exception e) {
		} finally {

		}
	}

	public DrawingObjectList(String objectString) {
		super();
		if (objectString == null || objectString.length() == 0)
			return;

		int index1 = 0;
		int index2 = objectString.indexOf('\n');
		while (index2 != -1) {
            AbstractDrawingObject o = AbstractDrawingObject.instantiate(objectString.substring(index1, index2));
			if (o != null)
				add(o);
			index1 = index2+1;
			index2 = objectString.indexOf('\n', index1);
			}
		}

	public String toString() {
		StringBuffer objectString = new StringBuffer();
		for (int i=0; i<size(); i++)
			objectString.append(get(i).getDescriptor()+"\n");
		return objectString.toString();
		}
	}
