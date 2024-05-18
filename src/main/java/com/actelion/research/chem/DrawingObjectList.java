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
			if (l != null)
				for (AbstractDrawingObject abstractDrawingObject : l)
					add(abstractDrawingObject.clone());
		} catch (Exception e) {
		} finally {

		}
	}

	public DrawingObjectList(String objectString) {
		super();
		if (objectString == null || objectString.isEmpty())
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
		StringBuilder objectString = new StringBuilder();
		for (AbstractDrawingObject abstractDrawingObject : this)
			objectString.append(abstractDrawingObject.getDescriptor()).append("\n");
		return objectString.toString();
		}
	}
