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
 * @author Thomas Sander
 */

package com.actelion.research.chem;

import java.io.Serializable;
import java.util.ArrayList;

public class UniqueStringList extends SortedStringList implements Serializable {
    static final long serialVersionUID = 0x20060720;

    private ArrayList<String> mOriginalOrder = new ArrayList<String>();
	private ArrayList<Integer> mIndexList = new ArrayList<Integer>();

    public int getSortedListIndex(String s) {
        return super.getListIndex(s);
        }

	public int getListIndex(String s) {
		int index = super.getListIndex(s);
		if (index == -1)
			return -1;

		return mIndexList.get(index).intValue();
		}

	public int addString(String theString) {
		int index = super.addString(theString);
		if (index == -1)
			return -1;

		int position = mOriginalOrder.size();

		mOriginalOrder.add(theString);
		mIndexList.add(index, new Integer(position));

		return position;
		}

    /**
     * @param i list index within list in original order
     * @return string at position i of list in order of the creation
     */
	public String getStringAt(int i) {
		return mOriginalOrder.get(i);
		}

	/**
	 * @param i list index within sorted list
	 * @return string at position i of sorted list
	 */
    public String getSortedStringAt(int i) {
        return super.getStringAt(i);
        }

	public String[] toArray() {
	    return mOriginalOrder.toArray(new String[0]);
	    }

    public String[] toSortedArray() {
        return super.toArray();
        }
    }
