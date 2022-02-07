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


public class SortedStringList implements Serializable {
    static final long serialVersionUID = 0x20060720;

    private ArrayList<String> mList = new ArrayList<String>();

	public boolean contains(String theString) {
		return getListIndex(theString) != -1;
		}

	public boolean equals(SortedStringList s) {
		if (mList.size() != s.mList.size())
			return false;

		for (int i=0; i<mList.size(); i++)
			if (!mList.get(i).equals(s.mList.get(i)))
				return false;

		return true;
		}

	public int getListIndex(String theString) {
		int vectorSize = mList.size();

		if (vectorSize == 0)
			return -1;

		int index = 1;
		while (2 * index <= vectorSize)
			index <<= 1;
		int increment = index;
		index--;

		while (increment != 0) {
			increment >>= 1;
			if (index >= vectorSize) {
				index -= increment;
				continue;
				}

			int comparison = theString.compareTo(mList.get(index));
			if (comparison == 0)
				return index;

			if (increment == 0)
				break;

			if (comparison < 0) {
				index -= increment;
				}
			else {
				index += increment;
				}
			}

		return -1;
		}

	public int addString(String theString) {
		int vectorSize = mList.size();

		if (vectorSize == 0) {
			mList.add(0, theString);
			return 0;
			}

		int index = 1;
		while (2 * index <= vectorSize)
			index <<= 1;
		int increment = index;
		index--;

		while (increment != 0) {
			increment >>= 1;
			if (index >= vectorSize) {
				index -= increment;
				continue;
				}

			int comparison = theString.compareTo(mList.get(index));
			if (comparison == 0)
				return -1;

			if (increment == 0)
				break;

			if (comparison < 0) {
				index -= increment;
				}
			else {
				index += increment;
				}
			}

		if ((index < vectorSize)
		 && (theString.compareTo(mList.get(index)) > 0))
			index++;

		mList.add(index,theString);
		return index;
		}

	public int getSize(){
		return mList.size();
		}

	public String getStringAt(int i) {
		return mList.get(i);
		}

	public String[] toArray() {
	    return mList.toArray(new String[0]);
	    }

	public void removeAllStrings() {
    	mList.clear();
		}
	}
