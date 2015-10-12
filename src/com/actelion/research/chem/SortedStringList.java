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
