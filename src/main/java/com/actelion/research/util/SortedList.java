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

package com.actelion.research.util;

import java.util.ArrayList;
import java.util.Comparator;


public class SortedList<T> {
    static final long serialVersionUID = 0x20160626;

    private ArrayList<T> mList = new ArrayList<T>();
	private Comparator mComparator;

	public SortedList() {
		this(null);
		}

	public SortedList(Comparator comparator) {
		mComparator = comparator;
		}

	public boolean contains(T object) {
		return getIndex(object) != -1;
		}

	public boolean equals(SortedList<T> s) {
		if (mList.size() != s.mList.size())
			return false;

		if (mComparator != null) {
			for (int i=0; i<mList.size(); i++)
				if (mComparator.compare(mList.get(i), s.mList.get(i)) != 0)
					return false;
			}
		else {
			for (int i=0; i<mList.size(); i++)
				if (!mList.get(i).equals(s.mList.get(i)))
					return false;
			}

		return true;
		}

	/**
	 * Returns the position index of object in the sorted list.
	 * If object is not in the list, -1 is returned.
	 * @param object
	 * @return
	 */
	public int getIndex(T object) {
		int index = getIndexOrInsertIndex(object);
		return (index < 0) ? -1 : index;
		}

	/**
	 * Returns the position index of object in the sorted list.
	 * If object is not in the list, -(insertIndex+1) is returned
	 * with insertIndex being that position where object would need
	 * to be inserted to keep a correct sort order.
	 * @param object
	 * @return
	 */
	public int getIndexOrInsertIndex(T object) {
		int vectorSize = mList.size();

		if (vectorSize == 0) {
			return -1;
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

			int comparison = compare(object, mList.get(index));
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

		if ((index < vectorSize)
	     && (compare(object, mList.get(index)) > 0))
			index++;

		return -(index+1);
		}

	private int compare(T o1, T o2) {
		return (mComparator != null) ? mComparator.compare(o1, o2) : ((Comparable<? super T>)o1).compareTo((T)o2);
		}

	/**
	 * If object is not member of this list, returns potential insert index.
	 * If object is member of this list, return list index of object.
	 * @param object
	 * @return list index or insert index
	 */
	public int getIndexBelowEqual(T object) {
       int index = getIndexOrInsertIndex(object);
       return (index < 0) ? -(index+1) : index;
       }

    /**
     * If object is not member of this list, returns potential insert index.
     * If object is member of this list, return list index of object incremented by one.
     * @param object
     * @return list index or insert index
     */
   public int getIndexAboveEqual(T object) {
       int index = getIndexOrInsertIndex(object);
       return (index < 0) ? -(index+1) : index+1;
       }

   /**
    * Adds object to the list provided that it doesn't contain
    * an object being considered equal by compareTo().
    * @param object
    * @return object's list index, no matter wether it was already in the list or after addition
    */
	public int add(T object) {
		int index = getIndexOrInsertIndex(object);
	    if (index < 0) {
	    	index = -(index+1);
	        mList.add(index, object);
	    	}

	    return index;
		}

	/**
	 * Adds object to the list provided that it doesn't contain
	 * an object being considered equal by compareTo().
	 * @param object
	 * @return whether the object was new and therefore added to the list
	 */
	public boolean addIfNew(T object) {
		int index = getIndexOrInsertIndex(object);
		if (index >= 0)
			return false;

		index = -(index+1);
		mList.add(index, object);
		return true;
		}

	public int size(){
		return mList.size();
		}

	/**
	 * Returns object at given index, or null if index==-1
	 * @param index existing index or < 0
	 * @return object or null
	 */
	public T get(int index) {
		return (index < 0) ? null : mList.get(index);
		}

	public T[] toArray(T[] e) {
	    return mList.toArray(e);
	    }

	public void remove(int index) {
	    mList.remove(index);
	    }

	public void removeAll() {
    	mList.clear();
		}
	}
