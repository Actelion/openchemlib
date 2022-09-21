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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;

public class UniqueList<T> extends SortedList<T> implements Serializable {
    static final long serialVersionUID = 0x20121016;

    private ArrayList<T> mOriginalOrder = new ArrayList<T>();
	private int[] mOriginalIndex;

	public UniqueList() {
		super();
	}

	public UniqueList(Comparator comparator) {
		super(comparator);
	}

	public int getSortedIndex(T s) {
        return super.getIndex(s);
        }

    @Override
	public boolean contains(T object) {
		return super.getIndex(object) != -1;	// uses the faster inherited getIndex() method
		}

    /**
     * When objects were added after the last getIndex() call,
     * then the original-index-map needs to be re-created.
     */
    @Override
	public int getIndex(T s) {
		int index = super.getIndex(s);
		if (index == -1)
			return -1;

		if (mOriginalIndex == null)
			createOriginalIndex();

		return mOriginalIndex[index];
		}

    @Override
	public int add(T s) {
    	int oldSize = size();
		int index = super.add(s);

		if (size() != oldSize) {
			mOriginalOrder.add(s);
			mOriginalIndex = null;
			return oldSize;
			}

		if (mOriginalIndex == null)
			createOriginalIndex();

		return mOriginalIndex[index];
		}

	public int add(int position, T s) {
		int index = super.add(s);
		if (index == -1)
			return -1;

		mOriginalOrder.add(position, s);
		mOriginalIndex = null;

		return position;
		}

	/**
     * @param i list index within list in original order
     * @return string at position i of list in order of the creation
     */
    @Override
	public T get(int i) {
		return mOriginalOrder.get(i);
		}

	/**
	 * @param i list index within sorted list
	 * @return string at position i of sorted list
	 */
    public T getSorted(int i) {
        return super.get(i);
        }

    @Override
	public T[] toArray(T[] e) {
	    return mOriginalOrder.toArray(e);
	    }

    public T[] toSortedArray(T[] e) {
        return super.toArray(e);
        }

    private void createOriginalIndex() {
    	mOriginalIndex = new int[size()];
    	int i=0;
    	for (T t:mOriginalOrder)
    		mOriginalIndex[super.getIndex(t)] = i++;
    	}
	}
