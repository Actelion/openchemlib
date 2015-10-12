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

import java.util.Arrays;

public class CanonizerBaseValue implements Comparable<CanonizerBaseValue> {
	private static final int BASE_VALUE_SIZE = 3;	// bond query features need third long
    public long[] mValue;
    private int mAtom;
    private int mIndex;
    private int mAvailableBits;

    public CanonizerBaseValue() {
        mValue = new long[BASE_VALUE_SIZE];
        }

    public void init(int atom) {
        mAtom = atom;
        mIndex = 0;
        mAvailableBits = 63;
        Arrays.fill(mValue, 0);
        }

    public void add(int data) {
        mValue[mIndex] += data;
        }

    public void add(int bits, int data) {
        if (mAvailableBits == 0) {
            mIndex++;
            mAvailableBits = 63;
            }
        if (mAvailableBits == 63) {
            mValue[mIndex] |= data;
            mAvailableBits -= bits;
            }
        else {
            if (mAvailableBits >= bits) {
                mValue[mIndex] <<= bits;
                mValue[mIndex] |= data;
                mAvailableBits -= bits;
                }
            else {
                mValue[mIndex] <<= mAvailableBits;
                mValue[mIndex] |= (data >> (bits - mAvailableBits));
                bits -= mAvailableBits;
                mIndex++;
                mAvailableBits = 63 - bits;
                mValue[mIndex] |= (data & ((1 << bits) - 1));
                }
            }
        }

    public int getAtom() {
        return mAtom;
        }

    public int compareTo(CanonizerBaseValue b) {
        if (mIndex != b.mIndex)
            return (mIndex < b.mIndex) ? -1 : 1;
        for (int i=0; i<mIndex; i++)
            if (mValue[i] != b.mValue[i])
                return (mValue[i] < b.mValue[i]) ? -1 : 1;
        return (mValue[mIndex] == b.mValue[mIndex]) ? 0
             : (mValue[mIndex] < b.mValue[mIndex]) ? -1 : 1;
        }
    }
