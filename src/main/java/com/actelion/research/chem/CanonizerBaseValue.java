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

import java.util.Arrays;

public class CanonizerBaseValue implements Comparable<CanonizerBaseValue> {
    public long[] mValue;
    private int mAtom;
    private int mIndex;
    private int mAvailableBits;

    /**
     * @param size depends on the maximum number of non-H neighbors
     *             and whether bond query features are present
     */
    public CanonizerBaseValue(int size) {
        mValue = new long[size];
        }

    public void init(int atom) {
        mAtom = atom;
        mIndex = 0;
        mAvailableBits = 63;
        Arrays.fill(mValue, 0);
        }

    public void add(long data) {
        mValue[mIndex] += data;
        }

    public void add(int bits, long data) {
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
        assert(mIndex == b.mIndex);
        for (int i=0; i<mIndex; i++)
            if (mValue[i] != b.mValue[i])
                return (mValue[i] < b.mValue[i]) ? -1 : 1;
        return (mValue[mIndex] == b.mValue[mIndex]) ? 0
             : (mValue[mIndex] < b.mValue[mIndex]) ? -1 : 1;
        }
    }
