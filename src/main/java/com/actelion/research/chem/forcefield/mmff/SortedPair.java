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

package com.actelion.research.chem.forcefield.mmff;

/**
 * SortedPair class is used in the Separation table. It holds a pair of
 * integers and ensures that they are sorted such that the first value is
 * larger than the second. The Sorted pair is used to ensure that
 * duplicate values are not stored in the separation table and that the
 * Hashtable can be keyed on both atom indices.
 */
public class SortedPair implements Comparable<SortedPair> {
    public final int a;
    public final int b;

    /**
     * Construct a new SortedPair. Make sure that A > B.
     */
    public SortedPair(int a, int b) {
        this.a = a>b ? a : b;
        this.b = a>b ? b : a;
    }

    /**
     * Returns the hash code of this sorted pair object. Defined as the
     * XOR of the hashes of A and B.
     *  @return The hash code.
     */
    @Override
    public int hashCode() {
        return new Integer(a).hashCode() ^ new Integer(b).hashCode();
    }

    /**
     * Returns true if this SortedPair is equal with another object.
     *  @param obj The object to compare with.
     *  @return True if SortedPair and object are equal.
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == this)
            return true;

        if (!(obj instanceof SortedPair))
            return false;

        SortedPair that = (SortedPair)obj;

        return a == that.a && b == that.b;
    }

    /**
     * Compares this object with another. This allows collections of
     * SortedPairs to be sorted using the java standard library sorting
     * algorithms.
     */
    @Override
    public int compareTo(SortedPair that) {
        if (a > that.a)
            return 1;
        if (a < that.a)
            return -1;
        if (b > that.b)
            return 1;
        if (b < that.b)
            return -1;
        return 0;
    }

    @Override
    public String toString() {
        return a+","+b;
    }
}
