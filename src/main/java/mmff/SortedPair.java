/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of ActelionMMFF94.
 * 
 * ActelionMMFF94 is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * ActelionMMFF94 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with ActelionMMFF94.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Paolo Tosco,Daniel Bergmann
 */

package mmff;

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
