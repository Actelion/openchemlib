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
 * The search class provides binary searching on arrays. This is used to
 * support searching arrays with duplicate values which the Java library
 * binary search does not support.
 */
public class Search {
    /**
     * Function used in swapping variables in one line. Swapping should be
     * done as follows:
     *      int a = 5;
     *      int b = 7;
     *      a = s(b, b = a);
     *      // now a = 7 and b = 5.
     * Generic function provided that the variables to be swapped are the
     * same type.
     */
    public static <T> T s(T a, T b) {
        return a;
    }

    /**
     * Binary searches for a value and returns the index to it. Works on
     * any given array that implements the ArraySearch interface. Note
     * that the array must be sorted. Duplicate values may be included,
     * in which case the 'findlow' parameter determines which end of the
     * continuous block of values to search. If the search value is unique
     * in the array then 'findlow' does not affect the returned index.
     *  @param col The column in the array to search, it should be an
     *      integer value.
     *  @param val The value to be searched for.
     *  @param lo The starting low index to begin searching from. This
     *      is included in the search.
     *  @param hi The ending high index to end searching in. This index is
     *      EXCLUDED in the search, with the value below it being searched.
     *  @param findlow True if the binary search should find the lowest value
     *      if there are duplicates, False if the binary search should find
     *      the highest value if there are duplicates.
     *  @param arr The array object to be searched. This array object must
     *      implement the ArraySearch interface.
     *  @return The lowest/highest index to the value.
     */
    public static int binary(int col, int val, int lo,
            int hi, boolean findlow, Searchable arr) {
        lo = Math.max(0,lo);
        hi = Math.min(arr.length(), hi);
        int at = (hi-lo)/2;
        int min = lo;
        int max = hi;

        // Catch cases where the value is obviously outside the array (it's
        // too large or small).
        if (arr.get(min, col) > val || arr.get(max-1, col) < val)
            return -1;

        while (hi >= lo) {
            at = lo + (hi-lo)/2;
            int vat = arr.get(at, col);

            if (vat == val) {
                if (findlow && at > min && arr.get(at-1, col) == val)
                    hi = at;
                else if (!findlow && at < max-1 && arr.get(at+1, col) == val)
                    lo = at;
                else
                    return at;
            } else if (vat > val)
                hi = at;
            else if (vat < val)
                lo = at;

            // we couldn't find the value
            if (hi-lo == 1 && arr.get(lo, col) < val && arr.get(hi, col) > val)
                break;
        }
        return -1;
    }

    /**
     * Helper function that acts like Search.binary(int[], int[],
     * Searchable) but takes a single int for cols and vals to search just
     * one column. Internally it wraps directly to the raw binary search
     * function.
     */
    public static int binary(int col, int val, Searchable data) {
        return binary(col, val, 0, data.length(), true, data);
    }

    /**
     * Binary search on a collection of values in columns in a table. This
     * does not assume any relation between the columns, only that columns
     * subsequently searched are sorted by value across the interval being
     * searched. In practice this means that the data should be sorted by
     * the first column, then the second column, then the third column etc.
     * This function can search an arbitrary number of columns, including
     * duplicate values. When searching the last column, if duplicate
     * values are found the binary search will return the index of the
     * lowest value.
     *  @param cols Array of column indicies to search in. cols[0] is
     *      searched first, then cols[1], then cols[2] etc.
     *  @param vals Array of values to be searched for in their
     *      corresponding columns at the indices in cols. This array
     *      should be the same length as cols.
     *  @param data An object implementing Searchable which holds the data
     *      to be searched.
     *  @return The index of the row that matches all the values in vals.
     *      Returns -1 if the value was not found or if cols and vals did
     *      not have equal lengths.
     */
    public static int binary(int[] cols, int[] vals, Searchable data) {
        if (cols.length != vals.length)
            return -1;

        int lo = 0;
        int hi = data.length();

        int i;
        for (i=0; i<cols.length-1; i++) {
            lo = binary(cols[i], vals[i], lo, hi+1, true, data);
            hi = binary(cols[i], vals[i], lo, hi+1, false, data);

            if (lo == -1 || hi == -1)
                return -1;
        }

        return binary(cols[i], vals[i], lo, hi+1, true, data);
    }
}
