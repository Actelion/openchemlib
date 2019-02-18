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

package mmff.table;

import mmff.Search;
import mmff.Searchable;
import mmff.Csv;
import mmff.Tables;

public class Charge {
    private final Object[][] pbci;
    private final Object[][] bci;

    public Charge(Tables t, String csv_pbci, String csv_bci) {
        pbci = Csv.readFile(csv_pbci);
        bci = Csv.readFile(csv_bci);
    }

    /**
     * Returns the formal charge adjustment for a given atom type.
     *  @param type The MMFF atom type for an atom.
     *  @return The formal charge adjustment as a double.
     */
    public double getFcadj(int type) {
        return ((Number) pbci[type-1][2]).doubleValue();
    }

    public double getPbci(int type) {
        return ((Number) pbci[type-1][1]).doubleValue();
    }

    /**
     * Gets the partial charge of a bond type and its two atom types.
     *  @param bondt The MMFF bond type.
     *  @param a1t The MMFF atom type of atom 1.
     *  @param a2t The MMFF atom type of atom 2.
     *  @return The partial charge.
     */
    public double getPartial(int bondt, int a1t, int a2t) {
        double sign = a1t > a2t ? 1.0 : -1.0;

        int atl = a1t > a2t ? a2t : a1t;
        int ath = a1t > a2t ? a1t : a2t;

        int a1lo = bci_binary_search(1, atl, 0, bci.length, true);
        int a1hi = bci_binary_search(1, atl, 0, bci.length, false);

        // Couldn't find the first atom type.
        if (a1lo == -1 || a1hi == -1)
            return getPbci(a1t) - getPbci(a2t);

        int a2lo = bci_binary_search(2, ath, a1lo, a1hi+1, true);
        int a2hi = bci_binary_search(2, ath, a1lo, a1hi+1, false);

        // couldn't find the second atom type.
        if (a2lo == -1 || a2hi == -1)
            return getPbci(a1t) - getPbci(a2t);

        if (bondt == 0 && get_bci_n(a2lo, 0) == 0)
            return sign * get_bci_f(a2lo, 3);
        else if (bondt == 1 && get_bci_n(a2hi, 0) == 1)
            return sign * get_bci_f(a2hi, 3);

        // Couldn't find the bond type.
        return getPbci(a1t) - getPbci(a2t);
    }

    /**
     * Binary search in the bci table. This is a wrapper function around
     * 'binary_search' which just passes a SearchArray object for bci.
     *  @param col The column in the array to search, it should be an
     *      integer value.
     *  @param val The value to be searched for.
     *  @param lo The starting low index to begin searching from. This
     *      index is included in the search.
     *  @param hi The ending high index to end searching in. This index is
     *      EXCLUDED in the search, with the value below it being searched.
     *  @param findlow True if the binary search should find the lowest
     *      value if there are duplicates, False if the binary search
     *      should find the highest value if there are duplicates.
     *  @return The lowest/highest index to the value.
     */
    public int bci_binary_search(int col, int val, int lo, int hi,
            boolean findlow) {
        return Search.binary(col, val, lo, hi, findlow, new Searchable(){
            public int get(int row, int col) {
                return ((Number) bci[row][col]).intValue();
            }

            public int length() {
                return bci.length;
            }
        });
    }

    public double get_bci_f(int row, int col) {
        return ((Number) bci[row][col]).doubleValue();
    }

    public int get_bci_n(int row, int col) {
        return ((Number) bci[row][col]).intValue();
    }

    public int get_bci_len() {
        return bci.length;
    }
}
