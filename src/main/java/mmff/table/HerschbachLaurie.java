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

import mmff.Tables;
import mmff.Csv;

/**
 * Table for Herschbach-Laurie version of Badger's rule. This table
 * provides parameters used in the empirical calculation of bond
 * stretching parameters.
 */
public final class HerschbachLaurie implements mmff.Searchable {
    private final Object[][] table;

    public HerschbachLaurie(Tables t, String csvpath) {
        table = Csv.readFile(csvpath);
    }

    @Override
    public int get(int row, int col) {
        return ((Number) table[row][col]).intValue();
    }

    @Override
    public int length() {
        return table.length;
    }

    /**
     * Gets a floating point value from a given row and column.
     *  @param row The row to fetch.
     *  @param col The column in the row to fetch.
     *  @return The floating point value at (row,col).
     */
    public double fget(int row, int col) {
        return ((Number) table[row][col]).doubleValue();
    }
}
