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
 * CovRad table, used in the bond stretch empirical calculations of 'kb'
 * and 'r0'.
 */
public final class CovRad implements mmff.Searchable {
    private final Object[][] table;

    public CovRad(Tables t, String csvpath) {
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
     * Returns the 'r0' (or covRad) parameter from the table given an
     * index.
     *  @param index The index in the table.
     *  @return The value of 'r0'.
     */
    public double r0(int index) {
        return ((Number) table[index][1]).doubleValue();
    }

    /**
     * Returns the 'chi' (or pauEle) parameter from the table at the given
     * index.
     *  @param index The index in the table.
     *  @return The value of 'chi'.
     */
    public double chi(int index) {
        return ((Number) table[index][2]).doubleValue();
    }
}
