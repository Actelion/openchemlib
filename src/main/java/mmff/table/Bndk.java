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
 * Bndk table, corresponds to the MMFFBNDK.PAR parameters table provided
 * in the MMFF literature. This table provides parameters used in the
 * empirical rules for bond-stretching force constants.
 */
public final class Bndk implements mmff.Searchable {
    private final Object[][] table;

    public Bndk(Tables t, String csvpath) {
        table = Csv.readFile(csvpath);
    }

    public int get(int row, int col) {
        return ((Number) table[row][col]).intValue();
    }

    public int length() {
        return table.length;
    }

    /**
     * Returns 'r0' the ideal bond length at a given index in the Bndk
     * table.
     *  @param index The index of the desired row.
     *  @return The value of 'r0' at index.
     */
    public double r0(int index) {
        return ((Number) table[index][2]).doubleValue();
    }

    /**
     * Returns 'kb' the force constant at a given index in the Bndk table.
     *  @param index The index of the desired row.
     *  @return The value of 'kb' at index.
     */
    public double kb(int index) {
        return ((Number) table[index][3]).doubleValue();
    }
}
