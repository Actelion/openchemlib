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

public final class VanDerWaals implements mmff.Searchable {
    public final double power = 0.25;
    public final double b = 0.2;
    public final double beta = 12.0;
    public final double darad = 0.8;
    public final double daeps = 0.5;
    // power      B       Beta     DARAD      DAEPS\n"
    // 0.25       0.2     12.      0.8        0.5\n"

    private final Object[][] table;

    public VanDerWaals(Tables t, String csvpath) {
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
     * Returns 'alpha-i' from the table.
     *  @param type The MMFF atom type.
     *  @return alpha-i.
     */
    public double alpha_i(int type) {
        return ((Number) table[index(type)][1]).doubleValue();
    }

    /**
     * Returns 'N-i' from the table.
     *  @param type The MMFF atom type.
     *  @return N-i.
     */
    public double n_i(int type) {
        return ((Number) table[index(type)][2]).doubleValue();
    }

    /**
     * Returns 'A-i' from the table.
     *  @param type The MMFF atom type.
     *  @return A-i.
     */
    public double a_i(int type) {
        return ((Number) table[index(type)][3]).doubleValue();
    }

    /**
     * Returns 'G-i' from the table.
     *  @param type The MMFF atom type.
     *  @return G-i.
     */
    public double g_i(int type) {
        return ((Number) table[index(type)][4]).doubleValue();
    }

    /**
     * Returns 'DA' from the table.
     *  @param type The MMFF atom type.
     *  @return DA.
     */
    public char da(int type) {
        return ((Character) table[index(type)][5]).charValue();
    }

    /**
     * Returns 'R*' which is derived from the table.
     *  @param type The MMFF atom type.
     *  @return R*.
     */
    public double r_star(int type) {
        int idx = index(type);
        return ((Number) table[idx][3]).doubleValue()
            * Math.pow(((Number) table[idx][1]).doubleValue(), power);
    }

    /**
     * Returns the index of a given atom type in the 'atprop' table.
     *  @param type The atom type.
     *  @return An index in the atprop table, or -1 if no index was found.
     */
    private int index(int type) {
        return mmff.Search.binary(0, type, this);
    }
}
