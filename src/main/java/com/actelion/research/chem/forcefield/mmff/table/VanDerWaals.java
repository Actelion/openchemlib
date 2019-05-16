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

package com.actelion.research.chem.forcefield.mmff.table;

import com.actelion.research.chem.forcefield.mmff.Csv;
import com.actelion.research.chem.forcefield.mmff.Tables;

public final class VanDerWaals implements com.actelion.research.chem.forcefield.mmff.Searchable {
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
        return com.actelion.research.chem.forcefield.mmff.Search.binary(0, type, this);
    }
}
