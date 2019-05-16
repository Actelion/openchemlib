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

/**
 * Bndk table, corresponds to the MMFFBNDK.PAR parameters table provided
 * in the MMFF literature. This table provides parameters used in the
 * empirical rules for bond-stretching force constants.
 */
public final class Bndk implements com.actelion.research.chem.forcefield.mmff.Searchable {
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
