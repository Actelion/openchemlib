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
 * Atom table, corresponds to the MMFFPROP.PAR parameters table provided
 * in the MMFF literature. This table specifies the chemical, topological
 * and geometrical properties associated with each of the MMFF atom types.
 */
public final class Atom implements com.actelion.research.chem.forcefield.mmff.Searchable {
    private final int[][] table;

    public Atom(Tables t, String csvpath) {
        table = Csv.readIntsFile(csvpath);
    }

    @Override
    public int get(int row, int col) {
        return table[row][col];
    }

    @Override
    public int length() {
        return table.length;
    }

    /**
     * Returns the ASPEC type of an atom given its MMFF type.
     *  @param type The MMFF atom type of the atom.
     *  @return The ASPEC value for the atom.
     */
    public int aspec(int type) {
        return table[index(type)][1];
    }

    /**
     * Returns the CRD type of an atom given its MMFF type.
     *  @param type The MMFF atom type of the atom.
     *  @return The CRD value for the atom.
     */
    public int crd(int type) {
        return table[index(type)][2];
    }

    /**
     * Returns the VAL type of an atom given its MMFF type.
     *  @param type The MMFF atom type of the atom.
     *  @return The VAL type for the atom.
     */
    public int val(int type) {
        return table[index(type)][3];
    }

    /**
     * Returns the PILP type of an atom given its MMFF type.
     *  @param type The MMFF atom type of the atom.
     *  @return The PILP type for the atom.
     */
    public int pilp(int type) {
        return table[index(type)][4];
    }

    /**
     * Returns the MLTB type of an atom given its MMFF type.
     *  @param type The MMFF atom type of the atom.
     *  @return The MLTB type for the atom.
     */
    public int mltb(int type) {
        return table[index(type)][5];
    }

    /**
     * Returns the Arom bool of an atom given its MMFF type.
     *  @param type The MMFF atom type of the atom.
     *  @return The Arom bool for the atom.
     */
    public boolean arom(int type) {
        return table[index(type)][6] > 0;
    }

    /**
     * Returns the linear bool of an atom given its MMFF type.
     *  @param type The MMFF atom type of the atom.
     *  @return The Linear bool for the atom.
     */
    public boolean linear(int type) {
        int idx = index(type);
        return idx >= 0 ? table[idx][7] > 0 : false;
    }

    /**
     * Returns the SBMB type of an atom given its MMFF type.
     *  @param type The MMFF atom type of the atom.
     *  @return The SBMB boolean for the atom.
     */
    public boolean sbmb(int type) {
        return table[index(type)][8] > 0;
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
