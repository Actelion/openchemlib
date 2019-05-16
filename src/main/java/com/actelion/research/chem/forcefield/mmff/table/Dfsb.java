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
import com.actelion.research.chem.forcefield.mmff.MMFFMolecule;
import com.actelion.research.chem.forcefield.mmff.PeriodicTable;
import com.actelion.research.chem.forcefield.mmff.Search;
import com.actelion.research.chem.forcefield.mmff.Searchable;
import com.actelion.research.chem.forcefield.mmff.Tables;

public final class Dfsb implements Searchable {
    private final Object[][] table;

    public Dfsb(Tables t, String csvpath) {
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
     * Returns the index that a given set of three atoms have in the Dfsb
     * table.
     *  @param mol The molecule containing the atoms.
     *  @param a1 Atom 1.
     *  @param a2 Atom 2.
     *  @param a3 Atom 3.
     *  @return The row index in the table.
     */
    public int index(MMFFMolecule mol, int a1, int a2, int a3) {
        int a1r = PeriodicTable.row(mol.getAtomicNo(a1));
        int a2r = PeriodicTable.row(mol.getAtomicNo(a2));
        int a3r = PeriodicTable.row(mol.getAtomicNo(a3));

        if (a1r > a3r)
            a1r = Search.s(a3r, a3r = a1r);

        return Search.binary(new int[]{1,0,2},
                new int[]{a2r,a1r,a3r}, this);
    }

    /**
     * Returns the equivalent 'kb' value for a given set of three atoms.
     *  @param mol The molecule containing the atoms.
     *  @param a1 Atom 1.
     *  @param a2 Atom 2.
     *  @param a3 Atom 3.
     *  @return Returns the 'kb' value or 0 if no entry was found.
     */
    public double kb(MMFFMolecule mol, int a1, int a2, int a3) {
        int a1r = PeriodicTable.row(mol.getAtomicNo(a1));
        int a3r = PeriodicTable.row(mol.getAtomicNo(a3));
        int at = a1r > a3r ? 4 : 3;
        int idx = index(mol, a1, a2, a3);

        return idx >= 0 ? ((Number) table[idx][at]).doubleValue() : 0.0;
    }
}
