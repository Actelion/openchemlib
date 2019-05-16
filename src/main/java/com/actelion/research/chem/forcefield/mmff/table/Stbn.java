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
import com.actelion.research.chem.forcefield.mmff.Search;
import com.actelion.research.chem.forcefield.mmff.Tables;

public final class Stbn implements com.actelion.research.chem.forcefield.mmff.Searchable {
    private final Object[][] table;
    private final Tables t;

    public Stbn(Tables t, String csvpath) {
        table = Csv.readFile(csvpath);
        this.t = t;
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
     * Returns the index of a row for a given molecule and three
     * connected atoms which form an angle.
     *  @param mol The molecule containing the atoms.
     *  @param a1 Atom 1 (atom i).
     *  @param a2 Atom 2, the central atom (atom j).
     *  @param a3 Atom 3 (atom k).
     *  @return The row index in the table.
     */
    public int index(MMFFMolecule mol, int a1, int a2, int a3) {
        int a1t = mol.getAtomType(a1);
        int a2t = mol.getAtomType(a2);
        int a3t = mol.getAtomType(a3);
        int angt = com.actelion.research.chem.forcefield.mmff.type.Angle.getStbnType(t, mol, a1, a2, a3);

        if (a1t > a3t)
            a1t = Search.s(a3t, a3t = a1t);

        return Search.binary(new int[]{2,1,3,0},
                new int[]{a2t,a1t,a3t,angt}, this);
    }

    /**
     * Returns 'kba' for a given index in the table.
     *  @param mol The molecule containing the atoms.
     *  @param a1 Atom 1 (atom i).
     *  @param a2 Atom 2, the central atom (atom j).
     *  @param a3 Atom 3 (atom k).
     *  @return The value of 'kba'.
     */
    public double kba(MMFFMolecule mol, int a1, int a2, int a3) {
        int a1t = mol.getAtomType(a1);
        int a3t = mol.getAtomType(a3);
        int b1t = com.actelion.research.chem.forcefield.mmff.type.Bond.getType(t, mol, a1, a2);
        int b2t = com.actelion.research.chem.forcefield.mmff.type.Bond.getType(t, mol, a2, a3);

        int idx = index(mol, a1, a2, a3);
        int at = a1t > a3t || a1t == a3t && b1t < b2t ? 1 : 0;

        if (idx >= 0)
            return ((Number) table[idx][4+at]).doubleValue();

        return t.dfsb.kb(mol, a1, a2, a3);
    }
}
