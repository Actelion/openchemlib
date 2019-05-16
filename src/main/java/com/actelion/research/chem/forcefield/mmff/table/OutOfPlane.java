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

import java.util.Arrays;

import com.actelion.research.chem.forcefield.mmff.Csv;
import com.actelion.research.chem.forcefield.mmff.MMFFMolecule;
import com.actelion.research.chem.forcefield.mmff.Search;
import com.actelion.research.chem.forcefield.mmff.Searchable;
import com.actelion.research.chem.forcefield.mmff.Tables;

public class OutOfPlane implements Searchable {
    private final Object[][] table;
    private final Tables t;

    /**
     *
     */
    public OutOfPlane(Tables t, String csvpath) {
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
     * Returns the Koop parameters from the OutOfPlane table.
     *  @param mol The molecule that the atoms are in.
     *  @param ac The central atom.
     *  @param a1 Atom neighbour 1.
     *  @param a2 Atom neighbour 2.
     *  @param a3 Atom neighbour 3.
     *  @return The Koop floating point value.
     */
    public double getKoop(MMFFMolecule mol, int ac, int a1, int a2,
                          int a3) {
        int act = mol.getAtomType(ac);
        int[] nbr = new int[]{
            mol.getAtomType(a1),
            mol.getAtomType(a2),
            mol.getAtomType(a3)};

        for (int i=0; i<4; i++) {
            int[] nbrf = new int[3];

            for (int n=0; n<3; n++)
                nbrf[n] = t.def.table[nbr[n]-1][i+1];
            Arrays.sort(nbrf);

            int index = Search.binary(new int[]{1,0,2,3},
                new int[]{act,nbrf[0],nbrf[1],nbrf[2]}, this);

            if (index >= 0)
                return ((Number) table[index][4]).doubleValue();
        }
        return 0.0;
    }

    public double getKoop(int index) {
        return ((Number) table[index][4]).doubleValue();
    }
}
