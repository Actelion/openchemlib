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
import com.actelion.research.chem.forcefield.mmff.Tables;

/**
 * Bond table, corresponds to the MMFFBOND.PAR parameters table provided
 * in the MMFF literature. This table holds parameters for bond stretching
 * interactions.
 */
public final class Bond implements com.actelion.research.chem.forcefield.mmff.Searchable {
    private final Object[][] table;
    private final Tables t;

    public Bond(Tables t, String csvpath) {
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
     * Returns the 'r0' value at the given index from the Bond table.
     *  @param index The table index for the row.
     *  @return 'r0' for that index.
     */
    public double r0(int index) {
        return ((Number) table[index][4]).doubleValue();
    }

    /**
     * Returns the 'kb' value at the given index from the Bond table.
     *  @param index The table index for the row.
     *  @return 'kb' for that index.
     */
    public double kb(int index) {
        return ((Number) table[index][3]).doubleValue();
    }

    /**
     * Returns 'r0' given a molecule and two atoms in that molecule that
     * form a bond.
     *  @param mol The molecule containing the bond.
     *  @param atom1 The first atom of the bond.
     *  @param atom2 The second atom of the bond.
     *  @return 'r0' for that bond.
     */
    public double r0(MMFFMolecule mol, int atom1, int atom2) {
        int bondt = com.actelion.research.chem.forcefield.mmff.type.Bond.getType(t, mol, mol.getBond(atom1, atom2));
        int a1t = mol.getAtomType(atom1);
        int a2t = mol.getAtomType(atom2);

        // Ensure that atom types are in ascending order.
        if (a1t > a2t)
            a1t = com.actelion.research.chem.forcefield.mmff.Search.s(a2t, a2t = a1t);

        // Search the bond tabke for the atom types and bond type.
        int index = com.actelion.research.chem.forcefield.mmff.Search.binary(new int[]{1,2,0},
                new int[]{a1t, a2t, bondt}, this);

        if (index >= 0)
            // Return r0 if it is in the table.
            return r0(index);
        else {
            int a1no = mol.getAtomicNo(atom1);
            int a2no = mol.getAtomicNo(atom2);

            if (a1no > a2no)
                a1no = com.actelion.research.chem.forcefield.mmff.Search.s(a2no, a2no = a1no);

            int a1cri = com.actelion.research.chem.forcefield.mmff.Search.binary(0, a1no, t.covrad);
            int a2cri = com.actelion.research.chem.forcefield.mmff.Search.binary(0, a2no, t.covrad);

            double r0_i = t.covrad.r0(a1cri);
            double r0_j = t.covrad.r0(a2cri);
            double chi_i = t.covrad.chi(a1cri);
            double chi_j = t.covrad.chi(a2cri);
            double c = a1no == 1 || a2no == 1 ? 0.050 : 0.085;

            return r0_i + r0_j - c*Math.pow(Math.abs(chi_i - chi_j), 1.4);
        }
    }

    /**
     * Returns 'kb' given a molecule and two atoms in that molecule that
     * form a bond.
     *  @param mol The molecule containing the bond.
     *  @param atom1 The first atom of the bond.
     *  @param atom2 The second atom of the bond.
     *  @return 'kb' for that bond.
     */
    public double kb(MMFFMolecule mol, int atom1,
                     int atom2) {
        int bondt = com.actelion.research.chem.forcefield.mmff.type.Bond.getType(t, mol, mol.getBond(atom1, atom2));
        int a1t = mol.getAtomType(atom1);
        int a2t = mol.getAtomType(atom2);

        // Ensure that atom types are in ascending order.
        if (a1t > a2t)
            a1t = com.actelion.research.chem.forcefield.mmff.Search.s(a2t, a2t = a1t);

        // Search the bond tabke for the atom types and bond type.
        int index = com.actelion.research.chem.forcefield.mmff.Search.binary(new int[]{1,2,0},
                new int[]{a1t, a2t, bondt}, this);

        if (index >= 0) {
            // Return r0 if it is in the table.
            return kb(index);
        } else {
            int a1no = mol.getAtomicNo(atom1);
            int a2no = mol.getAtomicNo(atom2);

            if (a1no > a2no)
                a1no = com.actelion.research.chem.forcefield.mmff.Search.s(a2no, a2no = a1no);

            double r0 = r0(mol, atom1, atom2);

            // Search Bndk.
            int bndk_idx = com.actelion.research.chem.forcefield.mmff.Search.binary(new int[]{0, 1},
                    new int[]{a1no, a2no}, t.bndk);

            if (bndk_idx >= 0) {
                // Found a value in Bndk.
                double coeff = Math.pow(t.bndk.r0(bndk_idx) / r0, 6);
                return t.bndk.kb(bndk_idx) * coeff;
            } else {
                // Fall back to HerschbachLaurie.
                int bidx = com.actelion.research.chem.forcefield.mmff.Search.binary(new int[]{0,1},
                        new int[]{PeriodicTable.rowtm(a1no),
                            PeriodicTable.rowtm(a2no)},
                        t.hblaurie);
                return Math.pow(10.0, -(r0 - t.hblaurie.fget(bidx, 2))
                        / t.hblaurie.fget(bidx, 3));
            }
        }
    }
}
