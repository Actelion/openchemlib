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

import mmff.Csv;
import mmff.PeriodicTable;
import mmff.Tables;
import mmff.MMFFMolecule;

/**
 * Bond table, corresponds to the MMFFBOND.PAR parameters table provided
 * in the MMFF literature. This table holds parameters for bond stretching
 * interactions.
 */
public final class Bond implements mmff.Searchable {
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
        int bondt = mmff.type.Bond.getType(t, mol, mol.getBond(atom1, atom2));
        int a1t = mol.getAtomType(atom1);
        int a2t = mol.getAtomType(atom2);

        // Ensure that atom types are in ascending order.
        if (a1t > a2t)
            a1t = mmff.Search.s(a2t, a2t = a1t);

        // Search the bond tabke for the atom types and bond type.
        int index = mmff.Search.binary(new int[]{1,2,0},
                new int[]{a1t, a2t, bondt}, this);

        if (index >= 0)
            // Return r0 if it is in the table.
            return r0(index);
        else {
            int a1no = mol.getAtomicNo(atom1);
            int a2no = mol.getAtomicNo(atom2);

            if (a1no > a2no)
                a1no = mmff.Search.s(a2no, a2no = a1no);

            int a1cri = mmff.Search.binary(0, a1no, t.covrad);
            int a2cri = mmff.Search.binary(0, a2no, t.covrad);

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
        int bondt = mmff.type.Bond.getType(t, mol, mol.getBond(atom1, atom2));
        int a1t = mol.getAtomType(atom1);
        int a2t = mol.getAtomType(atom2);

        // Ensure that atom types are in ascending order.
        if (a1t > a2t)
            a1t = mmff.Search.s(a2t, a2t = a1t);

        // Search the bond tabke for the atom types and bond type.
        int index = mmff.Search.binary(new int[]{1,2,0},
                new int[]{a1t, a2t, bondt}, this);

        if (index >= 0) {
            // Return r0 if it is in the table.
            return kb(index);
        } else {
            int a1no = mol.getAtomicNo(atom1);
            int a2no = mol.getAtomicNo(atom2);

            if (a1no > a2no)
                a1no = mmff.Search.s(a2no, a2no = a1no);

            double r0 = r0(mol, atom1, atom2);

            // Search Bndk.
            int bndk_idx = mmff.Search.binary(new int[]{0, 1},
                    new int[]{a1no, a2no}, t.bndk);

            if (bndk_idx >= 0) {
                // Found a value in Bndk.
                double coeff = Math.pow(t.bndk.r0(bndk_idx) / r0, 6);
                return t.bndk.kb(bndk_idx) * coeff;
            } else {
                // Fall back to HerschbachLaurie.
                int bidx = mmff.Search.binary(new int[]{0,1},
                        new int[]{PeriodicTable.rowtm(a1no),
                            PeriodicTable.rowtm(a2no)},
                        t.hblaurie);
                return Math.pow(10.0, -(r0 - t.hblaurie.fget(bidx, 2))
                        / t.hblaurie.fget(bidx, 3));
            }
        }
    }
}
