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
import mmff.Search;
import mmff.Searchable;
import mmff.Tables;
import mmff.MMFFMolecule;

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
