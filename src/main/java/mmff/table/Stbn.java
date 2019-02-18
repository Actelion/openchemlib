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

import mmff.Search;
import mmff.Tables;
import mmff.Csv;
import mmff.MMFFMolecule;

public final class Stbn implements mmff.Searchable {
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
        int angt = mmff.type.Angle.getStbnType(t, mol, a1, a2, a3);

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
        int b1t = mmff.type.Bond.getType(t, mol, a1, a2);
        int b2t = mmff.type.Bond.getType(t, mol, a2, a3);

        int idx = index(mol, a1, a2, a3);
        int at = a1t > a3t || a1t == a3t && b1t < b2t ? 1 : 0;

        if (idx >= 0)
            return ((Number) table[idx][4+at]).doubleValue();

        return t.dfsb.kb(mol, a1, a2, a3);
    }
}
