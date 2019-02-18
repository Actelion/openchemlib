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
import mmff.Tables;

/**
 * Atom table, corresponds to the MMFFPROP.PAR parameters table provided
 * in the MMFF literature. This table specifies the chemical, topological
 * and geometrical properties associated with each of the MMFF atom types.
 */
public final class Atom implements mmff.Searchable {
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
        return mmff.Search.binary(0, type, this);
    }
}
