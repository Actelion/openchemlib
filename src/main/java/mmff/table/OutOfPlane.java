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

import java.util.Arrays;

import mmff.Csv;
import mmff.Search;
import mmff.Searchable;
import mmff.Tables;
import mmff.MMFFMolecule;

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
