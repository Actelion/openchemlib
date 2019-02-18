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
package mmff;

public final class PeriodicTable {
    private PeriodicTable(){};

    /**
     * Returns which row in the periodic table an atom with given atomic
     * number is in.
     *  @param ano Atomic number of the given atom.
     *  @return The periodic table row of the given atom.
     */
    public static int row(int ano) {
        if (ano < 3)
            return 0;
        else if (ano < 11)
            return 1;
        else if (ano < 19)
            return 2;
        else if (ano < 37)
            return 3;
        else if (ano < 55)
            return 4;
        return 0;
    }

    /**
     * Returns which row in the periodic table an atom is in with transition
     * metals having their row number multiplied by 10. Hydrogen is counted in
     * row 0 while helium is counted in row 2.
     *  @param ano The atoms atomic number.
     *  @return The periodic table row.
     */
    public static int rowtm(int ano) {
        int row = 0;

        if (ano == 2)
            row = 1;
        else if (ano >= 3 && ano <= 10)
            row = 2;
        else if (ano >= 11 && ano <= 18)
            row = 3;
        else if (ano >= 19 && ano <= 36)
            row = 4;
        else if (ano >= 37 && ano <= 54)
            row = 5;

        if (ano >= 21 && ano <= 30 || ano >= 39 && ano <= 48)
            row *= 10;

        return row;
    }
}
