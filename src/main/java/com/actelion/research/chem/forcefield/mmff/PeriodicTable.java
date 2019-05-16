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
package com.actelion.research.chem.forcefield.mmff;

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
