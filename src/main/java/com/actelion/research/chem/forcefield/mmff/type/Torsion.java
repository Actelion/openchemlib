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

package com.actelion.research.chem.forcefield.mmff.type;

import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.forcefield.mmff.MMFFMolecule;
import com.actelion.research.chem.forcefield.mmff.Tables;

import java.util.HashSet;

/**
 * The torsion type class provides a static function for getting the
 * torsion type of a torsion angle.
 */
public final class Torsion {
    /**
     * Returns the torsion type of a torsion angle.
     *  @param mol The molecule that the atoms are in.
     *  @param a1 Atom 1.
     *  @param a2 Atom 2.
     *  @param a3 Atom 3.
     *  @param a4 Atom 4.
     *  @return The torsion type (combined tortype and secondary tortype).
     */
    public static int getType(Tables table, MMFFMolecule mol, int a1, int a2, int a3,
                              int a4) {
        int bond_jk = mol.getBond(a2, a3);

        int bijt = Bond.getType(table, mol, a1, a2);
        int bjkt = Bond.getType(table, mol, a2, a3);
        int bklt = Bond.getType(table, mol, a3, a4);
        int tort = bjkt;

        if (bjkt == 0 && mol.getBondOrder(bond_jk) == 1
                && (bijt == 1 || bklt == 1))
            tort = 2;

        int size = inRingOfSize(mol, a1, a2, a3, a4);

        if (size == 4 && mol.getBond(a1, a3) == -1
                && mol.getBond(a2, a4) == -1)
            return 40 + tort;

        if (size == 5 && (mol.getAtomType(a1) == 1
                    || mol.getAtomType(a2) == 1
                    || mol.getAtomType(a3) == 1
                    || mol.getAtomType(a4) == 1))
            return 50 + tort;

        return tort;
    }

    /**
     * Checks and returns the minimum ring size that a torsion is
     * contained in. Only checks rings of size 4 and 5.
     *  @param mol The molecule that the atoms are in.
     *  @param a1 Atom 1.
     *  @param a2 Atom 2.
     *  @param a3 Atom 3.
     *  @param a4 Atom 4.
     */
    public static int inRingOfSize(MMFFMolecule mol, int a1, int a2,
                                   int a3, int a4) {

        // If any of the four atoms don't have bonds between them, then they
        // don't form a torsion angle and there fore this invalid torsion
        // cannot be in a ring.
        if (mol.getBond(a1, a2) == -1 || mol.getBond(a2, a3) == -1
                || mol.getBond(a3, a4) == -1)
            return 0;

        // If atom 4 is connected to atom 1, then the torsion loops back on to
        // itself and is in a ring of size 4.
        if (mol.getBond(a4, a1) >= 0)
            return 4;


        RingCollection rings = mol.getRingSet();
        HashSet<Integer> tor = new HashSet<Integer>();
        tor.add(a1);
        tor.add(a2);
        tor.add(a3);
        tor.add(a4);

        // Loop over all the molecule rings, for all rings of size 5 check if
        // the set of ring atoms contains the torsion angle atoms with a set
        // intersection. If true then the ring contains the torsion.
        for (int r=0; r<rings.getSize(); r++)
            if (rings.getRingSize(r) == 5) {
                HashSet<Integer> ring = new HashSet<Integer>();

                for (int a : rings.getRingAtoms(r))
                    ring.add(a);

                if (ring.containsAll(tor))
                    return 5;
            }

        return 0;
    }
}
