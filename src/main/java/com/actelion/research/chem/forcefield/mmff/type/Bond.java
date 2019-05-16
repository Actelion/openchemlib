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

/**
 * The bond type class provides static functions for getting the MMFF bond
 * type (either by bond or by the atoms of a bond).
 */
public class Bond {
    /**
     * Returns the MMFF bond type of a bond on a molecule. This function
     * assumes that a valid bond is passed for the molecule.
     *  @param mol The molecule that contains the bond.
     *  @param bond The bond to be typed.
     *  @return The MMFF bond type.
     */
    public static int getType(Tables table, MMFFMolecule mol, int bond) {
        return getType(table, mol, mol.getBondAtom(0, bond),
                mol.getBondAtom(1, bond));
    }

    /**
     * Returns the MMFF bond type of a bond on a molecule. This function
     * assumes that both atoms passed are valid and have a valid bond
     * between them.
     *  @param mol The molecule that contains the bond.
     *  @param atom1 The first atom of the bond.
     *  @param atom2 The second atom of the bond.
     *  @return The MMFF bond type.
     */
    public static int getType(Tables table, MMFFMolecule mol, int atom1,
                              int atom2) {
        int a1t = mol.getAtomType(atom1);
        int a2t = mol.getAtomType(atom2);
        int bond = mol.getBond(atom1, atom2);

        boolean notInAromaticRing = true;

        RingCollection rings = mol.getRingSet();
        for (int r=0; r<rings.getSize() && notInAromaticRing; r++) {
            for (int b : rings.getRingBonds(r)) {
                if (b == bond && mol.ringIsMMFFAromatic(r)) {
                    notInAromaticRing = false;
                    break;
                }
            }
        }

        return (mol.getBondOrder(bond) == 1)
            && notInAromaticRing
            && (table.atom.arom(a1t) && table.atom.arom(a2t)
            || table.atom.sbmb(a1t) && table.atom.sbmb(a2t)) ? 1 : 0;
    }
}
