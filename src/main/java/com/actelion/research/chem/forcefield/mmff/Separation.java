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

import java.util.Hashtable;

/**
 * The Separation class is an efficient storage of an adjacency matrix
 * representing how many degrees of separation there are between any two
 * atoms. The Relations enum marks the relationship between two atoms. A
 * HashTable is used to store the relations between atoms. Atom pairs with
 * ONE_X relations are not stored (this is assumed to be the default case).
 */
public class Separation {
    /**
     * Relation class shows the relationship between two atoms A1 and A2.
     */
    public enum Relation {
        ONE_ONE,   // A1 is A2, an atom has a ONE_ONE relation with itself.
        ONE_TWO,   // A1 is directly bonded to A2.
        ONE_THREE, // A1 is indirectly bonded to A2 via one intermediary atom.
        ONE_FOUR,  // A1 is indirectly bonded to A2 via two intermediary atoms.
        ONE_X;     // (Default) A1 and A2 are separated by more than two atoms.
    }

    public Hashtable<SortedPair, Relation> table =
        new Hashtable<SortedPair, Relation>();

    /**
     * Constructs a new separation table for a molecule.
     *  @param mol The molecule that the table will describe.
     */
    public Separation(MMFFMolecule mol) {
        for (int atom=0; atom<mol.getAllAtoms(); atom++) {
            table.put(new SortedPair(atom, atom), Relation.ONE_ONE);

            for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                int nbr1 = mol.getConnAtom(atom, i);
                SortedPair keytwo = new SortedPair(atom, nbr1);

                table.put(keytwo, Relation.ONE_TWO);

                for (int j=0; j<mol.getAllConnAtoms(nbr1); j++) {
                    int nbr2 = mol.getConnAtom(nbr1, j);
                    SortedPair keythree = new SortedPair(atom, nbr2);

                    if (!table.containsKey(keythree)
                            || table.get(keythree) == Relation.ONE_FOUR)
                        table.put(keythree, Relation.ONE_THREE);

                    for (int k=0; k<mol.getAllConnAtoms(nbr2); k++) {
                        int nbr3 = mol.getConnAtom(nbr2, k);
                        SortedPair keyfour = new SortedPair(atom, nbr3);

                        if (!table.containsKey(keyfour))
                            table.put(keyfour, Relation.ONE_FOUR);
                    }
                }
            }
        }
    }

    /**
     * Returns the relation of a given pair of atoms. If no key was found
     * that is an indication that there is a ONE_X relationship between
     * the two atoms.
     *  @param key The sorted pair of atoms.
     *  @return The separation relationship between the two atoms.
     */
    public Relation get(SortedPair key) {
        return table.get(key) != null ? table.get(key) : Relation.ONE_X;
    }

    public Relation get(int a1, int a2) {
        return get(new SortedPair(a1, a2));
    }
}
