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
