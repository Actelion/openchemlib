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

package mmff.type;

import com.actelion.research.chem.RingCollection;
import java.util.HashSet;
import mmff.MMFFMolecule;
import mmff.Tables;

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
