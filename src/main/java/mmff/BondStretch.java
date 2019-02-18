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

import java.util.ArrayList;
import java.util.List;

/**
 * Bond stretching energy term class. This energy term represents the
 * bond stretching energy associated with two bonded atoms A1--A2.
 */
public class BondStretch implements EnergyTerm {
    public final int a1;
    public final int a2;
    public final double kb; // Force constant.
    public final double r0; // Ideal bond length.

    /**
     * Creates a new bond stretch given a force field and a bond. This is
     * a wrapper constructor.
     *  @param table The tables parameter object.
     *  @param mol The molecule.
     *  @param bond The bond index.
     */
    public BondStretch(Tables table, MMFFMolecule mol, int bond) {
        this(table, mol, mol.getBondAtom(0, bond), mol.getBondAtom(1, bond));
    }

    /**
     * Creates a new bond stretch given a force field and two bonded
     * atoms.
     *  @param table The tables parameter object.
     *  @param mol The molecule.
     *  @param a1 Atom 1 index.
     *  @param a2 Atom 2 index.
     */
    public BondStretch(Tables table, MMFFMolecule mol, int a1, int a2) {
        this.a1 = a1;
        this.a2 = a2;

        r0 = table.bond.r0(mol, a1, a2);
        kb = table.bond.kb(mol, a1, a2);
    }

    /**
     * Calculates the bond stretch energy.
     *  @param pos The atoms current positions array.
     *  @return The energy.
     */
    public double getEnergy(double[] pos) {
        final double c1 = 143.9325;
        final double cs = -2.0;
        final double c3 = 7.0 / 12.0;
        final double dist = new Vector3(pos,a1).distance(new Vector3(pos,a2));
        final double diff = (dist - r0)*(dist - r0);
        return (0.5*c1*kb*diff * (1.0 + cs*(dist - r0) + c3*cs*cs*diff));
    }

    /**
     * Calculates the gradient and adds it to the gradients array.
     *  @param pos The atoms current positions array.
     *  @param grad the atoms current gradients array.
     */
    public void getGradient(double[] pos, double[] grad) {
        final double cs = -2.0;
        final double c1 = Constants.MDYNE_A_TO_KCAL_MOL;
        final double c3 = 7.0/12.0;
        final double dist = new Vector3(pos,a1).distance(new Vector3(pos,a2));

        double distTerm = dist - r0;
        double dE_dr = c1*kb*distTerm *
            (1.0 + 1.5*cs*distTerm + 2.0*c3*cs*cs*distTerm*distTerm);

        if (dist > 0.0) {
            for (int i=0; i<3; i++) {
                grad[3*a1 + i] += dE_dr*(pos[3*a1 + i] - pos[3*a2 + i])/dist;
                grad[3*a2 + i] -= dE_dr*(pos[3*a1 + i] - pos[3*a2 + i])/dist;
            }
        }
    }

    /**
     * Finds all bond stretch energy terms in the current molecule.
     *  @param t The tables parameter object.
     *  @param mol The molecule.
     *  @return The bond stretch energy terms for this molecule.
     */
    public static List<BondStretch> findIn(Tables t, MMFFMolecule mol) {
        List<BondStretch> bstretches = new ArrayList<BondStretch>();

        for (int b=0; b<mol.getAllBonds(); b++)
            bstretches.add(new BondStretch(t, mol, b));

        return bstretches;
    }
}
