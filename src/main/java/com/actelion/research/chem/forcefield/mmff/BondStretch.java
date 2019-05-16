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
