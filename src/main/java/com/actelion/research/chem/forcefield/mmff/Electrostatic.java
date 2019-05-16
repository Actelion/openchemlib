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
 * Nonbonded electrostatic energy term class. This energy term represents
 * the energy from the electrostatic attraction/repulsion between any two
 * atoms A1..A2 which are in a 1,X relationship (X > 3; interactions where
 * X = 4 are scaled by 0.75). A cutoff (default: 100.0 angstrom) can be
 * set to skip computation of electrostatic interactions between atoms
 * separated by distances larger than the cutoff.
 * A constant or distance-dependent dielectric model can be used (default:
 * constant).
 */
public class Electrostatic implements EnergyTerm {
    public final MMFFMolecule mol;
    public final int a1;
    public final int a2;
    public final Separation.Relation rel;
    public final double charge_term;
    public final boolean distModel;

    /**
     * Construct a new electrostatic energy term.
     *  @param mol The molecule.
     *  @param a1 Index of atom 1 in mol.
     *  @param a2 Index of atom 2 in mol.
     *  @param rel The 1,X relationship between the two atoms (e.g., 1,2;
     *      1,3; 1,4; 1,X).
     *  @param chge1 The charge of atom 1.
     *  @param chge2 The charge of atom 2.
     *  @param distModel The distance model to use, true for "distance"
     *      and false for "constant".
     *  @param dielConst The dielectric constant.
     */
    public Electrostatic(MMFFMolecule mol, int a1, int a2,
                         Separation.Relation rel, double chge1, double chge2,
                         boolean distModel, double dielConst) {
        this.mol = mol;
        this.a1 = a1;
        this.a2 = a2;
        this.rel = rel;
        this.distModel = distModel;
        charge_term = chge1 * chge2 / dielConst;
    }

    /**
     * Constructor with default values for distModel and dielConst.
     *  @param mol The molecule.
     *  @param a1 Index of atom 1 in mol.
     *  @param a2 Index of atom 2 in mol.
     *  @param rel The separation relation between the two atoms (how many
     *      degrees of separation are between them).
     *  @param chge1 The charge of atom 1.
     *  @param chge2 The charge of atom 2.
     */
    public Electrostatic(MMFFMolecule mol, int a1, int a2,
                         Separation.Relation rel, double chge1, double chge2) {
        this(mol, a1, a2, rel, chge1, chge2, false, 1.0);
    }

    /**
     * Calculates the electrostatic energy.
     *  @return The energy.
     */
    @Override
    public double getEnergy(double[] pos) {
        double dist = new Vector3(pos, a1, a2).length();
        double corr_dist = dist + 0.05;
        double diel = 332.0716;

        if (distModel)
            corr_dist *= corr_dist;

        return diel * charge_term / corr_dist *
                (rel == Separation.Relation.ONE_FOUR ? 0.75 : 1.0);
    }

    /**
     * Calculates the gradient and adds it to the gradients array.
     *  @param pos The atoms current positions array.
     *  @param grad the atoms current gradients array.
     */
    @Override
    public void getGradient(double[] pos, double[] grad) {
        double dist = new Vector3(pos, a1, a2).length();
        double corr_dist = dist + 0.05;

        corr_dist *= (distModel ? corr_dist * corr_dist : corr_dist);

        double dE_dr = -332.0716 * (distModel ? 2.0 : 1.0)
            * charge_term / corr_dist
            * (rel == Separation.Relation.ONE_FOUR ? 0.75 : 1.0);

        for (int i=0; i<3; i++) {
            double dGrad = 0.02;
            if (dist > 0.0)
                dGrad = dE_dr * (pos[3*a1+i] - pos[3*a2+i]) / dist;

            grad[3*a1+i] += dGrad;
            grad[3*a2+i] -= dGrad;
        }
    }

    /**
     * Finds all Electrostatic energy terms in the current molecule.
     *  @param table The parameter tables to use.
     *  @param mol The molecule to search for Electrostatic forces.
     *  @param sep The separations table for molecule mol.
     *  @param nonbondedCutoff The nonbonded energy terms cutoff distance.
     *      Electrostatic energy terms between atoms further than this are
     *      not added.
     *  @param dielModel The dielectric model to use. True for 'distance'
     *      and False for 'linear'.
     *  @param dielConst The dielectric constant.
     *  @return The Electrostatic energy terms for this molecule.
     */
    public static List<Electrostatic> findIn(Tables table,
                                             MMFFMolecule mol, Separation sep, double nonbondedCutoff,
                                             boolean dielModel, double dielConst) {
        ArrayList<Electrostatic> eles = new ArrayList<Electrostatic>();
        double[] charges = com.actelion.research.chem.forcefield.mmff.type.Charge.getCharges(table, mol);

        for (int i=0; i<mol.getAllAtoms(); i++) {
            for (int j=0; j<i+1; j++) {

                Separation.Relation relation = sep.get(new SortedPair(i, j));
                if ((relation == Separation.Relation.ONE_FOUR
                        || relation == Separation.Relation.ONE_X)
                        && Math.abs(charges[i]) > 0.00001
                        && Math.abs(charges[j]) > 0.00001) {
                    if (new Vector3(mol, i, j).length() < nonbondedCutoff) {
                        eles.add(new Electrostatic(mol, i, j, relation,
                                    charges[i], charges[j], dielModel,
                                    dielConst));
                    }
                }
            }
        }

        return eles;
    }

    /**
     * Overloaded wrapper function for 'findIn' which sets default values
     * for the nonbonded cutoff, dielectric model and dielectric constant.
     *  @param table The parameter tables to use.
     *  @param mol The molecule to search for Electrostatic forces.
     *  @param sep The separations table for molecule mol.
     *  @return The Electrostatic energy terms for this molecule.
     */
    public static List<Electrostatic> findIn(Tables table,
                                             MMFFMolecule mol, Separation sep) {
        return findIn(table, mol, sep, 100.0, false, 1.0);
    }
}
