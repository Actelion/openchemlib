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
 * Nonbonded van der Waals energy term class. This energy term represents
 * the van der Waals interaction between two atoms A1..A2 which are in a
 * 1,X (X > 3) relationship. A cutoff (default: 100.0 angstrom) can be set
 * to skip computation of van der Waals interactions between atoms
 * separated by distances larger than the cutoff.
 */
public class VanDerWaals implements EnergyTerm {
    public final int a1t;
    public final int a2t;
    public final double rstar_ij;
    public final double well_depth;
    public final int a1;
    public final int a2;

    public final char da1;
    public final char da2;

    /**
     * Construct a new van der Waals energy term.
     *  @param table The tables parameter object.
     *  @param mol The molecule.
     *  @param a1 Index of atom 1 in mol.
     *  @param a2 Index of atom 2 (the central atom) in mol.
     */
    public VanDerWaals(Tables table, MMFFMolecule mol, int a1, int a2) {
        a1t = mol.getAtomType(a1);
        a2t = mol.getAtomType(a2);
        this.a1 = a1;
        this.a2 = a2;

        double rs = minimum(table);
        double wd = wellDepth(table, rs);
        da1 = table.vdws.da(a1t);
        da2 = table.vdws.da(a2t);

        if ((da1 == 'D' && da2 == 'A') || (da1 == 'A' && da2 == 'D')) {
            rs = rs * table.vdws.darad;
            wd *= table.vdws.daeps;
        }

        rstar_ij = rs;
        well_depth = wd;
    }

    private double minimum(Tables table) {
        double rs1 = table.vdws.r_star(a1t);
        double rs2 = table.vdws.r_star(a2t);
        char da1 = table.vdws.da(a1t);
        char da2 = table.vdws.da(a2t);
        double gamma_ij = (rs1 - rs2) / (rs1 + rs2);

        return (0.5*(rs1 + rs2)*(1.0 + (((da1 == 'D') || (da2 == 'D')) ? 0.0
            : table.vdws.b*(1.0 - Math.exp(-(table.vdws.beta)*gamma_ij*gamma_ij)))));
    }

    private double wellDepth(Tables table, double rs) {
        double gi1 = table.vdws.g_i(a1t);
        double gi2 = table.vdws.g_i(a2t);
        double alpha1 = table.vdws.alpha_i(a1t);
        double alpha2 = table.vdws.alpha_i(a2t);
        double ni1 = table.vdws.n_i(a1t);
        double ni2 = table.vdws.n_i(a2t);

        double rstar_ij2 = rs * rs;
        double c4 = 181.16;

        return (c4*gi1*gi2*alpha1*alpha2
            / ((Math.sqrt(alpha1 / ni1) + Math.sqrt(alpha2 / ni2))
            * rstar_ij2 * rstar_ij2 * rstar_ij2));
    }


    /**
     * Calculates the van der Waals energy.
     *  @param pos The atoms current positions array.
     *  @return The energy.
     */
    @Override
    public double getEnergy(double[] pos) {
        final double dist = new Vector3(pos, a1, a2).length();
        final double vdw1 = 1.07;
        final double vdw1m1 = vdw1 - 1.0;
        final double vdw2 = 1.12;
        final double vdw2m1 = vdw2 - 1.0;
        final double dist2 = dist * dist;
        final double dist7 = dist2 * dist2 * dist2 * dist;
        final double aTerm = vdw1 * rstar_ij / (dist + vdw1m1 * rstar_ij);
        final double aTerm2 = aTerm * aTerm;
        final double aTerm7 = aTerm2 * aTerm2 * aTerm2 * aTerm;
        final double rstar_ij2 = rstar_ij * rstar_ij;
        final double rstar_ij7 = rstar_ij2 * rstar_ij2 * rstar_ij2 * rstar_ij;
        final double bTerm = vdw2*rstar_ij7 / (dist7 + vdw2m1*rstar_ij7) - 2.0;
        return aTerm7 * bTerm * well_depth;
    }

    /**
     * Calculates the gradient and adds it to the gradients array.
     *  @param pos The atoms current positions array.
     *  @param grad the atoms current gradients array.
     */
    @Override
    public void getGradient(double[] pos, double[] grad) {
        final double vdw1 = 1.07;
        final double vdw1m1 = vdw1 - 1.0;
        final double vdw2 = 1.12;
        final double vdw2m1 = vdw2 - 1.0;
        final double vdw2t7 = vdw2 * 7.0;
        double dist = new Vector3(pos, a1, a2).length();
        double q = dist / rstar_ij;
        double q2 = q * q;
        double q6 = q2 * q2 * q2;
        double q7 = q6 * q;
        double q7pvdw2m1 = q7 + vdw2m1;
        double t = vdw1 / (q + vdw1 - 1.0);
        double t2 = t * t;
        double t7 = t2 * t2 * t2 * t;
        double dE_dr = well_depth / rstar_ij
                * t7 * (-vdw2t7 * q6 / (q7pvdw2m1 * q7pvdw2m1)
                + ((-vdw2t7 / q7pvdw2m1 + 14.0) / (q + vdw1m1)));

        for (int i=0; i<3; i++) {
            double dGrad = 0.01 * rstar_ij;
            if (dist > 0.0)
                dGrad = dE_dr*(pos[3*a1+i] - pos[3*a2+i])/dist;

            grad[3*a1+i] += dGrad;
            grad[3*a2+i] -= dGrad;
        }
    }

    /**
     * Finds all van der Waals energy terms in the current molecule.
     *  @param table The tables object.
     *  @param mol The molecule to search for van der Waals forces.
     *  @param sep The separations table for molecule mol.
     *  @return The van der Waals energy terms for this molecule.
     */
    public static List<VanDerWaals> findIn(Tables table, MMFFMolecule mol,
            Separation sep, double nonbondedCutoff) {
        ArrayList<VanDerWaals> vdws = new ArrayList<VanDerWaals>();

        for (int i=0; i<mol.getAllAtoms(); i++) {
            for (int j=i+1; j<mol.getAllAtoms(); j++) {
                SortedPair atompair = new SortedPair(i, j);

                Separation.Relation relation = sep.get(atompair);

                if ((relation == Separation.Relation.ONE_FOUR
                        || relation == Separation.Relation.ONE_X)
                        && new Vector3(mol, i, j).length() < nonbondedCutoff) {
                    vdws.add(new VanDerWaals(table, mol, i, j));
                }
            }
        }

        return vdws;
    }
}
