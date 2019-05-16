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
 * Angle bending energy term class. This energy term represents the
 * angle bending energy associated with three bonded atoms A1--A2--A3
 * with an angle at A2.
 */
public class AngleBend implements EnergyTerm {
    public final int a1;
    public final int a2; // Central atom.
    public final int a3;
    public final boolean isLinear;

    public final double ka;     // Force constant.
    public final double theta0; // Ideal angle.

    /**
     * Construct a new angle bend energy term.
     *  @param table The tables parameter object.
     *  @param mol The molecule.
     *  @param a1 Index of atom 1 in mol.
     *  @param a2 Index of atom 2 (the central atom) in mol.
     *  @param a3 Index of atom 3 in mol.
     */
    public AngleBend(Tables table, MMFFMolecule mol, int a1, int a2,
					 int a3) {
        this.a1 = a1;
        this.a2 = a2;
        this.a3 = a3;
        isLinear = table.atom.linear(mol.getAtomType(a2));

        theta0 = table.angle.theta(mol, a1, a2, a3);
        ka = table.angle.ka(mol, a1, a2, a3);
    }

    /**
     * Calculates the angle energy.
     *  @param pos The atoms current positions array.
     *  @return The energy.
     */
    @Override
    public double getEnergy(double[] pos) {
        double theta = new Vector3(pos, a2, a1).angle(new Vector3(pos, a2, a3));
        double angle = Math.toDegrees(theta) - theta0;

        final double cb = -0.006981317;
        final double c2 = Constants.MDYNE_A_TO_KCAL_MOL * Constants.DEG2RAD
            * Constants.DEG2RAD;

        // isLinear is a property of the central atom and can be found in the
        // prop table.
        if (isLinear)
            return Constants.MDYNE_A_TO_KCAL_MOL*ka*(1.0 + Math.cos(theta));

        return 0.5*c2*ka*angle*angle*(1.0 + cb*angle);
    }

    /**
     * Calculates the gradient and adds it to the gradients array.
     *  @param pos The atoms current positions array.
     *  @param grad the atoms current gradients array.
     */
    @Override
    public void getGradient(double[] pos, double[] grad) {
        Vector3 r0 = new Vector3(pos, a2, a1).normalise();
        Vector3 r1 = new Vector3(pos, a2, a3).normalise();

        double dist0 = new Vector3(pos, a2, a1).length();
        double dist1 = new Vector3(pos, a2, a3).length();

        double cosTheta = r0.cosAngle(r1);

        double sinThetaSq = 1.0 - cosTheta*cosTheta;
        double sinTheta = 1.0e-8;
        if (sinThetaSq > 0.0)
            sinTheta = Math.sqrt(sinThetaSq);

        double angleTerm = Constants.RAD2DEG * Math.acos(cosTheta) - theta0;
        double cb = -0.006981317;
        double c2 = Constants.MDYNE_A_TO_KCAL_MOL * Constants.DEG2RAD
            * Constants.DEG2RAD;

        double dE_dTheta = Constants.RAD2DEG*c2*ka*angleTerm
            * (1.0 + 1.5*cb*angleTerm);
        if (isLinear)
            dE_dTheta = -Constants.MDYNE_A_TO_KCAL_MOL * ka * sinTheta;

        double dCos_dS[] = new double[]{
            1.0/dist0*(r1.x - cosTheta*r0.x),
            1.0/dist0*(r1.y - cosTheta*r0.y),
            1.0/dist0*(r1.z - cosTheta*r0.z),
            1.0/dist1*(r0.x - cosTheta*r1.x),
            1.0/dist1*(r0.y - cosTheta*r1.y),
            1.0/dist1*(r0.z - cosTheta*r1.z)
        };

        grad[3*a1    ] += dE_dTheta*dCos_dS[0]/(-sinTheta);
        grad[3*a1 + 1] += dE_dTheta*dCos_dS[1]/(-sinTheta);
        grad[3*a1 + 2] += dE_dTheta*dCos_dS[2]/(-sinTheta);

        grad[3*a2    ] += dE_dTheta*(-dCos_dS[0] - dCos_dS[3])/(-sinTheta);
        grad[3*a2 + 1] += dE_dTheta*(-dCos_dS[1] - dCos_dS[4])/(-sinTheta);
        grad[3*a2 + 2] += dE_dTheta*(-dCos_dS[2] - dCos_dS[5])/(-sinTheta);

        grad[3*a3    ] += dE_dTheta*dCos_dS[3]/(-sinTheta);
        grad[3*a3 + 1] += dE_dTheta*dCos_dS[4]/(-sinTheta);
        grad[3*a3 + 2] += dE_dTheta*dCos_dS[5]/(-sinTheta);
    }

    /**
     * Helper function that builds a list of AngleBends for a molecule.
     *  @param t The tables object.
     *  @param mol The molecule to generate angles for.
     *  @return Am array of AngleBends.
     */
    public static List<AngleBend> findIn(Tables t, MMFFMolecule mol) {
        ArrayList<AngleBend> angles = new ArrayList<AngleBend>();

        for (int atom=0; atom<mol.getAllAtoms(); atom++) {
            if (mol.getAllConnAtoms(atom) > 1) {
                for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                    int nbr_i = mol.getConnAtom(atom, i);
                    for (int k=i+1; k<mol.getAllConnAtoms(atom); k++) {
                        int nbr_k = mol.getConnAtom(atom, k);
                        angles.add(new AngleBend(t, mol, nbr_i, atom, nbr_k));
                    }
                }
            }
        }

        return angles;
    }
}
