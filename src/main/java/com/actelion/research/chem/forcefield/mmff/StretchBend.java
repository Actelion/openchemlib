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
 * Stretch bending energy term class. This energy term represents the
 * bond stretching-angle bending energy associated with three bonded
 * atoms A1--A2--A3 with an angle at A2.
 */
public class StretchBend implements EnergyTerm {
    public final MMFFMolecule mol;
    public final int a1;
    public final int a2; // Central atom.
    public final int a3;

    public final double theta0;
    public final double kba_ijk;
    public final double kba_kji;
    public final double r0i;
    public final double r0k;

    /**
     * Constructs a new stretch bend object. Note that it only represents
     * the stretch bend energy in the IJK direction. A separate stretch
     * bend object should be created to represent the KJI direction with
     * a1 and a3 swapped.
     *  @param mol The molecule containing the atoms.
     *  @param a1 Atom 1 (atom i).
     *  @param a2 Atom 2, the central atom (atom j).
     *  @param a3 Atom 3 (atom k).
     */
    public StretchBend(Tables table, MMFFMolecule mol, int a1, int a2,
					   int a3) {
        this.mol = mol;
        this.a1 = a1;
        this.a2 = a2;
        this.a3 = a3;

        theta0 = table.angle.theta(mol, a1, a2, a3);
        r0i = table.bond.r0(mol, a1, a2);
        r0k = table.bond.r0(mol, a3, a2);
        kba_ijk = table.stbn.kba(mol, a1, a2, a3);
        kba_kji = table.stbn.kba(mol, a3, a2, a1);
    }

    /**
     * Calculates the stretch bend energy.
     *  @param pos The atoms current positions array.
     *  @return The energy.
     */
    @Override
    public double getEnergy(double[] pos) {
        double dist1 = new Vector3(pos, a2, a1).length();
        double dist2 = new Vector3(pos, a2, a3).length();
        double theta = new Vector3(pos,a2,a1).angle(new Vector3(pos,a2,a3));
        double factor = Constants.MDYNE_A_TO_KCAL_MOL * Constants.DEG2RAD
            * (Math.toDegrees(theta) - theta0);

        return factor*(dist1 - r0i)*kba_ijk + factor*(dist2 - r0k)*kba_kji;
    }

    /**
     * Calculates the gradient and adds it to the gradients array.
     *  @param pos The atoms current positions array.
     *  @param grad the atoms current gradients array.
     */
    @Override
    public void getGradient(double[] pos, double[] grad) {
        double dist1 = new Vector3(pos, a2, a1).length();
        double dist2 = new Vector3(pos, a2, a3).length();

        Vector3 p12 = new Vector3(pos, a2, a1).normalise();
        Vector3 p32 = new Vector3(pos, a2, a3).normalise();

        final double c5 = Constants.MDYNE_A_TO_KCAL_MOL * Constants.DEG2RAD;
        double cosTheta = p12.dot(p32);
        double sinThetaSq = 1.0 - cosTheta*cosTheta;
        double sinTheta = Math.max(sinThetaSq > 0.0
                ? Math.sqrt(sinThetaSq) : 0.0, 1.0e-8);
        double angleTerm = Constants.RAD2DEG * Math.acos(cosTheta) - theta0;
        double distTerm = Constants.RAD2DEG
                * (kba_ijk * (dist1 - r0i)
                +  kba_kji * (dist2 - r0k));

        double dCos_dS1 = 1.0 / dist1 * (p32.x - cosTheta * p12.x);
        double dCos_dS2 = 1.0 / dist1 * (p32.y - cosTheta * p12.y);
        double dCos_dS3 = 1.0 / dist1 * (p32.z - cosTheta * p12.z);

        double dCos_dS4 = 1.0 / dist2 * (p12.x - cosTheta * p32.x);
        double dCos_dS5 = 1.0 / dist2 * (p12.y - cosTheta * p32.y);
        double dCos_dS6 = 1.0 / dist2 * (p12.z - cosTheta * p32.z);

        grad[3*a1+0] += c5 * (p12.x * kba_ijk
            * angleTerm + dCos_dS1 / (-sinTheta) * distTerm);
        grad[3*a1+1] += c5 * (p12.y * kba_ijk
            * angleTerm + dCos_dS2 / (-sinTheta) * distTerm);
        grad[3*a1+2] += c5 * (p12.z * kba_ijk
            * angleTerm + dCos_dS3 / (-sinTheta) * distTerm);

        grad[3*a2+0] += c5 * ((-p12.x * kba_ijk
            - p32.x * kba_kji) * angleTerm
            + (-dCos_dS1 - dCos_dS4) / (-sinTheta) * distTerm);
        grad[3*a2+1] += c5 * ((-p12.y * kba_ijk
            - p32.y * kba_kji) * angleTerm
            + (-dCos_dS2 - dCos_dS5) / (-sinTheta) * distTerm);
        grad[3*a2+2] += c5 * ((-p12.z * kba_ijk
            - p32.z * kba_kji) * angleTerm
            + (-dCos_dS3 - dCos_dS6) / (-sinTheta) * distTerm);

        grad[3*a3+0] += c5 * (p32.x * kba_kji
            * angleTerm + dCos_dS4 / (-sinTheta) * distTerm);
        grad[3*a3+1] += c5 * (p32.y * kba_kji
            * angleTerm + dCos_dS5 / (-sinTheta) * distTerm);
        grad[3*a3+2] += c5 * (p32.z * kba_kji
            * angleTerm + dCos_dS6 / (-sinTheta) * distTerm);
    }

    /**
     * Checks if a Stretch Bend is valid
     */
    public static boolean valid(Tables table, MMFFMolecule mol, int a1,
								int a2, int a3) {
        int a2t = mol.getAtomType(a2);

        // Fail if the central atom is linear or the two neighbours are the
        // same atom.
        if (table.atom.linear(a2t) || a1 == a3)
            return false;

        // Fail if the three atoms do not form an angle.
        if (mol.getBond(a1, a2) == -1 || mol.getBond(a2, a3) == -1)
            return false;

        return true;
    }


    /**
     * Helper function that builds a list of StretchBends for a molecule.
     *  @param t The tables object.
     *  @param mol The molecule to generate angles for.
     *  @return An array of AngleBends.
     */
    public static List<StretchBend> findIn(Tables t, MMFFMolecule mol) {
        ArrayList<StretchBend> stbn = new ArrayList<StretchBend>();

        for (int atom=0; atom<mol.getAllAtoms(); atom++) {
            int atomt = mol.getAtomType(atom);

            if (mol.getAllConnAtoms(atom) <= 1 && t.atom.linear(atomt))
                continue;

            for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                int nbr_i = mol.getConnAtom(atom, i);
                for (int k=0; k<mol.getAllConnAtoms(atom); k++) {
                    int nbr_k = mol.getConnAtom(atom, k);

                    // Only add half of the bends.
                    if (nbr_i > nbr_k)
                        continue;

                    if (StretchBend.valid(t, mol, nbr_i, atom, nbr_k)) {
                        StretchBend bend = new StretchBend(t, mol, nbr_i,
                                atom, nbr_k);

                        // Don't bother adding any stretch bends whose kba
                        // would be 0.
                        if (Math.abs(bend.kba_ijk) > 0.0001
                                || Math.abs(bend.kba_kji) > 0.0001) {
                            stbn.add(bend);
                        }
                    }
                }
            }
        }

        return stbn;
    }
}
