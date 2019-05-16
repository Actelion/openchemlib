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
 * Torsional Angle energy term class. This energy term represents the
 * energy associated with the torsional angle formed by four atoms
 * A1..A4:
 *
 *      A1
 *        \
 *         A2--A3
 *               \
 *                A4
 *
 */
public class TorsionAngle implements EnergyTerm {
    public final int a1;
    public final int a2;
    public final int a3;
    public final int a4;

    public final int a1t;
    public final int a2t;
    public final int a3t;
    public final int a4t;

    public final double v1;
    public final double v2;
    public final double v3;

    /**
     * Construct a new torsion angle energy term.
     *  @param table The tables parameter object.
     *  @param mol The molecule.
     *  @param a1 Index of atom 1 in mol.
     *  @param a2 Index of atom 2 in mol.
     *  @param a3 Index of atom 3 in mol.
     *  @param a4 Index of atom 4 in mol.
     */
    public TorsionAngle(Tables table, MMFFMolecule mol, int a1, int a2,
						int a3, int a4) {
        this.a1 = a1;
        this.a2 = a2;
        this.a3 = a3;
        this.a4 = a4;

        a1t = mol.getAtomType(a1);
        a2t = mol.getAtomType(a2);
        a3t = mol.getAtomType(a3);
        a4t = mol.getAtomType(a4);

        com.actelion.research.chem.forcefield.mmff.table.Torsion.Kb kbs = table.torsion.getForceConstants(mol,
                a1, a2, a3, a4);
        v1 = kbs.v1;
        v2 = kbs.v2;
        v3 = kbs.v3;
    }

    /**
     * Calculates the torsional energy.
     *  @param pos The atoms current positions array.
     *  @return The energy.
     */
    @Override
    public double getEnergy(double[] pos) {
        Vector3 r1 = new Vector3(pos, a1, a2);
        Vector3 r2 = new Vector3(pos, a3, a2);
        Vector3 r3 = new Vector3(pos, a2, a3);
        Vector3 r4 = new Vector3(pos, a4, a3);

        Vector3 t1 = r1.cross(r2);
        Vector3 t2 = r3.cross(r4);
        double cosPhi = t1.cosAngle(t2);

        double cos2Phi = 2.0 * cosPhi * cosPhi - 1.0;
        double cos3Phi = cosPhi * (2.0 * cos2Phi - 1.0);

        return 0.5 * (v1*(1.0 + cosPhi)
                    + v2*(1.0 - cos2Phi)
                    + v3*(1.0 + cos3Phi));
    }

    /**
     * Calculates the gradient and adds it to the gradients array.
     *  @param pos The atoms current positions array.
     *  @param grad the atoms current gradients array.
     */
    @Override
    public void getGradient(double[] pos, double[] grad) {
        Vector3[] r = new Vector3[]{
            new Vector3(pos, a2, a1),
            new Vector3(pos, a2, a3),
            new Vector3(pos, a3, a2),
            new Vector3(pos, a3, a4)
        };

        Vector3[] t = new Vector3[]{
            r[0].cross(r[1]),
            r[2].cross(r[3])
        };

        double[] d = new double[]{
            t[0].length(),
            t[1].length()
        };

        if (Math.abs(d[0]) < 0.00001 || Math.abs(d[1]) < 0.00001)
            return;

        t[0] = t[0].normalise();
        t[1] = t[1].normalise();

        double cosPhi = t[0].dot(t[1]);
        double sinPhiSq = 1.0 - cosPhi * cosPhi;
        double sinPhi = ((sinPhiSq > 0.0) ? Math.sqrt(sinPhiSq) : 0.0);
        double sin2Phi = 2.0 * sinPhi * cosPhi;
        double sin3Phi = 3.0 * sinPhi - 4.0 * sinPhi * sinPhiSq;
        double dE_dPhi = 0.5 * (-(v1) * sinPhi + 2.0 * v2 * sin2Phi
                - 3.0 * v3 * sin3Phi);
        double sinTerm = -dE_dPhi * (Math.abs(sinPhi) < 0.00001
                ? (1.0 / cosPhi) : (1.0 / sinPhi));

        double[] dCos_dT = new double[]{
            1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
            1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
            1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
            1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
            1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
            1.0 / d[1] * (t[0].z - cosPhi * t[1].z)
        };

        grad[3*a1+0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
        grad[3*a1+1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
        grad[3*a1+2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);

        grad[3*a2+0] += sinTerm * (dCos_dT[1] * (r[1].z - r[0].z)
                + dCos_dT[2] * (r[0].y - r[1].y)
                + dCos_dT[4] * (-r[3].z)
                + dCos_dT[5] * (r[3].y));
        grad[3*a2+1] += sinTerm * (dCos_dT[0] * (r[0].z - r[1].z)
                + dCos_dT[2] * (r[1].x - r[0].x)
                + dCos_dT[3] * (r[3].z)
                + dCos_dT[5] * (-r[3].x));
        grad[3*a2+2] += sinTerm * (dCos_dT[0] * (r[1].y - r[0].y)
                + dCos_dT[1] * (r[0].x - r[1].x)
                + dCos_dT[3] * (-r[3].y)
                + dCos_dT[4] * (r[3].x));

        grad[3*a3+0] += sinTerm * (dCos_dT[1] * (r[0].z)
                + dCos_dT[2] * (-r[0].y)
                + dCos_dT[4] * (r[3].z - r[2].z)
                + dCos_dT[5] * (r[2].y - r[3].y));
        grad[3*a3+1] += sinTerm * (dCos_dT[0] * (-r[0].z)
                + dCos_dT[2] * (r[0].x)
                + dCos_dT[3] * (r[2].z - r[3].z)
                + dCos_dT[5] * (r[3].x - r[2].x));
        grad[3*a3+2] += sinTerm * (dCos_dT[0] * (r[0].y)
                + dCos_dT[1] * (-r[0].x)
                + dCos_dT[3] * (r[3].y - r[2].y)
                + dCos_dT[4] * (r[2].x - r[3].x));

        grad[3*a4+0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
        grad[3*a4+1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
        grad[3*a4+2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
    }

    /**
     * Checks that at least one of the constants is non-zero.
     *  @return True if any constant is non-zero, false otherwise.
     */
    public boolean nonZero() {
        return Math.abs(v1) > 0.001
            || Math.abs(v2) > 0.001
            || Math.abs(v3) > 0.001;
    }

    /**
     * Helper function that builds a list of TorsionAngles for a molecule.
     *  @param t The tables object.
     *  @param mol The molecule to generate torsions for.
     *  @return Am array of TorsionAngle.
     */
    public static List<TorsionAngle> findIn(Tables t, MMFFMolecule mol) {
        ArrayList<TorsionAngle> tors = new ArrayList<TorsionAngle>();

        for (int a1=0; a1<mol.getAllAtoms(); a1++) {
            for (int j=0; j<mol.getAllConnAtoms(a1); j++) {
                int a2 = mol.getConnAtom(a1, j);
                for (int k=0; k<mol.getAllConnAtoms(a2); k++) {
                    int a3 = mol.getConnAtom(a2, k);

                    if (a1 == a3)
                        continue;

                    for (int l=0; l<mol.getAllConnAtoms(a3); l++) {
                        int a4 = mol.getConnAtom(a3, l);

                        if (a2 == a4 || a1 == a4)
                            continue;

                        if (a4 > a1) {
                            TorsionAngle tor
                                = new TorsionAngle(t, mol, a1, a2, a3, a4);

                            if (tor.nonZero())
                                tors.add(tor);
                        }
                    }
                }
            }
        }

        return tors;
    }
}
