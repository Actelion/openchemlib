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

package com.actelion.research.chem.forcefield.mmff.table;

import com.actelion.research.chem.forcefield.mmff.Constants;
import com.actelion.research.chem.forcefield.mmff.Csv;
import com.actelion.research.chem.forcefield.mmff.MMFFMolecule;
import com.actelion.research.chem.forcefield.mmff.Search;
import com.actelion.research.chem.forcefield.mmff.Tables;

/**
 * Angle table, corresponds to the MMFFANG.PAR parameters table provided in
 * the MMFF literature. This table provides parameters for angle bending
 * interactions. This table is loaded from a CSV file as it is too large to
 * hard code in to the class file.
 */
public final class Angle implements com.actelion.research.chem.forcefield.mmff.Searchable {
    // Sort order: 4, 2, 1, 3
    // angle type, i type, j type, k type, ka, theta0
    private final Object[][] table;
    private final Tables t;

    public Angle(Tables t, String csvpath) {
        table = Csv.readFile(csvpath);
        this.t = t;
    }

    @Override
    public int get(int row, int col) {
        return ((Number) table[row][col]).intValue();
    }

    @Override
    public int length() {
        return table.length;
    }

    /**
     * Returns 'ka' the force constant for a given index in the table.
     *  @param index The row index in the table.
     *  @return The force constant.
     */
    public double ka(int index) {
        return ((Number) table[index][4]).doubleValue();
    }

    /**
     * Returns 'theta0' the equilibrium angle for a given index in the
     * table.
     *  @param index The row index in the table.
     *  @return The force constant.
     */
    public double theta(int index) {
        return ((Number) table[index][5]).doubleValue();
    }

    /**
     * Returns the index of a row for a given molecule and three connected
     * atoms which form an angle.
     *  @param mol The molecule that the atoms are in.
     *  @param a1 Atom 1.
     *  @param a2 Atom 2 (the central atom).
     *  @param a3 Atom 3.
     *  @return The index in the angle table, or -1 if no suitable entry was
     *      found.
     */
    public int index(MMFFMolecule mol, int a1, int a2,
                     int a3) {
        int a1t = mol.getAtomType(a1);
        int a2t = mol.getAtomType(a2);
        int a3t = mol.getAtomType(a3);
        int angt = com.actelion.research.chem.forcefield.mmff.type.Angle.getType(t, mol, a1, a2, a3);

        // Loop over the atom types, searching for an entry in the angles
        // table using equivalent types if no entry is found.
        int a1tf, a3tf, index = -1;
        for (int i=0; i<5 && index < 0; i++) {
            a1tf = t.def.table[a1t-1][i];
            a3tf = t.def.table[a3t-1][i];

            if (a1tf > a3tf)
                a1tf = Search.s(a3tf, a3tf = a1tf);

            index = Search.binary(new int[]{2,1,3,0},
                    new int[]{a2t,a1tf,a3tf,angt}, this);
        }
        return index;
    }

    /**
     * Returns 'theta0' the ideal angle given a molecule and three
     * connected atoms which form an angle.
     *  @param mol The molecule that the atoms are in.
     *  @param a1 Atom 1.
     *  @param a2 Atom 2 (the central atom).
     *  @param a3 Atom 3.
     *  @return The value of 'theta' from the angle table or calculated
     *      empirically.
     */
    public double theta(MMFFMolecule mol, int a1, int a2,
                        int a3) {
        int index = index(mol, a1, a2, a3);

        // If we didn't find an index.
        if (index < 0)
            return getEmpiricalTheta0(mol, a1, a2, a3);
        else
            // If we did find an index.
            return theta(index);
    }

    /**
     * Returns 'ka' given a molecule and three connected atoms which form
     * an angle.
     *  @param mol The molecule that the atoms are in.
     *  @param a1 Atom 1.
     *  @param a2 Atom 2 (the central atom).
     *  @param a3 Atom 3.
     *  @return The value of 'ka' from the angle table or calculated
     *      empirically.
     */
    public double ka(MMFFMolecule mol, int a1, int a2,
                     int a3) {
        int index = index(mol, a1, a2, a3);

        // If we didn't find an index.
        if (index < 0)
            return getEmpiricalKa(mol, a1, a2, a3,
                    getEmpiricalTheta0(mol, a1, a2, a3));
        else {
            // If we did find an index.
            if (Math.abs(ka(index)) < 0.001)
                return getEmpiricalKa(mol, a1, a2, a3, theta(index));
            else
                return ka(index);
        }
    }

    /**
     * Calculates an empirical value for theta0 (the ideal angle).
     *  @param mol The molecule that the atoms are in.
     *  @param a1 Atom 1.
     *  @param a2 Atom 2 (the central atom).
     *  @param a3 Atom 3.
     *  @return Theta0 in degrees.
     */
    private double getEmpiricalTheta0(MMFFMolecule mol, int a1, int a2, int a3) {
        if (com.actelion.research.chem.forcefield.mmff.type.Angle.inRingOfSize(mol, a1, a2, a3, 3))
            return 60.0;
        else if (com.actelion.research.chem.forcefield.mmff.type.Angle.inRingOfSize(mol, a1, a2, a3, 4))
            return 90.0;

        int a2t = mol.getAtomType(a2);
        switch (t.atom.crd(a2t)) {
            case 2:
                if (mol.getAtomicNo(a2) == 8)
                    return 105.0;
                else if (t.atom.linear(a2t))
                    return 180.0;
                break;
            case 3:
                if (t.atom.val(a2t) == 3 && t.atom.mltb(a2t) == 0) {
                    if (mol.getAtomicNo(a2) == 7)
                        return 107.0;
                    else
                        return 92.0;
                }
                break;
            case 4:
                return 109.45;
        }
        return 120.0;
    }

    /**
     * Calculates an empirical value for ka (the force constant).
     *  @param mol The molecule that the atoms are in.
     *  @param a1 Atom 1.
     *  @param a2 Atom 2 (the central atom).
     *  @param a3 Atom 3.
     *  @param theta0 The ideal angle theta0 as given either from the
     *      lookup table or calculated empirically.
     *  @return The empirical ka.
     */
    private double getEmpiricalKa(MMFFMolecule mol, int a1, int a2, int a3,
                                  double theta0) {
        double[] z = new double[]{0.0, 0.0, 0.0};
        double[] c = new double[]{0.0, 0.0, 0.0};
        int[] atno = new int[]{mol.getAtomicNo(a1), mol.getAtomicNo(a2),
            mol.getAtomicNo(a3)};
        double beta = 1.75;

        for (int i=0; i<3; i++) {
            switch (atno[i]) {
                // Hydrogen
                case 1:
                    z[i] = 1.395;
                    break;
                // Carbon
                case 6:
                    z[i] = 2.494;
                    c[i] = 1.016;
                    break;
                // Nitrogen
                case 7:
                    z[i] = 2.711;
                    c[i] = 1.113;
                    break;
                // Oxygen
                case 8:
                    z[i] = 3.045;
                    c[i] = 1.337;
                    break;
                // Fluorine
                case 9:
                    z[i] = 2.847;
                    break;
                // Silicon
                case 14:
                    z[i] = 2.350;
                    c[i] = 0.811;
                    break;
                // Phosphorus
                case 15:
                    z[i] = 2.350;
                    c[i] = 1.068;
                    break;
                // Sulfur
                case 16:
                    z[i] = 2.980;
                    c[i] = 1.249;
                    break;
                // Chlorine
                case 17:
                    z[i] = 2.909;
                    c[i] = 1.078;
                    break;
                // Bromine
                case 35:
                    z[i] = 3.017;
                    break;
                // Iodine
                case 53:
                    z[i] = 3.086;
                    break;
            }
        }

        double r0_ij = t.bond.r0(mol, a1, a2);
        double r0_jk = t.bond.r0(mol, a2, a3);
        double D = (r0_ij - r0_jk)*(r0_ij - r0_jk)
            /((r0_ij + r0_jk)*(r0_ij + r0_jk));
        double theta0_rad = Constants.DEG2RAD*theta0;

        if (com.actelion.research.chem.forcefield.mmff.type.Angle.inRingOfSize(mol, a1, a2, a3, 4))
            beta *= 0.85;
        else if (com.actelion.research.chem.forcefield.mmff.type.Angle.inRingOfSize(mol, a1, a2, a3, 3))
            beta *= 0.05;

        return beta * z[0] * c[1] * z[2]
        / ((r0_ij + r0_jk) * theta0_rad * theta0_rad * Math.exp(2.0 * D));
    }
}
