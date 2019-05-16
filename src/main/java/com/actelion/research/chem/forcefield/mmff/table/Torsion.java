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

import com.actelion.research.chem.forcefield.mmff.Csv;
import com.actelion.research.chem.forcefield.mmff.MMFFMolecule;
import com.actelion.research.chem.forcefield.mmff.PeriodicTable;
import com.actelion.research.chem.forcefield.mmff.Search;
import com.actelion.research.chem.forcefield.mmff.Tables;

public final class Torsion implements com.actelion.research.chem.forcefield.mmff.Searchable {
    // indices: 0, 1, 2, 3, 4
    // Sort order: 2, 3, 1, 4, 0
    // angle type, i type, j type, k type, l type, V1, V2, V3
    private final Object[][] table;
    private final Tables t;

    public Torsion(Tables t, String csvpath) {
        table = Csv.readFile(csvpath);
        this.t = t;
    }

    @Override
    public int get(int row, int col) {
        return ((Number) table[row][col]).intValue();
    }

    public static <T> T s(T a, T b) {
        return a;
    }

    @Override
    public int length() {
        return table.length;
    }


    public final class Kb {
        public double v1;
        public double v2;
        public double v3;

        public Kb(int index) {
            v1 = ((Number) table[index][5]).doubleValue();
            v2 = ((Number) table[index][6]).doubleValue();
            v3 = ((Number) table[index][7]).doubleValue();
        }

        public Kb() {
            v1 = 0.0;
            v2 = 0.0;
            v3 = 0.0;
        }

        public Kb(double v1, double v2, double v3) {
            this.v1 = v1;
            this.v2 = v2;
            this.v3 = v3;
        }

        public String toString() {
            return v1+","+v2+","+v3;
        }
    }

    public Kb get(int index) {
        return new Kb(index);
    }

    /**
     * Returns the index of a row for a given molecule and four connected
     * atoms which form a torsion angle.
     *  @param a1t Atom 1 type.
     *  @param a2t Atom 2 type.
     *  @param a3t Atom 3 type.
     *  @param a4t Atom 4 type.
     *  @param tort Torsion type.
     *  @return The index in the torsion table, or -1 if no suitable entry
     *      was found.
     */
    public int index(int a1t, int a2t, int a3t, int a4t, int tort) {
        return Search.binary(new int[]{2,3,1,4,0},
                new int[]{a2t,a3t,a1t,a4t,tort}, this);
    }

    /**
     * Gets the force constants associated with a torsion angle.
     */
    public Kb getForceConstants(MMFFMolecule mol, int a1, int a2, int a3, int a4) {
        int a1t = mol.getAtomType(a1);
        int a2t = mol.getAtomType(a2);
        int a3t = mol.getAtomType(a3);
        int a4t = mol.getAtomType(a4);
        int tort = com.actelion.research.chem.forcefield.mmff.type.Torsion.getType(t, mol, a1, a2, a3, a4);

        int tort1 = tort > 10 ? tort/10 : tort;
        int tort2 = tort > 10 ? tort - tort1*10 : 0;

        int idx = -1;

        int iter = 0;
        int iWildCard = 0;
        int lWildCard = 0;
        int canTorType = tort1;
        int maxIter = 5;

        while ((iter < maxIter && (idx == -1 || maxIter == 4))
                || (iter == 4 && tort1 == 5 && tort2 > 0)) {

            if (maxIter == 5 && iter == 4) {
                maxIter = 4;
                iter = 0;
                canTorType = tort2;
            }

            if (iter == 1) {
                iWildCard = 1;
                lWildCard = 3;
            } else if (iter == 2) {
                iWildCard = 3;
                lWildCard = 1;
            } else {
                iWildCard = iter;
                lWildCard = iter;
            }

            int canIAtomType = t.def.table[a1t-1][Math.min(iWildCard+1, 4)];
            int canJAtomType = a2t;
            int canKAtomType = a3t;
            int canLAtomType = t.def.table[a4t-1][Math.min(lWildCard+1, 4)];

            if (canJAtomType > canKAtomType) {
                canKAtomType = s(canJAtomType, canJAtomType = canKAtomType);
                canLAtomType = s(canIAtomType, canIAtomType = canLAtomType);
            } else if (canJAtomType == canKAtomType && canIAtomType > canLAtomType) {
                canLAtomType = s(canIAtomType, canIAtomType = canLAtomType);
            }

            idx = index(canIAtomType, canJAtomType, canKAtomType,
                    canLAtomType, canTorType);

            if (idx != -1 && maxIter == 4)
                break;

            iter++;
        }

        // A match was found in the torsion table.
        if (idx >= 0)
            return new Kb(idx);

        // Here on is the empirical rules.
        int bond = mol.getBond(a2, a3);

        double[] U = new double[]{0.0, 0.0};
        double[] V = new double[]{0.0, 0.0};
        double[] W = new double[]{0.0, 0.0};
        int[] atno = new int[]{mol.getAtomicNo(a2), mol.getAtomicNo(a3)};
        double N_jk = (t.atom.crd(a2t) - 1)*(t.atom.crd(a3t) - 1);

        for (int i=0; i<2; i++) {
            switch (atno[i]) {
                // carbon
                case 6:
                    U[i] = 2.0;
                    V[i] = 2.12;
                    break;

                // nitrogen
                case 7:
                    U[i] = 2.0;
                    V[i] = 1.5;
                    break;

                // oxygen
                case 8:
                    U[i] = 2.0;
                    V[i] = 0.2;
                    W[i] = 2.0;
                    break;

                // silicon
                case 14:
                    U[i] = 1.25;
                    V[i] = 1.22;
                    break;

                // phosphorus
                case 15:
                    U[i] = 1.25;
                    V[i] = 2.40;
                    break;

                // sulfur
                case 16:
                    U[i] = 1.25;
                    V[i] = 0.49;
                    W[i] = 8.0;
                    break;
            }
        }

        // --- Rule A ---
        if (t.atom.linear(a2t) || t.atom.linear(a3t))
            return new Kb(0.0, 0.0, 0.0);

        // --- Rule B ---
        if (t.atom.arom(a2t) && t.atom.arom(a3t) && mol.isAromaticBond(bond)) {

            double beta = (t.atom.val(a2t) == 3 && t.atom.val(a3t) == 4)
                || (t.atom.val(a2t) == 4 && t.atom.val(a3t) == 3) ? 3.0 : 6.0;
            double pi_jk = t.atom.pilp(a2t) == 0 && t.atom.pilp(a3t) == 0
                ? 0.5 : 0.3;

            return new Kb(0.0, beta * pi_jk * Math.sqrt(U[0] * U[1]), 0.0);
        }

        // --- Rule C ---
        if (mol.getBondOrder(bond) == 2) {
            double beta = 6.0;
            double pi_jk = t.atom.mltb(a2t) == 2 && t.atom.mltb(a3t) == 2
                ? 1.0 : 0.4;
            return new Kb(0.0, beta * pi_jk * Math.sqrt(U[0] * U[1]), 0.0);
        }

        // --- Rule D ---
        if (t.atom.crd(a2t) == 4 && t.atom.crd(a3t) == 4)
            return new Kb(0.0, 0.0, Math.sqrt(V[0] * V[1]) / N_jk);

        // --- Rule E ---
        if (t.atom.crd(a2t) == 4 && t.atom.crd(a3t) != 4) {
            if (((t.atom.crd(a3t) == 3) && (((t.atom.val(a3t) == 4)
                    || (t.atom.val(a3t) == 34)) || t.atom.mltb(a3t) > 0))
                    || ((t.atom.crd(a3t) == 2)
                    && ((t.atom.val(a3t) == 3) || t.atom.mltb(a3t) > 0)))
                return new Kb();
            else
                return new Kb(0.0, 0.0, Math.sqrt(V[0] * V[1]) / N_jk);
        }

        // --- Rule F ---
        if ((t.atom.crd(a3t) == 4) && (t.atom.crd(a2t) != 4)) {
            if (((t.atom.crd(a2t) == 3) && (((t.atom.val(a2t) == 4)
                    || (t.atom.val(a2t) == 34)) || t.atom.mltb(a2t) > 0))
                    || ((t.atom.crd(a2t) == 2)
                    && ((t.atom.val(a2t) == 3) || t.atom.mltb(a2t) > 0)))
                return new Kb();
            else
                return new Kb(0.0, 0.0, Math.sqrt(V[0] * V[1]) / N_jk);
        }

        // --- Rule G ---
        if (((mol.getBondOrder(bond) == 1)
                && t.atom.mltb(a2t) > 0 && t.atom.mltb(a3t) > 0)
                || (t.atom.mltb(a2t) > 0 && t.atom.pilp(a3t) > 0)
                || (t.atom.pilp(a2t) > 0 && t.atom.mltb(a3t) > 0)) {

            // --- Case 1 ---
            if (t.atom.pilp(a2t) > 0 && t.atom.pilp(a3t) > 0)
                return new Kb();

            // --- Case 2 ---
            if (t.atom.pilp(a2t) > 0 && t.atom.mltb(a3t) > 0) {
                double beta = 6.0;
                double pi_jk = 0.0;
                if (t.atom.mltb(a2t) == 1)
                    pi_jk = 0.5;
                else if (PeriodicTable.row(atno[0]) == 2 && PeriodicTable.row(atno[1]) == 2)
                    pi_jk = 0.3;
                else if (PeriodicTable.row(atno[0]) != 2 || PeriodicTable.row(atno[1]) != 2)
                    pi_jk = 0.15;

                return new Kb(0.0, beta * pi_jk * Math.sqrt(U[0] * U[1]), 0.0);
            }

            // --- Case 3 ---
            if (t.atom.pilp(a3t) > 0 && t.atom.mltb(a2t) > 0) {
                double beta = 6.0;
                double pi_jk = 0.0;
                if (t.atom.mltb(a3t) == 1)
                    pi_jk = 0.5;
                else if (PeriodicTable.row(atno[0]) == 2 && PeriodicTable.row(atno[1]) == 2)
                    pi_jk = 0.3;
                else if (PeriodicTable.row(atno[0]) != 2 || PeriodicTable.row(atno[1]) != 2)
                    pi_jk = 0.15;

                return new Kb(0.0, beta * pi_jk * Math.sqrt(U[0] * U[1]), 0.0);
            }

            // --- Case 4 ---
            if ((t.atom.mltb(a2t) == 1 || t.atom.mltb(a3t) == 1)
                    && (atno[0] != 6 || atno[1] != 6))
                return new Kb(0.0, 6.0 * 0.4 * Math.sqrt(U[0] * U[1]) , 0.0);

            // --- Case 5 ---
            return new Kb(0.0, 6.0 * 0.15 * Math.sqrt(U[0] * U[1]), 0.0);
        }

        // --- Rule H ---
        if ((atno[0] == 8 || atno[0] == 16) && (atno[1] == 8 || atno[1] == 16))
            return new Kb(0.0, -Math.sqrt(W[0] * W[1]), 0.0);

        return new Kb(0.0, 0.0, Math.sqrt(V[0] * V[1]) / N_jk);
    }
}
