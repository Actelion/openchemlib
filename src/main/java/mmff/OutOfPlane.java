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
 * Out of plane energy term class. This energy term represents the
 * out of plane energy associated with four atoms:
 *
 *             A1
 *            /
 *      A2--AC
 *            \
 *             A3
 *
 * Where AC is the central atom and A1, A2 and A3 are each connected to
 * AC.
 */
public class OutOfPlane implements EnergyTerm {
    int ac;
    int a1;
    int a2;
    int a3;
    double koop;

    /**
     * Construct a new out of plane energy term.
     *  @param table The tables parameter object.
     *  @param mol The molecule.
     *  @param a1 Index of atom 1 in mol.
     *  @param a2 Index of atom 2 (the central atom) in mol.
     *  @param a3 Index of atom 3 in mol.
     */
    public OutOfPlane(Tables table, MMFFMolecule mol, int ac, int a1,
					  int a2, int a3) {
        this.ac = ac; // j
        this.a1 = a1; // i
        this.a2 = a2; // k
        this.a3 = a3; // l

        koop = table.oop.getKoop(mol, ac, a1, a2, a3);
    }

    public double getKoop() {
        return koop;
    }

    /**
     * Calculates the out of plane energy.
     *  @param pos The atoms current positions array.
     *  @return The energy.
     */
    @Override
    public double getEnergy(double[] pos) {
        Vector3 rji = new Vector3(pos, ac, a1).normalise();
        Vector3 rjk = new Vector3(pos, ac, a2).normalise();
        Vector3 rjl = new Vector3(pos, ac, a3).normalise();
        Vector3 n = rji.cross(rjk).normalise();

        double chi = Constants.RAD2DEG * Math.asin(n.dot(rjl));
        double c2 = Constants.MDYNE_A_TO_KCAL_MOL * Constants.DEG2RAD
            * Constants.DEG2RAD;
        return 0.5 * c2 * koop * chi * chi;
    }

    /**
     * Calculates the gradient and adds it to the gradients array.
     *  @param pos The atoms current positions array.
     *  @param grad the atoms current gradients array.
     */
    @Override
    public void getGradient(double[] pos, double[] grad) {
        Vector3 rji = new Vector3(pos, ac, a1);
        Vector3 rjk = new Vector3(pos, ac, a2);
        Vector3 rjl = new Vector3(pos, ac, a3);

        final double dji = rji.length();
        final double djk = rjk.length();
        final double djl = rjl.length();

        rji = rji.normalise();
        rjk = rjk.normalise();
        rjl = rjl.normalise();
        Vector3 n = rji.negate().cross(rjk).normalise();

        final double c2 = Constants.MDYNE_A_TO_KCAL_MOL * Constants.DEG2RAD
            * Constants.DEG2RAD;

        double sinChi = rjl.dot(n);
        double cosChiSq = 1.0 - sinChi*sinChi;
        double cosChi = Math.max(cosChiSq > 0.0
                ? Math.sqrt(cosChiSq) : 0.0, 1.0e-8);
        double chi = Constants.RAD2DEG * Math.asin(sinChi);
        double cosTheta = rji.dot(rjk);
        double sinThetaSq = Math.max(1.0 - cosTheta * cosTheta, 1.0e-8);
        double sinTheta = Math.max(sinThetaSq > 0.0
                ? Math.sqrt(sinThetaSq) : 0.0, 1.0e-8);
        double dE_dChi = Constants.RAD2DEG * c2 * koop * chi;

        Vector3 t1 = rjl.cross(rjk);
        Vector3 t2 = rji.cross(rjl);
        Vector3 t3 = rjk.cross(rji);

        double term1 = cosChi * sinTheta;
        double term2 = sinChi / (cosChi * sinThetaSq);

        double[] tg1 = new double[]{
            (t1.x/term1 - (rji.x - rjk.x*cosTheta) * term2) / dji,
            (t1.y/term1 - (rji.y - rjk.y*cosTheta) * term2) / dji,
            (t1.z/term1 - (rji.z - rjk.z*cosTheta) * term2) / dji};
        double[] tg3 = new double[]{
            (t2.x/term1 - (rjk.x - rji.x*cosTheta) * term2) / djk,
            (t2.y/term1 - (rjk.y - rji.y*cosTheta) * term2) / djk,
            (t2.z/term1 - (rjk.z - rji.z*cosTheta) * term2) / djk};
        double[] tg4 = new double[]{
            (t3.x/term1 - rjl.x*sinChi/cosChi) / djl,
            (t3.y/term1 - rjl.y*sinChi/cosChi) / djl,
            (t3.z/term1 - rjl.z*sinChi/cosChi) / djl};

        for (int i=0; i<3; i++) {
            grad[3*a1 + i] +=  dE_dChi *  tg1[i];
            grad[3*ac + i] += -dE_dChi * (tg1[i] + tg3[i] + tg4[i]);
            grad[3*a2 + i] +=  dE_dChi *  tg3[i];
            grad[3*a3 + i] +=  dE_dChi *  tg4[i];
        }
    }


    /**
     * Checks if this OutOfPlane's atoms are the same as a given set of 4
     * atoms, it checks all possible permutations of a1, a2 and a3.
     *  @param ac The central atom.
     *  @param a1 Neighbouring atom 1.
     *  @param a2 Neighbouring atom 2.
     *  @param a3 Neighbouring atom 3.
     *  @return True if this OutOfPlane has the same permutation of
     *      neighbouring atoms as those provided.
     */
    public boolean equals(int ac, int a1, int a2, int a3) {
        return this.ac == ac
            && (this.a1 == a1 && this.a2 == a2 && this.a3 == a3
            ||  this.a1 == a1 && this.a2 == a3 && this.a3 == a2
            ||  this.a1 == a2 && this.a2 == a1 && this.a3 == a3
            ||  this.a1 == a2 && this.a2 == a3 && this.a3 == a1
            ||  this.a1 == a3 && this.a2 == a1 && this.a3 == a2
            ||  this.a1 == a3 && this.a2 == a2 && this.a3 == a1);
    }


    /**
     * Checks if this OutOfPlane term is exactly equal to a given set of
     * four atoms.
     *  @param ac The central atom.
     *  @param a1 Neighbouring atom 1.
     *  @param a2 Neighbouring atom 2.
     *  @param a3 Neighbouring atom 3.
     *  @return True if the out of plane atoms are the same as the
     *      provided atoms.
     */
    public boolean exactly(int ac, int a1, int a2, int a3) {
        return this.ac == ac
            && this.a1 == a1 && this.a2 == a2 && this.a3 == a3;
    }


    /**
     * Finds all out of plane angles in the current molecule.
     *  @param t The tables parameter object.
     *  @param mol The molecule to search in.
     *  @return An array of OutOfPlane angles.
     */
    public static List<OutOfPlane> findIn(Tables t, MMFFMolecule mol) {
        ArrayList<OutOfPlane> oops = new ArrayList<OutOfPlane>();

        for (int ac=0; ac<mol.getAllAtoms(); ac++) {
            if (mol.getAllConnAtoms(ac) != 3)
                continue;

            int a1 = mol.getConnAtom(ac, 0);
            int a2 = mol.getConnAtom(ac, 1);
            int a3 = mol.getConnAtom(ac, 2);

            oops.add(new OutOfPlane(t, mol, ac, a1, a2, a3));
            oops.add(new OutOfPlane(t, mol, ac, a1, a3, a2));
            oops.add(new OutOfPlane(t, mol, ac, a2, a3, a1));
        }

        return oops;
    }

    public static List<OutOfPlane> findIn(MMFFMolecule mol) {
        return findIn(ForceFieldMMFF94.table(ForceFieldMMFF94.MMFF94), mol);
    }
}
