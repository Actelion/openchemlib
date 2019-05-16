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

package com.actelion.research.chem.forcefield.mmff.type;

import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.forcefield.mmff.MMFFMolecule;
import com.actelion.research.chem.forcefield.mmff.RingBoolean;

import java.util.ArrayList;

/**
 * The Atom type class provides static functions to perform atom typing on
 * atoms in a Molecule. There are also several helper functions used in
 * Atom typing provided here.
 */
public class Atom {
    /**
     * Returns the MMFF type of an atom in a molecule.
     *  @param mol The molecule that the atom is in.
     *  @param atom The atom to assign a type for.
     *  @return The mmff atom type.
     */
    public static int getType(MMFFMolecule mol, int atom) {
        // Check if we have already cached the atom type, return it if we
        // have.
        if (mol.getAtomType(atom) > -1)
            return mol.getAtomType(atom);

        if (mol.getAtomicNo(atom) == 1)
            return getHydrogenType(mol, atom);
        else
            return getHeavyType(mol, atom);
    }

    /**
     * Returns the MMFF atom type for heavy (non-hydrogen) atoms.
     *  @param mol The molecule that the atom is in.
     *  @param atom The atom to assign a type for.
     *  @return The mmff atom type.
     */
    private static int getHeavyType(MMFFMolecule mol, int atom) {
        // Aromatic Atoms
        if (mol.isAromaticAtom(atom)) {
            // 5 member rings
            if (isInAromaticRingOfSize(mol, atom, 5)) {
                ArrayList<Integer> alphaHet = new ArrayList<Integer>(mol.getAllAtoms());
                ArrayList<Integer> betaHet = new ArrayList<Integer>(mol.getAllAtoms());

                boolean alphaBetaSameRing = false;
                boolean isAlphaOS = false;
                boolean isBetaOS = false;

                if (mol.getAtomicNo(atom) == 6 || mol.getAtomicNo(atom) == 7) {
                    for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                        int nbr = mol.getConnAtom(atom, i);

                        // Skip atoms not in the ring.
                        if (!isInAromaticRingOfSize(mol, nbr, 5))
                            continue;

                        if (inSameRing(mol, atom, nbr, 5)
                                && (mol.getAtomicNo(nbr) == 8
                                || mol.getAtomicNo(nbr) == 16
                                || mol.getAtomicNo(nbr) == 7
                                && degree(mol, nbr) == 3
                                && !isAtomNOxide(mol, nbr)))
                            alphaHet.add(nbr);

                        for (int j=0; j<mol.getAllConnAtoms(nbr); j++) {
                            int nbr2 = mol.getConnAtom(nbr, j);

                            if (nbr2 == atom)
                                continue;

                            if (!isInAromaticRingOfSize(mol, nbr2, 5))
                                continue;

                            if (inSameRing(mol, atom, nbr2, 5)
                                    && (mol.getAtomicNo(nbr2) == 8
                                    || mol.getAtomicNo(nbr2) == 16
                                    || mol.getAtomicNo(nbr2) == 7
                                    && degree(mol, nbr2) == 3
                                    && !isAtomNOxide(mol, nbr2)))
                                betaHet.add(nbr2);
                        }
                    }

                    for (int alpha : alphaHet)
                        if (mol.getAtomicNo(alpha) == 8 || mol.getAtomicNo(alpha) == 16) {
                            isAlphaOS = true;
                            break;
                        }

                    for (int beta : betaHet)
                        if (mol.getAtomicNo(beta) == 8 || mol.getAtomicNo(beta) == 16) {
                            isBetaOS = true;
                            break;
                        }

                    for (int i=0; i<alphaHet.size(); i++)
                        for (int j=i; j<betaHet.size(); j++)
                            if (inSameRing(mol, alphaHet.get(i), betaHet.get(j), 5)) {
                                alphaBetaSameRing = true;
                                break;
                            }
                }

                switch (mol.getAtomicNo(atom)) {
                    // ----- Carbon -----
                    case 6:
                        if (betaHet.isEmpty()) {
                            int nN = 0;
                            int nFormalCharge = 0;
                            int nIn5Ring = 0;
                            int nIn6Ring = 0;

                            for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                                int nbr = mol.getConnAtom(atom, i);

                                if (mol.getAtomicNo(nbr) == 7
                                        && degree(mol, nbr) == 3) {
                                    nN++;

                                    if (mol.getAtomCharge(nbr) > 0
                                            && !isAtomNOxide(mol, nbr))
                                        nFormalCharge++;

                                    if (isInAromaticRingOfSize(mol, nbr, 5))
                                        nIn5Ring++;

                                    if (isInAromaticRingOfSize(mol, nbr, 6))
                                        nIn6Ring++;
                                }
                            }

                            // CIM+ (Aromatic carbon between N's in
                            // imidazolium)
                            if ((nN == 2 && nIn5Ring > 0 || nN == 3 && nIn5Ring == 2)
                                    && nFormalCharge > 0 && nIn6Ring == 0)
                                return 80;
                        }

                        // If we have both or none of alpha/beta heteroatoms.
                        if (!(alphaHet.isEmpty() ^ betaHet.isEmpty())) {
                            boolean surroundedByBenzeneC = true;
                            boolean surroundedByArom = true;

                            for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                                int nbr = mol.getConnAtom(atom, i);

                                if (mol.getAtomicNo(nbr) != 6
                                        || !inRingOfSize(mol, nbr, 6))
                                    surroundedByBenzeneC = false;

                                if (inSameRing(mol, atom, nbr, 5)
                                        && !mol.isAromaticAtom(nbr))
                                    surroundedByArom = false;
                            }

                            // C5 (General Carbon in 5 membered heteroaromatic
                            // ring)
                            if (alphaHet.isEmpty() && betaHet.isEmpty()
                                    && !surroundedByBenzeneC && surroundedByArom)
                                return 78;

                            // C5 (General Carbon in 5 membered heteroaromatic
                            // ring)
                            if (!alphaHet.isEmpty() && !alphaHet.isEmpty()
                                    && (!alphaBetaSameRing || !isAlphaOS && !isBetaOS))
                                return 78;
                        }

                        // C5A (Aromatic 5-ring C, alpha to N:, O: or S:)
                        if (!alphaHet.isEmpty() && (betaHet.isEmpty() || isAlphaOS))
                            return 63;

                        // C5B (Aromatic 5-ring C, alpha to N:, O: or S:)
                        if (!betaHet.isEmpty() && (alphaHet.isEmpty() || isBetaOS))
                            return 64;
                        break;

                    // ----- Nitrogen -----
                    case 7:
                        // N5AX (N-oxide nitrogen in 5-ring alpha position)
                        // N5BX (N-oxide nitrogen in 5-ring beta position)
                        // N5OX (N-oxide nitrogen in other 5-ring position)
                        if (isAtomNOxide(mol, atom))
                            return 82;

                        if (alphaHet.isEmpty() && betaHet.isEmpty()) {
                            // NPYL (Aromatic 5-ring nitrogen with pi lone
                            // pair)
                            if (degree(mol, atom) == 3)
                                return 39;

                            // N5M (Nitrogen in 5-ring aromatic anion)
                            return 76;
                        }

                        // NIM+ (Aromatic nitrogen in imidazolium)
                        // N5A+ (Positive nitrogen in 5-ring alpha position)
                        // N5B+ (Positive nitrogen in 5-ring beta position)
                        // N5+ (Positive nitrogen in other 5-ring position)
                        if (degree(mol, atom) == 3)
                            if (alphaHet.isEmpty() ^ betaHet.isEmpty())
                                return 81;

                        // N5A (Aromatic 5-ring N, alpha to N:, O: or S:)
                        if (!alphaHet.isEmpty() && (betaHet.isEmpty() || isAlphaOS))
                            return 65;

                        // N5B (Aromatic 5-ring N, beta to N:, O: or S:)
                        if (!betaHet.isEmpty() && (alphaHet.isEmpty() || isBetaOS))
                            return 66;

                        // N5 (General Nitrogen in 5-membered heteroaromatic
                        // ring)
                        if (!alphaHet.isEmpty() && !betaHet.isEmpty())
                            return 79;
                        break;

                    // ----- Oxygen -----
                    case 8:
                        // OFUR (Aromatic 5-ring Oxygen with pi lone pair)
                        return 59;

                    // ----- Sulfur -----
                    case 16:
                        // STHI (Aromatic 5-ring Sulfur with pi lone pair)
                        return 44;
                }
            }

            // 6 member rings
            if (isInAromaticRingOfSize(mol, atom, 6)) {
                switch (mol.getAtomicNo(atom)) {
                    // ----- Carbon -----
                    case 6:
                        // CB (Aromatic Carbon (eg in benzene))
                        return 37;

                    // ----- Nitrogen -----
                    case 7:
                        // NPOX (Pyridinium N-oxide Nitrogen)
                        if (isAtomNOxide(mol, atom))
                            return 69;

                        // NPD+ (Aromatic Nitrogen in pyridinium)
                        if (degree(mol, atom) == 3)
                            return 58;

                        // NPYD (Aromatic Nitrogen with sigma lone pair)
                        return 38;
                }
            }
        }

        // Aliphatic Heavy Atoms
        switch (mol.getAtomicNo(atom)) {
            // ----- Lithium -----
            case 3:
                // LI+ (Lithium cation)
                if (mol.getConnAtoms(atom) == 0)
                    return 92;
                break;

            // ----- Carbon -----
            case 6:
                // 4 Neighbours
                if (degree(mol, atom) == 4) {
                    // CR3R (Aliphatic Carbon in 3-member ring)
                    if (mol.getAtomRingSize(atom) == 3)
                        return 22;

                    // CR4R (Aliphatic Carbon in 4-member ring)
                    if (mol.getAtomRingSize(atom) == 4)
                        return 20;

                    // CR (Alkyl Carbon)
                    return 1;
                }

                // 3 Neighbours
                if (degree(mol, atom) == 3) {
                    int nN2 = 0;
                    int nN3 = 0;
                    int nO  = 0;
                    int nS  = 0;
                    int doubleBondedElement = 0;

                    for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                        int nbr = mol.getConnAtom(atom, i);

                        // Find a double-bonded element
                        if (mol.getBondOrder(mol.getBond(atom, nbr)) == 2)
                            doubleBondedElement = mol.getAtomicNo(nbr);

                        // Count the terminal Oxygen/Sulfur atoms.
                        if (degree(mol, nbr) == 1) {
                            if (mol.getAtomicNo(nbr) == 8)
                                nO++;
                            if (mol.getAtomicNo(nbr) == 16)
                                nS++;
                        } else if (mol.getAtomicNo(nbr) == 7) {
                            if (degree(mol, nbr) == 3)
                                nN3++;
                            else if (degree(mol, nbr) == 2 &&
                                    mol.getBondOrder(mol.getBond(atom, nbr)) == 2)
                                nN2++;
                        }
                    }

                    // If there are two or more Nitrogens with 3 neighbours
                    // each, and there are no nitrogens with two neighbours
                    // only, and Carbon is double-bonded to Nitrogen.
                    if (nN3 >= 2 && nN2 == 0 && doubleBondedElement == 7) {
                        // CNN+ (Carbon in +N=C-N: resonance structures)
                        // CGD+ (Guanidinium Carbon)
                        return 57;
                    }

                    // If there are two terminal Oxygen/Sulfur atoms.
                    if (nO == 2 || nS == 2) {
                        // CO2M (Carbon in carboxylate anion)
                        // CS2M (Carbon in thiocarboxylate anion)
                        return 41;
                    }

                    // If the Carbon is in a 4-member ring and is
                    // double-bonded to another Carbon.
                    if (mol.getAtomRingSize(atom) == 4
                            && doubleBondedElement == 6) {
                        // CR4E (Olefinic Carbon in 4-member ring)
                        return 30;
                    }

                    // If the Carbon is double bonded to Nitrogen, Oxygen,
                    // Phosphorus or Sulfur.
                    if (doubleBondedElement == 7 || doubleBondedElement == 8
                            || doubleBondedElement == 15 || doubleBondedElement == 16) {
                        // C=N  (Imine-atomType Carbon)
                        // CGD  (Guanidine Carbon)
                        // C=O  (Generic carbonyl Carbon)
                        // C=OR (Ketone or aldehyde carbonyl Carbon)
                        // C=ON (Amide carbonyl Carbon)
                        // COO  (Carboxylic acid or ester carbonyl Carbon)
                        // COON (Carbamate carbonyl Carbon)
                        // COOO (Carbonic acid or ester carbonyl function)
                        // C=OS (Thioester carbonyl Carbon, double bonded to O)
                        // C=P  (Carbon doubly bonded to P)
                        // C=S  (Thioester Carbon, double bonded to S)
                        // C=SN (Thioamide Carbon, double bonded to S)
                        // CSO2 (Carbon in >C=SO2)
                        // CS=O (Sulfinyl Carbon in >C=S=O)
                        // CSS  (Thiocarboxylic acid or ester Carbon)
                        return 3;
                    }

                    // Generic SP2 Carbon
                    // C=C (Vinylic Carbon)
                    // CSP2 (Generic SP2 Carbon)
                    return 2;
                }

                // 2 Neighbours
                if (degree(mol, atom) == 2) {
                    // CSP (Acetylenic Carbon), =C= (Allenic Carbon)
                    return 4;
                }

                if (degree(mol, atom) == 1) {
                    // C%- (Isonitrile Carbon)
                    return 60;
                }
                break;

            // ----- Nitrogen -----
            case 7:
                // Count the number of terminal Oxygens as well as terminal
                // Oxygens bonded to neighbouring Phosphorus or Sulfur atoms.
                int nTermOtoN = 0;
                boolean isSulfonamideN = false;
                for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                    int nbr = mol.getConnAtom(atom, i);

                    if (mol.getAtomicNo(nbr) == 8 && degree(mol, nbr) == 1)
                        nTermOtoN++;

                    if (mol.getOccupiedValence(atom) + mol.getImplicitHydrogens(atom) >= 3
                            && (mol.getAtomicNo(nbr) == 15 || mol.getAtomicNo(nbr) == 16)) {
                        int nOtoSP = 0;

                        for (int j=0; j<mol.getAllConnAtoms(nbr); j++) {
                            int nbr2 = mol.getConnAtom(nbr, j);
                            if (mol.getAtomicNo(nbr2) == 8 && degree(mol, nbr2) == 1)
                                nOtoSP++;
                        }

                        if (!isSulfonamideN)
                            isSulfonamideN = nOtoSP >= 2;
                    }
                }

                // 4 neighbours.
                if (degree(mol, atom) == 4) {
                    // N3OX (sp3-hybridized N-oxide Nitrogen)
                    if (isAtomNOxide(mol, atom))
                        return 68;

                    // NR+ (Quaternary Nitrogen)
                    return 34;
                }

                // 3 neighbours.
                if (degree(mol, atom) == 3) {
                    // Total bond order >= 4
                    if (mol.getOccupiedValence(atom) + mol.getImplicitHydrogens(atom) >= 4) {
                        boolean doubleCN = false;

                        for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                            int nbr = mol.getConnAtom(atom, i);

                            if (mol.getBondOrder(mol.getBond(atom, nbr)) == 2) {
                                doubleCN = mol.getAtomicNo(nbr) == 7
                                    || mol.getAtomicNo(nbr) == 6;

                                if (mol.getAtomicNo(nbr) == 6) {
                                    for (int j=0; doubleCN && j<mol.getAllConnAtoms(nbr); j++) {
                                        int nbr2 = mol.getConnAtom(nbr, j);

                                        if (nbr2 == atom)
                                            continue;

                                        doubleCN = !(mol.getAtomicNo(nbr2) == 7
                                            && degree(mol, nbr2) == 3);
                                    }
                                }
                            }
                        }

                        // N20X (sp2-hybridized N-oxide Nitrogen)
                        if (nTermOtoN == 1)
                            return 67;

                        // NO2 (Nitrogen in nitro group)
                        // NO3 (Nitrogen in nitrate group)
                        if (nTermOtoN >= 2)
                            return 45;

                        // N+=C (Iminium Nitrogen)
                        // N+=N (Positively charged Nitrogen doubly bonded to
                        // N)
                        if (doubleCN)
                            return 54;
                    }

                    // Total bond order >= 3
                    if (mol.getOccupiedValence(atom) + mol.getImplicitHydrogens(atom) >= 3) {
                        boolean isNCOorNCS    = false;
                        boolean isNCNplus     = false;
                        boolean isNGDplus     = false;
                        boolean isNNNorNNC    = false;
                        boolean isNbrC        = false;
                        boolean isNbrBenzeneC = false;
                        int EdoubledC = 0;
                        int EtripledC = 0;
                        int N2toC   = 0;
                        int N3toC   = 0;
                        int OtoC    = 0;
                        int StoC    = 0;

                        for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                            int nbr = mol.getConnAtom(atom, i);

                            // Neighbour is Carbon
                            if (mol.getAtomicNo(nbr) == 6) {
                                isNbrC = true;

                                if (mol.isAromaticAtom(nbr)
                                        && mol.getAtomRingSize(nbr) == 6)
                                    isNbrBenzeneC = true;

                                N2toC = 0;
                                N3toC = 0;
                                OtoC = 0;
                                StoC = 0;
                                int nFormalCharge = 0;
                                int inAromatic6Ring = 0;

                                for (int j =0; j<mol.getAllConnAtoms(nbr); j++) {
                                    int nbr2 = mol.getConnAtom(nbr, j);
                                    int bond = mol.getBond(nbr, nbr2);

                                    // Check if we have an Oxygen or Sulfur
                                    // double bonded to this Carbon.
                                    if (mol.getBondOrder(bond) == 2
                                            && (mol.getAtomicNo(nbr2) == 8
                                            || mol.getAtomicNo(nbr2) == 16))
                                        isNCOorNCS = true;

                                    // Check and identify if an element is
                                    // double bonded to this Carbon.
                                    if (mol.getBondOrder(bond) == 2
                                            || (mol.isAromaticBond(bond)
                                            && (mol.getAtomicNo(nbr2) == 6
                                            || (mol.getAtomicNo(nbr2) == 7
                                            && inRings(mol, nbr2) == 1))))
                                        EdoubledC = mol.getAtomicNo(nbr2);

                                    // Check and identify if an element is
                                    // triple bonded to this Carbon.
                                    if (mol.getBondOrder(bond) == 3)
                                        EtripledC = mol.getAtomicNo(nbr2);

                                    // If this Carbon is bonded to a Nitrogen
                                    // with 3 neighbours.
                                    if (mol.getAtomicNo(nbr2) == 7
                                            && degree(mol, nbr2) == 3) {
                                        if (mol.getAtomCharge(nbr2) == 1)
                                            nFormalCharge++;

                                        if (isInAromaticRingOfSize(mol, nbr, 6))
                                            inAromatic6Ring++;

                                        int OtoN3 = 0;
                                        for (int k=0; k<mol.getAllConnAtoms(nbr2); k++) {
                                            int nbr3 = mol.getConnAtom(nbr2, k);
                                            if (mol.getAtomicNo(nbr3) == 8)
                                                OtoN3++;
                                        }

                                        if (OtoN3 < 2)
                                            N3toC++;
                                    }

                                    // If this Carbon is bonded to a Nitrogen
                                    // with 2 neighbours with a double or
                                    // aromatic bond.
                                    if (mol.getAtomicNo(nbr2) == 7
                                            && degree(mol, nbr2) == 2
                                            && (mol.getBondOrder(bond) == 2
                                             || mol.isAromaticBond(bond)))
                                        N2toC++;

                                    // If this Carbon is bonded to an aromatic
                                    // atom.
                                    if (mol.isAromaticAtom(nbr2)) {
                                        if (mol.getAtomicNo(nbr2) == 8)
                                            OtoC++;

                                        if (mol.getAtomicNo(nbr2) == 16)
                                            StoC++;
                                    }
                                }


                                if (EdoubledC == 7) {
                                    if (N3toC == 2 && N2toC == 0
                                            && nFormalCharge > 0
                                            && inAromatic6Ring == 0
                                            && degree(mol, nbr) < 4)
                                        isNCNplus = true;

                                    if (N3toC == 3)
                                        isNGDplus = true;
                                }

                            // Neighbour is Nitrogen
                            } else if (mol.getAtomicNo(nbr) == 7) {
                                int NtoN = 0;
                                int OtoN = 0;
                                int StoN = 0;

                                for (int j=0; j<mol.getAllConnAtoms(nbr); j++) {
                                    int nbr2 = mol.getConnAtom(nbr, j);
                                    int bond = mol.getBond(nbr, nbr2);

                                    if (mol.getBondOrder(bond) == 2) {
                                        if (mol.getAtomicNo(nbr2) == 6) {
                                            for (int k=0; k<mol.getAllConnAtoms(nbr2); k++) {
                                                int nbr3 = mol.getConnAtom(nbr2, k);

                                                if (nbr3 == atom)
                                                    continue;

                                                if (mol.getAtomicNo(nbr3) == 7)
                                                    NtoN++;
                                                else if (mol.getAtomicNo(nbr3) == 8)
                                                    OtoN++;
                                                else if (mol.getAtomicNo(nbr3) == 16)
                                                    StoN++;
                                            }

                                            if (NtoN == 0 && OtoN == 0 && StoN == 0
                                                    && !isNbrBenzeneC)
                                                isNNNorNNC = true;
                                        }

                                        if (mol.getAtomicNo(nbr2) == 7
                                                && !isNbrBenzeneC)
                                            isNNNorNNC = true;
                                    }
                                }
                            }
                        }

                        // If ipso Nitrogen is bonded to Carbon.
                        if (isNbrC) {
                            if (EtripledC == 7)
                                isSulfonamideN = true;

                            // NCN+ (Either Nitrogen in N+=C-N)
                            if (isNCNplus)
                                return 55;

                            // NGD+ (Guanidinium Nitrogen)
                            if (isNGDplus)
                                return 56;

                            // NC=C (Enamine or aniline nitrogen, deloc. lp)
                            // NC=N (Nitrogen in N-C=N with deloc. lp)
                            // NC=P (Nitrogen in N-C=P with deloc. lp)
                            // NC%C (Nitrogen attached to C-C triple bond)
                            if (!isNCOorNCS && !isSulfonamideN
                                    && (OtoC == 0 && StoC == 0 && isNbrBenzeneC
                                    || EdoubledC == 6  || EdoubledC == 7
                                    || EdoubledC == 15 || EtripledC == 6))
                                return 40;
                        }


                        // NC=O (Amide nitrogen)
                        // NC=S (Thioamide nitrogen)
                        // NN=C (Nitrogen in N-N=C moiety with deloc. lp)
                        // NN=N (Nitrogen in N-N=N moiety with deloc. lp)
                        if (!isSulfonamideN && (isNCOorNCS || isNNNorNNC))
                            return 10;
                    }
                }

                // 2 neighbours.
                if (degree(mol, atom) == 2) {
                    // Total bond order = 4
                    if (mol.getOccupiedValence(atom) + mol.getImplicitHydrogens(atom) == 4) {
                        boolean isIsonitrile = false;
                        for (int i=0; !isIsonitrile && i<mol.getAllConnAtoms(atom); i++) {
                            int nbr = mol.getConnAtom(atom, i);
                            isIsonitrile = mol.getBondOrder(mol.getBond(atom, nbr)) == 3;
                        }

                        // NR% (Isonitrile Nitrogen)
                        if (isIsonitrile)
                            return 61;

                        // =N= (Central Nitrogen in C=N=N or N=N=N)
                        return 53;
                    }

                    // Total bond order = 3
                    if (mol.getOccupiedValence(atom) + mol.getImplicitHydrogens(atom) == 3) {
                        boolean isNitroso = false;
                        boolean isImineOrAzo = false;

                        for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                            int nbr = mol.getConnAtom(atom, i);
                            if (mol.getBondOrder(mol.getBond(atom, nbr)) == 2) {
                                isNitroso = mol.getAtomicNo(nbr) == 8
                                    && nTermOtoN == 1;

                                isImineOrAzo = mol.getAtomicNo(nbr) == 6
                                    || mol.getAtomicNo(nbr) == 7;
                            }
                        }

                        // N=O (Nitrogen in nitroso group)
                        if (isNitroso && !isImineOrAzo)
                            return 46;

                        // N=C (Imine Nitrogen)
                        // N=N (Azo-group Nitrogen)
                        if (isImineOrAzo)
                            return 9;
                    }

                    // Total bond order >= 2
                    if (mol.getOccupiedValence(atom) + mol.getImplicitHydrogens(atom) >= 2) {
                        boolean isNSO = false;
                        for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                            int nbr = mol.getConnAtom(atom, i);

                            if (mol.getAtomicNo(nbr) == 16) {
                                int nTermOtoS = 0;

                                for (int j=0; j<mol.getAllConnAtoms(nbr); j++) {
                                    int nbr2 = mol.getConnAtom(nbr, j);
                                    if (mol.getAtomicNo(nbr2) == 8
                                            && degree(mol, nbr2) == 1)
                                        nTermOtoS++;
                                }

                                isNSO = nTermOtoS == 1;
                            }
                        }

                        // NSO (Divalent nitrogen replacing monovalent O in SO2 group)
                        if (isNSO)
                            return 48;

                        // NM (Anionic divalent Nitrogen)
                        if (!isSulfonamideN)
                            return 62;
                    }
                }

                // NSO2 (Sulfonamide nitrogen)
                // NSO3 (Sulfonamide nitrogen)
                // NC%N (Nitrogen attached to cyano group)
                if (isSulfonamideN)
                    return 43;

                // 1 neighbour
                if (degree(mol, atom) == 1) {
                    boolean isNSP = false;
                    boolean isNAZT = false;

                    for (int i=0; !isNSP && !isNAZT && i<mol.getAllConnAtoms(atom); i++) {
                        int nbr = mol.getConnAtom(atom, i);

                        isNSP = mol.getBondOrder(mol.getBond(atom, nbr)) == 3;

                        if (mol.getAtomicNo(nbr) == 7 && degree(mol, nbr) == 2) {
                            for (int j=0; !isNAZT && j<mol.getAllConnAtoms(nbr); j++) {
                                int nbr2 = mol.getConnAtom(nbr, j);

                                isNAZT = mol.getAtomicNo(nbr2) == 7
                                    && degree(mol, nbr2) == 2
                                    || mol.getAtomicNo(nbr2) == 6
                                    && degree(mol, nbr2) == 3;
                            }

                        }
                    }

                    // NSP (Triple bonded Nitrogen)
                    if (isNSP)
                        return 42;

                    // NAZT (Terminal nitrogen in azido or diazo group)
                    if (isNAZT)
                        return 47;
                }

                // If nothing was matched
                // NR (Amine Nitrogen)
                return 8;

            // ----- Oxygen -----
            case 8:
                // 3 neighbours. O+ (Oxonium Oxygen)
                if (degree(mol, atom) == 3)
                    return 49;
                // 2 neighbours.
                if (degree(mol, atom) == 2) {
                    // O=+ (Oxenium Oxygen)
                    if (mol.getOccupiedValence(atom) + mol.getImplicitHydrogens(atom) == 3)
                        return 51;

                    int nHtoO = 0;
                    for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                        int nbr = mol.getConnAtom(atom, i);

                        if (mol.getAtomicNo(nbr) == 1)
                            nHtoO++;
                    }

                    // OH2 (Oxygen in water)
                    if (nHtoO + mol.getImplicitHydrogens(atom) == 2)
                        return 70;

                    // Otherwise ipso must be one of:
                    // OC=O (Carboxylic acid or ester Oxygen)
                    // OC=C (Enolic or phenolic Oxygen)
                    // OC=N (Oxygen in -O-C=N- moiety)
                    // OC=S (Divalent Oxygen in thioacid or ester)
                    // ONO2 (Divalent nitrate "ether" Oxygen)
                    // ON=O (Divalent nitrate "ether" Oxygen)
                    // OSO3 (Divalent Oxygen in sulfate group)
                    // OSO2 (Divalent Oxygen in sulfite group)
                    // OSO  (One of two divalent Oxygens attached to Sulfur)
                    // OS=O (Divalent Oxygen in R(RO)S=O)
                    // -OS  (Other divalent Oxygen attached to Sulfur)
                    // OPO3 (Divalent Oxygen in phosphate group)
                    // OPO2 (Divalent Oxygen in phosphite group)
                    // OPO  (Divalent Oxygen, one of two Oxygens attached to P)
                    // -OP  (Other divalent Oxygen attached to Phosphorus)
                    return 6;
                }
                // 1 or 0 neighbours.
                if (mol.getConnAtoms(atom) <= 1) {
                    int nNtoS = 0;
                    int nOtoS = 0;
                    int nStoS = 0;

                    boolean isOxideOtoH = mol.getAllConnAtoms(atom)
                        - mol.getConnAtoms(atom)
                        + mol.getImplicitHydrogens(atom) > 0 ? true : false;
                    boolean isCarboxylateO            = false;
                    boolean isCarbonylO               = false;
                    boolean isOxideOtoC               = false;
                    boolean isNitrosoO                = false;
                    boolean isOxideOtoN               = false;
                    boolean isNOxideO                 = false;
                    boolean isNitroO                  = false;
                    boolean isThioSulfinateO          = false;
                    boolean isSulfateO                = false;
                    boolean isSulfoxideO              = false;
                    boolean isPhosphateOrPerchlorateO = false;

                    for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                        if (isOxideOtoH) break;
                        if (isCarboxylateO) break;
                        if (isCarbonylO) break;
                        if (isOxideOtoC) break;
                        if (isNitrosoO) break;
                        if (isOxideOtoN) break;
                        if (isNOxideO) break;
                        if (isNitroO) break;
                        if (isThioSulfinateO) break;
                        if (isSulfateO) break;
                        if (isSulfoxideO) break;
                        if (isPhosphateOrPerchlorateO) break;

                        int nbr = mol.getConnAtom(atom, i);

                        if (mol.getAtomicNo(nbr) == 6
                                || mol.getAtomicNo(nbr) == 7
                                || mol.getAtomicNo(nbr) == 16) {
                            for (int j=0; j<mol.getAllConnAtoms(nbr); j++) {
                                int nbr2 = mol.getConnAtom(nbr, j);

                                if (mol.getAtomicNo(nbr2) == 7
                                        && degree(mol, nbr2) == 2)
                                    nNtoS++;
                                if (mol.getAtomicNo(nbr2) == 8
                                        && degree(mol, nbr2) == 1)
                                    nOtoS++;
                                if (mol.getAtomicNo(nbr2) == 16
                                        && degree(mol, nbr2) == 1)
                                    nStoS++;
                            }
                        }

                        isOxideOtoH = mol.getAtomicNo(nbr) == 1;

                        // If ipso neighbour is a Carbon.
                        if (mol.getAtomicNo(nbr) == 6) {
                            isCarboxylateO = nOtoS == 2;

                            isCarbonylO = mol.getBondOrder(mol.getBond(atom, nbr)) == 2;

                            isOxideOtoC = mol.getBondOrder(mol.getBond(atom, nbr)) == 1
                                && nOtoS == 1;
                        }

                        // If ipso neighbour is a Nitrogen.
                        if (mol.getAtomicNo(nbr) == 7) {
                            isNitrosoO = mol.getBondOrder(mol.getBond(atom, nbr)) == 2;

                            if (mol.getBondOrder(mol.getBond(atom, nbr)) == 1 && nOtoS == 1) {
                                isOxideOtoN = degree(mol, nbr) == 2
                                    || mol.getOccupiedValence(nbr)
                                     + mol.getImplicitHydrogens(nbr) == 3;

                                isNOxideO = mol.getOccupiedValence(nbr)
                                     + mol.getImplicitHydrogens(nbr) == 4;
                            }

                            isNitroO = nOtoS >= 2;
                        }

                        // If ipso neighbour is a Sulfur.
                        if (mol.getAtomicNo(nbr) == 16) {
                            isThioSulfinateO = nStoS == 1;

                            isSulfateO = mol.getBondOrder(mol.getBond(atom, nbr)) == 1
                                || mol.getBondOrder(mol.getBond(atom, nbr)) == 2
                                && nOtoS + nNtoS > 1;

                            isSulfoxideO = mol.getBondOrder(mol.getBond(atom, nbr)) == 2
                                && nOtoS + nNtoS == 1;
                        }

                        isPhosphateOrPerchlorateO = mol.getAtomicNo(nbr) == 15
                            || mol.getAtomicNo(nbr) == 17;
                    }

                    // OM  (Oxide oxygen on sp3 carbon
                    // OM2 (Oxide oxygen on sp2 carbon
                    // OM  (Oxide oxygen on sp3 nitrogen (not in original MMFF.I Table III)
                    // OM2 (Oxide oxygen on sp2 nitrogen (not in original MMFF.I Table III)
                    if (isOxideOtoC || isOxideOtoN || isOxideOtoH)
                        return 35;

                    // O2CM (Oxygen in carboxylate group)
                    // ONX  (Oxygen in N-oxides)
                    // O2N  (Oxygen in nitro group)
                    // O2NO (Nitro-group oxygen in nitrate)
                    // O3N  (Nitrate anion oxygen)
                    // OSMS (Terminal oxygen in thiosulfinate anion)
                    // O-S  (Single terminal O on tetracoordinate sulfur)
                    // O2S  (One of 2 terminal O's on sulfur)
                    // O3S  (One of 3 terminal O's on sulfur)
                    // O4S  (Terminal O in sulfate anion)
                    // OP   (Oxygen in phosphine oxide)
                    // O2P  (One of 2 terminal O's on P)
                    // O3P  (One of 3 terminal O's on P)
                    // O4P  (One of 4 terminal O's on P)
                    // O4Cl (Oxygen in perchlorate anion)
                    if (isCarboxylateO || isNitroO || isNOxideO || isThioSulfinateO
                            || isSulfateO || isPhosphateOrPerchlorateO)
                        return 32;

                    // O=C  (Generic carbonyl oxygen)
                    // O=CN (Carbonyl oxygen in amides)
                    // O=CR (Carbonyl oxygen in aldehydes and ketones)
                    // O=CO (Carbonyl oxygen in acids and esters)
                    // O=N  (Nitroso oxygen)
                    // O=S  (Doubly bonded sulfoxide oxygen)
                    if (isCarbonylO || isNitrosoO || isSulfoxideO)
                        return 7;
                }
                break;

            // ----- Fluorine -----
            case 9:
                // F (Fluorine)
                if (mol.getConnAtoms(atom) == 1)
                    return 11;
                // F- (Fluoride anion)
                if (mol.getConnAtoms(atom) == 0)
                    return 89;
                break;

            // ----- Sodium -----
            case 11:
                // Na+ (Sodium cation)
                if (mol.getConnAtoms(atom) == 0)
                    return 93;
                break;

            // ----- Magnesium -----
            case 12:
                // Mg+2 (Dipositive Magnesium cation)
                if (mol.getConnAtoms(atom) == 0)
                    return 99;
                break;

            // ----- Silicon -----
            case 14:
                // Si (Silicon)
                return 19;

            // ----- Phosphorus -----
            case 15:
                // PO4 (Phosphate group phosphorus)
                // PO3 (Phosphorus with 3 attached oxygens)
                // PO2 (Phosphorus with 2 attached oxygens)
                // PO  (Phosphine oxide phosphorus)
                if (degree(mol, atom) == 4)
                    return 25;
                // P (Phosphorus in phosphines)
                if (degree(mol, atom) == 3)
                    return 26;
                // -P=C (Phosphorus doubly bonded to C)
                if (degree(mol, atom) == 2)
                    return 75;
                break;

            // ----- Sulfur -----
            case 16:
                // 3 or 4 neighbours.
                if (degree(mol, atom) == 3 || degree(mol, atom) == 4) {
                    int nONtoS = 0;
                    int nStoS = 0;
                    boolean isCdoubleS = false;

                    for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                        int nbr = mol.getConnAtom(atom, i);

                        // Check if ipso Sulfur is double bonded to Carbon.
                        if (mol.getAtomicNo(nbr) == 6
                                && mol.getBondOrder(mol.getBond(atom, nbr)) == 2)
                            isCdoubleS = true;

                        // If neighbour is terminal Oxygen/Sulfur or a
                        // secondary Nitrogen increment counter.
                        if ((mol.getConnAtoms(nbr) == 1 && mol.getAtomicNo(nbr) == 8)
                                || (degree(mol, nbr) == 2 && mol.getAtomicNo(nbr) == 7))
                            nONtoS++;

                        if (mol.getConnAtoms(nbr) == 1 && mol.getAtomicNo(nbr) == 16)
                            nStoS++;
                    }

                    // =SO2 (Sulfone Sulfur, doubly bonded to Carbon)
                    if (degree(mol, atom) == 3 && nONtoS == 2 && isCdoubleS
                            || degree(mol, atom) == 4)
                        return 18;

                    // SSOM (Tricoordinate Sulfur in anionic thiosulfinate
                    // group)
                    if ((nONtoS > 0 && nStoS > 0) || nONtoS == 2 && !isCdoubleS)
                        return 73;

                    // S=O (Sulfoxide Sulfur)
                    // S=N (Tricoordinate Sulfur doubly bonded to N)
                    return 17;
                }
                // 2 neighbours.
                if (degree(mol, atom) == 2) {
                    boolean isOdoubleS = false;

                    for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                        int nbr = mol.getConnAtom(atom, i);

                        // Check if ipso Sulfur is double bonded to Oxygen.
                        if (mol.getAtomicNo(nbr) == 8
                                && mol.getBondOrder(mol.getBond(atom, nbr)) == 2)
                            isOdoubleS = true;
                    }

                    // =S=O (Sulfinyl Sulfur)
                    if (isOdoubleS)
                        return 74;

                    // S (Thiol, Sulfide or Disulfide Sulfur)
                    return 15;
                }
                // 1 neighbour.
                if (degree(mol, atom) == 1) {
                    int nbrTermS = 0;
                    boolean isCdoubleS = false;

                    for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                        int nbr = mol.getConnAtom(atom, i);

                        // Count how many terminal Sulfurs are 2 jumps from
                        // atom.
                        for (int j=0; j<mol.getAllConnAtoms(nbr); j++) {
                            int nbr2 = mol.getConnAtom(nbr, j);

                            if (mol.getAtomicNo(nbr2) == 16
                                    && degree(mol, nbr2) == 1)
                                nbrTermS++;
                        }

                        // Check if ipso Sulfur is double bonded to Carbon.
                        if (mol.getAtomicNo(nbr) == 6
                                && mol.getBondOrder(mol.getBond(atom, nbr)) == 2)
                            isCdoubleS = true;
                    }

                    // If ipso Sulfur is double bonded to Carbon and the
                    // latter is not bonded to other terminal Sulfur atoms,
                    // then it is not a dithiocarboxylate, but a thioketone,
                    // etc.
                    // S=C (Sulfur doubly bonded to carbon)
                    if (isCdoubleS && nbrTermS != 2)
                        return 16;

                    // Ipso must be one of:
                    // S-P  (Terminal Sulfur bonded to P)
                    // SM   (Anionic terminal Sulfur)
                    // SSMO (Terminal Sulfur in thiosulfinate group)
                    return 72;
                }
                break;

            // ----- Chlorine -----
            case 17:
                // 4 Neighbours
                if (mol.getConnAtoms(atom) == 4) {
                    int nO = 0;
                    for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                        int nbr = mol.getConnAtom(atom, i);

                        if (mol.getAtomicNo(nbr) == 8)
                            nO++;
                    }

                    // CLO4 (Perchlorate anion Chlorine)
                    if (nO == 4)
                        return 77;
                }
                // Cl (Chloride)
                if (mol.getConnAtoms(atom) == 1)
                    return 12;
                // Cl- (Chloride anion)
                if (mol.getConnAtoms(atom) == 0)
                    return 90;
                break;

            // ----- Potassium -----
            case 19:
                // K+ (Potassium cation)
                if (mol.getConnAtoms(atom) == 0)
                    return 94;
                break;

            // ----- Calcium -----
            case 20:
                // Ca+2 (Dipositive Calcium cation)
                if (mol.getConnAtoms(atom) == 0)
                    return 96;
                break;

            // ----- Iron -----
            case 26:
                if (mol.getConnAtoms(atom) == 0) {
                    // Fe+2 (Dipositive Iron cation)
                    if (mol.getAtomCharge(atom) == 2)
                        return 87;

                    // Fe+3 (Tripositive Iron cation)
                    if (mol.getAtomCharge(atom) == 3)
                        return 88;
                }
                break;

            // ----- Copper -----
            case 29:
                if (mol.getConnAtoms(atom) == 0) {
                    // Cu+1 (Monopositive Copper cation)
                    if (mol.getAtomCharge(atom) == 1)
                        return 97;

                    // Cu+2 (Dipositive Copper cation)
                    if (mol.getAtomCharge(atom) == 2)
                        return 98;
                }
                break;

            // ----- Zinc -----
            case 30:
                // ZN+2 (Dipositive Zinc cation)
                if (mol.getConnAtoms(atom) == 0)
                    return 95;
                break;

            // ----- Bromine -----
            case 35:
                // Br (Bromine)
                if (mol.getConnAtoms(atom) == 1)
                    return 13;

                // Br- (Bromide anion)
                if (mol.getConnAtoms(atom) == 0)
                    return 91;
                break;

            // ----- Iodine -----
            case 53:
                // I (Iodine)
                if (mol.getConnAtoms(atom) == 1)
                    return 14;
                break;
        }

        // Fail with 0.
        return 0;
    }

    /**
     * Returns the MMFF atom type for hydrogen atoms.
     *  @param mol The molecule that the atom is in.
     *  @param atom The atom to assign a type for.
     *  @return The mmff atom type.
     */
    private static int getHydrogenType(MMFFMolecule mol, int atom) {
        for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
            int nbr = mol.getConnAtom(atom, i);

            switch (mol.getAtomicNo(nbr)) {
                // ----- Carbon -----
                case 6:
                    // HC  (Hydrogen attached to Carbon)
                    return 5;

                // ----- Nitrogen -----
                case 7:
                    switch (getHeavyType(mol, nbr)) {
                        case 8:
                            // HNR (Generic hydrogen on sp3 nitrogen, e.g. in amines)
                            // H3N (Hydrogen in ammonia)
                        case 39:
                            // HPYL (Hydrogen on nitrogen in pyrrole)
                        case 62:
                            // HNR (Generic hydrogen on sp3 nitrogen, e.g. in amines)
                        case 67:
                        case 68:
                            // HNOX (Hydrogen on N in a N-oxide)
                            return 23;

                        case 34:
                            // NR+ (Quaternary nitrogen)
                        case 54:
                            // N+=C (Iminium nitrogen)
                            // N+=N (Positively charged nitrogen doubly bonded to N)
                        case 55:
                            // HNN+ (Hydrogen on amidinium nitrogen)
                        case 56:
                            // HGD+ (Hydrogen on guanidinium nitrogen)
                        case 58:
                            // NPD+ (Aromatic nitrogen in pyridinium)
                        case 81:
                            // HIM+ (Hydrogen on imidazolium nitrogen)
                            return 36;

                        case 9:
                            // HN=N (Hydrogen on azo Nitrogen)
                            // HN=C (Hydrogen on imine Nitrogen)
                            return 27;

                        default:
                        // HNCC (Hydrogen on enamine nitrogen)
                        // HNCN (Hydrogen in H-N-C=N moiety)
                        // HNCO (Hydrogen on amide nitrogen)
                        // HNCS (Hydrogen on thioamide nitrogen)
                        // HNNC (Hydrogen in H-N-N=C moiety)
                        // HNNN (Hydrogen in H-N-N=N moiety)
                        // HNSO (Hydrogen on NSO, NSO2, or NSO3 nitrogen)
                        // HNC% (Hydrogen on N triply bonded to C)
                        // HSP2 (Generic hydrogen on sp2 nitrogen)
                        return 28;
                    }

                // ----- Oxygen -----
                case 8:
                    switch (getHeavyType(mol, nbr)) {
                        case 49:
                            // HO+ (Hydrogen on oxonium oxygen)
                            return 50;

                        case 51:
                            // HO=+ (Hydrogen on oxenium oxygen)
                            return 52;

                        case 70:
                            // HOH (Hydroxyl hydrogen in water)
                            return 31;

                        case 6:
                            boolean isHOCCorHOCN = false;
                            boolean isHOCO = false;
                            boolean isHOP = false;
                            boolean isHOS = false;

                            for (int j=0; j<mol.getAllConnAtoms(nbr); j++) {
                                int nbr2 = mol.getConnAtom(nbr, j);

                                if (mol.getAtomicNo(nbr2) == 6) {
                                    for (int k=0; k<mol.getAllConnAtoms(nbr2); k++) {
                                        int nbr3 = mol.getConnAtom(nbr2, k);

                                        if (nbr3 == nbr)
                                            continue;

                                        int bond = mol.getBond(nbr2, nbr3);

                                        if ((mol.getAtomicNo(nbr3) == 6
                                                || mol.getAtomicNo(nbr3) == 7)
                                                && (mol.getBondOrder(bond) == 2
                                                || mol.isAromaticBond(bond)))
                                            isHOCCorHOCN = true;

                                        if (mol.getAtomicNo(nbr3) == 8
                                                && mol.getBondOrder(bond) == 2)
                                            isHOCO = true;
                                    }
                                }

                                if (mol.getAtomicNo(nbr2) == 15)
                                    isHOP = true;

                                if (mol.getAtomicNo(nbr2) == 16)
                                    isHOS = true;
                            }

                            // HOCO (Hydroxyl hydrogen in carboxylic acids)
                            if (isHOCO || isHOP)
                                return 24;

                            // HOCC (Enolic or phenolic hydroxyl hydrogen)
                            // HOCN (Hydroxyl hydrogen in HO-C=N moiety)
                            if (isHOCCorHOCN)
                                return 29;

                            // HOS
                            // Hydrogen on oxygen attached to sulfur
                            if (isHOS)
                                return 33;

                        default:
                            // HO (Generic hydroxyl hydrogen)
                            // HOR (Hydroxyl hydrogen in alcohols)
                            return 21;
                    }

                // ----- Silicon -----
                case 14:
                    // HSI (Hydrogen attached to silicon)
                    return 5;

                // ----- Phosphorus -----
                case 15:
                    // HP (Hydrogen attached to phosphorus)
                    return 71;

                // ----- Sulfur -----
                case 16:
                    // HS (Hydrogen attached to sulfur)
                    // HS=N (Hydrogen attached to >S= sulfur doubly bonded to N)
                    return 71;
            }
        }

        // Fail with 0.
        return 0;
    }

    /**
     * Returns true if an atom is an N-oxide.
     *  @param mol The molecule that the atom is in.
     *  @param atom The atom to check.
     *  @return True if the atom is an N-oxide, false otherwise.
     */
    private static boolean isAtomNOxide(MMFFMolecule mol, int atom) {
        if (mol.getAtomicNo(atom) == 7 && degree(mol, atom) >= 3) {
            for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                int nbr = mol.getConnAtom(atom, i);
                if (mol.getAtomicNo(nbr) == 8 && degree(mol, nbr) == 1)
                    return true;
            }
        }
        return false;
    }

    /**
     * Returns the number of rings that an atom is a member of in a
     * molecule.
     *  @param mol The molecule that the atom is in.
     *  @param atom The atom to check.
     *  @return The nuber of rings that the atom is a member of.
     */
    public static int inRings(MMFFMolecule mol, int atom) {
        RingCollection rings = mol.getRingSet();
        int rc = 0;

        for (int r=0; r<rings.getSize(); r++)
            if (rings.getAtomIndex(r,atom) >= 0)
                rc++;

        return rc;
    }

    /**
     * Returns true if the atom is in an aromatic ring of given size.
     *  @param mol The molecule that the atom is in.
     *  @param atom The atom to check.
     *  @param size The size of aromatic rings to check.
     *  @return True if the atom is in an aromatic ring of given size,
     *      false otherwise.
     */
    public static boolean isInAromaticRingOfSize(MMFFMolecule mol, int atom,
                                                 int size) {
        if (mol.isAromaticAtom(atom)) {
            RingCollection rings = mol.getRingSet();

            for (int r=0; r<rings.getSize(); r++) {
                if (rings.getRingSize(r) != size
                        || !rings.isAtomMember(r,atom))
                    continue;

                //if (ringIsMMFFAromatic(mol, r) == RingBoolean.TRUE)
                if (mol.ringIsMMFFAromatic(r))
                    return true;
            }
        }
        return false;
    }

    /**
     * Determine if a ring is aromatic according to MMFF criteria. Only
     * designed to work with rings of size 5 and 6.
     *  @param mol The molecule that the ring is in.
     *  @param r The ring.
     *  @return True if the ring is aromatic, false otherwise.
     */
    public static RingBoolean ringIsMMFFAromatic(MMFFMolecule mol, int r) {
        RingCollection rings = mol.getRingSet();

        if (!rings.isAromatic(r))
            return RingBoolean.FALSE;

        if (rings.getRingSize(r) == 6) {
            for (int a : rings.getRingAtoms(r)) {
                if (mol.getOccupiedValence(a) + mol.getImplicitHydrogens(a)
                        != (degree(mol, a)+1))
                    return RingBoolean.FALSE;
            }
            for (int b : rings.getRingBonds(r)) {
                 int [] c = { -1, -1 };
                if (mol.getBondOrder(b) == 1) {
                    for (int j = 0; j <= 1; j++) {
                        int atom = mol.getBondAtom(j, b);
                        for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
                            int nbr = mol.getConnAtom(atom, i);
                            if (!rings.isAtomMember(r, nbr)
                                    && mol.getBondOrder(mol.getBond(atom, nbr)) == 2) {
                                c[j] = nbr;
                                break;
                            }
                        }
                    }
                    if (c[0] > -1 && c[1] > -1) {
                        for (int ri = 0; ri < rings.getSize(); ri++) {
                            if (rings.isAtomMember(ri, c[0])
                                    && rings.isAtomMember(ri, c[1])
                                    && !mol.isSetRingMMFFAromaticity(ri))
                                return RingBoolean.NOT_SET;

                            if (rings.isAtomMember(ri, c[0])
                                    && rings.isAtomMember(ri, c[1])
                                    && !mol.ringIsMMFFAromatic(ri))
                                return RingBoolean.FALSE;
                        }
                    }
                }
            }
        }

        if (rings.getRingSize(r) == 5) {
            int passes = 1;

            for (int a : rings.getRingAtoms(r)) {
                if (mol.getOccupiedValence(a) + mol.getImplicitHydrogens(a)
                        == degree(mol, a) && passes > 0) {
                    passes--;
                    continue;
                }
                if (mol.getOccupiedValence(a) + mol.getImplicitHydrogens(a)
                        != (degree(mol, a)+1)) {
                    return RingBoolean.FALSE;
                }
            }
        }

        return RingBoolean.TRUE;
    }

    /**
     * Returns the total number of atoms connected to an atom, including
     * implicit and explicit hydrogens.
     *  @param mol The molecule that the atom is in.
     *  @param atom The atom to count the neighbours of.
     *  @return The number of neighbours connected to this atom.
     */
    public static int degree(MMFFMolecule mol, int atom) {
        return mol.getAllConnAtoms(atom) + mol.getImplicitHydrogens(atom);
    }

    /**
     * Returns true if the atom is in any ring of a given size.
     *  @param mol The molecule that the atom is in.
     *  @param atom The atom to check.
     *  @param size The size of rings to check.
     *  @return True if the atom is in a ring of given size, false
     *      otherwise.
     */
    public static boolean inRingOfSize(MMFFMolecule mol, int atom, int size) {
        if (mol.isRingAtom(atom)) {
            RingCollection rings = mol.getRingSet();
            for (int r=0; r<rings.getSize(); r++)
                if (rings.getRingSize(r) == size && rings.isAtomMember(r,atom))
                    return true;
        }
        return false;
    }

    /**
     * Returns true if two atoms are in the same ring of a given size.
     *  @param mol The molecule that both atoms are in.
     *  @param a1 The first atom.
     *  @param a2 The second atom.
     *  @param size The size of the ring both atoms must be a member of.
     *  @return True if both atoms are in a ring of given size, false
     *      otherwise.
     */
    public static boolean inSameRing(MMFFMolecule mol, int a1, int a2, int size) {
        if (!mol.isRingAtom(a1) || !mol.isRingAtom(a2))
            return false;

        RingCollection rings = mol.getRingSet();

        for (int r=0; r<rings.getSize(); r++)
            if (rings.getRingSize(r) == size
                    && rings.isAtomMember(r, a1)
                    && rings.isAtomMember(r, a2))
                return true;

        return false;
    }
}
