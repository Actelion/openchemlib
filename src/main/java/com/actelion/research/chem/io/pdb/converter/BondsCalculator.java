package com.actelion.research.chem.io.pdb.converter;
/*
 * Copyright 2016 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of cif2molecule.
 *
 * cif2molecule is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * cif2molecule is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Antanas Vaitkus
 */
import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.BondLengthSet;
import com.actelion.research.util.IntQueue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import static com.actelion.research.chem.conf.VDWRadii.COVALENT_RADIUS;

/**
 * BondsCalculator is used to recreate the bonds and / or calculate the bonds orders
 * based on the 3D coordinates of the atoms
 */
public class BondsCalculator {
    private static final boolean DEBUG = false;


    private static final double BUMP_FACTOR = 0.5;
    private static final double SPECIAL_POSITION_CUTOFF = 0.01;



    /**
     * Calculates the bonds of a molecule by checking the distance between
     * all atoms. The bond order is not set with this function.
     *
     * Complexity O(nAtoms)
     * Memory O(nAtoms)
     *
     * @param mol
     *                  The molecule for which the bonds should be calculated.
     * @param lenient
     *                  Suppresses exceptions that would normally be raised
     *                  due to chemical discrepancies found in the molecule.
     * @throws Exception
     *                  Exception containing information about the exceeded
     *                  maximum valence of an atom.
     */
    public static void createBonds(Molecule3D mol, boolean lenient) throws Exception {

        if ( mol.getAllAtoms()==0 ) return;

        // Processing organic atoms first:
        // 1. Create a grid
        MoleculeGrid grid = new MoleculeGrid(mol);
        int[] neighborCount = new int[mol.getAllAtoms()];

        // 2. For each atom, check the neighbours and create
        //    a connection if the distance is close to the sum
        //    of covalent radii.
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            if (!mol.isOrganicAtom(i)) continue;

            // Get the neighbours
            Set<Integer> set = grid.getNeighbours(mol, i, 3.2);

            for (int j : set) {
                if (i > j) continue;
                if (!mol.isOrganicAtom(j)) continue;
                // Two hydrogen (H) atoms are very unlikely to form bonds in crystal structures
                if (mol.getAtomicNo(i) == 1 && mol.getAtomicNo(j) == 1) continue;

                double dist = getBondLength(mol, i, j);
                double idealDist = COVALENT_RADIUS[mol.getAtomicNo(i)] + COVALENT_RADIUS[mol.getAtomicNo(j)];

                if (dist > idealDist + .45) continue;

                if (neighborCount[i] >= maxNeighborCount(mol, i) ||
                        neighborCount[j] >= maxNeighborCount(mol, j)) {

                    int hypervalentAtom = (neighborCount[i] >= maxNeighborCount(mol, i)) ? i : j;

                    String error = "maximum valence of " + maxNeighborCount(mol, hypervalentAtom) +
                            " exceeded for atom '" + mol.getAtomName(hypervalentAtom) + "'";
                    if (!lenient) {
                        throw new BondCalculatorException(error);
                    } else {
                        System.err.println();
                    }
                    continue;
                }
                try {
                    mol.addBond(i, j, 1);
                    neighborCount[i]++;
                    neighborCount[j]++;
                } catch (Exception e) {
                    if (!lenient) throw e;
                }
            }
        }
        mol.ensureHelperArrays(Molecule.cHelperRings);

        // Deleting incorrectly assigned bonds
        for (int i = 0; i < mol.getRingSet().getSize(); i++) {
            int[] ringAtoms = mol.getRingSet().getRingAtoms(i);
            if (ringAtoms.length == 3) {
                if (markImprobableDiphosphorusBondForDeletionInThreeAtomRing(mol, i)) continue;
                markImprobableBondWithHydrogenForDeletionInThreeAtomRing(mol, i);
            } else if (ringAtoms.length == 4) {
                // Small rings sometimes get assigned incorrect interring bonds.
                // For example, a 4-membered ring with atom connectivity a1-a2,
                // a2-a3, a3-a4, a4-a1 sometimes get assigned an additional
                // a1-a3 (or a2-a4) bond.

                // The ring must be more or less flat to avoid deleting bonds
                // in molecules that have a closely packed tetrahedral geometry
                // (e. g. tetrahedrane).
                if (StrictMath.abs(getDihedral(mol, ringAtoms[0], ringAtoms[1],
                        ringAtoms[2], ringAtoms[3])) > 20 * StrictMath.PI / 180) continue;

                for (int ringAtom : ringAtoms) {
                    int inRingNeighbours = 0;
                    for (int k = 0; k < mol.getConnAtoms(ringAtom); k++) {
                        if (isRingAtom(mol, i, mol.getConnAtom(ringAtom, k))) {
                            inRingNeighbours++;
                        }
                    }
                    if (inRingNeighbours > 2) {
                        int longestBondIndex = -1;
                        double longestBondLength = -Double.MAX_VALUE;
                        for (int k = 0; k < mol.getConnAtoms(ringAtom); k++) {
                            if (!isRingAtom(mol, i, mol.getConnAtom(ringAtom, k))) continue;
                            if (longestBondLength < getBondLength(mol, mol.getConnBond(ringAtom, k))) {
                                longestBondIndex = k;
                                longestBondLength = getBondLength(mol, mol.getConnBond(ringAtom, k));
                            }
                        }
                        //  if (isPlanar(mol, ringAtom, mol.getConnAtom(ringAtom, longestBondIndex))) {
                        mol.markBondForDeletion(mol.getConnBond(ringAtom, longestBondIndex));
                        //   }
                    }
                }
            }
        }
        mol.deleteMarkedAtomsAndBonds();

        // Processing non-organic atoms
        grid = new MoleculeGrid(mol);
        for( int i = 0; i < mol.getAllAtoms(); i++ ) {
            if (mol.isOrganicAtom(i)) continue;
            // TODO: since these are not covalent bonds, maybe none of the metals should be excluded?
            // Not Lithium or Magnesium
            // TODO: also Be in organometallic compounds?
            // TODO: allow Alkali metals to coordinate water?
     //       if ( (isAlkaliMetalAtom(mol, i) && mol.getAtomicNo(i) != 3 ) ||
     //             (isAlkalineEarthMetalAtom(mol, i) &&  mol.getAtomicNo(i) != 12 ) ) continue;
            Set<Integer> set = grid.getNeighbours(mol, i, 3.2);

            for (int j : set) {
                double dist = getBondLength( mol, i, j );
                double idealDist = COVALENT_RADIUS[mol.getAtomicNo(i)] + COVALENT_RADIUS[mol.getAtomicNo(j)];

                if (dist > idealDist + .45) continue;
                mol.addBond(i, j, Molecule3D.cBondTypeMetalLigand);
            }
        }

        int ringNumber = mol.getRingSet().getSize();
        Coordinates[] ringCenters = new Coordinates[ringNumber];
        for (int i = 0; i < ringNumber; i++) {
            ringCenters[i] = getAtomRingCenter(mol, i);
        }
        ArrayList<Integer>[] atomsToRings = getAtomToRings(mol);
        for ( int i = 0; i < mol.getAllAtoms(); i++ ) {
            if ( !mol.isMetalAtom(i) ) continue;
            boolean[] ringProcessed  = new boolean[ringNumber];
            boolean[] ringCoordinated = new boolean[ringNumber];
            for ( int j = mol.getAllConnAtoms(i); j < mol.getAllConnAtomsPlusMetalBonds(i); j++ ) {

                int ligand = mol.getConnAtom(i, j);
                if ( mol.isMetalAtom(ligand) ) continue;
                if ( atomsToRings[ligand] == null ) continue;

                int closestRing = -1;
                double closestDist = 0;
                for ( int ringIndex: atomsToRings[ligand] ) {
                    if ( closestRing < 0 ||
                            ringCenters[ringIndex].distance(mol.getCoordinates(i)) < closestDist ) {
                        closestDist = ringCenters[ringIndex].distance(mol.getCoordinates(i));
                        closestRing = ringIndex;
                    }
                }
                if ( ringProcessed[closestRing] ) continue;

                int bondCount = 0;
                for ( int k = 0; k < mol.getAllConnAtomsPlusMetalBonds(i); k++ ) {
                    if ( isRingAtom(mol, closestRing, mol.getConnAtom(i, k) ) ) bondCount++;
                }
                // TODO: check is ring is (potentially) aromatic?
                ringProcessed[closestRing] = true;
                if ( bondCount > mol.getRingSet().getRingSize(closestRing) / 2 ) {
                    ringCoordinated[closestRing] = true;
                }
            }
            for ( int j = 0; j < ringNumber; j++ ) {
                if ( ringCoordinated[j] ) {
                    int[] ringAtoms = mol.getRingSet().getRingAtoms(j);
                    for (int ringAtom : ringAtoms) {
                        mol.addBond(i, ringAtom, Molecule.cBondTypeMetalLigand);
                    }
                }
            }
        }
        mol.deleteMarkedAtomsAndBonds();

        mol.ensureHelperArrays(Molecule.cHelperNeighbours);

        // TODO: ask about the applicability of this part
        // remove erroneously assigned metal-ligand bonds
        ArrayList<Integer> bondsToDelete = new ArrayList<>();
        for (int i = 0; i < mol.getAllBonds(); i++) {
            if ( mol.getBondType(i) != Molecule3D.cBondTypeMetalLigand ) continue;

            int metal;
            int ligand;
            if ( mol.isMetalAtom( mol.getBondAtom(0, i) ) ) {
                metal  = mol.getBondAtom(0, i);
                ligand = mol.getBondAtom(1, i);
            } else {
                metal  = mol.getBondAtom(1, i);
                ligand = mol.getBondAtom(0, i);
            }

            for (int j = 0; j < mol.getAllConnAtomsPlusMetalBonds(ligand); j++ ) {
                int ligandNeighbour = mol.getConnAtom(ligand, j);
                if ( ligandNeighbour == metal ) continue;
                double sigma = 0.2;
                if ( getBondLength( mol, metal, ligand ) - getBondLength(mol, metal, ligandNeighbour) > sigma) {
                    bondsToDelete.add(i);
                }
            }
        }
        bondsToDelete.forEach(mol::markBondForDeletion);
        mol.deleteMarkedAtomsAndBonds();
    }

    /**
     * Evaluates if the given three-membered ring contains an unlikely bond
     * between two phosphorus atoms. The bond is considered unlikely if it
     * occurs in a three-membered ring and is above the specified bond length
     * threshold.
     *
     * Such incorrect bonds are likely to be placed in molecules where two
     * phosphorus (P) atoms are bridged by a third atom forming an acute angle.
     * For example, COD entry 4071806 contains the PCP pattern that is incorrectly
     * fused into a P1CP1 ring due to the distance between the two P atoms (2.426 A)
     * falling into the range of previously observed long P-P bonds. It should be noted,
     * however, that not all diphosporus bonds appearing in three-membered rings are
     * incorrect as illustrated by COD entry 4062100. As such, additional bond length
     * threshold is used to determine is the bond should be marked for deletion.
     *
     * @param mol
     *          Molecule that contains the ring.
     * @param ringIndex
     *          Index of the ring in the molecule.
     * @return
     *          Boolean value denoting if an improbable diphosphorus
     *          bond was identified and marked for deletion.
     */
    private static boolean markImprobableDiphosphorusBondForDeletionInThreeAtomRing(Molecule3D mol, int ringIndex) {
        if (mol.getRingSet().getRingSize(ringIndex) != 3) return false;

        int[] ringBonds = mol.getRingSet().getRingBonds(ringIndex);
        int longestDiphosphoroBondIndex = -1;
        double longestDiphosphoroBondLength = 0;
        for (int ringBond : ringBonds) {
            int a1 = mol.getBondAtom(0, ringBond);
            if (mol.getAtomicNo(a1) != 15) continue;
            int a2 = mol.getBondAtom(1, ringBond);
            if (mol.getAtomicNo(a2) != 15) continue;

            if (getBondLength(mol, ringBond) > longestDiphosphoroBondLength) {
                longestDiphosphoroBondIndex = ringBond;
                longestDiphosphoroBondLength = getBondLength(mol, ringBond);
            }
        }

        if (longestDiphosphoroBondIndex == -1) return false;
        // TODO: eventually replace the hardcoded value with one based on statistics
        // Explanation: the current cutoff values of 2.42 A was selected
        // after analysing about a hundred structures that were assigned
        // extremely long diphosphorus bonds. The longest correct
        // disphophorus bond in a three-membered ring with the bond
        // length of 2.409 was observed in COD entry 7107333.
        if (longestDiphosphoroBondLength < 2.42) return false;

        mol.markBondForDeletion(longestDiphosphoroBondIndex);

        return true;
    }

    /**
     * Evaluates if the given three-membered ring contains an unlikely bond
     * that involves a hydrogen atom.
     *
     * Such incorrect bonds are placed in molecules where a hydrogen atom
     * is positioned withing a bonding distance of two atoms. For example,
     * COD entry 7027071 contains a hydrogen atom ('H4C') that is within
     * a 0.8 A distance of the nitrogen atom ('N4') and a 1.46 A distance
     * of a carbon atom ('C8'). As a result, an incorrect bond between
     * the 'H4C' and 'C8' atoms gets placed forming a three-membered ring
     * (C1HN1).
     *
     * @param mol
     *          Molecule that contains the ring.
     * @param ringIndex
     *          Index of the ring in the molecule.
     * @return
     *          Boolean value denoting if an improbable bond that involves
     *          a hydrogen atom was identified and marked for deletion.
     */
    private static boolean markImprobableBondWithHydrogenForDeletionInThreeAtomRing(Molecule3D mol, int ringIndex){
        if (mol.getRingSet().getRingSize(ringIndex) != 3) return false;

        // TODO: allow two boron atoms to share a hydrogen
        boolean bondMarkedForDeletion = false;
        int[] ringAtoms = mol.getRingSet().getRingAtoms(ringIndex);
        for (int j = 0; j < ringAtoms.length; j++) {
            if (mol.getAtomicNo(ringAtoms[j]) != 1) continue;
            int n1 = ringAtoms[((j - 1) + ringAtoms.length) % ringAtoms.length];
            int n2 = ringAtoms[((j + 1) + ringAtoms.length) % ringAtoms.length];
            int boundAtom = n1;
            int unboundAtom = n2;
            if (getBondLength(mol, ringAtoms[j], n1) > getBondLength(mol, ringAtoms[j], n2)) {
                boundAtom = n2;
                unboundAtom = n1;
            }

            mol.markBondForDeletion(mol.getBond(ringAtoms[j], unboundAtom));
            System.err.printf(
                    "WARNING",
                    "a bond between atoms '" + mol.getAtomName(ringAtoms[j]) + "' and '" +
                            mol.getAtomName(unboundAtom) + "' was initially assigned based on the interatomic " +
                            "distance (" + "%.3f" + " nm), but was " +
                            "subsequently removed based on the molecular geometry -- atoms '" +
                            mol.getAtomName(ringAtoms[j]) + "', '" + mol.getAtomName(boundAtom) + "' and '" +
                            mol.getAtomName(unboundAtom) + "' are very unlikely to form a three-membered atom ring" + "\n",
                    getBondLength(mol, ringAtoms[j], unboundAtom)/10
            );
            bondMarkedForDeletion = true;
        }

        return bondMarkedForDeletion;
    }

    /**
     * Calculates the geometrical center of the atom ring.
     * @param mol
     *          The molecule that contains the ring.
     * @param ringIndex
     *          The ring index in the molecule.
     * @return
     *          Coordinates of the geometrical center of the atom ring.
     */
    private static Coordinates getAtomRingCenter(Molecule3D mol, int ringIndex) {
        Coordinates c = new Coordinates();
        for ( int i: mol.getRingSet().getRingAtoms(ringIndex) ) {
            c.x += mol.getAtomX(i);
            c.y += mol.getAtomY(i);
            c.z += mol.getAtomZ(i);
        }
        int ringAtomCount = mol.getRingSet().getRingAtoms(ringIndex).length;
        c.x /= ringAtomCount;
        c.y /= ringAtomCount;
        c.z /= ringAtomCount;

        return c;
    }

    /**
     * Calculates for the organic subset of atoms the maximum number of neighbors
     * @param mol
     *          The molecule that contains the atom.
     * @param atom
     *          The index of the atom in the molecule.
     * @return
     *          The maximum neighbour count for the atom.
     */
    private static int maxNeighborCount(Molecule3D mol, int atom) {
        int atomicNo = mol.getAtomicNo(atom);
        return atomicNo == 5 ? 6    // B
                : atomicNo <= 7 ? 4    // N, C
                : atomicNo == 8 ? 3    // O
                : atomicNo == 9 ? 1    // F
                : atomicNo == 14 ? 6   // Si
                : atomicNo == 15 ? 6   // P
                : atomicNo == 16 ? 6   // S
                : atomicNo == 17 ? 4   // Cl
                : atomicNo == 33 ? 6   // As
                : atomicNo == 34 ? 6   // Se
                : atomicNo == 35 ? 6   // Br
                : atomicNo == 52 ? 6
                : atomicNo == 53 ? 6 : 8;
    }

    /**
     * Calculate the bond orders of the molecule (without knowing the hydrogens).
     * The calculation is based on the bond distance between each atoms.
     *
     * http://www.ccp14.ac.uk/ccp/web-mirrors/i_d_brown/valence.txt
     * s = exp((Ro - R)/B)
     *
     * The implementation of this method is also most similar to the algorithm
     * described in:
     * http://www.daylight.com/meetings/mug01/Sayle/m4xbondage.html
     *
     * @param mol
     *              The molecule for which the bond orders should be calculated.
     */
     public static void calculateBondOrders(Molecule3D mol) {

        mol.ensureHelperArrays(Molecule.cHelperRings);
        /*
         * Pass 1: assign higher order bonds in unambiguous locations
         */
        assignObviousHigherOrderBonds(mol);

        // Place disilene bonds
        assignDisileneBonds(mol);

        //////////////////////////////////
        // Functional Group Recognition
        // Hybridization State Determination
        resolveKnownIons(mol);

        int[] spOrder = calculateHybridizationStates(mol);
        recogniseFunctionalGroups(mol, spOrder);
        // FIXME: this method should eventually be moved to the start of the method
        resolveSharedHydrogenAtoms(mol);
        // hybridization states need to be recalculated since
        // the resolveSharedHydrogenAtoms() method might affect
        // the order of atoms in the array
        spOrder = calculateHybridizationStates(mol);

        boolean[] visited = new boolean[mol.getAllBonds()];
        if (DEBUG) {
            System.err.println("> DEBUG: entering preliminary pass phase");
        }
        //Preliminary pass: process obvious bonds outside rings
        for (int i = 0; i < mol.getBonds(); i++) {
            if (mol.getBondType(i) == Molecule.cBondTypeMetalLigand) continue;
            if (mol.isRingBond(i)) continue;
            int a1 = mol.getBondAtom(0, i);
            int a2 = mol.getBondAtom(1, i);
            if (mol.getLowestFreeValence(a1) < 2) continue;
            if (mol.getLowestFreeValence(a2) < 2) continue;

            // if (!isPlanar(mol, a1, a2)) continue;
            if (spOrder[a1] != 1 || spOrder[a2] != 1) continue;
            int index = BondLengthSet.getBondIndex(3, false, false, mol.getAtomicNo(a1), mol.getAtomicNo(a2), 2, 2);
            if (index == -1) continue;

            double bondLength = BondLengthSet.getBondLength(index);
            double bondLengthStd = BondLengthSet.getBondStdDev(index);
            // Used to be 3 order > TRIPLE_BOND_CUTOFF
            if (bondLength + bondLengthStd > getBondLength(mol, i)) {
                mol.setBondOrder(i, 3);
                visited[i] = true;
                if (DEBUG) {
                    System.err.println(">> DEBUG: placing a triple bond between atoms '" +
                            mol.getAtomName(mol.getBondAtom(0, i)) + "' and '" +
                            mol.getAtomName(mol.getBondAtom(1, i)) + "' based on " +
                            "the geometry of the molecule and bond length statistics.");
                }
            }
        }
        if (DEBUG) {
            System.err.println("> DEBUG: exiting preliminary pass phase");
        }

        identifySquarates(mol, spOrder);

        if (DEBUG) {
            System.err.println("> DEBUG: entering exocyclic double bond detection phase.");
        }
        // go through all exocyclic atoms and set the most likely double bonds

        boolean[] isNonAromatic = assignExocyclicDoubleBonds(mol, spOrder);
        // check exocyclic atoms that are part of non-aromatic rings
        for (int i = 0; i < mol.getRingSet().getSize(); i++) {
            if (isNonAromatic[i]) continue;
            int[] ringAtoms = mol.getRingSet().getRingAtoms(i);
            for (int ringAtom : ringAtoms) {
                if (mol.getLowestFreeValence(ringAtom) < 1) continue;
                for (int j = 0; j < mol.getConnAtoms(ringAtom); j++) {
                    // check if atom is exocyclic and capable of a double bond
                    int exoAtom = mol.getConnAtom(ringAtom, j);
                    if (mol.getAtomCharge(exoAtom) != 0) continue;
                    if (mol.getLowestFreeValence(exoAtom) < 1) continue;
                    if (spOrder[exoAtom] != 2) continue;

                    //  if ( mol.getAtomicNo(exoAtom) == 8 ) continue;
                    // Check if the neighbour atom doesn't belong to the same ring
                    if (isRingAtom(mol, i, exoAtom)) continue;
                    boolean allowsDoubleBonds = true;
                    if (mol.isRingAtom(exoAtom)) {
                        for (int k = 0; k < mol.getRingSet().getSize(); k++) {
                            if (isRingAtom(mol, k, exoAtom) && !isNonAromatic[k]) {
                                allowsDoubleBonds = false;
                                break;
                            }
                        }
                    }

                    if (!allowsDoubleBonds) continue;
                    if (!isPlanar(mol, ringAtom, exoAtom)) continue;

                    int minPiCount = (mol.getConnAtoms(exoAtom) > 1) ? 1 : 0;

                    double averageSingleBondLength = getIdealisedBondLength(1, mol.getAtomicNo(ringAtom), mol.getAtomicNo(exoAtom), 1, minPiCount);
                    double averageDoubleBondLength = getIdealisedBondLength(2, mol.getAtomicNo(ringAtom), mol.getAtomicNo(exoAtom), 1, 1);
                    double averageDoubleBondStd = BondLengthSet.getBondStdDev(BondLengthSet.getBondIndex(2, false, false, mol.getAtomicNo(ringAtom), mol.getAtomicNo(exoAtom), 1, 1));
                    double calculatedLength = getBondLength(mol, ringAtom, exoAtom);
                    double doubleBondDiff = StrictMath.abs(averageDoubleBondLength - calculatedLength);

                    if (averageSingleBondLength < averageDoubleBondLength) {
                        if (DEBUG) {
                            System.err.println(">>> DEBUG: single '" + mol.getAtomicNo(ringAtom) + "-" +
                                    mol.getAtomicNo(exoAtom) + "' bond seems to be shorter than a the " +
                                    "double bond of the same type -- skipping the estimation of this bond.");
                        }
                        continue;
                    }

                    if (doubleBondDiff < 0 || doubleBondDiff <= averageDoubleBondStd * 1.5) {
                        mol.setBondOrder(mol.getBond(ringAtom, exoAtom), 2);
                        if (DEBUG) {
                            System.err.println(">>> DEBUG: placing a double bond between atoms '" +
                                    mol.getAtomName(ringAtom) + "' and '" + mol.getAtomName(exoAtom) + "'.");
                        }
                    }
                }
            }
        }
        if (DEBUG) {
            System.err.println("> DEBUG: exiting exocyclic ring atom double bond detection phase.");
        }

        /////////////////////////////////////////////////////////
        // Aromatic Ring Perception
        // This procedure calculates a normal to the ring and check that all
        // atoms in the ring and their neighbours are within the plane defined by the normal
        if ( DEBUG ) {
            System.err.println("> DEBUG: entering aromatic ring perception phase");
        }
        RingCollection ringSet = mol.getRingSet();
        boolean[] aromaticRing = findAromaticRings(mol, spOrder);
        for (int i = 0; i < ringSet.getSize(); i++) {
            if ( aromaticRing[i] ) {
                boolean isResolved = false;
                if (resolveTetraatomicAromaticRingIons(mol,i)) {
                    isResolved = true;
                // Cyclothiaselenazenium
                } else if ( ringSet.getRingSize(i) == 7 ) {
                    int NSBondCount = 0;
                    int SSBondCount = 0;
                    for (int bond : ringSet.getRingBonds(i)) {
                        if ( mol.getAtomicNo(mol.getBondAtom(0, bond) ) == 16 &&
                             mol.getAtomicNo(mol.getBondAtom(1, bond) ) == 16 ) {
                            SSBondCount++;
                        } else if (
                            ( mol.getAtomicNo(mol.getBondAtom(0, bond)) == 7 &&
                              mol.getAtomicNo(mol.getBondAtom(1, bond)) == 16 ) ||
                            ( mol.getAtomicNo(mol.getBondAtom(0, bond)) == 16 &&
                              mol.getAtomicNo(mol.getBondAtom(1, bond)) == 7 ) ) {
                            NSBondCount++;
                        } else {
                            break;
                        }
                    }
                    if ( NSBondCount == 6 && SSBondCount == 1 ) {
                        int SSBond = -1;
                        int ringSize = ringSet.getRingSize(i);
                        for (int j = 0; j < ringSet.getRingBonds(i).length; j++) {
                            int bond = ringSet.getRingBonds(i)[j];
                            if ( mol.getAtomicNo(mol.getBondAtom(0, bond) ) == 16 &&
                                 mol.getAtomicNo(mol.getBondAtom(1, bond) ) == 16 ) {
                                SSBond = j;
                                break;
                            }
                        }
                        if (SSBond == -1) {
                            throw new RuntimeException("the S-S bond could not be found even though " +
                                "it was previously detected.");
                        } else {
                            int doubleBond1 = ringSet.getRingBonds(i)[(SSBond - 2 + ringSize) % ringSize];
                            int doubleBond2 = ringSet.getRingBonds(i)[(SSBond + 2 + ringSize) % ringSize];
                            int chargedNBond = ringSet.getRingBonds(i)[(SSBond + 3 + ringSize) % ringSize];

                            mol.setBondOrder(doubleBond1, 2);
                            int sulphurAtom = mol.getAtomicNo(mol.getBondAtom(0,doubleBond1)) == 16 ? 0 : 1;
                            mol.setAtomCharge(mol.getBondAtom(sulphurAtom, doubleBond1), 1);
                            mol.setBondOrder(doubleBond2, 2);
                            sulphurAtom = mol.getAtomicNo(mol.getBondAtom(0,doubleBond2)) == 16 ? 0 : 1;
                            mol.setAtomCharge(mol.getBondAtom(sulphurAtom, doubleBond2), 1);
                            int nitrogenAtom = mol.getAtomicNo(mol.getBondAtom(0,doubleBond1)) == 7 ? 0 : 1;
                            mol.setAtomCharge(mol.getBondAtom(nitrogenAtom, chargedNBond), 1);
                        }
                        isResolved = true;
                    }
                }

                if (!isResolved) {
                    for (int bond : ringSet.getRingBonds(i)) {
                        mol.setBondType(bond, Molecule3D.cBondTypeDelocalized);
                    }
                }
            }
        }
        //Aromatizer
        AromaticityResolver resolver = new AromaticityResolver(mol);
        resolver.locateDelocalizedDoubleBonds(null, true, true);
        mol.ensureHelperArrays(Molecule3D.cHelperRings);
        if (DEBUG) {
            for (int i = 0; i < ringSet.getSize(); i++) {
                if ( aromaticRing[i] ) {
                    System.err.println(">> DEBUG: entering ring number " + i );
                    for (int bond : ringSet.getRingBonds(i)) {
                        System.err.println(">>> DEBUG: aromatic bond between atoms " +
                                mol.getAtomName(mol.getBondAtom(0, bond) )+ "' and '" +
                                mol.getAtomName(mol.getBondAtom(1, bond) ) + "' was marked as a " +
                                "bond of order " + mol.getBondOrder(bond));
                    }
                }
            }
            System.err.println("> DEBUG: exiting aromatic ring perception phase");
        }
        assignObviousHigherOrderBonds(mol);

        // Find pi bonds
        if ( DEBUG ) {
            System.err.println("> DEBUG: entering pi bond detection phase");
        }
        boolean[] visitedAtoms = new boolean[mol.getAtoms()];
        List<Integer> piBonds = new LinkedList<>();
        for (int i = 0; i < mol.getAtoms(); i++) {
            // skip atoms that cannot accommodate an additional double bond
            if (mol.getConnAtoms(i) < 2) continue;
            if (spOrder[i] != 2) continue;
            if (mol.getAtomPi(i) != 0) continue;
            if (mol.isAromaticAtom(i)) continue;
            if (mol.getAtomCharge(i) != 0) continue;
            if (mol.getLowestFreeValence(i) < 1) continue;

            visitedAtoms[i] = true;
            for (int j = 0; j < mol.getConnAtoms(i); j++) {
                int a1 = mol.getConnAtom(i, (j + 1) % mol.getConnAtoms(i));
                int a2 = mol.getConnAtom(i, j);

                if (visitedAtoms[a2]) continue;
                // skip atoms that cannot accommodate an additional double bond
                if (spOrder[a2] != 2) continue;
                if (mol.getAtomPi(a2) != 0) continue;
                if (mol.getAtomCharge(a2) != 0) continue;
                if (mol.getLowestFreeValence(a2) < 1) continue;

                if (!fitsAromaticLengthConstraints(mol, i, a2)) continue;

                Coordinates origin = mol.getCoordinates(i);
                Coordinates normal = calculatePlaneNormal(origin, mol.getCoordinates(a1), mol.getCoordinates(a2));
                if (normal.distSq() == 0) continue;

                boolean isFlat = true;
                // TODO: join them in a list?
                for (int k = 0; k < mol.getConnAtoms(i); k++) {
                    if (calculatePointToPlaneDistance(normal, origin, mol.getCoordinates(mol.getConnAtom(i, k))) > 0.4) {
                        isFlat = false;
                        break;
                    }
                }

                for (int k = 0; k < mol.getConnAtoms(a2); k++) {
                    if (calculatePointToPlaneDistance(normal, origin, mol.getCoordinates(mol.getConnAtom(a2, k))) > 0.4) {
                        isFlat = false;
                        break;
                    }
                }

                if (isFlat) {
                    piBonds.add(mol.getConnBond(i, j));
                    if (DEBUG) {
                        System.err.println(">> DEBUG: bond between atoms '" + mol.getAtomName(i) +
                                "' and '" + mol.getAtomName(mol.getConnAtom(i, j)) + "' has been marked " +
                                "as delocalized");
                    }
                }
            }
        }
        if ( DEBUG ) {
            System.err.println("> DEBUG: exiting pi bond detection phase");
        }

        for (int piBond : piBonds) {
            mol.setBondType(piBond, Molecule3D.cBondTypeDelocalized);
        }
        mol.ensureHelperArrays(Molecule3D.cHelperRings);
        List<List<Integer>> fragments = getPiBondFragments(mol);
        for (List<Integer> fragment : fragments) {
            if (!(mol.getLowestFreeValence(fragment.get(0)) > 0 &&
                    mol.getLowestFreeValence(fragment.get(fragment.size() - 1)) > 0)) {
                continue;
            }

            // if the fragment is capable of interchanging single-double bonds
            // then the more likely pattern should be selected
            double[] singleBonds = new double[fragment.size() - 1];
            double[] doubleBonds = new double[fragment.size() - 1];
            for (int j = 0; j < fragment.size() - 1; j++) {
                singleBonds[j] = estimatePiBondFitness(mol, fragment.get(j), fragment.get(j + 1), 1);
                doubleBonds[j] = estimatePiBondFitness(mol, fragment.get(j), fragment.get(j + 1), 2);
            }

            if (fragmentCanResonate(mol, fragment)) {
                double leftToRight = 0;
                double rightToLeft = 0;
                for (int j = 0; j < fragment.size() - 1; j++) {
                    // start with a single bond
                    leftToRight += (j % 2) == 0 ? singleBonds[j] : doubleBonds[j];
                    // start with a double bond
                    rightToLeft += (j % 2) == 0 ? doubleBonds[j] : singleBonds[j];
                }

                if (fragment.size() % 2 == 1) {
                    for (int j = 0; j < fragment.size() - 1; j++) {
                        int bondOrder;
                        if (leftToRight < rightToLeft) {
                            bondOrder = (j % 2 == 0 ? 1 : 2);
                        } else {
                            bondOrder = (j % 2 == 0 ? 2 : 1);
                        }
                        mol.setBondOrder(mol.getBond(fragment.get(j), fragment.get(j + 1)), bondOrder);
                        if (DEBUG) {
                            System.err.println(">>> DEBUG: a bond of order " + bondOrder + " was placed between " +
                                    "atoms '" + mol.getAtomName(fragment.get(j)) + "' and '" +
                                    mol.getAtomName(fragment.get(j + 1)) + "' outside of the aromatizer");
                        }
                    }
                }
            } else {
                for (int j = 0; j < fragment.size() - 1; j++) {
                    if (doubleBonds[j] < 0.03) {
                        mol.setBondOrder(mol.getBond(fragment.get(j), fragment.get(j + 1)), 2);
                    }
                }
            }
        }
        mol.ensureHelperArrays(Molecule3D.cHelperRings);
        resolver = new AromaticityResolver(mol);
        resolver.locateDelocalizedDoubleBonds(null, true, true);
        mol.ensureHelperArrays(Molecule3D.cHelperRings);

        if ( DEBUG ) {
            for (int piBond : piBonds) {
                int a1 = mol.getBondAtom(0, piBond);
                int a2 = mol.getBondAtom(1, piBond);
                System.err.println(">>>> DEBUG: a bond of order " + mol.getBondOrder(piBond) + " was " +
                        "placed between atoms '" + mol.getAtomName(a1) + "' and '" + mol.getAtomName(a2) +
                        "'.");
            }
        }
        // Try to recognise the N-C-N pattern and set the charges accordingly
        int index = BondLengthSet.getBondIndex(1, true, false, 6, 7, 1 , 1);
        double bondLength = BondLengthSet.getBondLength(index);
        double bondStd    = BondLengthSet.getBondStdDev(index);
        for ( int i = 0; i < mol.getAtoms(); i++ ) {
            if ( mol.getAtomicNo(i) == 6 && spOrder[i] == 2 && mol.getConnAtoms(i) > 1 ) {
                for (int j = 0; j < mol.getConnAtoms(i); j++) {
                    int a1 = mol.getConnAtom(i, j);
                    int a2 = mol.getConnAtom(i, (j + 1) % mol.getConnAtoms(i));

                    if (spOrder[a1] != 2 || spOrder[a2] != 2) continue;
                    if (mol.getAtomicNo(a1) != 7 || mol.getAtomicNo(a2) != 7) continue;

                    double a1Distance = getBondLength(mol, i, a1);
                    double a2Distance = getBondLength(mol, i, a2);
                    if ( ( a1Distance > bondLength + 3 * bondStd ) ||
                           mol.getBondOrder(mol.getBond(i, a1)) == 2 ) continue;
                    if ( ( a2Distance > bondLength + 3 * bondStd ) ||
                           mol.getBondOrder(mol.getBond(i, a1)) == 2 ) continue;
                    int shorterBondNitrogen = a1Distance < a2Distance ? a1 : a2;

                    // TODO: add an angle constraint
                    if ( mol.getAtomPi(i) == 0 ) {
                        if ( mol.getAtomCharge(a1) != 0 || mol.getAtomCharge(a2) != 0 ) continue;
                        if ( mol.getAtomPi(a1) != 0 || mol.getAtomPi(a2) != 0 ) continue;

                        if ( mol.getBondOrder(mol.getBond(i, a1)) < 2 &&
                             mol.getBondOrder(mol.getBond(i, a2)) < 2) {
                            mol.setBondOrder(mol.getBond(i, shorterBondNitrogen), 2);
                            mol.setAtomCharge(shorterBondNitrogen, +1);
                            mol.ensureHelperArrays(Molecule.cHelperRings);
                        }
                    } else {
                    // Switch N-C=N to N=C-N in if the latter resonance form seems more likely

                        // exocyclic aromatic resonance should not be considered
                        if ( mol.isAromaticAtom(a1) != mol.isAromaticAtom(a2) ) continue;

                        int longerBondNitrogen = ( shorterBondNitrogen == a1 ) ? a2 : a1;
                        if ( mol.getBondOrder(mol.getBond(i, longerBondNitrogen)) == 2 &&
                             mol.getAtomCharge(longerBondNitrogen) == 1 &&
                             mol.getAtomPi(shorterBondNitrogen) == 0 ) {
                            mol.setBondOrder(mol.getBond(i, longerBondNitrogen), 1);
                            mol.setBondOrder(mol.getBond(i, shorterBondNitrogen), 2);
                            mol.setAtomCharge(longerBondNitrogen, 0);
                            mol.setAtomCharge(longerBondNitrogen, +1);
                            mol.ensureHelperArrays(Molecule.cHelperRings);
                        }
                    }
                }
            }
        }
        mol.ensureHelperArrays(Molecule.cHelperRings);
        /*
         * 2nd pass: find obvious double bonds on sp and sp2 B, C and N outside aromatic rings
         */
        assignObviousDoubleBonds(mol, spOrder);
        /*
         * 3rd pass: double and triple bonds that correspond to the bond length
         * (but might deviate from the expected geometry)
         */
        IntQueue queue = new IntQueue();
        for (int i = 0; i < mol.getAllBonds(); i++) {
            if ( mol.getBondType(i) == Molecule3D.cBondTypeMetalLigand ) continue;
            if (!visited[i] && !mol.isAromaticBond(i)) queue.push(i);
            while (!queue.isEmpty()) {
                int bnd = queue.pop();
                if (visited[bnd]) continue;
                visited[bnd] = true;
                int a1 = mol.getBondAtom(0, bnd);
                int a2 = mol.getBondAtom(1, bnd);

                if ( a1 >= spOrder.length || a2 >= spOrder.length) continue;

                // Push the neighbour bonds into the queue, but
                // avoid processing hydrogen and metal atoms
                for (int j = 0; j < mol.getConnAtoms(a1); j++) {
                    queue.push(mol.getConnBond(a1, j));
                }
                for (int j = 0; j < mol.getConnAtoms(a2); j++) {
                    queue.push(mol.getConnBond(a2, j));
                }

                //Compute the free valence and increase the bond order if needed
                double realBondLength = getBondLength(mol, a1, a2);
                int doubleBondIndex = BondLengthSet.getBondIndex(2, false, false, mol.getAtomicNo(a1), mol.getAtomicNo(a2), 1, 1);
                if ( doubleBondIndex == -1 ) continue;
                double idealisedDoubleBondLength = BondLengthSet.getBondLength(doubleBondIndex);
                double idealisedDoubleBondStdDev = BondLengthSet.getBondStdDev(doubleBondIndex);
                if ( idealisedDoubleBondStdDev == BondLengthSet.DEFAULT_BOND_STDDEV ) {
                    idealisedDoubleBondStdDev = 0.1;
                }

                if ( realBondLength < idealisedDoubleBondLength + 3 * idealisedDoubleBondStdDev &&
                         !mol.isAromaticAtom(a1) && !mol.isAromaticAtom(a2) ) {
                    //Special case CS
                    if (mol.getAtomicNo(a1) == 16 && mol.getAllConnAtoms(a1) <= 2) continue;
                    if (mol.getAtomicNo(a2) == 16 && mol.getAllConnAtoms(a2) <= 2) continue;

                    if ( mol.getAtomCharge(a1) != 0 || mol.getAtomCharge(a2) != 0 ) continue;

                    int freeValence1 = StrictMath.max( mol.getLowestFreeValence(a1),
                                            StrictMath.max( ( mol.getAtomicNo(a1) == 7 && spOrder[a1] < 3 ) ? 1 : 0,
                                                            ( mol.getAtomicNo(a1) == 5 && spOrder[a1] < 3 ) ? 1 : 0 ) );
                    int freeValence2 = StrictMath.max( mol.getLowestFreeValence(a2),
                                            StrictMath.max( ( mol.getAtomicNo(a2) == 7 && spOrder[a2] < 3 ) ? 1 : 0,
                                                            ( mol.getAtomicNo(a2) == 5 && spOrder[a2] < 2 ) ? 1 : 0 ) );

                    if ( freeValence1 < 1 || freeValence2 < 1 ) continue;

                    boolean aligned = ( spOrder[a1] == 1 && spOrder[a2] == 1 );

                    boolean fitsTripleBondLength;
                    int tripleBondIndex = BondLengthSet.getBondIndex(3, false, false, mol.getAtomicNo(a1), mol.getAtomicNo(a2), 2, 2);

                    if ( tripleBondIndex != -1 ) {
                        double idealisedTripleBondLength = BondLengthSet.getBondLength(tripleBondIndex);
                        double idealisedTripleBondStdDev = BondLengthSet.getBondStdDev(tripleBondIndex);
                        if (idealisedTripleBondStdDev == BondLengthSet.DEFAULT_BOND_STDDEV) {
                            idealisedTripleBondStdDev = 0.1;
                        }
                        fitsTripleBondLength = realBondLength > idealisedTripleBondLength + 3 * idealisedTripleBondStdDev;
                    } else {
                        fitsTripleBondLength = realBondLength < idealisedDoubleBondLength;
                    }

                    if ( aligned && mol.getAtomPi(a1) < 2 && mol.getAtomPi(a2) < 2 && fitsTripleBondLength ) {
                        mol.setBondOrder(bnd, 3);
                    } else if ( mol.getAtomPi(a1) < ( spOrder[a1] == 1 ? 2 : 1 ) &&
                                mol.getAtomPi(a2) < ( spOrder[a2] == 1 ? 2 : 1 ) ) {
                        mol.setBondOrder(bnd, 2);
                    }
                }
            }
        }
        mol.ensureHelperArrays(Molecule3D.cHelperRings);
        /*
         * 4th pass: increase bond order if the geometry requires it
         * (even though the bond length might seem too great)
         */
        int[] localAtomPi = getAtomPiCount(mol);
        for (int i = 0; i < mol.getAtoms(); i++) {
            if (!canAccommodateDoubleBond(mol, i, localAtomPi[i], spOrder[i])) continue;
            // sp3 atoms smaller than Si should not have double bonds
            if (spOrder[i] == 3 && mol.getAtomicNo(i) < 14) continue;
            for (int j = 0; j < mol.getConnAtoms(i); j++) {
                int a1 = mol.getConnAtom(i, j);
                if (canAccommodateDoubleBond(mol, a1, localAtomPi[a1], spOrder[a1])) {
                    mol.setBondOrder(mol.getBond(i, a1), 2);
                    localAtomPi[i] += 1;
                    localAtomPi[a1] += 1;
                }
                if (localAtomPi[i] >= (spOrder[i] == 1 ? 2 : 1)) break;
            }
        }

        // 5th pass: redistributing charges
        /*
         * v0 is the start atom
         * v1 is the 1st order neighbour of atom v0
         * v2 is the 2nd order neighbour of atom v0
         * v3 is the 3rd order neighbour of atom v0
         */
        localAtomPi = getAtomPiCount(mol);
        for (int v0 = 0; v0 < mol.getAllAtoms(); v0++) {
            if ( mol.getAtomicNo(v0) == 6 ) {
                // C[-]-
                if ( spOrder[v0] != 2 ) continue;
                if ( mol.getLowestFreeValence(v0) != 1 ) continue;
                if ( localAtomPi[v0] != 0 ) continue;
                if ( hasMetalLigandBonds(mol, v0) ) continue;
                // O=C-C(-)-C=O -> O=C-C=C-O(-)
                List<Integer> terminalOxygens = new LinkedList<>();
                // C(-)-N -> C->N(+)
                List<Integer> sp2NitrogenNeighbours = new LinkedList<>();
                for (int j = 0; j < mol.getConnAtoms(v0); j++) {
                    int v1 = mol.getConnAtom(v0, j);
                    if ( mol.getAtomicNo(v1) == 6 ) {
                        // C[-]-C(R)=
                        if ( spOrder[v1] != 2 ) continue;
                        if ( localAtomPi[v1] == 0 ) continue;

                        for (int k = 0; k < mol.getConnAtoms(v1); k++) {
                            int v2 = mol.getConnAtom(v1, k);
                            // select the second order neighbour atom that is double bonded
                            // to the first order neighbour atom
                            if ( v2 == v0 ) continue;
                            if ( mol.getBondOrder(mol.getConnBond(v1, k)) != 2 ) continue;
                            if ( mol.getAtomicNo(v2) == 8 &&
                                 mol.getConnAtoms(v2) == 1 ) {
                                terminalOxygens.add(v2);
                            } else if ( mol.getAtomicNo(v2) == 6 ) {
                                // C[-]-C(R)=C(R)(R)
                                for ( int l = 0; l < mol.getConnAtoms(v2); l++ ) {
                                    if ( localAtomPi[v0] > 0 ) break;
                                    int v3 = mol.getConnAtom(v2, l);
                                    if ( v1 == v3 ) continue;
                                    if ( mol.getAtomicNo(v3) != 6 ) continue;
                                    if ( spOrder[v3] != 2 ) continue;
                                    if ( mol.getLowestFreeValence(v3) != 1 ) continue;
                                    if ( hasMetalLigandBonds(mol, v3) )  continue;

                                    mol.setBondOrder(mol.getBond(v0, v1), 2);
                                    mol.setBondOrder(mol.getBond(v1, v2), 1);
                                    mol.setBondOrder(mol.getBond(v2, v3), 2);
                                    localAtomPi[v0]++;
                                    localAtomPi[v3]++;
                                    break;
                                }
                            }
                        }
                    } else if (mol.getAtomicNo(v1) == 7 && spOrder[v1] == 2 &&
                            mol.getAtomPi(v1) == 0 && mol.getImplicitHydrogens(v1) == 0) {
                        sp2NitrogenNeighbours.add(j);
                    }
                }

                // FIXME: the valency of associated C atoms might change due to previous loop
                // This needs serious refactoring
                // moving charges from carbon atoms to oxygen atoms where possible
                if ( mol.getLowestFreeValence(v0) > 0 ) {
                    if (terminalOxygens.size() > 0 ) {
                        mol.setBondOrder(mol.getConnBond(terminalOxygens.get(0), 0), 1);
                        mol.setBondOrder(mol.getBond(mol.getConnAtom(terminalOxygens.get(0), 0), v0), 2);
                    } else if (sp2NitrogenNeighbours.size() > 0) {
                        int[] sortedAtom = sortAtomsByBondLength(mol, v0, sp2NitrogenNeighbours);
                        mol.setBondOrder(mol.getConnBond(v0, sortedAtom[0]), 2);
                        // nitrogen atoms may contain a negative charge that was assigned during previous processing
                        int nAtom = mol.getConnAtom(v0, sortedAtom[0]);
                        if (mol.getAtomCharge(nAtom) < 0) {
                            mol.setAtomCharge(nAtom, mol.getAtomCharge(nAtom) + 1);
                        }
                    }
                }
            } else if ( mol.getAtomicNo(v0) == 16 ) {
                // Sp2 sulphur with +1 charge with no metal bonds
                if ( spOrder[v0] != 2 ) continue;
                if ( mol.getAtomPi(v0) == 0 ) continue;
                if ( mol.getAtomCharge(v0) != 1 ) continue;
                if ( hasMetalLigandBonds(mol, v0) ) continue;
                // (S+)=C-C=(S+) -> S-C=C-S
                // (S+)=C-C -> (S+)-C=C
                for (int j = 0; j < mol.getConnAtoms(v0); j++) {
                    // Sp2 carbon that is double bonded to the v0 sulphur
                    int v1 = mol.getConnAtom(v0, j);
                    if ( mol.getAtomicNo(v1) != 6 ) continue;
                    if ( mol.getBondOrder( mol.getBond(v0, v1) ) != 2 ) continue;

                    for (int k = 0; k < mol.getConnAtoms(v1); k++) {
                        // Sp2 carbon that is single bonded to the v1 carbon
                        int v2 = mol.getConnAtom(v1, k);
                        if ( mol.getAtomicNo(v2) != 6 ) continue;
                        if ( mol.getBondOrder(mol.getBond(v1,v2)) != 1 ) continue;

                        if ( mol.getLowestFreeValence(v2) == 1 ) {
                            mol.setBondOrder(mol.getBond(v1, v2), 2);
                            mol.setBondOrder(mol.getBond(v0, v1), 1);
                            mol.setAtomCharge(v0, 0);
                        } else if (mol.getAtomPi(v2) != 0){
                            for (int l = 0; l < mol.getConnAtoms(v2); l++) {
                                int v3 = mol.getConnAtom(v2, l);
                                if ( mol.getAtomicNo(v3) != 16 ) continue;
                                if ( mol.getAtomCharge(v3) != 1 ) continue;
                                if ( mol.getBondOrder(mol.getBond(v2, v3)) != 2 ) continue;

                                mol.setBondOrder(mol.getBond(v0, v1), 1);
                                mol.setBondOrder(mol.getBond(v1, v2), 2);
                                mol.setBondOrder(mol.getBond(v2, v3), 1);
                                mol.setAtomCharge(v0, 0);
                                mol.setAtomCharge(v3, 0);
                            }
                        }
                    }
                }
            }
        }
        mol.ensureHelperArrays(Molecule.cHelperRings);

        // Set charges for atoms other than metals
        assignAtomChargesBasedOnValence(mol, spOrder);

        /*
         * Transform the molecule to acquire a more reasonable resonance structure
         */
        // R1(+)=X-R2(-) -> R1-X=R2
        for (int i = 0; i < mol.getAtoms(); i++) {
            if ( mol.getAtomCharge(i) != 1 ) continue;
            if ( mol.getAtomPi(i) == 0 ) continue;

            for ( int j = 0; j < mol.getConnAtoms(i); j++ ) {
                int a1 = mol.getConnAtom(i, j);
                if ( mol.getAtomCharge(a1) == 0 && mol.getConnBondOrder(i, j) == 2 ) {
                    for (int k = 0; k < mol.getConnAtoms(a1); k++) {
                        int a2 = mol.getConnAtom(a1, k);
                        double R1XBondLengthDiff = BondLengthSet.getBondLength(
                                BondLengthSet.getBondIndex(2, false, false,
                                        mol.getAtomicNo(i),
                                        mol.getAtomicNo(a1), 1, 1) ) - getBondLength(mol, i, a1);
                        if ( mol.getAtomCharge(i) < 0 && spOrder[a2] == 2 &&
                                mol.getAtomPi(a2) == 0 &&
                                ( mol.getAtomCharge(i) + mol.getAtomCharge(a2) ) == 0 ) {
                            double R2XBondLengthDiff = BondLengthSet.getBondLength(
                                    BondLengthSet.getBondIndex(2, false, false,
                                            mol.getAtomicNo(a1),
                                            mol.getAtomicNo(a2), 1, 1)) - getBondLength(mol, a1, a2);
                            if ( R2XBondLengthDiff < R1XBondLengthDiff ) {
                                mol.setBondOrder(mol.getBond(a1, a2), 2);
                                mol.setBondOrder(mol.getBond(i, a1), 1);
                                mol.setAtomCharge(i, 0);
                                mol.setAtomCharge(a2, 0);
                                break;
                            }
                        }
                    }
                }
            }
        }

        processKnownIonsInvolvingMetals(mol);

        assignMetalCharges(mol);

        mol.ensureHelperArrays(Molecule.cHelperNeighbours);
    }

    /**
     * Gets the count of used pi electrons for each atom in the molecule.
     * @param mol
     *          The molecule that contains the bonded atoms.
     * @return
     *          Counts of used pi electron.
     */
    private static int[] getAtomPiCount(Molecule3D mol) {
        int[] atomPiCount = new int[mol.getAllAtoms()];
        for ( int i = 0; i < atomPiCount.length; i++ ) {
            if (mol.getAtomicNo(i) == 1) {
                atomPiCount[i] = 0;
            } else {
                atomPiCount[i] = mol.getAtomPi(i);
            }
        }
        return atomPiCount;
    }

    /**
     * Assigns obvious higher order bonds where the geometry explicitly requires it.
     *
     * @param mol
     *          Molecule that contains the bonded atoms.
     */
    private static void assignObviousHigherOrderBonds(Molecule3D mol) {
        // Place obvious triple bonds for sp atoms
        boolean changeOccurred = true;
        int tripleBondIterationNo = 1;
        while (changeOccurred) {
            if (DEBUG) {
                System.err.println(">> DEBUG: entered obvious triple bond assignment iteration " +
                        tripleBondIterationNo + ".");
            }
            changeOccurred = assignObviousTripleBonds(mol);
            if (DEBUG) {
                System.err.println(">> DEBUG: exited obvious triple bond assignment iteration " +
                        tripleBondIterationNo + ".");
            }
            tripleBondIterationNo++;
        }
    }

    /**
     * Assigns obvious triple bonds where the geometry explicitly requires it.
     *
     * @param mol
     *          Molecule that contains the bonded atoms.
     * @return
     *          Boolean value denoting if any triple bonds were assigned.
     */
    private static boolean assignObviousTripleBonds(Molecule3D mol) {
        // the local pi count array is used to avoid calling ensureHelperArrays()
        // each time a triple bond is assigned
        boolean changeOccured = false;
        int[] localAtomPi = getAtomPiCount(mol);
        for (int i = 0; i < mol.getAtoms(); i++) {
            if (mol.isMetalAtom(i)) continue;
            if (localAtomPi[i] != 0) continue;
            if (calculateHybridizationState(mol, i) != 1) continue;
            if (mol.getAllConnAtoms(i) != 2) continue; // skip terminal atoms

            boolean allowNeighbourCharge = true;
            if (!canAccommodateTripleBond(mol, i, false)) {
                if (!canAccommodateTripleBond(mol, i, true)) continue;
                allowNeighbourCharge = false;
            }

            int a1 = mol.getConnAtom(i, 0);
            int a2 = mol.getConnAtom(i, 1);

            // skip if both atoms are capable of at least a single double bond
            boolean a1CanCarryDoubleBond = canAccommodateDoubleBond(mol, a1, localAtomPi[a1], calculateHybridizationState(mol, a1));
            boolean a2CanCarryDoubleBond = canAccommodateDoubleBond(mol, a2, localAtomPi[a2], calculateHybridizationState(mol, a2));
            boolean a1CanCarryTripleBond = canAccommodateTripleBond(mol, a1, allowNeighbourCharge);
            boolean a2CanCarryTripleBond = canAccommodateTripleBond(mol, a2, allowNeighbourCharge);

            if ((a1CanCarryDoubleBond || a1CanCarryTripleBond) &&
                    (a2CanCarryDoubleBond || a2CanCarryTripleBond)) {
                continue;
            }

            if (!a1CanCarryTripleBond && !a2CanCarryTripleBond) continue;
            int tripleBondAtom = a1CanCarryTripleBond ? a1 : a2;

            mol.setBondOrder(mol.getBond(i, tripleBondAtom), 3);
            localAtomPi[i] = 2;
            localAtomPi[tripleBondAtom] = 2;
            changeOccured = true;
            if (DEBUG) {
                System.err.println(">>> DEBUG: placed a triple bond between atoms '" +
                        mol.getAtomName(i) + "' and '" + mol.getAtomName(tripleBondAtom) +
                        "' without relying on the bond length statistics.");
            }
        }
        mol.ensureHelperArrays(Molecule.cHelperRings);

        return changeOccured;
    }

    /**
     * Evaluates is the given atom is capable of carrying a triple bond.
     *
     * @param mol
     *          Molecule that contains the atom.
     * @param atom
     *          Index of the atom that should be evaluated.
     * @param allowCharge
     *          Boolean value denoting if the atom is allowed to gain
     *          a formal charge in order to accommodate the triple bond
     *          (i.e. N+, B-).
     * @return
     *          Boolean value denoting if the atom is capable of carrying a triple bond.
     */
    private static boolean canAccommodateTripleBond(Molecule3D mol, int atom, boolean allowCharge) {
        if (mol.getAtomicNo(atom) == 1) return false;
        if (mol.isHalogene(atom)) return false;

        int requiredValence = 2;
        if (allowCharge) {
            // boron (B) with the charge of -1 may carry a triple bond
            if (mol.getAtomicNo(atom) == 5 ) requiredValence = 1;
            // nitrogen (N) with the charge of +1 may carry a triple bond
            if (mol.getAtomicNo(atom) == 7 ) requiredValence = 1;
        }
        if (mol.getLowestFreeValence(atom) < requiredValence) return false;

        int spOrder = calculateHybridizationState(mol, atom);
        return spOrder == 1;
    }

    /**
     * Assigns obvious double bonds to boron, nitrogen and carbon atoms.
     *
     * @param mol
     *          The molecule that contains the bonded atoms.
     * @param spOrder
     *          The sp hybridization type for each atom.
     */
    private static void assignObviousDoubleBonds(Molecule3D mol, int[] spOrder) {
        // using the local pi count array to avoid calling ensureHelperArrays
        // each time a double bond is assigned
        int[] atomPi = getAtomPiCount(mol);
        for (int i = 0; i < mol.getAtoms(); i++) {
            if ( mol.isAromaticAtom(i) ) continue;
            if ( mol.getAtomCharge(i) != 0 ) continue;
            // B, N, C with valence lower than 4
            if ( mol.getAtomicNo(i) > 3 && mol.getAtomicNo(i) < 8 && mol.getOccupiedValence(i) < 4 ) {
                if ( spOrder[i] == 1 ) {
                    for (int j = 0; j < mol.getConnAtoms(i); j++) {
                        int a1 = mol.getConnAtom(i, j);
                        if ( mol.getLowestFreeValence(a1) > 0 &&
                                mol.getAtomCharge(a1) == 0 &&
                                spOrder[a1] < 3 && atomPi[a1] <
                                ( spOrder[a1] == 1 ? 2 :
                                        spOrder[a1] == 2 ? 1 : 0 ) ) {
                            mol.setBondOrder(mol.getConnBond(i, j), 2);
                            atomPi[i]++;
                            atomPi[a1]++;
                        }
                    }
                } else if (spOrder[i] == 2) {
                    if (atomPi[i] != 0) continue;
                    List<Integer> potentialDoubleBondNeighbours = new LinkedList<>();
                    // nitrogen atoms with filled valencies are handled a little differently
                    boolean isFilledNitrogenAtom =  mol.getAtomicNo(i) == 7 && mol.getLowestFreeValence(i) == 0;
                    for (int j = 0; j < mol.getConnAtoms(i); j++) {
                        int neighbour = mol.getConnAtom(i, j);
                        if ( mol.getAtomCharge(neighbour) != 0 ) continue;
                        int freeValence = StrictMath.max( mol.getLowestFreeValence(neighbour),
                                mol.getAtomicNo(neighbour) == 5 && spOrder[neighbour] < 3  ? 1 : 0 );
                        if ( freeValence < 1 ) continue;
                        if ( isFilledNitrogenAtom && mol.getAtomicNo(neighbour) == 7 && spOrder[neighbour] == 2 ) continue;
                        if ( atomPi[neighbour] < ( spOrder[neighbour] == 1 ? 2 :
                                spOrder[neighbour] == 2 ? 1 : 0 ) ||
                                mol.getAtomicNo(neighbour) == 15  )  {
                            potentialDoubleBondNeighbours.add(j);
                        }
                    }

                    if ( potentialDoubleBondNeighbours.size() > 1 ) {
                        double minDiff = -Double.MAX_VALUE;
                        int minDiffNeighbour = -1;
                        for (Integer neighbourIndex : potentialDoubleBondNeighbours) {
                            int neighbour = mol.getConnAtom(i, neighbourIndex);
                            double singleBondLength =BondLengthSet.getBondLength(
                                    BondLengthSet.getBondIndex(1, false, false, mol.getAtomicNo(i), mol.getAtomicNo(neighbour), 1, 1));
                            double doubleBondLength = BondLengthSet.getBondLength(
                                    BondLengthSet.getBondIndex(2, false, false, mol.getAtomicNo(i), mol.getAtomicNo(neighbour), 1, 1));
                            // TODO: think of a more legitimate way of detecting statistical anomalies
                            // safeguard against statistical anomalies
                            if ( singleBondLength < doubleBondLength || singleBondLength - doubleBondLength < 0.1 ) {
                                continue;
                            }
                            double bondLengthDifference = doubleBondLength - getBondLength(mol, i, neighbour);
                            if (bondLengthDifference > minDiff) {
                                minDiff = bondLengthDifference;
                                minDiffNeighbour = neighbour;
                            }
                        }
                        if (minDiffNeighbour > -1) {
                            mol.setBondOrder(mol.getBond(i, minDiffNeighbour), 2);
                            atomPi[i]++;
                            atomPi[minDiffNeighbour]++;
                        }
                    }
                }
            }
        }
        mol.ensureHelperArrays(Molecule.cHelperRings);
    }

    /**
     * Assigns obvious disilene bonds.
     *
     * @param mol
     *          The molecule that contains the bonded atoms.
     */
    private static void assignDisileneBonds(Molecule3D mol) {

        double doubleBondCutoff = 2.29;
        // FIXME: properly check if a normal index is returned
        int bondIndex = BondLengthSet.getBondIndex(2, false, false, 14, 14, 1, 1 );
        double siliconDoubleBondLength = BondLengthSet.getBondLength( bondIndex );
        double siliconDoubleBondStdDev = BondLengthSet.getBondStdDev( bondIndex );
        if ( siliconDoubleBondLength != -1 && siliconDoubleBondStdDev != -1 ) {
            if ( ( siliconDoubleBondLength + 3 * siliconDoubleBondStdDev ) < doubleBondCutoff ) {
                doubleBondCutoff = siliconDoubleBondLength + 3 * siliconDoubleBondStdDev;
            }
        }

        int[] atomPi = getAtomPiCount(mol);
        for (int i = 0; i < mol.getAtoms(); i++) {
            if ( mol.getAtomicNo(i) != 14 ) continue;
            if ( atomPi[i] != 0 ) continue;
            if ( mol.getOccupiedValence(i) > 3 ) continue;
            if ( mol.isAromaticAtom(i) ) continue;
            for (int j = 0; j < mol.getConnAtoms(i); j++) {
                int a1 = mol.getConnAtom(i, j);
                // ring bonds are handled in aromaticity perception stage
                if ( mol.isRingBond( mol.getBond( i, a1 ) ) ) continue;
                if ( mol.getAtomicNo(a1) != 14 ) continue;
                if ( mol.getOccupiedValence(a1) > 3 ) continue;
                if ( atomPi[a1] != 0 ) continue;
                if ( mol.isAromaticAtom(a1) ) continue;
                if ( getBondLength( mol, i, a1 ) < doubleBondCutoff ) {
                    mol.setBondOrder( mol.getBond(i, a1), 2 );
                    atomPi[i]++;
                    atomPi[a1]++;
                    break;
                }
            }
        }
    }

    /**
     * Detects and assigns exocyclic double bonds that are more likely than
     * double bonds inside the rings.
     * @param mol
     *          Molecule that contains the atoms.
     * @param spOrder
     *          The sp hybridization type for each atom.
     * @return
     *          Boolean values denoting is an exocyclic bond has been placed
     *          for the corresponding atom ring.
     */
    private static boolean[] assignExocyclicDoubleBonds(Molecule3D mol, int[] spOrder) {
        boolean[] hasExocyclicDoubleBond = new boolean[mol.getRingSet().getSize()];
        int[] localAtomPi = getAtomPiCount(mol);
        for (int i = 0; i < mol.getRingSet().getSize(); i++) {
            int[] ringAtoms = mol.getRingSet().getRingAtoms(i);
            for (int ringAtom : ringAtoms) {
                // skip the atom if the only exocyclic neighbours are hydrogen atoms
                if (mol.getConnAtoms(ringAtom) == 2) continue;

                // additional exclusion rules for certain elements
                if (mol.getAtomicNo(ringAtom) == 14 ||      // Silicon (Si)
                        mol.getAtomicNo(ringAtom) == 15 ||  // Phosphorus (P)
                        mol.getAtomicNo(ringAtom) == 32) {  // Germanium (Ge)
                    // skip atoms with full valency
                    if (mol.getFreeValence(ringAtom) < 1) continue;
                    // skip atoms that already contain a double bond
                    if (localAtomPi[ringAtom] > 0) continue;
                }

                // sp3d2 P and Te
                if (spOrder[ringAtom] == 5 &&
                        (mol.getAtomicNo(ringAtom) == 15 || mol.getAtomicNo(ringAtom) == 52)) continue;
                // skip the atom if it is not capable of a double bond
                if (mol.getLowestFreeValence(ringAtom) <= 0 && mol.getAtomicNo(ringAtom) < 14) {
                    continue;
                }
                // find the best fit for an aromatic double bond inside the ring
                double inRingMinDeviation = Double.MAX_VALUE;
                for (int j = 0; j < mol.getConnAtoms(ringAtom); j++) {
                    int ringAtomNeighbour = mol.getConnAtom(ringAtom, j);
                    if (isRingAtom(mol, i, ringAtomNeighbour)) {
                        int bondTypeIndex = BondLengthSet.getBondIndex(2, true, false, mol.getAtomicNo(ringAtom),
                                mol.getAtomicNo(ringAtomNeighbour), 1, 1);
                        /*
                         * FIXME: the logic used for exocyclic bonds involving silicon (Si) should be generalised.
                         * In the future the logic should be generalised for all elements to initially try
                         * to retrieve the aromatic double bond distance and if it is not available fall back
                         * to the simple double bond one.
                         */
                        if (mol.getAtomicNo(ringAtom) == 14 || mol.getAtomicNo(ringAtomNeighbour) == 14) {
                            if (bondTypeIndex == -1) {
                                bondTypeIndex = BondLengthSet.getBondIndex(2, false, false, mol.getAtomicNo(ringAtom),
                                        mol.getAtomicNo(ringAtomNeighbour), 1, 1);
                            }
                        }
                        double currentDeviation = determineExpectedBondDeviation(mol,
                                mol.getBond(ringAtom, ringAtomNeighbour),
                                bondTypeIndex);
                        if (currentDeviation < inRingMinDeviation) {
                            inRingMinDeviation = currentDeviation;
                        }
                    }
                }

                for (int k = 0; k < mol.getConnAtoms(ringAtom); k++) {
                    // check if the ring atom is still capable of
                    // carrying an additional double bond
                    if (mol.getAtomicNo(ringAtom) == 15) { // phosphorus (P)
                        if (localAtomPi[ringAtom] > 0) break;
                    }
                    // check if atom is exocyclic and capable of a double bond
                    int exoAtom = mol.getConnAtom(ringAtom, k);

                    if (mol.isRingBond(mol.getBond(ringAtom, exoAtom))) continue;
                    if (mol.getAtomCharge(exoAtom) != 0) continue;
                    if (mol.getLowestFreeValence(exoAtom) < 1) continue;

                    // additional exclusion rules for certain elements
                    if (mol.getAtomicNo(exoAtom) == 14 ||      // Silicon (Si)
                            mol.getAtomicNo(exoAtom) == 15 ||  // Phosphorus (P)
                            mol.getAtomicNo(exoAtom) == 32) {  // Germanium (Ge)
                        // skip atoms with full valency
                        if (mol.getFreeValence(exoAtom) < 1) continue;
                        // skip atoms that already contains a double bond
                        if (localAtomPi[exoAtom] > 0) continue;
                    }
                    // exocyclic double-bonded atoms should normally be planar,
                    // however, some large atoms may disrupt the planarity
                    if (!isPlanar(mol, ringAtom, exoAtom) &&
                            mol.getAtomicNo(ringAtom) < 14 &&
                            mol.getAtomicNo(exoAtom) < 14) {
                        continue;
                    }
                    // FIXME: pentavalent phosphorus (P) should be considered conjungated,
                    // but not strictly aromatic as stated in 10.1016/j.ccr.2015.02.010
                    if ((spOrder[exoAtom] == 2 || getAllConnAtoms(mol, exoAtom, true) == 1)) {
                        // check exocyclic atom bond length
                        int minPiCount = (getAllConnAtoms(mol, exoAtom, true) > 1) ? 1 : 0;
                        double averageSingleBondLength = getIdealisedBondLength(1, mol.getAtomicNo(ringAtom), mol.getAtomicNo(exoAtom), 1, minPiCount);
                        double averageSingleBondStd = BondLengthSet.getBondStdDev(BondLengthSet.getBondIndex(1, false, false, mol.getAtomicNo(ringAtom), mol.getAtomicNo(exoAtom), 1, minPiCount));
                        double averageDoubleBondLength = getIdealisedBondLength(2, mol.getAtomicNo(ringAtom), mol.getAtomicNo(exoAtom), 1, 1);
                        double averageDoubleBondStd = BondLengthSet.getBondStdDev(BondLengthSet.getBondIndex(2, false, false, mol.getAtomicNo(ringAtom), mol.getAtomicNo(exoAtom), 1, 1));
                        double calculatedLength = getBondLength(mol, ringAtom, exoAtom);
                        double singleBondDiff = calculatedLength - averageSingleBondLength;
                        double doubleBondDiff = calculatedLength - averageDoubleBondLength;

                        if (averageSingleBondLength < averageDoubleBondLength) {
                            if (DEBUG) {
                                System.err.println(">>> DEBUG: single '" + mol.getAtomicNo(ringAtom) + "-" +
                                        mol.getAtomicNo(exoAtom) + "' bond seems to be shorter than a the " +
                                        "double bond of the same type -- skipping the estimation of this bond.");
                            }
                            continue;
                        }

                        // add an exocyclic double bond is the exocyclic atom
                        // has no alternative way to fill its valency
                        if (mol.getConnAtoms(exoAtom) > 1 && !hasMetalLigandBonds(mol, exoAtom)) {
                            int doubleBondNeighbours = 0;
                            for (int j = 0; j < mol.getConnAtoms(exoAtom); j++) {
                                if (mol.getLowestFreeValence(mol.getConnAtom(exoAtom, j)) != 0) {
                                    doubleBondNeighbours++;
                                }
                            }
                            if (doubleBondNeighbours == 1 && doubleBondDiff <= 3 * averageDoubleBondStd) {
                                mol.setBondOrder(mol.getBond(ringAtom, exoAtom), 2);
                                localAtomPi[ringAtom] += 1;
                                localAtomPi[exoAtom] += 1;
                                continue;
                            }
                        }

                        // skip the atom if we are more likely to find a double bond inside the ring
                        if (doubleBondDiff > inRingMinDeviation) continue;

                        if (doubleBondDiff < 0 || doubleBondDiff <= 1.5 * averageDoubleBondStd
                                || StrictMath.abs(singleBondDiff) > 2 * averageSingleBondStd) {
                            mol.setBondOrder(mol.getBond(ringAtom, exoAtom), 2);
                            localAtomPi[ringAtom] += 1;
                            localAtomPi[exoAtom] += 1;
                            hasExocyclicDoubleBond[i] = true;
                            if (DEBUG) {
                                System.err.println(">>> DEBUG: placing a double bond between atoms '" +
                                        mol.getAtomName(ringAtom) + "' and '" + mol.getAtomName(exoAtom) + "'.");
                            }
                        }
                    }
                }
            }
        }

        mol.ensureHelperArrays(Molecule3D.cHelperRings);
        return hasExocyclicDoubleBond;
    }

    /**
     * Assigns positive charges to metal atoms in well-known ions.
     * @param mol
     *          The molecule that contains the metal atoms.
     */
    private static void processKnownIonsInvolvingMetals(Molecule3D mol) {
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            if (!mol.isMetalAtom(i)) continue;

            if (resolveMetalHexafluoride(mol, i)) continue;
            if (resolveGoldInLinearCoordination(mol, i)) continue;

            // Aluminium (Al)
            if (mol.getAtomicNo(i) == 13) {
                // AlCl4(-)
                if (mol.getAllConnAtomsPlusMetalBonds(i) == 4) {
                    boolean isIon = true;
                    for (int j = 0; j < mol.getAllConnAtomsPlusMetalBonds(i); j++) {
                        int connAtom = mol.getConnAtom(i, j);
                        // is terminal Cl atom
                        if (!(mol.getAtomicNo(connAtom) == 17 && mol.getAllConnAtoms(connAtom) == 0)) {
                            isIon = false;
                            break;
                        }
                    }
                    if (isIon) {
                        mol.setAtomCharge(i, 3);
                    }
                }
            }
        }
    }

    /**
     * Evaluates if the given atom serves as the central metal atom of a hexafluoride complex
     * and assigns appropriate charges if needed. Only hexafluorides that can be resolved
     * unambiguously are modified.
     *
     * @param mol
     *          Molecule that contains the metal atom.
     * @param metalAtom
     *          Metal atom to be evaluated.
     * @return
     *          Boolean value denoting if the charge of the metal atom was modified.
     */
    private static boolean resolveMetalHexafluoride(Molecule3D mol, int metalAtom) {
        int charge;
        switch (mol.getAtomicNo(metalAtom)) {
            case 13: // AlF6(3-)
            case 26: // FeF6(3-)
                charge = +3;
                break;
            case 51: // SbF6(-1)
            case 79: // AuF6(-1)
                charge = +5;
                break;
            // -1 is treated as a special value
            default:
                charge = -1;
        }
        if (charge == -1) return false;

        if (!isCentralHexafluoriteIonAtom(mol, metalAtom)) return false;

        mol.setAtomCharge(metalAtom, charge);
        return true;
    }

    /**
     * Evaluates if the given atom is a gold atom in a linear coordination
     * and assigns appropriate charges if needed.
     *
     * NOTE: gold atoms appearing in linear coordination most likely have
     * the oxidation state of +1.
     *
     * @param mol
     *          Molecule that contains the metal atom.
     * @param metalAtom
     *          Metal atom to be evaluated.
     * @return
     *          Boolean value denoting if the charge of the metal atom was modified.
     */
    private static boolean resolveGoldInLinearCoordination(Molecule3D mol, int metalAtom) {
        if (mol.getAtomicNo(metalAtom) != 79) return false; // Gold (Au)
        if (mol.getAllConnAtomsPlusMetalBonds(metalAtom) != 2) return false;

        int n1 = mol.getConnAtom(metalAtom, 0);
        int n2 = mol.getConnAtom(metalAtom, 1);
        // 2.80 rad ~ 160 deg
        if (getAngle(mol, n1, metalAtom, n2) < 2.80) return false;
        mol.setAtomCharge(metalAtom, +1);

        return true;
    }

    /**
     * Evaluates if the given atom serves as the central atom of a hexafluoride complex.
     * @param mol
     *          Molecule that contains the atom.
     * @param atom
     *          Atom to be evaluated.
     * @return
     *          Boolean value denoting if given atom serves as the central atom of
     *          a hexafluoride complex
     */
    static boolean isCentralHexafluoriteIonAtom(Molecule3D mol, int atom) {
        boolean isIon = true;
        if (mol.getAllConnAtomsPlusMetalBonds(atom) != 6) return false;
        for (int j = 0; j < mol.getAllConnAtomsPlusMetalBonds(atom); j++) {
            int connAtom = mol.getConnAtom(atom, j);
            // is terminal F atom
            if (!(mol.getAtomicNo(connAtom) == 9 && mol.getAllConnAtoms(connAtom) == 0)) {
                isIon = false;
                break;
            }
        }

        return isIon;
    }

    private static void assignAtomChargesBasedOnValence(Molecule3D mol, int[] spOrder) {
        // Set charges for non-metal atoms
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            if (mol.isMetalAtom(i)) continue;

            // FIXME: check and properly test if the ( mol.getAtomCharge(i) != 0 )
            // can be refactored to the start of this loop

            // only apply this set of rules to atoms that do not yet have a charge
            if (mol.getAtomCharge(i) == 0) {
                if (mol.getAtomicNo(i) == 5) { // boron (B)
                    if (mol.getOccupiedValence(i) == 4) {
                        mol.setAtomCharge(i, -1);
                    }
                    // } else if ( getAllConnAtoms(mol, i, true) == 6 ) {
                    //       mol.setAtomAbnormalValence(i, 6);
                    //   }
                }

                // Phoshorus (P)
                if (mol.getAtomicNo(i) == 15) {
                    // Phoshorus (P) atoms with 4 neighbours and no higher order bonds have a +1 formal charge
                    if (getAllConnAtoms(mol, i, true) == 4 && mol.getAtomPi(i) == 0) {
                        mol.setAtomCharge(i, +1);
                    }
                    // Phosphorus (P) atoms with the valence of 6 have a -1 formal charge
                    if (mol.getOccupiedValence(i) == 6) {
                        mol.setAtomCharge(i, -1);
                    }
                }

                // Sulphur (S)
                if (mol.getAtomicNo(i) == 16) {
                    // Sulphur (S) sp3 atoms with 3 neighbours have a +1 formal charge
                    if (spOrder[i] == 3 && getAllConnAtoms(mol, i, true) == 3) {
                        mol.setAtomCharge(i, +1);
                    }
                }

                // Arsernic (As) atoms with 4 neighbours and no double bonds have a +1 formal charge
                if (mol.getAtomicNo(i) == 33) {
                    if (getAllConnAtoms(mol, i, true) == 4 && mol.getAtomPi(i) == 0) {
                        mol.setAtomCharge(i, +1);
                    }
                }
            }

            // Hydrogen (H) ions
            if (mol.getAtomicNo(i) == 1 && mol.getAllConnAtoms(i) == 0) {
                if (mol.getAllConnAtomsPlusMetalBonds(i) > 0 &&
                        mol.isMetalAtom(mol.getConnAtom(i, 0))) {
                    // Hydride anion H- (usually bonded to metals)
                    mol.setAtomCharge(i, -1);
                } else {
                    // Hydron cation H+
                    mol.setAtomCharge(i, +1);
                }
            }

            // Removing implicit hydrogens by setting a negative charge.
            // '_atom_site_attached_hydrogens' data item values are taken
            // into consideration
            if (mol.getImplicitHydrogens(i) > 0 || mol.getAttachedHydrogenCount(i) > 0) {
                // sp2 carbon with three non-hydrogen neighbours should have a positive charge
                if (mol.getAtomicNo(i) == 6 && spOrder[i] == 2 &&
                        mol.getConnAtoms(i) == 3 && mol.getAtomPi(i) == 0
                        && mol.getAtomCharge(i) == 0 && !mol.isAromaticAtom(i)) {
                    mol.setAtomCharge(i, +1);
                } else {
                    mol.setAtomCharge(i, mol.getAtomCharge(i) - mol.getImplicitHydrogens(i)
                            + mol.getAttachedHydrogenCount(i));
                }
            }

            // Adding a positive charge to the elements that have a high valency
            if (mol.getOccupiedValence(i) > mol.getMaxValence(i) && mol.getAtomCharge(i) == 0) {
                mol.setAtomCharge(i, mol.getAtomCharge(i) + mol.getOccupiedValence(i) - mol.getMaxValence(i));
            }
        }
    }

    /**
     * Assigns positive charges to metal atoms if needed to balance out the overall
     * charge of the molecule. It is assumed that the metal atoms were not assigned
     * a charge up to this point.
     * @param mol
     *          The molecule that contains the metal atoms.
     */
    private static void assignMetalCharges(Molecule3D mol) {

        int netCharge = 0;
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            netCharge += mol.getAtomCharge(i);
        }
        if ( netCharge >= 0 ) {
            return;
        }

        List<Integer> metalAtoms = new ArrayList<>();
        for (int i = 0; i < mol.getAtoms(); i++) {
            if ( !mol.isMetalAtom(i) ) continue;
            if ( mol.getAtomCharge(i) != 0 ) continue;
            metalAtoms.add(i);
        }
        if ( metalAtoms.size() == 0 ) {
            return;
        }

        // Put all of the remaining negative charge on the single metal atom
        if ( metalAtoms.size() == 1 ) {
            mol.setAtomCharge(metalAtoms.get(0), -Math.max(netCharge, -8));
        } else {
        /*
        // Process the metal coordination sphere surroundings first
        // All negatively charged ligands are assumed to increase the
        // oxidation state by one despite their actual charge
        if ( netCharge < 0 ) {
            for (Integer metal : metalAtoms) {
                if (maximumIonChargeReached(mol, metal)) continue;
                if (!hasMetalLigandBonds(mol, metal))  continue;

                for (int k = mol.getAllConnAtoms(metal); k < mol.getAllConnAtomsPlusMetalBonds(metal); k++) {
                    int ligandCharge = mol.getAtomCharge(mol.getConnAtom(metal, k)) < 0 ? 1 : 0;
                    if (ligandCharge < 0) {
                        int oldMetalCharge = mol.getAtomCharge(metal);
                        int newMetalCharge =
                                StrictMath.min( oldMetalCharge + StrictMath.abs(ligandCharge),
                                        getMaximumIonCharge(mol, metal) );
                        mol.setAtomCharge(metal, newMetalCharge);
                        netCharge += newMetalCharge - oldMetalCharge;
                    }
                }
            }
        }
        */

            // TODO: select only atoms that do not have maximum charge

            // If all metals are of the same type, valence and charge
            // and the net charge is divisible by the number of metals
            // atoms, distribute the net charge equally
            boolean similarMetals = false;
            if ( StrictMath.abs(netCharge) % metalAtoms.size() == 0 ) {
                similarMetals = true;
                for (Integer metalAtom : metalAtoms) {
                    similarMetals = similarMetals &&
                            mol.getAtomicNo(metalAtoms.get(0)) == mol.getAtomicNo(metalAtom) &&
                            mol.getAtomCharge(metalAtoms.get(0)) == mol.getAtomCharge(metalAtom) &&
                            mol.getAllConnAtoms(metalAtoms.get(0)) == mol.getAllConnAtoms(metalAtom);
                }
            }

            if (similarMetals) {
                int increment = Math.min( StrictMath.abs(netCharge) / metalAtoms.size(), 8 );
                for (Integer metalAtom : metalAtoms) {
                    mol.setAtomCharge(metalAtom, increment);
                }
            } else {
                // Otherwise try to distribute it over all of the metal atoms
                boolean changeOccurred = true;
                while (netCharge < 0 && changeOccurred) {
                    changeOccurred = false;
                    for (Integer metal : metalAtoms) {
                        int newMetalCharge = mol.getAtomCharge(metal) + 1;

                        if (newMetalCharge <= getMaximumIonCharge(mol, metal)) {
                            netCharge += 1;
                            mol.setAtomCharge(metal, newMetalCharge);
                            changeOccurred = true;
                        }
                        if (netCharge == 0) break;
                    }
                }
            }
        }

        // Set abnormal valence to metals that exhibit uncommon
        // charges in order to generate a proper id codes
        for (Integer metal : metalAtoms) {
            // TODO: extend the logic to other metals as well
            // Currently, the logic is only applied to Aluminium (Al).
            // The code should work is other metals as well, but
            // more extensive testing is needed
            if (mol.getAtomicNo(metal) != 13) continue;
            if (Molecule.cAtomValence[mol.getAtomicNo(metal)] == null) continue;
            if (mol.getMaxValence(metal) > 0) {
                mol.setAtomAbnormalValence(metal, mol.getDefaultMaxValenceUncharged(metal) - mol.getMaxValence(metal));
            }
        }
    }

    /**
     * Creates a list of atom ring indexes that an atom is a
     * part of for each atom of the molecule.
     *
     * @param mol
     *          The molecule that contains the atoms and rings.
     * @return
     *          Array of atom ring index lists. The array indexes
     *          match the indexes of the molecule atoms. If an
     *          atom is not a part of any rings, null value is
     *          placed in the array instead of an empty list.
     */
    private static ArrayList<Integer>[] getAtomToRings(Molecule3D mol) {
        ArrayList<Integer>[] atomToRings = new ArrayList[mol.getAllAtoms()];

        RingCollection ringSet = mol.getRingSet();
        for ( int r = 0; r < ringSet.getSize(); r++ ) {
            int[] ringAtom = ringSet.getRingAtoms(r);
            for (int atom: ringAtom) {
                if (atomToRings[atom] == null) {
                    atomToRings[atom] = new ArrayList<>();
                }
                atomToRings[atom].add(r);
            }
        }
        return atomToRings;
    }

    /**
     * Determines which atom rings in the molecule are aromatic.
     * The aromaticity is determined by checking the ring geometry
     * against various constraints like atom chemical type, atom
     * planarity, atom bond lengths, atom bond angles.
     *
     * @param mol
     *          The molecule that contains the atom rings.
     * @param spOrder
     *          The sp hybridization type for each atom.
     * @return
     *          An array of boolean values marking which atom rings could
     *          potentially accommodate aromaticity. The indexes of the
     *          returned array correspond to the atom ring indexes in
     *          the molecule.
     */
    private static boolean[] findAromaticRings(Molecule3D mol, int[] spOrder) {
        RingCollection ringSet = mol.getRingSet();
        boolean[] looksLikeAromaticRing = new boolean[ringSet.getSize()];

        for (int ringNo = 0; ringNo < ringSet.getSize(); ringNo++) {
            looksLikeAromaticRing[ringNo] = looksLikeAromaticRing(mol, ringNo, spOrder);
        }

        return looksLikeAromaticRing;
    }

    /**
     * Determines is the given atom ring could potentially accommodate aromaticity.
     * False positive values may be returned.
     * @param mol
     *          Molecule that contains the atom ring.
     * @param ringNo
     *          Index of the ring in the molecule.
     * @param spOrder
     *          The sp hybridization type for each atom.
     * @return
     *          Boolean values denoting is the ring could potentially accommodate aromaticity.
     */
    private static boolean looksLikeAromaticRing(Molecule3D mol, int ringNo, int[] spOrder) {
        if (mol.getRingSet().getRingSize(ringNo) < 3) return false;
        if (mol.getRingSet().getRingSize(ringNo) > 7) return false;

        int[] ringAtoms = mol.getRingSet().getRingAtoms(ringNo);
        // Rings that are potentially interacting with a metal using
        // pi bonds should be allowed to be a little more deformed
        boolean isCoordinatedRing = isRingCoordinated(mol, ringNo);
        for (int i = 0; i < ringAtoms.length; i++) {
            int ringAtom = ringAtoms[i];
            if (!looksLikeAromaticAtom(mol, ringAtom, spOrder[ringAtom])) return false;

            int a0 = ringAtoms[(i - 1 + ringAtoms.length) % ringAtoms.length];
            int a1 = ringAtoms[(i) % ringAtoms.length];
            int a2 = ringAtoms[(i + 1) % ringAtoms.length];
            int a3 = ringAtoms[(i + 2) % ringAtoms.length];

            // Make sure that all boron, carbon, nitrogen and silicon
            // atoms in the ring are planar, but allow some deviation
            // due to neighbouring sulphur atoms
            if ((mol.getAtomicNo(a1) == 5 ||
                    mol.getAtomicNo(a1) == 6 ||
                    mol.getAtomicNo(a1) == 7 ||
                    mol.getAtomicNo(a1) == 14) &&
                    spOrder[ringAtom] != 2) {
                if (!((mol.getAtomicNo(a0) == 16 || mol.getAtomicNo(a2) == 16) &&
                        (getAllConnAtoms(mol, a1, true)) < 4)) {
                    return false;
                }
            }

            if (!(fitsAromaticLengthConstraints(mol, a1, a2) ||
                    hasMetalLigandBonds(mol, a1) ||
                    hasMetalLigandBonds(mol, a2) ||
                    isCoordinatedRing)) {
                if (DEBUG) {
                    System.err.println(">>> DEBUG: bond length between atoms '" + mol.getAtomName(a1) + "' and '" +
                            mol.getAtomName(a2) + "' caused ring " + ringNo + " to be marked as non aromatic.");
                }
                return false;
            }

            // check if all of the ring atoms lie in a plane
            double dihedral = getDihedral(mol, a0, a1, a2, a3);
            if (Math.abs(Math.toDegrees(dihedral)) > 20) {
                return false;
            }

            if (mol.getAtomicNo(a1) == 5 || mol.getAtomicNo(a1) == 6 || mol.getAtomicNo(a1) == 7) {
                if (isCoordinatedRing) continue;
                for (int j = 0; j < mol.getAllConnAtoms(a1); j++) {
                    if (isRingAtom(mol, ringNo, mol.getConnAtom(a1, j))) continue;
                    dihedral = getDihedral(mol, mol.getConnAtom(a1, j), a1, a2, a3);
                    if (StrictMath.min(Math.abs(Math.toDegrees(dihedral)), 180 - Math.abs(Math.toDegrees(dihedral))) > 25) {
                        if (DEBUG) {
                            System.err.println(">>> DEBUG: the torsion angle involving atoms '"
                                    + mol.getAtomName(mol.getConnAtom(a1, j)) + "', '"
                                    + mol.getAtomName(a1) + "' , "
                                    + mol.getAtomName(a2) + "' and '"
                                    + mol.getAtomName(a3) + "' caused ring " + ringNo +
                                    " to be marked as non aromatic.");
                        }
                        return false;
                    }
                }
            }
        }

        // excluding 3-membered rings that are part of a tetrahedron
        if (ringAtoms.length == 3) {
            int ringAtom = ringAtoms[0];
            int exocyclicAtom = -1;
            if (mol.getConnAtoms(ringAtom) == 4) {
                for (int j = 0; j < mol.getConnAtoms(ringAtom); j++) {
                    if (!isRingAtom(mol, ringNo, mol.getConnAtom(ringAtom, j))) {
                        exocyclicAtom = mol.getConnAtom(ringAtom, j);
                        break;
                    }
                }
            }
            return (mol.getBond(ringAtoms[1], exocyclicAtom) != -1) &&
                    (mol.getBond(ringAtoms[2], exocyclicAtom) != -1);
        }

        return true;
    }

    /**
     * Evaluates if the given atom could potentially participate in aromaticity.
     * False positive values may be returned.
     *
     * @param mol
     *          Molecule that contains the atom.
     * @param atom
     *          Index of the atom.
     * @param spOrder
     *          The sp hybridization type for each atom.
     * @return
     *          Boolean value denoting if the atom could participate in aromaticity.
     */
    private static boolean looksLikeAromaticAtom(Molecule3D mol, int atom, int spOrder) {
        // boron (B), carbon (C), nitrogen (N), oxygen (O), silicon (Si)
        if (mol.getAtomicNo(atom) < 15) {
            if (getAllConnAtoms(mol, atom, true) > 3) return false;
        }
        // germanium (Ge)
        if (mol.getAtomicNo(atom) == 32) {
            if (getAllConnAtoms(mol, atom, true) > 3) return false;
            if (spOrder != 2) return false;
        }

        if (getAllConnAtoms(mol, atom, true) > 4) return false;

        // Atoms that already have double bonds cannot participate in aromaticity
        // TODO: check if it is only true in some cases, sulphur (S) for example
        if (mol.getAtomPi(atom) != 0) return false;

        return true;
    }

    /**
     * Calculates the hybridization state for each atom in a molecule.
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @return
     *          Array of integers with values {1,2,3,4,5} that correspond to
     *          the the {sp1, sp2, sp3, sp3d1, sp3d2} hybridization states
     *          accordingly.
     */
    private static int[] calculateHybridizationStates(Molecule3D mol) {
        int[] spOrder = new int[mol.getAllAtoms()];
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            spOrder[i] = calculateHybridizationState(mol, i);
            if (DEBUG) {
                System.err.println("Atom '" + mol.getAtomName(i) + "' was assigned hybridization state " + spOrder[i]);
            }
        }
        return spOrder;
    }

    /**
     * Calculates the hybridization state for the given atom.
     *
     * @param mol
     *          Molecule that contains the atom.
     * @return
     *          Integers from the set of {1,2,3,4,5} that correspond to
     *          the the {sp1, sp2, sp3, sp3d1, sp3d2} hybridization states
     *          accordingly.
     */
    private static int calculateHybridizationState(Molecule3D mol, int i) {
        int spOrder = 0; // special value denoting an unrecognised hybridization state
        if (getAllConnAtoms(mol, i, true) <= 1) {
            // TODO: ask if this can always be assumed
            if (mol.getAtomicNo(i) == 8) {
                spOrder = 2;
            } else {
                spOrder = 1;
            }
        } else if (getAllConnAtoms(mol, i, true) == 2) {
            if (getAllConnAtoms(mol, i, false) <= 1) {
                spOrder = 2;
            } else {
                double angle = getAngle(mol, mol.getConnAtom(i, 0), i, mol.getConnAtom(i, 1));
                if (Math.abs(angle - Math.PI) < Math.PI / 6) {
                    spOrder = 1;
                } else {
                    spOrder = 2;
                }
            }
        } else if (getAllConnAtoms(mol, i, true) == 3) {
            if (getAllConnAtoms(mol, i, false) < 2) {
                spOrder = 0; // TODO: to be changed to 2 or 3
            } else if (getAllConnAtoms(mol, i, false) < 3) {
                spOrder = 2;
            } else {
                Coordinates origin = mol.getCoordinates(i);

                Coordinates normal = calculatePlaneNormal(origin,
                        mol.getCoordinates(mol.getConnAtom(i, 0)),
                        mol.getCoordinates(mol.getConnAtom(i, 1)));

                if (normal.distSq() > 0) {
                    if (calculatePointToPlaneDistance(normal, origin, mol.getCoordinates(mol.getConnAtom(i, 2))) < 0.4) {
                        spOrder = 2;
                    } else {
                        spOrder = 3;
                    }
                } else {
                    spOrder = 3;
                }
            }
        } else if (getAllConnAtoms(mol, i, true) == 4) {
            spOrder = 3;
        } else if (getAllConnAtoms(mol, i, true) == 5) {
            spOrder = 4;
        } else if (getAllConnAtoms(mol, i, true) == 6) {
            spOrder = 5;
        }

        return spOrder;
    }

    /**
     * Assigns bond orders and charges to hydrogen (H) atoms that
     * are bonded to two non-metal(lic) atoms.
     *
     * @param mol
     *      The molecule that contains the shared hydrogen atoms.
     */
    private static void resolveSharedHydrogenAtoms(Molecule3D mol) {
        // Shared hydrogen atoms are handled in one of the two ways
        // depending on the bond lengths to the neighbouring atoms:
        // - in case one of the bonds is significantly shorter,
        //   the longer bond gets removed;
        // - in case neither of the bonds is significantly shorter,
        //   both bonds are marked as metal-ligand (coordinate) bonds
        //   and the hydrogen atom is assigned a +1 charge.
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            // R-H-R bridging hydrogen atom
            if (mol.getAtomicNo(i) != 1) continue;
            if (mol.getConnAtoms(i) != 2) continue;

            // boron compounds must be handled separately
            if (mol.getAtomicNo(mol.getConnAtom(i, 0)) == 5) continue;
            if (mol.getAtomicNo(mol.getConnAtom(i, 1)) == 5) continue;

            int bond1 = mol.getConnBond(i, 0);
            int bond2 = mol.getConnBond(i, 1);
            double bondLength1 = getBondLength(mol, bond1);
            double bondLength2 = getBondLength(mol, bond2);
            if (StrictMath.abs(bondLength1 - bondLength2) < 0.1) {
                mol.setBondType(bond1, Molecule3D.cBondTypeMetalLigand);
                mol.setBondType(bond2, Molecule3D.cBondTypeMetalLigand);
                mol.setAtomCharge(i, +1);
            } else {
                // delete the longer single bond
                if (bondLength1 > bondLength2) {
                //    mol.setBondType(bond1, Molecule3D.cBondTypeMetalLigand);
                    mol.markBondForDeletion(bond1);
                } else {
                //    mol.setBondType(bond2, Molecule3D.cBondTypeMetalLigand);
                    mol.markBondForDeletion(bond2);
                }
            }
        }
        mol.deleteMarkedAtomsAndBonds();
        mol.ensureHelperArrays(Molecule.cHelperNeighbours);
    }

    /**
     * Recognises known ions in a molecule and assigns the appropriate
     * bond orders and atom charges if needed.
     *
     * @param mol
     *          Molecule that contains the ions.
     */
    private static void resolveKnownIons(Molecule3D mol) {
        for (int i = 0; i < mol.getAtoms(); i++) {
            if (mol.getAtomicNo(i) != 6) continue;
            if (calculateHybridizationState(mol, i) != 2) continue;
            if (mol.getAllConnAtoms(i) != 3) continue;
            if (mol.getAttachedHydrogenCount(i) != 0) continue;

            resolveTricyanomethanideIon(mol, i);
        }
    }

    /**
     * Evaluates if the given atom acts as the central atom of the tricyanomethanide
     * ion ([C-](C#N)(C#N)(C#N)) and assigns appropriate bond orders and atom charges
     * if needed.
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param cCentralAtom
     *          Index of the central carbon atom to be evaluated.
     * @return
     *          Boolean value denoting if any bond order or atom charge
     *          modifications were applied.
     */
    private static boolean resolveTricyanomethanideIon(Molecule3D mol, int cCentralAtom) {
        if (mol.getAtomicNo(cCentralAtom) != 6) return false;
        if (calculateHybridizationState(mol, cCentralAtom) != 2) return false;
        if (mol.getAllConnAtoms(cCentralAtom) != 3) return false;
        if (mol.getAttachedHydrogenCount(cCentralAtom) != 0) return false;

        int[] cAtoms = new int[3];
        int[] nAtoms = new int[3];
        for (int i = 0; i < mol.getAllConnAtoms(cCentralAtom); i++) {
            int cAtom = mol.getConnAtom(cCentralAtom, i);
            if (mol.getAllConnAtoms(cAtom) != 2) return false;
            for (int j = 0; j < mol.getAllConnAtoms(cAtom); j++) {
                int nAtom = mol.getConnAtom(cAtom, j);
                if (cCentralAtom == nAtom) continue;
                if (!looksLikeCyanoGroup(mol, cAtom, nAtom)) return false;
                cAtoms[i] = cAtom;
                nAtoms[i] = nAtom;
            }
        }

        // [C-](C#N)(C#N)(C#N)
        for (int j = 0; j < 3; j++) {
            mol.setBondOrder(mol.getBond(cAtoms[j], nAtoms[j]), 3);
            mol.setAtomCharge(cCentralAtom, -1);
        }
        mol.ensureHelperArrays(Molecule3D.cHelperNeighbours);
        return true;
    }

    /**
     * Evaluates if the given atom ring looks like a tetraatomic aromatic ring ion
     * and assigns appropriate bond orders and atom charges if needed.
     *
     * @param mol
     *          Molecule that contains the atoms ring.
     * @param ringNo
     *          Index of the ring in the molecule.
     * @return
     *          Boolean value denoting if any bond order or atom charge
     *          modifications were applied.
     */
    private static boolean resolveTetraatomicAromaticRingIons(Molecule3D mol, int ringNo) {
        if (!looksLikeTetraatomicAromaticRingIon(mol, ringNo)) return false;

        boolean isChalcogen = true;
        int atomicNo = mol.getAtomicNo(mol.getRingSet().getRingAtoms(ringNo)[0]);
        switch (atomicNo) {
            case 7: // nitrogen (N)
            case 15: // phosphorus (P);
            case 33: // arsenic (As);
                isChalcogen = false;
                break;
            case 8: // oxygen (O)
            case 16: // sulphur (S)
            case 34: // selenium (Se); tetraselenium ion ([Se]1[Se+]=[Se+][Se]1)
            case 52: // tellurium (Te)
                break;
            default:
                return false;
        }

        int[] ringBonds = mol.getRingSet().getRingBonds(ringNo);
        int shortestBondIndex = 0;
        double shortestBondLength = getBondLength(mol, ringBonds[shortestBondIndex]);
        for (int i = 0; i < ringBonds.length; i++) {
            if (getBondLength(mol, ringBonds[i]) < shortestBondLength) {
                shortestBondLength = getBondLength(mol, ringBonds[i]);
                shortestBondIndex = i;
            }
        }

        mol.setBondOrder(ringBonds[shortestBondIndex], 2);
        int dbAtom1 = mol.getBondAtom(0, ringBonds[shortestBondIndex]);
        int dbAtom2 = mol.getBondAtom(1, ringBonds[shortestBondIndex]);
        if (isChalcogen) {
            // chalcogen ring, i.e. [Se]1[Se][Se+]=[Se+]1
            mol.setAtomCharge(dbAtom1, +1);
            mol.setAtomCharge(dbAtom2, +1);
        } else {
            // pnictogen ring, i.e. [P-]1[P-][P]=[P]1
            for (int atom : mol.getRingSet().getRingAtoms(ringNo)) {
                if (atom != dbAtom1 && atom != dbAtom2) {
                    mol.setAtomCharge(atom, -1);
                }
            }
        }

        return true;
    }

    /**
     * Evaluates if the given atom ring looks like a tetraatomic aromatic ring ion,
     * i.e. the tetraselenium ion ([Se]1[Se+]=[Se+][Se]1).
     *
     * @param mol
     *          Molecule that contains the atoms ring.
     * @param ringNo
     *          Index of the ring in the molecule.
     * @return
     *          Boolean value denoting if the ring looks like a tetraatomic aromatic
     *          ring ion.
     */
    private static boolean looksLikeTetraatomicAromaticRingIon(Molecule3D mol, int ringNo) {
        int[] ringAtoms = mol.getRingSet().getRingAtoms(ringNo);
        if (ringAtoms.length != 4) return false;

        // rings must be homocyclic
        int atomicNo = mol.getAtomicNo(ringAtoms[0]);
        for (int atom : ringAtoms) {
            if (mol.getAtomicNo(atom) != atomicNo) return false;
            if (getAllConnAtoms(mol, atom, true) != 2) return false;
            if (mol.getAtomRingCount(atom, 7) != 1) return false;
        }

        return true;
    }

    /**
     * Recognises functional groups in a molecule and sets the bonds and charges accordingly.
     * @param mol
     *          The molecule that contains the functional groups.
     * @param spOrder
     *          Array that contains the orbital hybridization type for each atom.
     */
    private static void recogniseFunctionalGroups(Molecule3D mol, int[] spOrder) {

        // TODO: consider the atom type of terminal bonds when sorting
        // TODO: ask about the best way to keep track of new double bonds
        // calling molecule helpers to ensure pi count after each double
        // bond assignment is inefficient so a local array is used instead
        int[] localAtomPiCount = getAtomPiCount(mol);
        for(int i = 0; i < mol.getAtoms(); i++) {
            if ( mol.getAtomCharge(i) != 0 ) continue;
            if ( mol.getAtomRadical(i) != Molecule.cAtomRadicalStateNone ) continue;
            if ( spOrder[i] == 2 && localAtomPiCount[i] != 0 && mol.getAtomicNo(i) < 7 ) continue;
            // Terminal oxygen atoms capable of double bonds
            List<Integer> terminalOxygenAtoms = new LinkedList<>();
            // All terminal atoms capable of double bonds (including the oxygen atoms)
            List<Integer> terminalDoubleBondAtoms = new LinkedList<>();
            int oxygenCount = 0;
            int terminalAtomCount = 0;
            for (int j = 0; j < mol.getConnAtoms(i); j++ ) {
                int neighbour = mol.getConnAtom(i, j);
                boolean isTerminal = ( getAllConnAtoms(mol, neighbour, true) == 1 );
                boolean isOxygen   = ( mol.getAtomicNo(neighbour) == 8 );
                if (isOxygen) oxygenCount++;
                if ( isTerminal && isOxygen ) {
                    terminalOxygenAtoms.add(j);
                }

                if ( mol.getAtomicNo(neighbour) == 1 ) continue;
                if ( isTerminal ) {
                    terminalAtomCount++;
                    // for now only consider terminal N, O, P and S for double bonds
                    if ( localAtomPiCount[neighbour] == 0 &&
                         ( mol.getAtomicNo(neighbour) == 7 ||
                           mol.getAtomicNo(neighbour) == 8 ||
                           mol.getAtomicNo(neighbour) == 15 ||
                           mol.getAtomicNo(neighbour) == 16 ) ) {
                        terminalDoubleBondAtoms.add(j);
                    }
                }
            }

            if ( getAllConnAtoms(mol, i, true) == 1 ) {
                if (mol.getConnAtoms(i) == 0) continue;
                int a1 = mol.getConnAtom(i, 0);

                if (resolveDiatomicPnictogenIons(mol, i)) continue;
                if (resolveDiatomicChalcogenIons(mol, i)) continue;

                if ( mol.getAtomicNo(i) == 6 ) {
                    // R-N(+)=-=C(-)
                    if (mol.getAtomicNo(a1) == 7 && spOrder[a1] == 1) {
                        mol.setBondOrder(mol.getBond(i, a1), 3);
                        mol.setAtomCharge(i, -1);
                        mol.setAtomCharge(a1, +1);
                        localAtomPiCount[i]  += 2;
                        localAtomPiCount[a1] += 2;
                    }
                } else if ( mol.getAtomicNo(i) == 7 ) {
                    if (mol.getAtomicNo(a1) == 16 ) {
                        if ( getAllConnAtoms(mol, a1, true) == 1 ) {
                            if ( hasMetalLigandBonds(mol, i) ) {
                                mol.setBondOrder(mol.getBond(i, a1), 2);
                                localAtomPiCount[i]++;
                                localAtomPiCount[a1]++;
                                mol.setAtomCharge(i, -1);
                            } else {
                                mol.setBondOrder(mol.getBond(i, a1), 3);
                                localAtomPiCount[i] = localAtomPiCount[i] + 2;
                                localAtomPiCount[a1] = localAtomPiCount[a1] + 2;
                                mol.setAtomCharge(a1, 1);
                            }
                        } else if (spOrder[a1] == 1) {
                            mol.setBondOrder(mol.getBond(i, a1), 2);
                            localAtomPiCount[i]++;
                            localAtomPiCount[a1]++;
                        } else if (mol.getAtomicNo(a1) == 16 && spOrder[a1] == 3) {
                            // R-S#N
                            // FIXME: eventually replace the hardcoded S#N bond length value once enough statistics
                            // are gathered
                            if (getBondLength(mol, mol.getBond(i, a1)) < 1.47) {
                                mol.setBondOrder(mol.getBond(i, a1), 3);
                                localAtomPiCount[i] = localAtomPiCount[i] + 2;
                                localAtomPiCount[a1] = localAtomPiCount[a1] + 2;
                            }
                        }
                    }
                } else if ( mol.getAtomicNo(i) == 8 ) {
                    // TODO: change O-C distances with proper values
                    if ( mol.getAtomicNo(a1) == 6 ) {
                        // C-=-O
                        if ( getAllConnAtoms(mol, a1, true) == 1 ) {
                            mol.setBondOrder(mol.getBond(i, a1), 3);
                            mol.setAtomCharge(a1, -1);
                            mol.setAtomCharge(i, +1);
                            localAtomPiCount[i]  += 2;
                            localAtomPiCount[a1] += 2;
                        }
                    } else if ( mol.getAtomicNo(a1) == 7 ) {
                    // TODO: add the correct processing of the nitrosonium ion (N-=-O+)
                    // N=O ion
                        if ( getAllConnAtoms(mol, a1, true) == 1 ) {
                            mol.setBondOrder(mol.getBond(i, a1), 2);
                            mol.setAtomRadical(a1, Molecule3D.cAtomRadicalStateD);
                            localAtomPiCount[i]++;
                            localAtomPiCount[a1]++;
                        }
                    } else if ( mol.getAtomicNo(a1) == 8 ) {
                    // O-O ions
                        double sigma = 0.05;
                        if ( getAllConnAtoms(mol, a1, true) == 1 ) {
                            if (getBondLength(mol, i, a1) < 1.21 + sigma) {
                                mol.setBondOrder(mol.getBond(i, a1), 2);
                                localAtomPiCount[i]++;
                                localAtomPiCount[a1]++;
                            } else if (getBondLength(mol, i, a1) < 1.28 + sigma) {
                                mol.setBondOrder(mol.getConnBond(i, 0), 1);
                                mol.setAtomCharge(i, -1);
                                mol.setAtomRadical(a1, Molecule3D.cAtomRadicalStateD);
                            } else {
                                mol.setBondOrder(mol.getConnBond(i, 0), 1);
                                mol.setAtomCharge(i, -1);
                                mol.setAtomCharge(a1, -1);
                            }
                        }
                    }
                }
            } else if ( mol.getAllConnAtoms(i) == 2 ) {
                if ( localAtomPiCount[i] == 2 ) continue;
                // assigning neighbour atoms ordered by chemical number
                int a1 = mol.getAtomicNo(mol.getConnAtom(i, 0)) <=
                         mol.getAtomicNo(mol.getConnAtom(i, 1)) ?
                         mol.getConnAtom(i, 0) : mol.getConnAtom(i, 1);
                int a2 = a1 != mol.getConnAtom(i, 0) ?
                         mol.getConnAtom(i, 0) : mol.getConnAtom(i, 1);

                if ( mol.getAtomicNo(i) == 5 ) {
                    if ( spOrder[i] == 1 ) {
                        if ( terminalDoubleBondAtoms.size() == 2 ) {
                            mol.setBondOrder(mol.getBond(i, a1), 2);
                            mol.setBondOrder(mol.getBond(i, a2), 2);
                            localAtomPiCount[i]  += 2;
                            localAtomPiCount[a1] += 1;
                            localAtomPiCount[a2] += 1;
                            mol.setAtomCharge(i, -1);
                        }
                    }
                } else if ( mol.getAtomicNo(i) == 6 ) {
                    if ( spOrder[i] == 1 ) {
                        if (terminalDoubleBondAtoms.size() == 1) {
                            int terminalAtom = mol.getConnAtom(i, terminalDoubleBondAtoms.get(0));
                            int nonTerminalAtom = (terminalAtom != a1) ? a1 : a2;

                            // cyano group (R-C#N)
                            if (resolveCyanoGroup(mol, i, terminalAtom)) {
                                localAtomPiCount[i] += 2;
                                localAtomPiCount[terminalAtom] += 2;
                                continue;
                            }

                            // isocyanate group (R-N=C=O)
                            if (resolveIsocyanateGroup(mol, i, nonTerminalAtom, terminalAtom)) {
                                localAtomPiCount[i] += 2;
                                localAtomPiCount[terminalAtom] += 1;
                                localAtomPiCount[nonTerminalAtom] += 1;
                                continue;
                            }

                            // isothiocyanate group (R-N=C=S)
                            if (resolveIsothiocyanateGroup(mol, i, nonTerminalAtom, terminalAtom)) {
                                localAtomPiCount[i] += 2;
                                localAtomPiCount[terminalAtom] += 1;
                                localAtomPiCount[nonTerminalAtom] += 1;
                                continue;
                            }
                        } else if (terminalDoubleBondAtoms.size() == 2) {
                            // thiocyanate ion (S=C=N[-])
                            if (resolverThiocyanateIon(mol, i, a2, a1)) {
                                localAtomPiCount[i] += 2;
                                localAtomPiCount[a1] += 1;
                                localAtomPiCount[a2] += 1;
                                continue;
                            }

                            // atoms are of the same species
                            if (mol.getAtomicNo(a1) == mol.getAtomicNo(a2)) {
                                mol.setBondOrder(mol.getBond(i, a1), 2);
                                mol.setBondOrder(mol.getBond(i, a2), 2);
                                localAtomPiCount[i] += 2;
                                localAtomPiCount[a1] += 1;
                                localAtomPiCount[a2] += 1;
                            }
                        }
                    }
                } else if ( mol.getAtomicNo(i) == 7 ) {
                    if ( spOrder[i] == 1 ) {
                        if (terminalDoubleBondAtoms.size() > 0) {
                            // azide anion ([N-]=[N+]=[N-])
                            if (resolveAzideIon(mol, i)) {
                                localAtomPiCount[i] += 2;
                                localAtomPiCount[a1] += 1;
                                localAtomPiCount[a2] += 1;
                                continue;
                            }
                            // azide group (R[N]=[N+]=[N-])
                            if (resolveAzideGroup(mol, i)) {
                                localAtomPiCount[i] = mol.getAtomPi(i);
                                localAtomPiCount[a1] = mol.getAtomPi(a1);
                                localAtomPiCount[a2] = mol.getAtomPi(a2);
                                continue;
                            }

                            if ((mol.getAtomicNo(a1) == 6 && mol.getAtomicNo(a2) == 7) ||
                                    (mol.getAtomicNo(a1) == 7 && mol.getAtomicNo(a2) == 6)) {
                                int carbonAtom = (mol.getAtomicNo(a1) == 6) ? a1 : a2;
                                int nitrogenAtom = (mol.getAtomicNo(a1) == 7) ? a1 : a2;
                                if ( getAllConnAtoms(mol, nitrogenAtom, true) == 1 &&
                                     spOrder[carbonAtom] == 2) {
                                    double doubleNNBond = getIdealisedBondLength(2, 7, 7, 1, 2);
                                    double tripleNNBond = getIdealisedBondLength(3, 7, 7, 2, 2);
                                    double singleCNBond = getIdealisedBondLength(1, 6, 7, 1, 2);
                                    double doubleCNBond = getIdealisedBondLength(2, 6, 7, 1, 2);
                                    double NNBondLength = getBondLength(mol, i, nitrogenAtom);
                                    double CNBondLength = getBondLength(mol, i, carbonAtom);

                                    double C2N2N = StrictMath.max( doubleNNBond - NNBondLength, 0 ) +
                                                   StrictMath.max( doubleCNBond - CNBondLength, 0 );
                                    double C1N3N = StrictMath.max( tripleNNBond - NNBondLength, 0 ) +
                                                   StrictMath.max( singleCNBond - CNBondLength, 0 );
                                    // RC-[N+]#N
                                    if ( C1N3N < C2N2N )  {
                                        mol.setBondOrder(mol.getBond(i, nitrogenAtom), 3);
                                        mol.setAtomCharge(i, +1);
                                        localAtomPiCount[i]  += 2;
                                        localAtomPiCount[nitrogenAtom] += 1;
                                    } else {
                                        // RC=[N+]=[N-]
                                        mol.setBondOrder(mol.getBond(i, nitrogenAtom), 2);
                                        mol.setBondOrder(mol.getBond(i, carbonAtom), 2);
                                        mol.setAtomCharge(i, +1);
                                        mol.setAtomCharge(nitrogenAtom, -1);
                                        localAtomPiCount[i]  += 2;
                                        localAtomPiCount[nitrogenAtom] += 1;
                                        localAtomPiCount[carbonAtom] += 1;
                                    }
                                }
                            // O=N=O, S=N=S
                            } else if ( terminalOxygenAtoms.size() == 2 ||
                                    mol.getAtomicNo(a1) == 16 && mol.getAtomicNo(a2) == 16 ) {
                                mol.setAtomCharge(i, +1);
                                for (Integer terminalDoubleBondAtom : terminalDoubleBondAtoms) {
                                    mol.setBondOrder(mol.getConnBond(i, terminalDoubleBondAtom), 2);
                                    localAtomPiCount[i] += 1;
                                    localAtomPiCount[mol.getConnAtom(i, terminalDoubleBondAtom)] += 1;
                                }
                            }
                        } else {
                            // R-P=[N+]=P-R
                            if ( mol.getAtomicNo(a1) == 15 && spOrder[a1] == 3 &&
                                 mol.getAtomicNo(a2) == 15 && spOrder[a2] == 3 ) {
                                // cyclophosphaneses
                                if ( mol.isRingAtom(i) ) {
                                    double NPBondLength1 = getBondLength(mol, i, mol.getConnAtom(i, 0));
                                    double NPBondLength2 = getBondLength(mol, i, mol.getConnAtom(i, 1));
                                    if ( NPBondLength1 < NPBondLength2 && localAtomPiCount[mol.getConnAtom(i, 0)] == 0) {
                                        mol.setBondOrder(mol.getConnBond(i, 0), 2);
                                        localAtomPiCount[i] += 1;
                                        localAtomPiCount[mol.getConnAtom(i, 0)] += 1;
                                    } else if ( localAtomPiCount[mol.getConnAtom(i, 1)] == 0 ) {
                                        mol.setBondOrder(mol.getConnBond(i, 1), 2);
                                        localAtomPiCount[i] += 1;
                                        localAtomPiCount[mol.getConnAtom(i, 1)] += 1;
                                    }
                                // R-P=N(+)=P-R (bridging)
                                } else {
                                    mol.setBondOrder(mol.getConnBond(i, 0), 2);
                                    mol.setBondOrder(mol.getConnBond(i, 1), 2);
                                    mol.setAtomCharge(i, +1);
                                    localAtomPiCount[i] += 2;
                                    localAtomPiCount[mol.getConnAtom(i, 0)] += 1;
                                    localAtomPiCount[mol.getConnAtom(i, 1)] += 1;
                                }
                            }
                        }
                    } else if ( spOrder[i] == 2 ) {
                        if (terminalOxygenAtoms.size() == 1) {
                            if (resolveNitrosoGroup(mol, i)) {
                                localAtomPiCount[i] += 1;
                                localAtomPiCount[mol.getConnAtom(i, terminalOxygenAtoms.get(0))] += 1;
                            }
                        } else {
                            // R-P=N(+)-P-R
                            if (mol.getAtomicNo(a1) != 15) continue;
                            if (spOrder[a1] != 3) continue;
                            if (mol.getAtomicNo(a2) != 15) continue;
                            if (spOrder[a2] != 3) continue;

                            // skip if both phosphorus atoms already have a double bond
                            if (localAtomPiCount[a1] != 0 && localAtomPiCount[a2] != 0) continue;

                            // select the shortest bond out of the available ones
                            int doubleBondAtom = a1;
                            if (localAtomPiCount[a1] != 0) {
                                doubleBondAtom = a2;
                            } else if (localAtomPiCount[a2] == 0) {
                                if (getBondLength(mol, mol.getBond(i, a1)) > getBondLength(mol, mol.getBond(i, a2))) {
                                    doubleBondAtom = a2;
                                }
                            }
                            mol.setBondOrder(mol.getBond(i, doubleBondAtom), 2);
                            localAtomPiCount[i] += 1;
                            localAtomPiCount[doubleBondAtom] += 1;
                        }
                    }
                // S
                } else if ( mol.getAtomicNo(i) == 16 ) {
                    if (spOrder[i] == 2) {
                        // SO2
                        if (terminalOxygenAtoms.size() > 0) {
                            int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalOxygenAtoms);
                            setTerminalDoubleBonds(mol, i, sortedAtoms, 1);
                            mol.setAtomCharge(i, +1);
                        }
                    }
                // F, Cl, Br, I
                } else if ( mol.getAtomicNo(i) ==  7 ||
                            mol.getAtomicNo(i) == 17 ||
                            mol.getAtomicNo(i) == 35 ||
                            mol.getAtomicNo(i) == 53 ) {
                    if (terminalDoubleBondAtoms.size() > 0) {
                        // F is not capable of a double bond
                        if (mol.getAtomicNo(i) == 7) {
                            mol.setAtomCharge(i, +1);
                        } else {
                            int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalDoubleBondAtoms);
                            setTerminalDoubleBonds(mol, i, sortedAtoms, 1);
                        }
                    } else {
                        // R-X-R (linear)
                        if (spOrder[i] == 1) {
                            mol.setAtomCharge(i, -1);
                        // R-X-R (bent)
                        } else if (spOrder[i] == 2) {
                            mol.setAtomCharge(i, +1);
                        }
                    }
                }
            } else if ( mol.getAllConnAtoms(i) == 3 ) {

                int a, b;
                int a1 = mol.getConnAtom(i, 0);
                int a2 = mol.getConnAtom(i, 1);
                int a3 = mol.getConnAtom(i, 2);

                if (mol.getAtomicNo(i) == 6 && spOrder[i] == 2) {
                    if (mol.getAtomRingSize(i) > 0) {
                        if (terminalDoubleBondAtoms.isEmpty()) continue;
                        // thioketone in TTF: C(R)(R)(=S)
                        if (mol.getAtomRingSize(i) != 5) continue;
                        if (mol.getAtomicNo(a1) != 16) continue;
                        if (mol.getAtomicNo(a2) != 16) continue;
                        if (mol.getAtomicNo(a3) != 16) continue;

                        mol.setBondOrder(mol.getConnBond(i, terminalDoubleBondAtoms.get(0)), 2);
                    } else {
                        //C(N)(N)(=O)

                        //C(R)(O)(=O)
                        if (terminalOxygenAtoms.size() >= 2) {
                            int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalOxygenAtoms);
                            setTerminalDoubleBonds(mol, i, sortedAtoms, 1);
                            continue;
                        }

                        //C(R)(OR)(=O)
                        a = connectedAtom(mol, i, 8, 2, 0, 0);
                        b = connectedBond(mol, i, 8, 1);
                        if (a >= 0 && b >= 0) {
                            mol.setBondOrder(b, 2);
                            continue;
                        }

                        //C(R)(SR)(=O)
                        a = connectedAtom(mol, i, 16, 2, 0, 0);
                        b = connectedBond(mol, i, 8, 1);
                        if (a >= 0 && b >= 0) {
                            mol.setBondOrder(b, 2);
                            continue;
                        }

                        //C(R)(NR)(=O)
                        a = connectedAtom(mol, i, 7, 3, 0, 0);
                        b = connectedBond(mol, i, 8, 1);
                        if (a >= 0 && b >= 0) {
                            mol.setBondOrder(b, 2);
                            continue;
                        }

                        //C(R)(SR)(=S)
                        a = connectedAtom(mol, i, 16, 2, 0, 0);
                        b = connectedBond(mol, i, 16, 1);
                        if (a >= 0 && b >= 0) {
                            mol.setBondOrder(b, 2);
                            continue;
                        }

                        //C(R)(NR)(=S)
                        a = connectedAtom(mol, i, 7, 3, 0, 0);
                        b = connectedBond(mol, i, 16, 1);
                        if (a >= 0 && b >= 0) {
                            mol.setBondOrder(b, 2);
                            continue;
                        }

                        //C(CR)(N)(=N)
                        if (terminalDoubleBondAtoms.size() == 2 &&
                                mol.getAtomicNo(terminalDoubleBondAtoms.get(0)) == 7 &&
                                mol.getAtomicNo(terminalDoubleBondAtoms.get(1)) == 7 &&
                                (mol.getAtomicNo(a1) == 6 ||
                                        mol.getAtomicNo(a2) == 6 ||
                                        mol.getAtomicNo(a3) == 6)) {
                            int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalDoubleBondAtoms);
                            setTerminalDoubleBonds(mol, i, sortedAtoms, 1);
                        }
                    }
                } else if( mol.getAtomicNo(i) == 7 ) {
                    // TODO: evaluate the probability of a certain SP order from the geometry
                    // slightly disordered nitrate ion
                    if ( spOrder[i] == 2 ) {
                        //N(R)(R)C=O -> Amide
                        a = connectedAtom(mol, i, 6, 2, 8, 1);
                        b = connectedBond(mol, a, 8, 1);
                        if (a >= 0 && b >= 0) {
                            mol.setBondOrder(b, 2);
                            localAtomPiCount[i] += 1;
                            localAtomPiCount[b] += 1;
                        } else
                        // N(O)(=O) -> Nitro
                        // nitrate ion
                        // NF2O(1+)
                        if ( terminalAtomCount > 1 && terminalOxygenAtoms.size() > 0) {
                            int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalOxygenAtoms);
                            setTerminalDoubleBonds(mol, i, sortedAtoms, 1);
                            mol.setAtomCharge(i, +1);
                        }
                    }
                } else if (mol.getAtomicNo(i) == 16) {
                    if (spOrder[i] == 2) {
                        // sulphur trioxide
                        if (terminalOxygenAtoms.size() > 0) {
                            int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalOxygenAtoms);
                            setTerminalDoubleBonds(mol, i, sortedAtoms, 3);
                        }
                    } else if (spOrder[i] == 3) {
                        // sulphite ion
                        if (terminalOxygenAtoms.size() > 0) {
                            int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalOxygenAtoms);
                            setTerminalDoubleBonds(mol, i, sortedAtoms, 1);
                        }
                    }
                // Cl, Br, I
                } else if ( mol.getAtomicNo(i) == 17 ||
                            mol.getAtomicNo(i) == 35 ||
                            mol.getAtomicNo(i) == 53 ) {
                    if ( spOrder[i] == 3 && terminalDoubleBondAtoms.size() > 0 ) {
                        int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalDoubleBondAtoms);
                        setTerminalDoubleBonds(mol, i, sortedAtoms, 2);
                    }
                } else if ( mol.getAtomicNo(i) == 34 ) {
                    if (oxygenCount == 3) {
                        if (spOrder[i] == 2) {
                            // selenium trioxide
                            int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalOxygenAtoms);
                            setTerminalDoubleBonds(mol, i, sortedAtoms, 2);
                        } else if (spOrder[i] == 3) {
                            // selenite ion
                            int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalOxygenAtoms);
                            setTerminalDoubleBonds(mol, i, sortedAtoms, 1);
                        }
                    } else {
                        mol.setAtomCharge(i, +1);
                    }
                // Te
                } else if( mol.getAtomicNo(i) == 52 ) {
                    mol.setAtomCharge(i, +1);
                }
            } else if( mol.getAllConnAtoms(i) == 4 ) {
                // phosphorus (P) and arsenic (As)
                if( mol.getAtomicNo(i) == 15 || mol.getAtomicNo(i) == 33 ) {
                    if (terminalDoubleBondAtoms.isEmpty()) continue;
                    // skip is the atom already has a double bond
                    if (localAtomPiCount[i] != 0) continue;
                    int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalDoubleBondAtoms);
                    setTerminalDoubleBonds(mol, i, sortedAtoms, 1);
                    localAtomPiCount[i] += 1;
                    localAtomPiCount[mol.getConnAtom(i, sortedAtoms[0])] += 1;
                } else if( mol.getAtomicNo(i)==16 ) {
                    int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalDoubleBondAtoms);
                    setTerminalDoubleBonds(mol, i, sortedAtoms, 2);
                // Cl, Br, I
                } else if( mol.getAtomicNo(i) == 17 ||
                           mol.getAtomicNo(i) == 35 ||
                           mol.getAtomicNo(i) == 53 ) {
                    if ( terminalDoubleBondAtoms.size() > 0 ) {
                        int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalDoubleBondAtoms);
                        setTerminalDoubleBonds(mol, i, sortedAtoms, 3);
                    } else {
                        // Set a negative charge to hypervalent atom
                        mol.setAtomCharge(i, -1);
                    }
                } else if( mol.getAtomicNo(i) == 34 ) {
                    if ( oxygenCount == 4 && terminalOxygenAtoms.size() > 1 ) {
                        int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalOxygenAtoms);
                        setTerminalDoubleBonds(mol, i, sortedAtoms, 2);
                    }
                }
            } else if( mol.getAllConnAtoms(i) == 5 ) {
                if ( mol.getAtomicNo(i) == 14 ) {        // Silicon (Si)
                    mol.setAtomCharge(i, -1);
                } else if ( mol.getAtomicNo(i) == 53 ) { // Iodine (I)
                    if ( terminalOxygenAtoms.size() > 0 ) {
                        int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalOxygenAtoms);
                        setTerminalDoubleBonds(mol, i, sortedAtoms, 2);
                    }
                }
            } else if ( mol.getAllConnAtoms(i) == 6) {
                if (mol.getAtomicNo(i) == 14) {        // Silicon (Si)
                    mol.setAtomCharge(i, -2);
                } else if (mol.getAtomicNo(i) == 32) { // Germanium (Ge)
                    mol.setAtomCharge(i, -2);
                } else if (mol.getAtomicNo(i) == 33) { // Arsenic (As)
                    mol.setAtomCharge(i, -1);
                } else if (mol.getAtomicNo(i) == 34) { // Selenium (Se)
                    if ( terminalDoubleBondAtoms.size() > 0 ) {
                        int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalDoubleBondAtoms);
                        setTerminalDoubleBonds(mol, i, sortedAtoms, 1);
                    } else {
                        mol.setAtomCharge(i, -2);
                    }
                } else if( mol.getAtomicNo(i) == 52 ) { // Tellurium (Te)
                    // telluric acid has a net charge of 0
                    // even though it has 6 neighbours
                    if ( terminalDoubleBondAtoms.size() == 0 && oxygenCount == 0 ) {
                        mol.setAtomCharge(i, -2);
                    }
                } else if( mol.getAtomicNo(i) == 53 ) { // Iodine (I)
                    if ( terminalDoubleBondAtoms.size() > 0 ) {
                        int[] sortedAtoms = sortAtomsByBondLength(mol, i, terminalDoubleBondAtoms);
                        setTerminalDoubleBonds(mol, i, sortedAtoms, 1);
                    } else {
                        mol.setAtomCharge(i, +1);
                    }
                }
            }
        }
        mol.ensureHelperArrays(Molecule3D.cHelperRings);
    }

    /**
     * Evaluates if the given nitrogen atom serves as a part of a nitroso (R-N=O)
     * group and assigns appropriate bond orders and atom charges if needed.
     *
     * Current implementation uses the following assumptions that are not explicitly checked:
     *  - the provided nitrogen atom has only two neighbour atoms;
     *  - one of the neighbour atoms is an oxygen atom capable of a double bond;
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param nAtom
     *          Index of the nitrogen atom to be evaluated.
     * @return
     *          Boolean value denoting if any bond order or atom charge
     *          modifications were applied.
     */
    private static boolean resolveNitrosoGroup(Molecule3D mol, int nAtom) {
        int oAtom;
        int rAtom;
        if (mol.getAtomicNo(mol.getConnAtom(nAtom, 0)) == 8) {
            oAtom = mol.getConnAtom(nAtom, 0);
            rAtom = mol.getConnAtom(nAtom, 1);
        } else {
            oAtom = mol.getConnAtom(nAtom, 1);
            rAtom = mol.getConnAtom(nAtom, 0);
        }

        // c-nitroso group (RC-N=O)
        if (mol.getAtomicNo(rAtom) == 6) {
            mol.setBondOrder(mol.getBond(nAtom, oAtom), 2);
            return true;
        }

        return false;
    }

    /**
     * Evaluates if the given atoms could potentially comprise a cyano group (R-C#N).
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param cAtom
     *          Index of the carbon atom to be evaluated.
     * @param nAtom
     *          Index of the nitrogen atom to be evaluated.
     * @return
     *          Boolean value denoting if the provided atoms could potentially
     *          comprise the group.
     */
    private static boolean looksLikeCyanoGroup(Molecule3D mol, int cAtom, int nAtom) {
        if (mol.getAtomicNo(cAtom) != 6) return false;
        if (calculateHybridizationState(mol, cAtom) != 1) return false;

        if (mol.getAtomicNo(nAtom) != 7) return false;
        if (mol.getAllConnAtoms(nAtom) != 1) return false;
        if (mol.getAttachedHydrogenCount(nAtom) != 0) return false;

        return true;
    }

    /**
     * Evaluates if the given atoms comprise a cyano group (R-C#N)
     * and assigns appropriate bond orders and atom charges if needed.
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param cAtom
     *          Index of the carbon atom to be evaluated.
     * @param nAtom
     *          Index of the nitrogen atom to be evaluated.
     * @return
     *          Boolean value denoting if any bond order or atom charge
     *          modifications were applied.
     */
    private static boolean resolveCyanoGroup(Molecule3D mol, int cAtom, int nAtom) {
        if (!looksLikeCyanoGroup(mol, cAtom, nAtom)) return false;
        // cyano group (R-C#N)
        mol.setBondOrder(mol.getBond(cAtom, nAtom), 3);

        return true;
    }

    /**
     * Evaluates if the given atoms comprise an isocyanate group (R-N=C=O)
     * and assigns appropriate bond orders and atom charges if needed.
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param cAtom
     *          Index of the carbon atom to be evaluated.
     * @param nAtom
     *          Index of the nitrogen atom to be evaluated.
     * @param oAtom
     *          Index of the oxygen atom to be evaluated.
     * @return
     *          Boolean value denoting if any bond order or atom charge
     *          modifications were applied.
     */
    private static boolean resolveIsocyanateGroup(Molecule3D mol, int cAtom, int nAtom, int oAtom) {
        if (mol.getAtomicNo(nAtom) != 7) return false;

        if (mol.getAtomicNo(cAtom) != 6) return false;
        if (calculateHybridizationState(mol, cAtom) != 1) return false;

        if (mol.getAtomicNo(oAtom) != 8) return false;
        if (mol.getAllConnAtoms(oAtom) != 1) return false;
        if (mol.getAttachedHydrogenCount(oAtom) != 0) return false;

        // isocyanate group (N=C=O)
        mol.setBondOrder(mol.getBond(nAtom, cAtom), 2);
        mol.setBondOrder(mol.getBond(cAtom, oAtom), 2);

        return true;
    }

    /**
     * Evaluates if the given atoms comprise an isothiocyanate group (R-N=C=S)
     * and assigns appropriate bond orders and atom charges if needed.
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param cAtom
     *          Index of the carbon atom to be evaluated.
     * @param nAtom
     *          Index of the nitrogen atom to be evaluated.
     * @param sAtom
     *          Index of the sulphur atom to be evaluated.
     * @return
     *          Boolean value denoting if any bond order or atom charge
     *          modifications were applied.
     */
    private static boolean resolveIsothiocyanateGroup(Molecule3D mol, int cAtom, int nAtom, int sAtom) {
        if (mol.getAtomicNo(nAtom) != 7) return false;

        if (mol.getAtomicNo(cAtom) != 6) return false;
        if (calculateHybridizationState(mol, cAtom) != 1) return false;

        if (mol.getAtomicNo(sAtom) != 16) return false;
        if (mol.getAllConnAtoms(sAtom) != 1) return false;
        if (mol.getAttachedHydrogenCount(sAtom) != 0) return false;

        // isothiocyanate group (N=C=S)
        mol.setBondOrder(mol.getBond(nAtom, cAtom), 2);
        mol.setBondOrder(mol.getBond(cAtom, sAtom), 2);

        return true;
    }

    /**
     * Evaluates if the given atoms comprise a thiocyanate ion (S=C=[N-])
     * and assigns appropriate bond orders and atom charges if needed.
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param cAtom
     *          Index of the carbon atom to be evaluated.
     * @param sAtom
     *          Index of the sulphur atom to be evaluated.
     * @param nAtom
     *          Index of the nitrogen atom to be evaluated.
     * @return
     *          Boolean value denoting if any bond order or atom charge
     *          modifications were applied.
     */
    private static boolean resolverThiocyanateIon(Molecule3D mol, int cAtom, int sAtom, int nAtom) {
        if (mol.getAtomicNo(cAtom) != 6) return false;
        if (calculateHybridizationState(mol, cAtom) != 1) return false;
        if (mol.getAtomicNo(nAtom) != 7) return false;
        if (mol.getAtomicNo(sAtom) != 16) return false;

        // TODO: consider the resonance structure of SCN ([S-]-C=N <=> S=C=N[-])?
        // thiocyanate ion (S=C=N[-])
        mol.setBondOrder(mol.getBond(sAtom, cAtom), 2);
        mol.setBondOrder(mol.getBond(cAtom, nAtom), 2);
        mol.setAtomCharge(nAtom, -1);

        return true;
    }

    /**
     * Evaluate if the given atom could potentially act as the central
     * atom of a azide group (R[N]=[N+]=[N-]).
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param nCentralAtom
     *          Index of the nitrogen atom to be evaluated.
     * @return
     *          Boolean value denoting if the given atom could potentially
     *          act as the central atom of an azide group.
     */
    private static boolean looksLikeUnresolvedAzideGroup(Molecule3D mol, int nCentralAtom) {
        if (mol.getAtomicNo(nCentralAtom) != 7) return false;
        if (calculateHybridizationState(mol, nCentralAtom) != 1) return false;
        if (mol.getConnAtoms(nCentralAtom) != 2) return false;
        if (getAllConnAtoms(mol, nCentralAtom, true) != 2) return false;

        int nAtom1 = mol.getConnAtom(nCentralAtom, 0);
        int nAtom2 = mol.getConnAtom(nCentralAtom, 1);

        if (mol.getAtomicNo(nAtom1) != 7) return false;
        if (getAllConnAtoms(mol, nAtom1, true) > 2) return false;

        if (mol.getAtomicNo(nAtom2) != 7) return false;
        if (getAllConnAtoms(mol, nAtom2, true) > 2) return false;

        // check if there is at least one terminal atom
        if (getAllConnAtoms(mol, nAtom1, true) != 1 &&
                getAllConnAtoms(mol, nAtom2, true) != 1) {
            return false;
        }

        return true;
    }

    /**
     * Evaluate if the given atom could potentially act as the central
     * atom of an azide group (R[N]=[N+]=[N-]) and assigns appropriate bond
     * orders and atom charges for the entire group if needed.
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param nCentralAtom
     *          Index of the nitrogen atom to be evaluated.
     * @return
     *          Boolean value denoting if any bond order or atom charge
     *          modifications were applied.
     */
    private static boolean resolveAzideGroup(Molecule3D mol, int nCentralAtom) {
        if (!looksLikeUnresolvedAzideGroup(mol, nCentralAtom)) return false;

        int nAtom1 = mol.getConnAtom(nCentralAtom, 0);
        int nAtom2 = mol.getConnAtom(nCentralAtom, 1);
        int shorterBondAtom;
        int longerBondAtom;
        if (getBondLength(mol, nCentralAtom, nAtom1) < getBondLength(mol, nCentralAtom, nAtom2)) {
            shorterBondAtom = nAtom1;
            longerBondAtom = nAtom2;
        } else {
            shorterBondAtom = nAtom2;
            longerBondAtom = nAtom1;
        }

        double shorterBondLength = getBondLength(mol, nCentralAtom, shorterBondAtom);
        double longerBondLength = getBondLength(mol, nCentralAtom, longerBondAtom);
        boolean accommodatesTripleBond = calculateHybridizationState(mol, shorterBondAtom) == 1 &&
                StrictMath.abs(shorterBondAtom - longerBondAtom) > 0.05;
        if (accommodatesTripleBond) {
            double NNSingleBondLength = getIdealisedBondLength(2, 7, 7, 1, 2);
            double NNDoubleBondLength = getIdealisedBondLength(2, 7, 7, 1, 2);
            double NNTripleBondLength = getIdealisedBondLength(3, 7, 7, 2, 2);

            double N2N2N = StrictMath.abs(NNDoubleBondLength - shorterBondLength) +
                    StrictMath.abs(NNDoubleBondLength - longerBondLength);
            double N1N3N = StrictMath.abs(NNTripleBondLength - shorterBondLength) +
                    StrictMath.abs(NNSingleBondLength - longerBondLength);

            accommodatesTripleBond = N1N3N < N2N2N;
        }

        // choose between the resonance structures
        // RN=[N+]=[N-] and R[N-]-[N+]#N
        if (accommodatesTripleBond) {
            mol.setBondOrder(mol.getBond(nCentralAtom, shorterBondAtom), 3);
            mol.setBondOrder(mol.getBond(nCentralAtom, longerBondAtom), 1);
            mol.setAtomCharge(nCentralAtom, +1);
            // in case of NH2-[N+]#N
            // TODO: check if this conditional is needed
            if (mol.getFreeValence(longerBondAtom) > 0) {
                mol.setAtomCharge(longerBondAtom, -1);
            }
        } else {
            mol.setBondOrder(mol.getBond(nCentralAtom, shorterBondAtom), 2);
            mol.setBondOrder(mol.getBond(nCentralAtom, longerBondAtom), 2);
            mol.setAtomCharge(nCentralAtom, +1);
            if (mol.getFreeValence(shorterBondAtom) > 0 || calculateHybridizationState(mol, shorterBondAtom) == 1) {
                mol.setAtomCharge(shorterBondAtom, -1);
            }
            if (mol.getFreeValence(longerBondAtom) > 0 || calculateHybridizationState(mol, longerBondAtom) == 1) {
                mol.setAtomCharge(longerBondAtom, -1);
            }
        }
        mol.ensureHelperArrays(Molecule.cHelperRings);

        return true;
    }

    /**
     * Evaluate if the given atom could potentially act as the central
     * atom of a azide ion ([N-]=[N+]=[N-]).
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param nCentralAtom
     *          Index of the nitrogen atom to be evaluated.
     * @return
     *          Boolean value denoting if the given atom could potentially
     *          act as the central atom of an azide ion.
     */
    private static boolean looksLikeUnresolvedAzideIon(Molecule3D mol, int nCentralAtom) {
        if (!looksLikeUnresolvedAzideGroup(mol, nCentralAtom)) return false;

        int nAtom1 = mol.getConnAtom(nCentralAtom, 0);
        int nAtom2 = mol.getConnAtom(nCentralAtom, 1);

        if (getAllConnAtoms(mol, nAtom1, true) != 1) return false;
        if (getAllConnAtoms(mol, nAtom2, true) != 1) return false;

        return true;
    }

    /**
     * Evaluate if the given atom could potentially act as the central
     * atom of an azide ion ([N-]=[N+]=[N-]) and assigns appropriate bond
     * orders and atom charges for the entire ion if needed.
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param nCentralAtom
     *          Index of the nitrogen atom to be evaluated.
     * @return
     *          Boolean value denoting if any bond order or atom charge
     *          modifications were applied.
     */
    private static boolean resolveAzideIon(Molecule3D mol, int nCentralAtom) {
        if (!looksLikeUnresolvedAzideIon(mol, nCentralAtom)) return false;

        int a1 = mol.getConnAtom(nCentralAtom, 0);
        int a2 = mol.getConnAtom(nCentralAtom, 1);
        mol.setBondOrder(mol.getBond(nCentralAtom, a1), 2);
        mol.setBondOrder(mol.getBond(nCentralAtom, a2), 2);
        mol.setAtomCharge(nCentralAtom, +1);
        mol.setAtomCharge(a1, -1);
        mol.setAtomCharge(a2, -1);
        mol.ensureHelperArrays(Molecule.cHelperRings);

        return true;
    }

    /**
     * Evaluates if the given atom is part of a diatomic pnictogen ion of
     * the forms [E-2][E-2], [E-1]=[E-1] where E = P, As (i.e. [P-2][P-2],
     * [P-1]=[P-1]) and assigns appropriate bond orders and atom charges
     * for the entire ion if needed.
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param atom
     *          Index of the atom to be evaluated.
     * @return
     *          Boolean value denoting if any atom charge modifications
     *          were applied.
     */
    private static boolean resolveDiatomicPnictogenIons(Molecule3D mol, int atom) {
        if (!mol.isNitrogenFamily(atom)) return false;
        // skip nitrogen (N) for now
        if (mol.getAtomicNo(atom) == 7) return false;
        if (getAllConnAtoms(mol, atom, true) != 1) return false;
        if (mol.getConnAtoms(atom) != 1) return false;

        int nAtom = mol.getConnAtom(atom, 0);
        if (!mol.isNitrogenFamily(nAtom)) return false;
        // skip nitrogen (N) for now
        if (mol.getAtomicNo(nAtom) == 7) return false;
        if (mol.getConnAtoms(nAtom) != 1) return false;

        // double bonded form ([E-]=[E-])
        if (mol.getAtomicNo(atom) == 15 && mol.getAtomicNo(nAtom) == 15) {
            // TODO: eventually replace the hardcoded value with one based on statistics
            // Explanation: the provided cutoff value was chosen based
            // on bond distances from 6 crystal structures that were
            // show to contain the P=P ligand. The hardcoded value may
            // eventually be removed if sufficiently granular and accurate
            // statistics become available.
            double phosphorusDoubleBondCuttoff = 2.05;
            if (phosphorusDoubleBondCuttoff > getBondLength(mol, atom, nAtom)) {
                mol.setBondOrder(mol.getBond(atom, nAtom), 2);
                mol.setAtomCharge(atom, -1);
                mol.setAtomCharge(nAtom, -1);
                return true;
            }
        }
        // TODO: investigate if a cutoff value for As is also needed

        // single bonded form ([E-2][E-2])
        mol.setAtomCharge(atom, -2);
        mol.setAtomCharge(nAtom, -2);

        return true;
    }

    /**
     * Evaluates if the given atom is part of a diatomic chalcogen ion of
     * the form [E-][E-] where E = S, Se, Te (i.e. [S-][S-]) and assigns
     * appropriate atom charges for the entire ion if needed.
     *
     * @param mol
     *          Molecule that contains the atoms.
     * @param atom
     *          Index of the atom to be evaluated.
     * @return
     *          Boolean value denoting if any atom charge modifications
     *          were applied.
     */
    private static boolean resolveDiatomicChalcogenIons(Molecule3D mol, int atom) {
        if (!mol.isChalcogene(atom)) return false;
        // skip oxygen (O) for now
        if (mol.getAtomicNo(atom) == 8) return false;
        if (getAllConnAtoms(mol, atom, true) != 1) return false;
        if (mol.getConnAtoms(atom) != 1) return false;

        int nAtom = mol.getConnAtom(atom, 0);
        if (!mol.isChalcogene(nAtom)) return false;
        // skip oxygen (O) for now
        if (mol.getAtomicNo(nAtom) == 8) return false;
        if (mol.getConnAtoms(nAtom) != 1) return false;

        mol.setAtomCharge(atom, -1);
        mol.setAtomCharge(nAtom, -1);

        return true;
    }

    static void coordinationBondsToCovalent(Molecule3D mol) {
        for (int i = 0; i < mol.getRingSet().getSize(); i++ ) {
            if ( isRingCoordinated(mol, i) ) {
                int ringAtom = mol.getRingSet().getRingAtoms(i)[0];
                int metalAtom = -1;
                for ( int j = mol.getAllConnAtoms(ringAtom); j < mol.getAllConnAtomsPlusMetalBonds(ringAtom); j++) {
                    if ( mol.isMetalAtom( mol.getConnAtom(ringAtom, j) ) ) {
                        metalAtom = mol.getConnAtom(ringAtom, j);
                        break;
                    }
                }
                for (int j = 0; j < mol.getRingSet().getRingSize(i); j++) {
                    mol.markBondForDeletion(mol.getBond(mol.getRingSet().getRingAtoms(i)[j], metalAtom));
                }
            }
        }
        mol.deleteMarkedAtomsAndBonds();

        mol.ensureHelperArrays(Molecule3D.cHelperNeighbours);
        for (int i = 0; i < mol.getAllBonds(); i++ ) {
            if ( mol.getBondType(i) == Molecule3D.cBondTypeMetalLigand ) {
                mol.setBondOrder(i, 1);
                int a1 = mol.getBondAtom(0, i);
                int a2 = mol.getBondAtom(1, i);
                if ( mol.isMetalAtom(a1) && mol.isMetalAtom(a2) ) continue;
                if ( mol.isMetalAtom(a1) || mol.getAtomicNo(a1) == 1 ) {
                    mol.setAtomCharge( a1, mol.getAtomCharge(a1) - 1 );
                    mol.setAtomCharge( a2, mol.getAtomCharge(a2) + 1 );
                } else if ( mol.isMetalAtom(a2) || mol.getAtomicNo(a2) == 1  ) {
                    mol.setAtomCharge( a1, mol.getAtomCharge(a1) + 1 );
                    mol.setAtomCharge( a2, mol.getAtomCharge(a2) - 1 );
                } else {
                    throw new RuntimeException( "a metal-ligand bond was assigned to a non-metal atoms " +
                                                "'" + mol.getAtomName(a1) + "' and '" + mol.getAtomName(a2) + "'" );
                }
            }
        }
    }

    private static boolean isPlanar(Molecule3D mol, int a1, int a2) {
        Coordinates ci = mol.getCoordinates(a1);
        Coordinates u = null, v = null;

        for (int i = 0; v == null && i < mol.getAllConnAtoms(a1); i++) {
            if ( u == null ) {
                u = mol.getCoordinates(mol.getConnAtom(a1, i)).subC(ci);
            } else {
                v = mol.getCoordinates(mol.getConnAtom(a1, i)).subC(ci);
            }
        }
        for (int i = 0; v == null && i < mol.getAllConnAtoms(a2); i++) {
            if ( u == null ) {
                u = mol.getCoordinates(mol.getConnAtom(a2, i)).subC(ci);
            } else {
                v = mol.getCoordinates(mol.getConnAtom(a2, i)).subC(ci);
            }
        }

        if ( u == null ) return false;

        Coordinates normal = u.cross(v);
        if ( normal.distSq() == 0 ) return false; //what to do?
        normal = normal.unitC();

        Coordinates cj = mol.getCoordinates(a2);
        for(int k = 0; k < mol.getAllConnAtoms(a2); k++) {
            Coordinates ck = mol.getCoordinates(mol.getConnAtom(a2, k));
            if(Math.abs(ck.subC(cj).dot(normal))>0.25) return false;
        }
        for(int k=0; k<mol.getAllConnAtoms(a1); k++) {
            Coordinates ck = mol.getCoordinates(mol.getConnAtom(a1, k));
            if(Math.abs(ck.subC(cj).dot(normal))>0.25) return false;
        }
        return true;
    }

    /**
     * Util function for substructure searches. It searches for the substructure
     * of form a-b-c where a, b, c are consecutively connected atoms and a != c.
     * If a negative value is provided as the atomic number or valency of atom c
     * this atom will be ignored and only the a-b substructure search will be
     * carried out.
     *
     * @param mol
     *          The molecule that contains the atoms.
     * @param aIndex
     *          The index of the atom 'a' in the molecule. The atom 'a' is used
     *          as the starting position of the substructure search.
     * @param bAtomicNo
     *          The atomic number of the substructure atom 'b' that should be
     *          connected directly to the atom 'a'.
     * @param bValence
     *          The valence of the substructure atom 'b' that should be connected
     *          directly to the atom 'a'.
     * @param cAtomicNo
     *          The atomic number of the substructure atom 'c' that should be
     *          connected directly to the atom 'b'.
     * @param cValence
     *          The valence of the substructure atom 'c' that should be connected
     *          directly to the atom 'b'.
     * @return
     *          Index of the atom 'b' in the molecule is the search is successful
     *          and -1 otherwise.
     */
    private static int connectedAtom( Molecule3D mol, int aIndex, int bAtomicNo,
                                      int bValence, int cAtomicNo, int cValence ) {
        loop: for(int i=0; i<mol.getAllConnAtoms(aIndex); i++) {
            int atm = mol.getConnAtom(aIndex, i);

            // Checking if the atom 'b' that is connected to atom 'a'
            // is of the the requested chemical type and valence
            if ( bAtomicNo > 0 && mol.getAtomicNo(atm)     != bAtomicNo ) continue;
            if ( bValence  > 0 && mol.getAllConnAtoms(atm) != bValence  ) continue;

            // Checking if the atom 'c' that is connected to atom 'b'
            // is of the the requested chemical type and valence
            if ( cAtomicNo > 0 || cValence > 0 ) {
                for ( int j = 0; j < mol.getAllConnAtoms(atm); j++ ) {
                    int otherAtm = mol.getConnAtom(atm, j);
                    if(otherAtm==aIndex) continue loop;
                    if ( cAtomicNo > 0 && mol.getAtomicNo(otherAtm)     != cAtomicNo ) continue loop;
                    if ( cValence  > 0 && mol.getAllConnAtoms(otherAtm) != cValence  ) continue loop;
                }
            }

            return atm;
        }
        return -1;
    }

    private static double getBondLength(Molecule3D mol, int a1, int a2) {
        return mol.getCoordinates(a1).distance(mol.getCoordinates(a2));
    }

    private static double getBondLength(Molecule3D mol, int bond) {
        return getBondLength(mol, mol.getBondAtom(0, bond),  mol.getBondAtom(1, bond));
    }

    private static int connectedBond(Molecule3D mol, int a, int atomicNo, int valence) {
        if(a<0) return -1;
        for(int i=0; i< mol.getAllConnAtoms(a); i++) {
            int atm = mol.getConnAtom(a, i);
            if(atomicNo>0 && mol.getAtomicNo(atm)!=atomicNo) continue;
            if(valence>0 && mol.getAllConnAtoms(atm)!=valence) continue;
            return mol.getConnBond(a, i);
        }
        return -1;
    }

    /**
     * Returns the average bond length value based on bond order
     * and the type and pi electron count of the bonded atoms.
     * The values are taken from the statistical model
     * constructed from previously observed bond parameter
     * values.
     * @param bondOrder
     *          Bond order expressed as a positive integer.
     * @param atomicNo1
     *          Atomic number of the first atom.
     * @param atomicNo2
     *          Atomic number of the second atom.
     * @param atomPi1
     *          Pi electron count of the first atom.
     * @param atomPi2
     *          Pi electron count of the second atom.
     * @return
     *          The average bond length.
     */
    private static double getIdealisedBondLength(int bondOrder, int atomicNo1, int atomicNo2, int atomPi1, int atomPi2) {
        return BondLengthSet.getBondLength(BondLengthSet.getBondIndex(bondOrder, false, false, atomicNo1, atomicNo2, atomPi1, atomPi2));
    }

    /**
     * Produce a logical estimate of the bond being aromatic based
     * strictly on length constraints. The bond length is compared
     * to the average single bond.
     * @param mol
     *          The molecule that contains the bond.
     * @param a1
     *          Index of the first atom that forms a bond.
     * @param a2
     *          Index of the second atom that forms a bond.
     * @return
     *          Estimate deciding if the bond could be aromatic
     *          based strictly on the bond length.
     */
    private static boolean fitsAromaticLengthConstraints(Molecule3D mol, int a1, int a2) {
        int atomType1 = Math.min(mol.getAtomicNo(a1), mol.getAtomicNo(a2));
        int atomType2 = Math.max(mol.getAtomicNo(a1), mol.getAtomicNo(a2));

        int pi1 = BondLengthSet.isPiConsidered(atomType1) ? 1 : -1;
        int pi2 = BondLengthSet.isPiConsidered(atomType2) ? 1 : -1;

        // Try to get statistics for the aromatic single bond
        boolean aromaticConstraintsUsed = true;
        int index = BondLengthSet.getBondIndex(1, true, false, atomType1, atomType2, pi1, pi2);
        // If no such bond exist, try to get statistics for the
        // regular single bond
        if ( index == -1 ) {
            aromaticConstraintsUsed = false;
            index = BondLengthSet.getBondIndex(1, false, false, atomType1, atomType2, pi1, pi2);
            if (index == -1) return false;
        }
        double bondLength = BondLengthSet.getBondLength(index);
        double sigma = BondLengthSet.getBondStdDev(index);
        if ( sigma < 0 ) sigma = 0.1;
        double distance = getBondLength(mol, a1, a2);

        // allow 3 standard deviations for aromatic bonds and 1 for non-aromatic
        return ( distance < bondLength + ( aromaticConstraintsUsed ? 3 : 1 ) * sigma);
    }

    /**
     * Produce a numerical estimate of how well the bond length fits the statistical model
     * constructed from previously observed bond parameter values. The estimate is provided
     * as a deviation value from the idealised value (smaller values mean better fitness).
     * Currently it only compares the bond length with the average bond length of the same
     * type, but a more intricate model might be implemented in the future.
     * @param mol
     *          The molecule that contains the bond.
     * @param a1
     *          Index of the first atom that forms a bond.
     * @param a2
     *          Index of the second atom that forms a bond.
     * @param bondOrder
     *          Order of the bond. Currently only values '1' (single bond)
     *          and '2' (double bond) are supported.
     * @return
     *          Bond length deviation from the previously observed bond parameter value.
     */
    private static double estimatePiBondFitness(Molecule3D mol, int a1, int a2, int bondOrder) {

        int atomType1 = Math.min(mol.getAtomicNo(a1), mol.getAtomicNo(a2));
        int atomType2 = Math.max(mol.getAtomicNo(a1), mol.getAtomicNo(a2));

        double singleBond = getIdealisedBondLength(1, atomType1, atomType2, 1, 1);
        double doubleBond = getIdealisedBondLength(2, atomType1, atomType2, 1, 1);
        boolean statisticsAvailable = singleBond != -1 || doubleBond != -1;

        if (!statisticsAvailable) {
            throw new RuntimeException("bond length statistics for bond '" + atomType1
                    + "-" + atomType2 + "' is not available");
        }

        double idealLength = ( bondOrder == 1 ) ? singleBond : doubleBond;

        if ( bondOrder == 2 && getBondLength(mol, a1, a2) < idealLength ) {
            return 0;
        }
        return StrictMath.abs( getBondLength(mol, a1, a2) - idealLength );
    }

    private static boolean isAlkaliMetalAtom (Molecule3D mol, int atom) {
        return ( mol.getAtomicNo(atom) ==  3 ) ||
               ( mol.getAtomicNo(atom) == 11 ) ||
               ( mol.getAtomicNo(atom) == 19 ) ||
               ( mol.getAtomicNo(atom) == 37 ) ||
               ( mol.getAtomicNo(atom) == 55 ) ||
               ( mol.getAtomicNo(atom) == 87 );
    }

    private static boolean isAlkalineEarthMetalAtom(Molecule3D mol, int atom) {
        return ( mol.getAtomicNo(atom) ==  4 ) ||
               ( mol.getAtomicNo(atom) == 12 ) ||
               ( mol.getAtomicNo(atom) == 20 ) ||
               ( mol.getAtomicNo(atom) == 38 ) ||
               ( mol.getAtomicNo(atom) == 56 ) ||
               ( mol.getAtomicNo(atom) == 88 );
    }

    List<String> checkMoleculeForBumps(Molecule3D mol) {
        List<String> bumps = new LinkedList<>();
        double distance;
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            for (int j = i+1; j < mol.getAllAtoms(); j++) {
                if ( i < j ) {
                    distance = mol.getCoordinates(i).distance(mol.getCoordinates(j));
                    if ( distance < SPECIAL_POSITION_CUTOFF ||
                            distance < BUMP_FACTOR * ( COVALENT_RADIUS[mol.getAtomicNo(i)] +
                                                       COVALENT_RADIUS[mol.getAtomicNo(j)] ) ) {
                        bumps.add( "atoms '" + mol.getAtomName(i) + "' and '" + mol.getAtomName(j) + "' " +
                                "are too close (distance = " + String.format("%6.4f", distance) + ") and " +
                                "are considered a bump" );
                    }
                }
            }
        }
        return bumps;
    }

    private static boolean hasMetalLigandBonds(Molecule3D mol, int atom) {
        return (mol.getMetalBondedConnAtoms(atom) > 0);
    }

    private static boolean isCommonIonCharge(int atomicNo, int charge) {

        for ( int i = 0; i < Molecule.cCommonOxidationState[atomicNo].length; i++ ) {
            if ( Molecule.cCommonOxidationState[atomicNo][i] == charge ) {
                return true;
            } else if ( Molecule.cCommonOxidationState[atomicNo][i] > charge ) {
                return false;
            }
        }

        return false;
    }

    private static boolean isCommonIonCharge(Molecule3D mol, int atomicNo, int charge) {
        return isCommonIonCharge(mol.getAtomicNo(atomicNo), charge);
    }

    public static int getMaximumIonCharge(Molecule3D mol, int atom) {
        if ( Molecule.cCommonOxidationState[mol.getAtomicNo(atom)] == null ) {
            return 0;
        } else {
            return Molecule.cCommonOxidationState[mol.getAtomicNo(atom)]
                    [Molecule.cCommonOxidationState[mol.getAtomicNo(atom)].length - 1];
        }
    }

    private static boolean maximumIonChargeReached(Molecule3D mol, int atom) {
        return mol.getAtomCharge(atom) >= getMaximumIonCharge(mol, atom);
    }

    private static int getDelocalizedBondCount(Molecule3D mol, int atom) {
        return getDelocalizedBonds(mol, atom).size();
    }

    private static List<Integer> getDelocalizedBonds(Molecule3D mol, int atom) {
        List<Integer> delocalizedBonds = new LinkedList<>();
        for( int i = 0; i < mol.getAllConnAtoms(atom); i++) {
            if ( mol.getBondType(mol.getConnBond(atom, i)) == Molecule3D.cBondTypeDelocalized ) {
                delocalizedBonds.add(mol.getConnBond(atom, i));
            }
        }
        return delocalizedBonds;
    }

    private static List<List<Integer>> getPiBondFragments (Molecule3D mol) {

        List<List<Integer>> fragments = new LinkedList<>();
        mol.ensureHelperArrays(Molecule3D.cHelperRings);
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            if (getDelocalizedBondCount(mol, i) == 1) {
                // Check if the selected atom is not already at the end of analyzed fragments
                boolean used = false;
                for (List<Integer> fragment1 : fragments) {
                    if (fragment1.get(fragment1.size() - 1) == i) {
                        used = true;
                        break;
                    }
                }
                if ( used ) { continue; }

                List<Integer> fragment = new LinkedList<>();
                fragment.add(i);
                int prevBond = getDelocalizedBonds(mol, i).get(0);
                int atom;
                if ( mol.getBondAtom( 0, prevBond ) != i ) {
                    atom = mol.getBondAtom( 0, prevBond );
                } else {
                    atom = mol.getBondAtom( 1, prevBond );
                }

                boolean finished = false;
                while (!finished) {
                    // Currently only process linear pi-bond systems
                    // TODO: extend the logic to process branched systems
                    if (getDelocalizedBondCount(mol, atom) > 2) {
                        break;
                    }

                    fragment.add(atom);
                    if (getDelocalizedBondCount(mol, atom) == 1) {
                        fragments.add(fragment);
                        finished = true;
                    } else {
                        int currentBond = ( getDelocalizedBonds(mol, atom).get(0) != prevBond ) ?
                                getDelocalizedBonds(mol, atom).get(0) :
                                getDelocalizedBonds(mol, atom).get(1);
                        if ( mol.getBondAtom( 0, currentBond ) != atom ) {
                            atom = mol.getBondAtom( 0, currentBond );
                        } else {
                            atom = mol.getBondAtom( 1, currentBond );
                        }
                        prevBond = currentBond;
                    }
                }
            }
        }

        return fragments;
    }

    private static boolean fragmentCanResonate(Molecule3D mol, List<Integer> fragment) {
        for (Integer atom : fragment) {
            if (mol.getFreeValence(atom) - 1 < 0) {
                return false;
            }
        }
        return true;
    }

    private static int[] sortAtomsByBondLength(Molecule3D mol, int a, List<Integer> bondAtoms) {
        double[] bondLength = new double[bondAtoms.size()];
        for (int i = 0; i < bondLength.length; i++ ) {
            bondLength[i] = getBondLength(mol, a, mol.getConnAtom(a, bondAtoms.get(i)));
        }
        // using primitive bubble sort for now to sort the indexes
        int[] indexes = new int[bondAtoms.size()];
        for (int i = 0; i < bondLength.length; i++) {
            indexes[i] = i;
        }

        for (int i = 0; i < bondLength.length; i++) {
            for (int j = i+1; j < bondLength.length; j++ ) {
                if( bondLength[indexes[i]] > bondLength[indexes[j]] ) {
                    int tmp = indexes[i];
                    indexes[i] = indexes[j];
                    indexes[j] = tmp;
                }
            }
        }
        // replacing sort indexes with atom indexes
        for (int i = 0; i < bondLength.length; i++) {
            indexes[i] = bondAtoms.get(indexes[i]);
        }

        return indexes;
    }

    private static Coordinates calculatePlaneNormal (Coordinates p1, Coordinates p2, Coordinates p3) {
        Coordinates u = p1.subC(p2);
        Coordinates v = p1.subC(p3);
        return u.cross(v);
    }

    private static double calculatePointToPlaneDistance ( Coordinates normal, Coordinates p0, Coordinates p1 ) {
        Coordinates v1 = p0.subC(p1);
        return StrictMath.abs(normal.unitC().dot(v1) / v1.dist());
    }

    private static void setTerminalDoubleBonds( Molecule3D mol, int a, int[] sortedAtoms,
                                                int maxDoubleBondCount ) {
        int doubleBondCount = 0;
        for (int sortedAtom : sortedAtoms) {
            if (doubleBondCount < maxDoubleBondCount) {
                mol.setBondOrder(mol.getConnBond(a, sortedAtom), 2);
                doubleBondCount++;
            } else {
                mol.setBondOrder(mol.getConnBond(a, sortedAtom), 1);
                mol.setAtomCharge(mol.getConnAtom(a, sortedAtom), -1);
            }
        }
    }

    // TODO: ask to add to the molecule
    private static boolean isRingAtom(Molecule3D mol, int ring, int a) {
        int[] ringAtoms = mol.getRingSet().getRingAtoms(ring);
        for (int ringAtom : ringAtoms) {
            if (ringAtom == a) return true;
        }
        return false;
    }

    private static void identifySquarates(Molecule3D mol, int[] spOrder) {
        for (int i = 0; i < mol.getAtoms(); i++) {
            if (mol.getAtomRingSize(i) != 4) continue;
            if (getRingCount(mol, i) != 1) continue;

            int ring = getSmallestRing(mol, i);
            int[] ringAtoms = mol.getRingSet().getRingAtoms(ring);
            boolean[] allowsExocyclicDoubleBond = new boolean[ringAtoms.length];
            int[] exocyclicAtoms = new int[ringAtoms.length];
            boolean fitsConstraints = true;
            for (int j = 0; j < ringAtoms.length; j++) {
                int a = ringAtoms[j];
                if (mol.getAtomicNo(a) == 6 && mol.getConnAtoms(a) == 3 &&
                        spOrder[a] == 2 && mol.getLowestFreeValence(a) > 0) {
                    int a1 = ringAtoms[(j - 1 + ringAtoms.length) % ringAtoms.length];
                    int a2 = ringAtoms[(j + 1 + ringAtoms.length) % ringAtoms.length];

                    if (StrictMath.abs(90 - StrictMath.toDegrees(
                            getAngle(mol, a1, a, a2))) > 5) {
                        fitsConstraints = false;
                        break;
                    }
                    for (int k = 0; k < mol.getConnAtoms(a); k++) {
                        if (!isRingAtom(mol, ring, mol.getConnAtom(a, k))) {
                            allowsExocyclicDoubleBond[j] = (mol.getLowestFreeValence(mol.getConnAtom(a, k)) > 0);
                            exocyclicAtoms[j] = mol.getConnAtom(a, k);
                            break;
                        }
                    }
                } else {
                    fitsConstraints = false;
                    break;
                }
            }

            if (fitsConstraints) {
                for (int a1 = 0; a1 < ringAtoms.length; a1++) {
                    int a2 = (a1 + 1) % ringAtoms.length;
                    int a3 = (a1 + 2) % ringAtoms.length;
                    int a4 = (a1 + 3) % ringAtoms.length;
                    if (allowsExocyclicDoubleBond[a1] && allowsExocyclicDoubleBond[a2]) {
                        mol.setBondOrder(mol.getBond(ringAtoms[a1], exocyclicAtoms[a1]), 2);
                        mol.setBondOrder(mol.getBond(ringAtoms[a2], exocyclicAtoms[a2]), 2);
                        mol.setBondOrder(mol.getBond(ringAtoms[a3], ringAtoms[a4]), 2);
                        if (allowsExocyclicDoubleBond[a3]) mol.setAtomCharge(exocyclicAtoms[a3], -1);
                        if (allowsExocyclicDoubleBond[a4]) mol.setAtomCharge(exocyclicAtoms[a4], -1);
                        break;
                    }
                }
            }
        }
    }

    private static int getRingCount(Molecule3D mol, int a) {
        int ringCount = 0;
        RingCollection rings = mol.getRingSet();
        for (int i = 0; i < mol.getRingSet().getSize(); i++) {
            int[] ringAtoms = rings.getRingAtoms(i);
            for (int ringAtom : ringAtoms) {
                if (ringAtom == a) {
                    ringCount++;
                    break;
                }
            }
        }
        return ringCount;
    }

    private static int getSmallestRing(Molecule3D mol, int a) {
        int smallestRing = -1;
        int ringSize = -1;
        RingCollection rings = mol.getRingSet();
        for (int i = 0; i < mol.getRingSet().getSize(); i++) {
            int[] ringAtoms = rings.getRingAtoms(i);
            for (int ringAtom : ringAtoms) {
                if (ringAtom == a) {
                    if (ringSize == -1 || ringAtoms.length < ringSize) {
                        ringSize = ringAtoms.length;
                        smallestRing = i;
                        break;
                    }
                }
            }
        }
        return smallestRing;
    }

    private static int getAllConnAtoms(Molecule3D mol, int a, boolean includeAttachedHydrogens) {
        if ( includeAttachedHydrogens ) {
            return mol.getAllConnAtoms(a) + mol.getAttachedHydrogenCount(a);
        } else {
            return mol.getAllConnAtoms(a);
        }
    }

    private static boolean isRingCoordinated(Molecule3D mol, int ringNo) {
        int[] ringAtoms = mol.getRingSet().getRingAtoms(ringNo);
        boolean atomsAreCoordinated = true;
        for (int ringAtom : ringAtoms) {
            atomsAreCoordinated = atomsAreCoordinated && hasMetalLigandBonds(mol, ringAtom);
        }
        return atomsAreCoordinated;
    }

    /**
     * Determines if the given atom is capable of accommodating an additional double bond.
     * @param mol
     *          Molecule that contains the atom.
     * @param atom
     *          Index of the atom inside the molecule data structure.
     * @param atomPi
     *          The number of pi electrons that the atom is currently using.
     * @param spOrder
     *          The sp order (1, 2, 3) of the atom.
     * @return
     *          Logical value denoting is an atom can accommodate an additional double bond.
     */
    private static boolean canAccommodateDoubleBond(Molecule3D mol, int atom, int atomPi, int spOrder) {

        // Tetravalent silicon (Si) cannot accommodate any double bonds
        if (mol.getAtomicNo(atom) == 14 && mol.getOccupiedValence(atom) > 3) return false;
        // Germanium (Ge)
        if (mol.getAtomicNo(atom) == 32) {
            // only sp3 germanium (Ge) can potentially accommodate a double bond
            if (mol.getOccupiedValence(atom) > 3) return false;
            if (spOrder != 2) return false;
        }

        // Phosphorus (P)
        if (mol.getAtomicNo(atom) == 15) {
            if (atomPi > 0 && spOrder != 1) return false;
        }

        if ( mol.getLowestFreeValence(atom) < 1 ) return false;
        if ( mol.getAtomCharge(atom) != 0 ) return false;
        // Sp1-like atoms are not allowed to have more than 2 pi electrons
        if ( spOrder == 1 && atomPi > 1 ) return false;
        // Sp2-like atoms are not allowed to have more than 1 pi electron
        if ( spOrder == 2 && atomPi > 0 ) return false;
        // Atoms smaller than phosphorus are not allowed a double bond in Sp3-like configuration
        if ( spOrder == 3 && mol.getAtomicNo(atom) < 14 ) return false;

        return true;
    }

    /**
     * Calculates the bond length deviation from the statistical
     * average for the bond assuming it is of the given bond type.
     * @param mol
     *          Molecule that contains the bond.
     * @param bond
     *          The bond that is being evaluated.
     * @param bondTypeIndex
     *          The bond type index that is used to retrieve
     *          the average bond length and standard deviation
     *          for the given bond type.
     * @return
     *          The difference between the actual and the average
     *          bond lengths.
     */
    private static double determineExpectedBondDeviation(Molecule3D mol, int bond, int bondTypeIndex) {
        return ( getBondLength(mol, bond ) - BondLengthSet.getBondLength(bondTypeIndex));
    }

    /**
     * Calculates the angle a1-a2-a3 between 3 atoms.
     * @param mol
     *          The molecule that contains the atoms.
     * @param a1
     *          Index of the first atom in the molecule.
     * @param a2
     *          Index of the second atom in the molecule.
     * @param a3
     *          Index of the third atom in the molecule.
     * @return
     *          The angle a1-a2-a3 in radians.
     */
    private static double getAngle(Molecule mol, int a1, int a2, int a3) {
        Coordinates c1 = mol.getCoordinates(a1);
        Coordinates c2 = mol.getCoordinates(a2);
        Coordinates c3 = mol.getCoordinates(a3);

        return c1.subC(c2).getAngle(c3.subC(c2));
    }

    /**
     * Calculates the dihedral angle a1-a2-a3-a4 between 4 atoms.
     * @param mol
     *         The molecule that contains the atoms.
     * @param a1
     *         Index of the first atom in the molecule.
     * @param a2
     *         Index of the second atom in the molecule.
     * @param a3
     *         Index of the third atom in the molecule.
     * @param a4
     *         Index of the fourth atom in the molecule.
     * @return
     *         The dihedral angle a1-a2-a3-a4 in radians.
     */
    private static double getDihedral(Molecule mol, int a1, int a2, int a3, int a4) {
        Coordinates c1 = mol.getCoordinates(a1);
        Coordinates c2 = mol.getCoordinates(a2);
        Coordinates c3 = mol.getCoordinates(a3);
        Coordinates c4 = mol.getCoordinates(a4);
        return c1.getDihedral(c2, c3, c4);
    }
    
	public static int connected(StereoMolecule mol, int a, int atomicNo, int bondOrder) {
		for(int i=0; i<mol.getAllConnAtoms(a); i++) {
			int atm = mol.getConnAtom(a, i);
			if(atomicNo>=0 && mol.getAtomicNo(atm)!=atomicNo) continue;
			if(bondOrder>0 && mol.getConnBondOrder(a, i)!=bondOrder) continue;
			return atm;
		}
		return -1;
	}
    
	public static boolean aromatize(Molecule3D mol, Set<Integer> aromaticAtoms, Set<Integer> aromaticBonds) {
		ArrayList<Integer>[] atomToRings = getAtomToRings(mol);
		RingCollection ringSet = mol.getRingSetSimple();
		return aromatize(mol,atomToRings,ringSet,aromaticAtoms,aromaticBonds);
	}

	public static boolean aromatize(Molecule3D mol, ArrayList<Integer>[] atomToRings, RingCollection ringSet, Set<Integer> aromaticAtoms, Set<Integer> aromaticBonds) {
		//RingCollection ringSet = mol.getRingSetSimple();

		//Flag the aromatic rings
		boolean[] aromaticRings = new boolean[ringSet.getSize()];
		for (int i = 0; i < ringSet.getSize(); i++) {
			//Is is an aromatic ring
			boolean isAromatic = true;
			int ringSize = ringSet.getRingSize(i);
			for (int j = -1; j < ringSize; j++) {
				int[] ringAtom = ringSet.getRingAtoms(i);
				int a1 = j==-1? ringAtom[ringSize-1]: ringAtom[j];
				int a2 = j==ringSize-1? ringAtom[0]: ringAtom[j+1];
				
				int b = mol.getBond(a1, a2);
				if(!aromaticBonds.contains(b)) {
					isAromatic = false;
				}
			}
				
			aromaticRings[i] = isAromatic;
		}
		Set<Integer> nonAromaticAtoms = new HashSet<Integer>();
		for (int i=0;i<mol.getAllAtoms();i++) nonAromaticAtoms.add(i);
		nonAromaticAtoms.removeAll(aromaticAtoms);		
		
		//Launch the aromatizer
		boolean ok = true;
		for (int i = 0; i < aromaticRings.length; i++) {
			if(aromaticRings[i]) {
				boolean success = aromatize(mol, atomToRings, ringSet, i, aromaticRings, nonAromaticAtoms, new boolean[mol.getAllAtoms()], 0, ringSet.getRingSize(i)%2, new ArrayList<Integer>(), true);
				if(!success) success = aromatize(mol, atomToRings, ringSet, i, aromaticRings, nonAromaticAtoms, new boolean[mol.getAllAtoms()], 0, ringSet.getRingSize(i)%2, new ArrayList<Integer>(), false);
				if(!success) {
					System.out.println("Could not aromatize ring "+i);
					aromaticRings[i] = false;
					ok = false;
				}				
			}
		}
		return ok; 
	}
	
	private static boolean aromatize(Molecule3D mol, ArrayList<Integer>[] atomToRings, RingCollection ringSet, int index, boolean[] aromatic, Set<Integer> nonAromaticAtoms, boolean[] visited, int seed, int left, List<Integer> bondsMade, boolean easy) {
		//RingCollection ringSet = mol.getRingSetSimple();

		//Ends if the ring has been fully visited
		int[] ring = (int[]) ringSet.getRingAtoms(index);
		boolean allVisited = true;
		int bestSeed = -1;
		for (int i = 0; i < ring.length; i++) {
			if(!visited[ring[(seed + i)%ring.length]]) {				 
				if(bestSeed<0) bestSeed = seed + i; 
				allVisited = false;				
			}
		}
		if(allVisited) {
			return true;
		} else  {
			seed = bestSeed;
		}
		
		
		int a = ring[seed%ring.length];
		int ap = ring[(seed+1)%ring.length];
		if(visited[a]) { //already treated
			
			return aromatize(mol, atomToRings, ringSet, index, aromatic, nonAromaticAtoms, visited, seed+1, left, bondsMade, easy);
			
		} else {
		
			//Try to create a double bond from the atom a to a connected one
			for (int j = -1; j < mol.getAllConnAtoms(a); j++) {
				int a2 = j==-1? ap: mol.getConnAtom(a, j);
				
				if(visited[a2]) continue;
				if(j>=0 && a2==ap) continue;
				
				if(nonAromaticAtoms.contains(a2)) continue;
				if(nonAromaticAtoms.contains(a)) continue;
				
				if(mol.getAtomicNo(a)==8 || mol.getAtomicNo(a)==16) continue;
				if(mol.getAtomicNo(a2)==8 || mol.getAtomicNo(a2)==16) continue;
				if(easy && mol.getFreeValence(a) <= 0) continue;
				if(easy && mol.getFreeValence(a2) <= 0) continue;
				if(connected(mol, a, -1, 2)>=0) continue;
				if(connected(mol, a2, -1, 2)>=0) continue;
				
				visited[a] = visited[a2] = true;
				int b = mol.getBond(a, a2);
				mol.setBondOrder(b, 2);
				
				//Test whole ring 
				List<Integer> trackBondsMade = new ArrayList<Integer>();
				boolean success = aromatize(mol, atomToRings, ringSet, index, aromatic, nonAromaticAtoms, visited, seed+1, left, trackBondsMade, easy);

				//Test connecting rings
				if(success) {
					List<Integer> rings = atomToRings[a2];
					if(rings.size()>1) {
						for (int r : rings) {
							if(r!=index && r>=0 && r<aromatic.length && aromatic[r]) {
								int newSeed;
								for(newSeed=0; ringSet.getRingAtoms(r)[newSeed]!=a2; newSeed++) {}
							
								//System.out.println("try connected ring "+r+" / "+newSeed);
								success = aromatize(mol, atomToRings, ringSet, r, aromatic, nonAromaticAtoms, visited, newSeed, ringSet.getRingSize(r)%2, trackBondsMade, easy);
								//System.out.println("try connected ring "+r +" -> " +success+" "+trackBondsMade.size());
								if(!success) break;
							}
						}
					}
				}
				
				if(success) {
					//It works!!!
					bondsMade.add(b);
					bondsMade.addAll(trackBondsMade);
					return true;
				} else {
					//Backtrack changes
					visited[a] = visited[a2] = false;
					mol.setBondOrder(b, 1);
					for (int b2 : trackBondsMade) {
						//System.out.println("retrack "+mol.getBondAtom(0, b2)+"-"+mol.getBondAtom(1, b2));
						mol.setBondOrder(b2, 1);
						visited[mol.getBondAtom(0, b2)] = visited[mol.getBondAtom(1, b2)] = false;
					}
				} 
			}
			
			//Try to skip this atom
			if(left>0 && (mol.getAtomicNo(a)!=6)) {
				visited[a] = true;
				List<Integer> trackBondsMade = new ArrayList<Integer>();
				boolean success = aromatize(mol, atomToRings, ringSet, index, aromatic, nonAromaticAtoms, visited, seed+1, left-1, trackBondsMade, easy);
				if(success) {
					bondsMade.addAll(trackBondsMade);
					return true;
				} else {
					visited[a] = false;
					for (int b2 : trackBondsMade) {
						mol.setBondOrder(b2, 1);
						visited[mol.getBondAtom(0, b2)] = visited[mol.getBondAtom(1, b2)] = false;
					}
				}
			}
			
			return false;
		}
	}
}