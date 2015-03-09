/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Christian Rufener
 */

package com.actelion.research.chem;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.BitSet;
import java.util.Enumeration;
import java.util.Hashtable;


/**
 * Generator of a path-based Fingerprint
 * Not thread safe!
 */
public class FingerPrintGenerator
{

    private static final int MAX_BITS = 512;
    private static final int MAX_DEPTH = 6;
    private static final boolean DEBUG = false;
    private static int debugCounter = 0;
    private Hashtable paths;


    /**
     * Generates a fingerprint of the default size for the given ExtendedMolecule
     *
     * @param mol The ExtendedMolecule for which a Fingerprint is generated
     * @return The Fingerprint (A one-dimensional bit array)
     */
    public BitSet getFingerprint(StereoMolecule mol)
    {
        return getFingerprint(mol, MAX_BITS);
    }


    /**
     * Generates the fingerprint of a molecule
     *
     * @param mol  The ExtendedMolecule for which a Fingerprint is generated
     * @param size The desired size of the fingerprint
     * @return The Fingerprint (A one-dimensional bit array)
     */
    public BitSet getFingerprint(StereoMolecule mol, int size)
    {
        String path = null;
        int index;

        mol.ensureHelperArrays(Molecule.cHelperRings);
        findAllPaths(mol);
        BitSet bs = new BitSet(size);
        for (Enumeration e = paths.elements(); e.hasMoreElements(); ) {
            path = (String) e.nextElement();
            index = new java.util.Random(path.hashCode()).nextInt(size);
            if (DEBUG)
                System.out.print(this.toString().substring(40) + " ");
//                System.out.println("Setting bit " + position + " for " + path);
            bs.set(index);
        }
        return bs;
    }


    /**
     * Checks whether query is contained  (i.e. subset of a possible sub-structure) of reference
     *
     * @param reference Reference BitSet
     * @param query     Query BitSet
     * @return true, if query is a subset of reference
     * @keyword substructure search
     */
    private static boolean matches(BitSet reference, BitSet query)
    {
        BitSet cp = (BitSet) reference.clone();
        cp.and(query);
        if (cp.equals(query)) {
            return true;
        }
        return false;
    }

    public static byte[] getBitSetBits(java.util.BitSet bs)
    {
        int size = bs.size();
        int t = size % 8 == 0 ? size / 8 : size / 8 + 1;
        byte r[] = new byte[t];
        int bite;
        byte bit;
        if (DEBUG)
            System.out.println();
        for (int i = 0; i < size; i++) {
            if (bs.get(i)) {
                bit = (byte) (1 << (7 - (i % 8)));
                bite = i / 8;
                r[bite] |= bit;
                if (DEBUG)
                    System.out.print("1");
            } else {
                if (DEBUG)
                    System.out.print("0");
            }
        }
        if (DEBUG)
            System.out.println();

        return r;
    }

    /**
     * Find all paths up to MAX_DEPTH length.
     * The paths are acquired by a number of depth first searches, one for each atom.
     *
     * @param mol The Molecule which is to be searched.
     */
    private void findAllPaths(StereoMolecule mol)
    {
        paths = new Hashtable();
        debugCounter = 0;
        int atoms = mol.getAllAtoms();
        String s;
        boolean[] flags = new boolean[atoms];
        for (int atom = 0; atom < atoms; atom++) {
            s = mol.getAtomLabel(atom);
            for (int i = 0; i < atoms; i++)
                flags[i] = false;
            addPath(s);
            if (DEBUG)
                System.out.println("\t***Starting at atom " + atom + " with symbol " + mol.getAtomLabel(atom));
            traverseDFS(mol, -1, atom, s, 0, flags);
        }
    }


    /**
     * Performs a recursive depth first search
     *
     * @param mol          The Molecule to be searched
     * @param lastAtom     The Atom we came from
     * @param rootAtom     The Atom to start the search at
     * @param path         The Path that has been generated so far
     * @param depth        The current depth in this recursive search
     * @param flags        Helper flags
     */
    private void traverseDFS(StereoMolecule mol, int lastAtom, int rootAtom, String path, int depth, boolean flags[])
    {
        int connAtoms = mol.getConnAtoms(rootAtom);
        int nextAtom = 0, bond = 0;
        StringBuilder newPath = new StringBuilder();

        depth++;
        // Check all connected atoms
        for (int i = 0; i < connAtoms; i++) {
            // Flag starting point
            flags[rootAtom] = true;
            nextAtom = mol.getConnAtom(rootAtom, i);
            if (nextAtom == lastAtom) {
                continue;
            }

            // Not yet seen
            if (!flags[nextAtom]) {
                bond = mol.getConnBond(rootAtom, i);
                newPath.setLength(0);
                newPath.append(path);
                // Construct a string with Daylight kind of bond encoding
                if (mol.isDelocalizedBond(bond) || mol.isAromaticBond(bond)) {
                    newPath.append(":");
                } else if (mol.getBondOrder(bond) == 1) {
                    newPath.append("-");
                } else if (mol.getBondOrder(bond) == 2) {
                    newPath.append("=");
                } else if (mol.getBondOrder(bond) == 3) {
                    newPath.append("#");
                } else
                    System.out.println(
                        "FingerPrintGenerator.depthFirstSearch() " +
                            "Error: Invalid Bond order! " + mol.getBondOrder(bond));
                newPath.append(mol.getAtomLabel(nextAtom));

                // Mark nextAtom as visited
                flags[nextAtom] = true;
                addPath(newPath.toString());
                if (depth == MAX_DEPTH) {
                    // Unmark nextAtoms so we might visit it again during another run
                    flags[nextAtom] = false;
                    // System.out.println("Reached end of depth!!! " + newPath);
                }
                if (depth < MAX_DEPTH) {
                    // System.out.println("Visiting " + rootAtom + " + " + nextAtom + " " + newPath);
                    traverseDFS(mol, rootAtom, nextAtom, newPath.toString(), depth, flags);
                    // Unmark nextAtoms so we might visit it again during another run
                    flags[nextAtom] = false;
                }
            }
        }
    }

    private boolean addPath(String newPath)
    {
        String storePath = newPath;
        String reversePath = new StringBuilder(storePath).reverse().toString();
        boolean ok = false;
        // Find the "smaller" version
        if (reversePath.compareTo(newPath) < 0) {
            /* reversePath is smaller than newPath so we keep reversePath */
            storePath = reversePath;
        }
        if (DEBUG)
            System.out.println("Checking for existence of Path " + storePath);
        if (!paths.containsKey(storePath)) {
            paths.put(storePath, storePath);
            ok = true;
            if (DEBUG) {
                debugCounter++;
                System.out.println("Storing path no. " + debugCounter + ": " + storePath + ", Hash: " + storePath.hashCode());
            }
        }
        return ok;

    }


    public static void main(String args[])
    {
        String query;
        String test = "sFp@DiTt@@@@ S~x>xixix>";
        IDCodeParser p = new IDCodeParser(false);
        StereoMolecule m = new StereoMolecule();
        BitSet referencebs = null, querybs = null;
        BufferedReader r = new BufferedReader(new InputStreamReader(System.in));
        FingerPrintGenerator fp = null;
        try {
            if (args.length > 0) {
                System.out.println("Query set");
                test = args[0];
                p.parse(m, test);
                m.ensureHelperArrays(Molecule.cHelperRings);
                fp = new FingerPrintGenerator();
                querybs = fp.getFingerprint(m);
                while (true) {
                    System.out.print("Molecule:");
                    query = r.readLine();
                    if (query != null && query.trim().length() > 0) {
                        p.parse(m, query);
                        m.ensureHelperArrays(Molecule.cHelperRings);
                        fp = new FingerPrintGenerator();
                        referencebs = fp.getFingerprint(m);
                        if (FingerPrintGenerator.matches(referencebs, querybs))
                            System.out.println(query.trim() + "\tOK");
                        else
                            System.out.println(query.trim() + "\tNOT OK");
                    }
                }
            } else {
                test = r.readLine();
                p.parse(m, test);
                m.ensureHelperArrays(Molecule.cHelperRings);
                fp = new FingerPrintGenerator();
                querybs = fp.getFingerprint(m);
                int count = 0;
                int notok = 0;

                while (true) {
                    query = r.readLine();
                    if (query != null && query.trim().length() > 0) {
                        System.out.println("*********************");
                        query = query.trim();
                        p.parse(m, query);
                        m.ensureHelperArrays(Molecule.cHelperRings);
                        fp = new FingerPrintGenerator();
                        referencebs = fp.getFingerprint(m);
                        if (!FingerPrintGenerator.matches(referencebs, querybs))
                            ++notok;
                        count++;
                    } else
                        break;
                }
                System.out.println("Total/Not OK " + count + " " + notok);
            }
        } catch (Exception e) {
            System.err.println("Server Exception: " + e);
            e.printStackTrace();
        }

    }
}
