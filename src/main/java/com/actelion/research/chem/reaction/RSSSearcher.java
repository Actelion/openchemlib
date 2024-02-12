/*
 * Project: DD_core
 * @(#)ReactionRSS.java
 *
 * Copyright (c) 1997- 2014
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * Author: Christian Rufener
 */

package com.actelion.research.chem.reaction;

import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.ArrayUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Project:
 * User: rufenec
 * Date: 9/30/2014
 * Time: 1:38 PM
 */
public class RSSSearcher
{

    private static final boolean debug = false;
    private static final boolean debug2 = false;

    private static void debug(String format,Object ...args)
    {
        if (debug)
            System.out.printf(format,args);
    }

    private static void debug2(String format,Object ...args)
    {
        if (debug2)
            System.out.printf(format,args);
    }

    private static class MapsList extends ArrayList<int[]>
    {
        public String toString()
        {
            String s = new String();
            for (int[] a : this) {
                s += "\n {";
                for (int i : a) {
                    s += " " + i;
                }
                s += "}";
            }
            return s;
        }
    }


    /**
     * Try to match a query reaction in a target reaction
     * Algo:
     * First we check that there are more product reactants than query reactants
     * Idem for products
     *
     * @param queryRxn
     * @param targetRxn
     * @return
     */
    public static boolean match(Reaction queryRxn, Reaction targetRxn)
    {

        List<int[]> reactantMatchList = new ArrayList<int[]>();
        List<int[]> productMatchList = new ArrayList<int[]>();

        int numOfQueryReactants = queryRxn.getReactants();
        int numOfTargetReactants = targetRxn.getReactants();
        int numOfQueryProducts = queryRxn.getProducts();
        int numOfTargetProducts = targetRxn.getProducts();

        if (numOfQueryReactants > numOfTargetReactants || numOfQueryProducts > numOfTargetProducts)
            return false;

        debug("MATCHING Start\n");

        if (numOfQueryReactants <= numOfTargetReactants) {
            debug("MATCHING Reactants\n");
            List<StereoMolecule> reactantList = getReactants(targetRxn);
            ListPermutator<StereoMolecule> permutator = new ListPermutator<StereoMolecule>(reactantList);

            while (permutator.hasNext()) {
                List<StereoMolecule> targetReactants = permutator.next();
                List<int[]> matchList = new ArrayList<int[]>();
                for (int i = 0; i < numOfQueryReactants; i++) {
                    StereoMolecule queryMol = queryRxn.getReactant(i);
                    // Since the targetReactant list is a permutation,
                    // we use the same index as for the query
                    StereoMolecule targetMol = targetReactants.get(i);
//                    debug("Attempt to match\n%s : %s\n%s : %s\n",
//                        queryMol.getIDCode(), queryMol,
//                        targetMol.getIDCode(), targetMol);
                    debug2("Matching reactant\n");
                    MapsList matched = findMatchingMaps(queryMol, targetMol);
                    // molecules did not even match on SSS, so try another permutation
                    if (matched == null) {
                        debug("Did not match SSS\n");
                        break;
                    }
                    if (matched.size() > 0) {
                        debug("Matched with matchlist\n");
                        matchList = add(matchList, matched);
                    } else {
                        debug("NO match with matchlist\n");
                        //break; // this seemed to be wrong
                    }
                }
                // Choose the best match:
                // We consider longer match list as a better solution
                if (getMaxSize(matchList) > getMaxSize(reactantMatchList)) {
                    reactantMatchList = matchList;
                }
            }
            debug("Reactant matchlist %d\n", reactantMatchList.size());
        }


        if (numOfQueryProducts <= numOfTargetProducts) {
            debug("MATCHING PRODUCTS\n");
            List<StereoMolecule> productList = getProducts(targetRxn);
            ListPermutator<StereoMolecule> permute = new ListPermutator<StereoMolecule>(productList);
            while (permute.hasNext()) {
//                debug("Next product permutation\n");
                List<StereoMolecule> targetProducts = permute.next();
                List<int[]> matchList = new ArrayList<int[]>();
                for (int i = 0; i < numOfQueryProducts; i++) {
                    StereoMolecule queryMol = queryRxn.getProduct(i);
                    StereoMolecule target = targetProducts.get(i);
//                    debug("Attempt to match\n%s : %s\n%s : %s\n",
//                        queryMol.getIDCode(), queryMol,
//                        target.getIDCode(), target);
                    debug2("Matching product\n");
                    List<int[]> matched = findMatchingMaps(queryMol, target);
                    // molecules did not even match on SSS, so try another permutation
                    if (matched == null) {
                        debug("Did not match SSS\n");
                        break;
                    }
                    if (matched.size() > 0) {
                        debug("Matched with matchlist\n");
                        matchList = add(matchList, matched);
                    } else {
                        debug("NO match with matchlist\n");
                        //break; // this seemed to be wrong
                       // break;
                    }
                }
                if (getMaxSize(matchList) > getMaxSize(productMatchList)) {
                    // Found a better match, so use this?
                    productMatchList = matchList;
                }
            }
            debug("Product matchlist ist " + productMatchList.size());
        }
        boolean ok = false;

        // Open issues
        // What if  reactant Matchlist is empty
//        for (int[] rs : reactantMatchList) {
////            Arrays.sort(rs);
//            debug("Reactant List: ");
//            for (int j : rs) {
//                debug("%d,", j);
//            }
//        }
//         debug("\n");

//        for (int[] ps : productMatchList) {
////            Arrays.sort(ps);
//            debug("Product List: ");
//            for (int j : ps) {
//                debug("%d,", j);
//            }
//        }
//        debug("\n");

        boolean sort = true;
        if (sort) {
            debug2("Sorting\n");
            for (int[] rs : reactantMatchList) {
                Arrays.sort(rs);
            }
            for (int[] ps : productMatchList) {
                Arrays.sort(ps);
            }
        }

        for (int[] rs : reactantMatchList) {
            for (int[] ps : productMatchList) {
                if (Arrays.equals(rs, ps)) {
                    ok = true;
                    return ok;
                }
            }
        }
//        debug("Query did not match!");
        return ok;
    }

    /**
     * Returns a List of sorted arrays containing the mapping numbers of the target molecule
     * which have been matched by the query sub-structure
     * Please note the list contains only valid mapping numbers. Unmapped atoms
     * which have been matched as well are ignored
     * If the structures did not match SSS wise null will be returned
     * @param query
     * @param target
     * @return null if the simple SSS failed
     */
    private static MapsList findMatchingMaps(StereoMolecule query, StereoMolecule target)
    {
        MapsList ret = new MapsList();
        // First Performance check: Don't SSS for bigger query mols
        if (query.getAllAtoms() <= target.getAllAtoms()) {
            // Get the list of map numbers of the query
            int[] queryMaps = getMapList(query);
            boolean found = false;
            for (int i : queryMaps) {
                if (i != 0) {
                    found = true;
                    break;
                }
            }
            if (found) {
                // Leave if query has no maps
                SSSearcher searcher = new SSSearcher();
                boolean fragment = query.isFragment();
//                debug("Query check on target \n%s\n%s\n",query.getIDCode(),target.getIDCode());
                query.setFragment(true);
                searcher.setMol(query, target);
                int count = 0;

                // SSS first
                if ((count = searcher.findFragmentInMolecule()) > 0) {
                    // so we found the query in the target {count} times
                    // Get list of the matched indizes of the target molecule
                    List<int[]> sssMatchList = searcher.getMatchList();
                    for (int i = 0; i < count; i++) {
                        int[] mapList = new int[0];
                        // these are the indizes of the atoms found in the target
                        int[] matchedSSSAtoms = sssMatchList.get(i);
                        // Query Atom[n] matches matchedSSSAtoms[n] => Target Atom Index

                        // queryMaps[0] = 5 ; means mapping number of query atom 0 is 5
                        // matchedSSSAtoms[0] = 4 ; means query Atom 0 matched on Target Atom 4
                        // targetMaps[0] = target.getAtomMapNo(matchedSSSAtoms[0]) => x;

                        // Get the corresponding mapping numbers of these matched target atoms
                        int[] targetMaps = getMapList(target, matchedSSSAtoms);

                        int index = mapList.length;
                        // Make room for more mapping numbers in the maplist
                        mapList = copyOf(mapList, mapList.length + matchedSSSAtoms.length);
                        // And append the mapping numbers
                        debug2("Query Map Arr\t: %s\n", ArrayUtils.toString(queryMaps));
                        debug2("Target Map Arr\t: %s\n", ArrayUtils.toString(targetMaps));

                        for (int k = 0; k < matchedSSSAtoms.length; k++) {
//                            System.out.printf("targetMap[%d] = %d queryMap[%d] = %d\n",k,targetMaps[k],k,queryMaps[k]);
                            if (targetMaps[k] != 0 && queryMaps[k] != 0) {
                                mapList[index++] = targetMaps[k];
                            }
                        }
                        // Remove the unmapped entries
                        mapList = removeZeros(mapList);
                        debug2("Query Map List\t: %s\n", ArrayUtils.toString(removeZeros(queryMaps)));
                        debug2("Target Map List\t: %s\n", ArrayUtils.toString(mapList));
                        ret.add(mapList);
                    }
//                    debug("Matched!\n");
                } else {
                    ret = null; // signal not found!
//                    debug("Did not match!\n");
                }
                query.setFragment(fragment);
            }
        }
        return ret;
    }



    public static boolean matchKeys(byte[] tK,byte[] qK)
    {
        if (qK == null || tK == null || qK.length != tK.length)
            return false;
        for (int i = 0; i < qK.length; i++) {
            if (qK[i] > tK[i])
                return false;
        }
        return true;
    }


    private static int getMaxSize(List<int[]> foo)
    {
        int size = 0;
        for (int[] k : foo) {
            size = Math.max(size,k.length);
        }
        return size;
    }



    private static List<StereoMolecule> getProducts(Reaction r)
    {
        List<StereoMolecule> list = new ArrayList<StereoMolecule>();
        for (int i = 0; i < r.getProducts(); i++) {
            list.add(r.getProduct(i));
        }
        return list;
    }

    private static List<StereoMolecule> getReactants(Reaction r)
    {
        List<StereoMolecule> list = new ArrayList<StereoMolecule>();
        for (int i = 0; i < r.getReactants(); i++) {
            list.add(r.getReactant(i));
        }
        return list;
    }

    /*
        Before adding
    +++++++++++++++++++++++++++++++++++
      sourceList             listToAdd
      ___                    ___   ___
      |3|                    |7|   |4|
      |2|                    |8|   |3|
      |1|                    |9|   |4|
      ___                    |4|   |5|
                             |5|   |2|
                             ___   ___

       After adding
    +++++++++++++++++++++++++++++++++++
      sourceList             listToAdd
      ___   ___
      |3    |3|
      |2|   |2|
      |1|   |1|
      |7|   |4|
      |8|   |3|
      |9|   |4|
      |4|   |5|
      |5|   |2|
      ___   ___
     */

    /**
     * Append target list(s) to source list(s)
     * if # of targets > 1 the source lists needs to be cloned n-1 times
     * so effectively do the cross product
     * Maybe we need a simpler solution
     * However the final check for equality of reactant and product mapping is simple then
     * @param sourceList
     * @param listToAdd
     * @return
     */
    private static List<int[]> add(List<int[]> sourceList, List<int[]> listToAdd)
    {
        int sizeofListToAdd = listToAdd.size();
        int originalSourceListSize = sourceList.size();
        if (sizeofListToAdd > 1) {
            int sourceListSize = sourceList.size();
            if (sourceListSize == 0) {
                for (int[] t : listToAdd) {
                    sourceList.add(t);
                }
            } else {
                // OK there are multiple lists to add
                // clone the source list n-1 times
                for (int i = 1; i < sizeofListToAdd; i++) {
                    for (int j = 0; j < sourceListSize; j++) {
                        int[] s = sourceList.get(j);
                        sourceList.add(s.clone());
                    }
                }
                // Add the elements of each list to add at the end of the sourcelist
                for (int i = 0; i < sizeofListToAdd; i++) {
                    int[] t = listToAdd.get(i);
                    for (int j = 0; j < originalSourceListSize; j++) {
                        int index = i * originalSourceListSize + j;
                        int[] s = sourceList.get(index);
                        // Create a new array from the current array to hold t.length more elements
                        // and copy the s.length no of elements into it
                        int[] q = copyOf(s, s.length + t.length);
                        // Append the t array at the end
                        System.arraycopy(t, 0, q, s.length, t.length);
                        // Replace the original element in the list
                        sourceList.set(index, q);
                    }
                }
            }
        } else if (sizeofListToAdd == 1) {
            if (sourceList.size() == 0) {
                sourceList.add(listToAdd.get(0));
            } else {
                int[] t = listToAdd.get(0);
                for (int i = 0; i < sourceList.size(); i++) {
                    int[] s = sourceList.get(i);
                    int[] q = copyOf(s, s.length + t.length);
                    System.arraycopy(t, 0, q, s.length, t.length);
                    sourceList.set(i, q);
                }
            }
        }
        return sourceList;
    }

    private static int[] removeZeros(int[] array)
    {
        int count = 0;
        int[] t = new int[array.length];
        for (int i = 0; i < array.length; i++) {
            if (array[i] != 0) {
                t[count++] = array[i];
            }
        }
        return copyOf(t,count);
    }

//    private static int[] sortAndRemoveZeros(int[] array)
//    {
//        Arrays.sort(array);
//        int index = 0;
//        for (int i = 0; i < array.length; i++) {
//            if (array[i] == 0) {
//                index++;
//            }
//        }
//        return copyOfRange(array, index, array.length);
//    }


    /**
     * Algorithm:
     * (This is not 100% correct, but for now lets go with it)
     * If query molecules are SSS of target molecules
     * Find the (combined) mapping numbers for each side of the target reaction by using the SSS match list from each Q / T comparison
     * If the mapping numbers on both sides of the target are equal then we have a match
     *
     * let ML={}
     * for each query reactant:
     *  if query reactant matches SSS in target reactant
     *      let QL = matching atoms in target (matchlist)
     *      let AM = list of Atom Maps of QL
     *      let ML += AM
     *  end if
     * end for
     *
     * let MP={}
     * for each query product
     *  if query product matches SSS in target product
     *      let QL = matching atoms in target (matchlist)
     *      let AM = list of Atom Maps of QL
     *      let MP += AM
     *  end if
     * end for
     *
     * ML = eliminate 0 map nos from ML and sort
     * MP = eliminate 0 map nos from MP and sort
     *
     * if (ML == MP)
     *   => MATCH
     * else
     *   -> NO MATCH
     */


    static int[] getMapList(StereoMolecule m, int atoms[])
    {
        int[] ret = new int[atoms.length];
        for (int i = 0; i < atoms.length; i++) {
            ret[i] = m.getAtomMapNo(atoms[i]);
        }
        return ret;
    }

    /**
     * Returns an array of mapping number for this molecule,
     * the index into this array corresponds to the atom index in the molecule
     * Note: Unmapped mapping numbers are included
     *
     * @param m
     * @return
     */
    static int[] getMapList(StereoMolecule m)
    {
        int atoms = m.getAllAtoms();
        int[] ret = new int[atoms];
        for (int i = 0; i < atoms; i++) {
            ret[i] = m.getAtomMapNo(i);
        }
        return ret;
    }


//    static int getAtomByMap(StereoMolecule m, int mapNo)
//    {
//        int atoms = m.getAllAtoms();
//        for (int i = 0; i < atoms; i++) {
//            if (m.getAtomMapNo(i) == mapNo) {
//                return i;
//            }
//        }
//        return -1;
//    }



    private static class ListPermutator<T> //implements Iterator?
    {

        int total;
        int index = 0;
        int count = 0;
        List<T> list;

        ListPermutator(List<T> list)
        {
            this.list = new ArrayList<T>(list);
            total = fac(list.size());
        }

        boolean hasNext()
        {
            return count < total;
        }

        List<T> next()
        {
            if (count == 0) {
                count++;
                return list;
            }
            permute(list, index);
            count++;
            index = (index + 1) % (list.size() - 1);
            return list;
        }

        private void permute(List<T> arr, int k)
        {
            java.util.Collections.swap(arr, k, k + 1);
        }

        private int fac(int c)
        {
            if (c <= 1) {
                return 1;
            }
            return c * fac(c - 1);
        }


    }


    /** We need these since we need 1.5 compliance for ORACLE db i.e Cartridge */
    private static int[] copyOf(int[] original, int newLength) {
        int[] copy = new int[newLength];
        System.arraycopy(original, 0, copy, 0,
                         Math.min(original.length, newLength));
        return copy;
    }


    private static int[] copyOfRange(int[] original, int from, int to) {
        int newLength = to - from;
        if (newLength < 0)
            throw new IllegalArgumentException(from + " > " + to);
        int[] copy = new int[newLength];
        System.arraycopy(original, from, copy, 0,
                         Math.min(original.length - from, newLength));
        return copy;
    }


}
