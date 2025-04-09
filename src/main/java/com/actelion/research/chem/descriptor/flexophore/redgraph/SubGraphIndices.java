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
 * @author Modest v. Korff
 */

package com.actelion.research.chem.descriptor.flexophore.redgraph;

import com.actelion.research.chem.ExtendedMoleculeFunctions;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.hash.HashSetInt;

import java.util.*;

/**
 * SubGraphIndices
 *
 * A class to handle indices in a molecular graph.
 *
 * It is just a hash set for integer with some methods.
 *
 * Created by korffmo1 on 23.02.16.
 */
public class SubGraphIndices {

    private HashSetInt hsIndexAtom;


    public SubGraphIndices() {
        hsIndexAtom = new HashSetInt();
    }

    public SubGraphIndices(int[] arrIndexAtom) {
        hsIndexAtom = new HashSetInt(arrIndexAtom);
    }

    public int getNumIndices() {
        return hsIndexAtom.size();
    }

    public void clear() {
        hsIndexAtom.clear();
    }

    public void addIndex(int indexAtom) {
        hsIndexAtom.add(indexAtom);
    }

    public void addIndex(int[] arrIndexAtom) {
        hsIndexAtom.add(arrIndexAtom);
    }

    public int[] getAtomIndices() {
        return hsIndexAtom.getValues();
    }

    public boolean contains(int indexAtom) {
        return hsIndexAtom.contains(indexAtom);
    }
    public boolean equalIndices(int [] arrIndexAtom) {
        if(hsIndexAtom.size()!=arrIndexAtom.length)
            return false;

        boolean eq=true;

        for (int index : arrIndexAtom) {
            if(!hsIndexAtom.contains(index)){
                eq=false;
                break;
            }
        }

        return eq;
    }

    @Override
    public boolean equals(Object obj) {
        if(!(obj instanceof SubGraphIndices)){
            return false;
        }

        SubGraphIndices sg = (SubGraphIndices) obj;
        if(hsIndexAtom.size() != sg.getNumIndices()){
            return false;
        }

        int [] a0 = hsIndexAtom.getValues();
        int [] a1 = sg.getAtomIndices();
        Arrays.sort(a0);
        Arrays.sort(a1);
        boolean eq = true;
        for (int i = 0; i < a0.length; i++) {
            if(a0[i] != a1[i]){
                eq = false;

            }
        }

        return eq;
    }

    public boolean isOverlap(SubGraphIndices frag) {

        boolean overlap = false;
        int[] arrFragIndexAt = frag.hsIndexAtom.getValues();
        for (int indexAtomFrag : arrFragIndexAt) {
            if (hsIndexAtom.contains(indexAtomFrag)) {
                overlap = true;
                break;
            }
        }
        return overlap;

    }

    public int getNumOverlappingIndices(SubGraphIndices frag) {

        int nOverlapping = 0;

        int[] arrFragIndexAt = frag.hsIndexAtom.getValues();
        for (int indexAtomFrag : arrFragIndexAt) {
            if (hsIndexAtom.contains(indexAtomFrag)) {
                nOverlapping++;
            }
        }
        return nOverlapping;
    }

    public void merge(SubGraphIndices frag) {
        int[] arrFragIndexAt = frag.hsIndexAtom.getValues();
        for (int indexAtomFrag : arrFragIndexAt) {
            hsIndexAtom.add(indexAtomFrag);
        }
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        int [] a = hsIndexAtom.getValues();
        Arrays.sort(a);
        for (int i = 0; i < a.length; i++) {
            sb.append(a[i]);
            if(i<a.length-1){
                sb.append(",");
            }
        }

        return sb.toString();
    }


    public static boolean isOnlyCarbon(StereoMolecule mol, SubGraphIndices sgi){
        boolean carbon=true;
        for (int atomIndex : sgi.getAtomIndices()) {
            if(mol.getAtomicNo(atomIndex)!=6){
                carbon=false;
                break;
            }
        }
        return carbon;
    }
    public static boolean isCharged(StereoMolecule mol, SubGraphIndices sgi){
        int charge=0;
        // A nitro group compensates to 0.
        for (int atomIndex : sgi.getAtomIndices()) {
            charge+=mol.getAtomCharge(atomIndex);
        }
        return (charge!=0);
    }

    /**
     * Merges fragments containing a common atom index. The fragments in the list are merged and the list is
     * shrinked.
     *
     * @param liFragment
     */
    public static void merge(List<SubGraphIndices> liFragment) {
        boolean merged = true;
        while (merged) {
            merged = false;
            for (int i = 0; i < liFragment.size(); i++) {
                SubGraphIndices frag1 = liFragment.get(i);
                for (int j = liFragment.size() - 1; j > i; j--) {
                    SubGraphIndices frag2 = liFragment.get(j);
                    if (frag1.isOverlap(frag2)) {
                        frag1.merge(frag2);
                        liFragment.remove(j);
                        merged = true;
                    }
                }
            }
        }
    }

    /**
     * Merges overlapping fragments. They may share a minimum number of overlapping indices.
     * @param liFragment
     * @param minNumIndicesOverlapping
     * @return
     */
    public static List<SubGraphIndices> mergeOverlapping(List<SubGraphIndices> liFragment, int minNumIndicesOverlapping) {

        List<OverlappingFragments> liOverlappingFragments = new ArrayList<>();
        for (int i = 0; i < liFragment.size(); i++) {
            SubGraphIndices frag1 = liFragment.get(i);
            OverlappingFragments overlappingFragments = new OverlappingFragments(frag1);
            liOverlappingFragments.add(overlappingFragments);
        }

        boolean merged = true;
        while (merged) {
            merged = false;
            for (int i = 0; i < liOverlappingFragments.size(); i++) {
                OverlappingFragments of1 = liOverlappingFragments.get(i);
                for (int j = liOverlappingFragments.size() - 1; j > i; j--) {
                    OverlappingFragments of2 = liOverlappingFragments.get(j);
                    if (of1.getNumLargestOverlap(of2) >= minNumIndicesOverlapping) {
                        of1.add(of2);
                        liOverlappingFragments.remove(j);
                    }
                }
            }
        }

        List<SubGraphIndices> liFragmentMerged = new ArrayList<>();
        for (OverlappingFragments olf : liOverlappingFragments) {
            liFragmentMerged.add(olf.getSubGraphIndices());
        }

        return liFragmentMerged;
    }

    /**
     * Adds the atom indices to the hash set.
     * @param hs
     * @param liFragment
     */
    protected static void addAtomIndices(HashSetInt hs, List<SubGraphIndices> liFragment) {
        for (SubGraphIndices fragment : liFragment) {
            hs.add(fragment.getAtomIndices());
        }
    }

    /**
     *
     * @param mol
     * @param liSGI
     * @param indexSGI index for the {@link SubGraphIndices} under consideration.
     * @return true if the {@link SubGraphIndices} under consideration connects two or more other {@link SubGraphIndices}.
     */
    public static boolean isLinker(StereoMolecule mol, List<SubGraphIndices> liSGI, int indexSGI){
        int atoms = mol.getAtoms();

        SubGraphIndices sgi = liSGI.get(indexSGI);

        boolean [] arrAtomIndicesUsedMap = new boolean[atoms];

        for (int indexAtmUsed : sgi.getAtomIndices()) {
            arrAtomIndicesUsedMap[indexAtmUsed] = true;
        }

        boolean [] mapSGI = new boolean[liSGI.size()];

        for (int indAtmStart : sgi.getAtomIndices()) {
            LinkedList<Integer> liIndAtm = new LinkedList<>();
            liIndAtm.add(indAtmStart);

            out:
            while (!liIndAtm.isEmpty()) {
                int indAtm = liIndAtm.poll();
                int nConnAtms = mol.getConnAtoms(indAtm);
                for (int i = 0; i < nConnAtms; i++) {
                    int indAtmConn = mol.getConnAtom(indAtm, i);
                    if (arrAtomIndicesUsedMap[indAtmConn]) {
                        continue;
                    }
                    arrAtomIndicesUsedMap[indAtmConn] = true;
                    // Is in neighbour sgi?
                    for (int j = 0; j < liSGI.size(); j++) {
                        if(indexSGI==j)
                            continue;
                        SubGraphIndices sgi2Check = liSGI.get(j);
                        if(sgi2Check.contains(indAtmConn)){
                            mapSGI[j]=true;
                            break out;
                        }
                    }
                    liIndAtm.add(indAtmConn);
                }
            }
        }

        int ccLinked=0;
        for (boolean b : mapSGI) {
            if(b)ccLinked++;
        }
        return (ccLinked>1);
    }


    public static Comparator<SubGraphIndices> getComparatorNumIndices() {

        return new Comparator<SubGraphIndices>() {
            @Override
            public int compare(SubGraphIndices f1, SubGraphIndices f2) {
                int cmp = 0;
                if (f1.hsIndexAtom.size() > f2.hsIndexAtom.size()) {
                    cmp = 1;
                } else if (f1.hsIndexAtom.size() < f2.hsIndexAtom.size()) {
                    cmp = -1;
                }
                return cmp;
            }
        };
    }

    private static class OverlappingFragments {
        List<SubGraphIndices> liSubGraphIndices;
        public OverlappingFragments(SubGraphIndices sg) {
            liSubGraphIndices = new ArrayList<>();
            liSubGraphIndices.add(sg);
        }

       public OverlappingFragments(SubGraphIndices sg1, SubGraphIndices sg2) {
            liSubGraphIndices = new ArrayList<>();
            liSubGraphIndices.add(sg1);
            liSubGraphIndices.add(sg2);
        }


        public void add(SubGraphIndices sg){
            liSubGraphIndices.add(sg);
        }

        public void add(OverlappingFragments of){
            liSubGraphIndices.addAll(of.liSubGraphIndices);
        }


        public int getNumLargestOverlap(OverlappingFragments of){
            int nOverlapMax = 0;
            for (SubGraphIndices sgi0 : liSubGraphIndices) {
                for (SubGraphIndices sgi1 : of.liSubGraphIndices) {
                    int n = sgi0.getNumOverlappingIndices(sgi1);
                    if(n > nOverlapMax){
                        nOverlapMax = n;
                    }
                }
            }
            return nOverlapMax;
        }

        public boolean containsSubGraph(OverlappingFragments of) {
            boolean contains = false;
            containsB:
            for (SubGraphIndices sgi1 : liSubGraphIndices) {
                for (SubGraphIndices sgi2 : of.liSubGraphIndices) {
                    if(sgi1.equals(sgi2)){
                        contains = true;
                        break containsB;
                    }
                }
            }
            return contains;
        }

        public SubGraphIndices getSubGraphIndices(){
            SubGraphIndices.merge(liSubGraphIndices);
            if(liSubGraphIndices.size() != 1) {
                throw new RuntimeException("This should not happen.");
            }
            return liSubGraphIndices.get(0);
        }
    }
}
