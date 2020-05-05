package com.actelion.research.chem.descriptor.flexophore.redgraph;

import com.actelion.research.util.hash.HashSetInt;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * SubGraphIndices
 *
 * A class to handle indices in a molecular graph.
 *
 * It is just a hash set for integer with some methods.
 *
 * <p>Copyright: Actelion Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 23.02.16.
 */
public class SubGraphIndices {

    private HashSetInt hsIndexAtom;

    public SubGraphIndices() {
        hsIndexAtom = new HashSetInt();
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

    protected static void addAtomIndices(HashSetInt hs, List<SubGraphIndices> liFragment) {

        for (SubGraphIndices fragment : liFragment) {

            hs.add(fragment.getAtomIndices());

        }
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
