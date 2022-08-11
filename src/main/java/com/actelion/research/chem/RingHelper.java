package com.actelion.research.chem;

import com.actelion.research.chem.ExtendedMoleculeFunctions;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;

import java.util.*;

/**
 * Modest v. Korff
 * Idorsia Pharmaceuticals Ltd.
 * 07.07.2022 Start implementation
 **/
public class RingHelper {

    public static final int MAX_RING_SIZE = 100;

    private int [][] arrTopoDist;

    private RingCollection ringCollection;
    private List<int []> liRing;

    public RingHelper(StereoMolecule mol) {
        ringCollection = new RingCollection(mol,RingCollection.MODE_SMALL_AND_LARGE_RINGS, MAX_RING_SIZE);
        arrTopoDist = ExtendedMoleculeFunctions.getTopologicalDistanceMatrix(mol);
        liRing = new ArrayList<>();
    }

    public RingCollection getRingCollection() {
        return ringCollection;
    }

    public int getLargestRingSize(){
        int maxRingSize = 0;

        liRing.clear();
        int nRings = ringCollection.getSize();
        for (int i = 0; i < nRings; i++) {
            liRing.add(ringCollection.getRingAtoms(i));
        }

        Collections.sort(liRing, new Comparator<int[]>() {
            @Override
            public int compare(int[] o1, int[] o2) {
                int cmp=0;
                if(o1.length>o2.length){
                    cmp=1;
                }else if(o1.length<o2.length){
                    cmp=-1;
                }
                return cmp;
            }
        });

        Collections.reverse(liRing);

        for (int [] arrRing : liRing){

            if(!isEnclosingRing(arrRing)){
                maxRingSize = arrRing.length;
                // System.out.println(Arrays.toString(arrRing));
                break;
            }
        }

        return maxRingSize;
    }

    public boolean isEnclosingRing(int [] arrIndexRingAtomsRing2Check){

        boolean enclosingRing = false;

        int offset=1;

        int index2Width = arrIndexRingAtomsRing2Check.length / 2 + 1;

        leave:
        for (int i = 0; i < index2Width; i++) {

            int indAt1 = arrIndexRingAtomsRing2Check[i];

            for (int j = offset; j < index2Width; j++) {

                int index2 = i+j;

                if(index2 >= arrIndexRingAtomsRing2Check.length){
                    index2 = index2-arrIndexRingAtomsRing2Check.length;
                }

                int topoDistRingPath = index2-i;
                if(topoDistRingPath > index2Width){
                    topoDistRingPath = i-index2-arrIndexRingAtomsRing2Check.length;
                }

                int indAt2 = arrIndexRingAtomsRing2Check[index2];

                int minTopoDist = arrTopoDist[indAt1][indAt2];

                if(minTopoDist<topoDistRingPath){
                    enclosingRing = true;
                    break leave;
                }
            }
        }

        return enclosingRing;
    }

}
