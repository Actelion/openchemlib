package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;
import com.actelion.research.util.StringFunctions;

/**
 * SlidingWindowDistHist
 * Created by korffmo1 on 01.03.16.
 */
public class SlidingWindowDistHist {


    private double [] arrFilter;
    private byte [] arrTmp;

    private int lenFilHalf;

    public SlidingWindowDistHist(double[] arrFilter) {
        this.arrFilter = arrFilter;
        if(arrFilter.length % 2 == 0){
            throw new RuntimeException("Odd number of filter values needed!");
        }
        lenFilHalf = arrFilter.length / 2;
        int n = ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE;
        int maxNumDistHist = ((n*n)-n)/2;
        arrTmp = new byte[ConstantsFlexophoreGenerator.BINS_HISTOGRAM * maxNumDistHist];
    }

    public void apply(DistHist distHist){
        int nNodes = distHist.getNumPPNodes();
        int n = ConstantsFlexophoreGenerator.BINS_HISTOGRAM;
        int end = n - lenFilHalf;

        // System.out.println(StringFunctions.toStringShort(distHist.getDistHists()));

        for (int i = 0; i < nNodes; i++) {
            for (int j = i+1; j < nNodes; j++) {
                int indexStart = distHist.getIndexPosStartForDistHist(i,j);
                for (int k = lenFilHalf; k < end; k++) {
                    double v = 0;
                    for (int l = 0; l < arrFilter.length; l++) {
                        int ind = k - lenFilHalf + l + indexStart;
                        v +=  distHist.getValueAtAbsolutePosition(ind) * arrFilter[l];
                    }
                    arrTmp[indexStart+k] = (byte)(v+0.5);
                }
            }
        }


        System.arraycopy(arrTmp, 0, distHist.arrDistHists, 0, distHist.arrDistHists.length);

//        System.out.println("\n\n");
//        System.out.println(StringFunctions.toStringShort(distHist.arrDistHists));
    }

    public byte [] getBlurred(byte [] hist){

        byte [] histBlurred = new byte[ConstantsFlexophoreGenerator.BINS_HISTOGRAM];

        int n = ConstantsFlexophoreGenerator.BINS_HISTOGRAM;
        int end = n - lenFilHalf;

        for (int k = lenFilHalf; k < end; k++) {
            double v = 0;
            for (int l = 0; l < arrFilter.length; l++) {
                int ind = k - lenFilHalf + l;
                v +=  hist[ind] * arrFilter[l];
            }
            histBlurred[k] = (byte)(v+0.5);
        }
        return histBlurred;
    }
}
