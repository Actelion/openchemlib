package com.actelion.research.chem.descriptor;

import com.actelion.research.calc.Logarithm;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.datamodel.IntVec;

/**
 * DescriptorHandlerBinarySkelSpheres
 *
 * This descriptor is a binary representation of the Skeleton Spheres descriptor. The correlation with the
 * SkeletonSpheres descriptor was calculated with R²=0.9891. Basis for the correlation were one million similartiy values
 * calculated from 1000 RND samples from
 * /home/korffmo1/Projects/Software/Development/VirtualScreening/data/CalibrationDataset/gpcr_ligandsAllDescriptors.dwar
 * the fraction of 0.25 of the highest similarity values was used to calculate the Pearson correlation coefficient.
 *
 *
 * <p>Copyright: Actelion Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 09.05.17.
 *
 */
public class DescriptorHandlerBinarySkelSpheres extends AbstractDescriptorHandlerFP<StereoMolecule> {



    public static final String VERSION = DescriptorConstants.DESCRIPTOR_BINARY_SKELETONSPHERES.version;


    //
    // Parameter for standardization.
    //
    // Widens the distribution
    private static final double FAC_WIDTH = 8;

    // Moves the distribution along the axis.
    private static final double OFFSET_MOVE = 0.5;

    // Centers the distribution around 0.
    private static final double OFFSET_CENTER = 0.2;

    private static DescriptorHandlerBinarySkelSpheres DEFAULT_INSTANCE;


    private DescriptorHandlerSkeletonSpheres dhSkeletonSpheres;

    public DescriptorHandlerBinarySkelSpheres() {

        dhSkeletonSpheres = new DescriptorHandlerSkeletonSpheres();
    }



    @Override
    public DescriptorInfo getInfo() {
        return DescriptorConstants.DESCRIPTOR_BINARY_SKELETONSPHERES;
    }

    @Override
    public String getVersion() {
        return VERSION;
    }

    @Override
    public int [] createDescriptor(StereoMolecule mol) {

        byte [] arrSkelSpheres = dhSkeletonSpheres.createDescriptor(mol);

        return createDescriptorFromSkelSpheresArrayCompressed(arrSkelSpheres);
    }


    /**
     * Calculates a binary vector with half number of bits as fields in the original SkeletonSpheres descriptor.
     * Th length is 512 bits.
     * Two fields of the original SkeletonSpheres descriptor are summarized.
     *
     * @param arrSkelSpheres
     * @return
     */
    public static int [] createDescriptorFromSkelSpheresArrayCompressed(byte [] arrSkelSpheres) {

        final int thresh = 8;

        if ((arrSkelSpheres == null) || (arrSkelSpheres.length==0)){
            return null;
        }

        final int numInteger = arrSkelSpheres.length / (Integer.SIZE * 2);
        final int lenInBit = arrSkelSpheres.length / 2;

        IntVec iv = new IntVec(numInteger);

        for (int i = 0; i < arrSkelSpheres.length; i+=2) {

            int freq = arrSkelSpheres[i] + arrSkelSpheres[i+1];

            if(freq>0) {

                // int bitsSet = Logarithm.log2(freq)+1;
                // int bitsSet = freq;
                int bitsSet = 0;
                if(freq > thresh){
                    bitsSet = thresh+Logarithm.log2(freq-thresh)+1;
                } else {
                    bitsSet = freq;
                }

                for (int j = 0; j < bitsSet; j++) {

                    int index= (i / 2) + j;

                    // int indexCircular = getCircularIndexOffset(index, lenInBit);
                    int indexCircular = getCircularIndexOffset(index, lenInBit);

                    iv.setBit(indexCircular);
                }
            }
        }

        return iv.get();
    }

    private static int getCircularIndexOffset(int index, int lenDescriptor){

        int indexNew = index;

        if(indexNew < lenDescriptor){
            return  indexNew;
        }

        indexNew = index-lenDescriptor;

        return indexNew;
    }


    @Override
    public DescriptorHandler getThreadSafeCopy() {

        return new DescriptorHandlerBinarySkelSpheres();
    }

    /**
     * Calculates the similarity by the number of common bits devided by the total number of bts.
     * @param a1
     * @param a2
     * @return normalized similarity.
     */
    @Override
    public float getSimilarity(int [] a1, int [] a2){

        int bcAND=0;
        int bcOR=0;

        for (int i = 0; i < a1.length; i++) {

            final int v1 = a1[i];
            final int v2 = a2[i];

            bcAND += Integer.bitCount(v1 & v2);

            bcOR += Integer.bitCount(v1 | v2);

        }

//        System.out.println("bcAND " + bcAND);
//        System.out.println("bcOR " + bcOR);


        double score = (double)bcAND/((double)bcOR);
        // double score = (double)bcAND/((bcAND * 0.5) + (bcOR * 0.5));

        // return (float)score;

        return (float)correctionTS(score);
    }

    public static DescriptorHandlerBinarySkelSpheres getDefaultInstance() {

        synchronized(DescriptorHandlerBinarySkelSpheres.class) {
            if (DEFAULT_INSTANCE == null)
                DEFAULT_INSTANCE = new DescriptorHandlerBinarySkelSpheres();
        }

        return DEFAULT_INSTANCE;
    }

    /**
     * The parameter were derived from one million similarity scores with
     * com.actelion.research.chem.descriptor.util.SimilarityCalibration
     * and
     * com.actelion.research.chem.descriptor.util.CalculateStandardsationFactorFromSimilarityMatrix
     *
     * @param s
     * @return
     */
    public static double standardize(double s){

        if(s <= 0){
            return 0;
        } else if(s>=1.0){
            return 1.0;
        }

        final double sc = s-OFFSET_CENTER;

        double v =  Math.exp(-FAC_WIDTH * (sc)) - OFFSET_MOVE;

        return v;
    }

    public static double correctionTS(double s){

//        double t1 = 0.615;
//        double t2 = 0.615;
        double t1 = 0.7;
        double t2 = 0.7;

        double v = 1.0-Math.pow(1-Math.pow(s, t1) ,1.0/t2);

        if(v <= 0){
            return 0;
        } else if(v>=1.0){
            return 1.0;
        }

        return v;

    }
//    public static double correctionTS(double s){
//
//        double t1 = 0.5;
//        double t2 = 0.8;
//
//        double v = 1.0-Math.pow(1-Math.pow(s-0.2, t1) ,1.0/t2) + 0.4;
//
//        if(v <= 0){
//            return 0;
//        } else if(v>=1.0){
//            return 1.0;
//        }
//
//        return v;
//
//    }

}
