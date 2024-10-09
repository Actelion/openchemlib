package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.chem.descriptor.flexophore.PPNode;

public class ConstantsFlexophoreGenerator {

    public static final double SUM_VAL_HIST = 100.0;
    /**
     * Defines the resolution for the range. This is also the length for the DistHist vector.
     */
    public static final int BINS_HISTOGRAM = 80;

    public static final int INTERACTION_TYPE_NONE = -1;


    /**
     * Range histogram in Angstrom. As long as the molecule coordinates are in Angstrom.
     */
    public static final int RANGE_HISTOGRAM = 40;

    public static final int MAX_VAL_INTERACTION_TYPE = (int)Math.pow(2, Byte.SIZE * PPNode.NUM_BYTES_INTERACTION_TYPE)-1;


    public static final boolean OPTIMIZE_RIGID_FRAGS = false;


    /**
     * Filter for sliding window to blurr distance histograms
     * 07.04.2020
     */
    // public static final double [] FILTER = FILTER05;
    public static final double [] FILTER07_ = {0.125, 0.125, 0.125,0.25,0.125, 0.125, 0.125};
    public static final double [] FILTER05_ = {0.125, 0.25,0.25,0.25, 0.125};

    public static final double [] FILTER = FILTER07_;

    public static double getResolution(){
        return RANGE_HISTOGRAM / (double) BINS_HISTOGRAM;
    }
}
