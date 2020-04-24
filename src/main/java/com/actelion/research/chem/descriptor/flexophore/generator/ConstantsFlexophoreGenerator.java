package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.chem.descriptor.flexophore.PPNode;

public class ConstantsFlexophoreGenerator {


    /**
     * Defines the resolution for the range.
     */
    public static final int BINS_HISTOGRAM = 80;
    /**
     * Range histogram in Angstrom.
     */
    public static final int RANGE_HISTOGRAM = 40;

    public static final int MAX_VAL_INTERACTION_TYPE = (int)Math.pow(2, Byte.SIZE * PPNode.NUM_BYTES_INTERACTION_TYPE)-1;


    public static final boolean OPTIMIZE_RIGID_FRAGS = false;


    public static double getResolution(){
        return RANGE_HISTOGRAM / (double) BINS_HISTOGRAM;
    }
}
