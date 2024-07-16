package com.actelion.research.chem.descriptor.flexophore;

public class ConstantsFlexophore {

    public static final byte MODE_SOFT_PPPOINTS = 0;
    public static final byte MODE_HARD_PPPOINTS = 1;

    public static final int FLAG_CENTER_ATOM = 1<<4;

    public static final int MAX_NUM_NODES_FLEXOPHORE = 64;


    public static final String TAG_FLEXOPHORE_OBJECT =  "FlexDecoded";



    public static final int LABEL_LOW = 0;
    public static final int LABEL_NORMAL = 1;
    public static final int LABEL_MANDATORY = 2;
    public static final double VAL_WEIGHT_LOW = 0.5;
    public static final double VAL_WEIGHT_NORMAL = 1;
    public static final double VAL_WEIGHT_MANDATORY = 2;
}
