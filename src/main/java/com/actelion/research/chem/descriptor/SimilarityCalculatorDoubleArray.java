package com.actelion.research.chem.descriptor;

import com.actelion.research.calc.distance.DistanceMetrics;

/**
 * SimilarityCalculatorDoubleArray
 * Created by korffmo1 on 12.07.18.
 */
public class SimilarityCalculatorDoubleArray implements ISimilarityCalculator<double []>{

    public static final String NAME = "SimilarityCalculatorDoubleArray";
    public static final String SHORT_NAME = "SimCalcDblArray";


    @Override
    public float getSimilarity(double[] d1, double[] d2) {

        if (d1 == null || d2 == null)
            return Float.NaN;

        return (float) DistanceMetrics.getCosine(d1, d2);

    }

    @Override
    public SimilarityCalculatorInfo getInfo() {
        return new SimilarityCalculatorInfo(NAME, SHORT_NAME);
    }

    @Override
    public ISimilarityCalculator<double[]> getThreadSafeCopy() {
        return new SimilarityCalculatorDoubleArray();
    }
}
