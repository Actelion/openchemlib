package com.actelion.research.calc.distance;

/**
 * DistanceMetrics
 * <p>Copyright: Idorsia Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Idorsia Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 13.09.17.
 */
public class DistanceMetrics {


    public static final double getJaccard(int [] a1, int [] a2) {

        double sumMax = 0;
        double sumMin = 0;

        for (int i = 0; i < a1.length; i++) {
            sumMax += Math.max(a1[i], a2[i]);
            sumMin += Math.min(a1[i], a2[i]);
        }

        return sumMin / sumMax;
    }

    public static double getCosine(double[] vectorA, double[] vectorB) {
        double dotProduct = 0.0;
        double normA = 0.0;
        double normB = 0.0;
        for (int i = 0; i < vectorA.length; i++) {
            dotProduct += vectorA[i] * vectorB[i];
            normA += vectorA[i] * vectorA[i];
            normB += vectorB[i] * vectorB[i];
        }
        return dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
    }
}
