package com.actelion.research.calc;

import com.actelion.research.util.DoubleVec;

public class VectorSimilarity {


     public static double getTanimotoSimilarity(double[] d1, double[] d2) {

        double sum = 0;
        double dAtB = mult(d1, d2);
        double dAtA = mult(d1, d1);
        double dBtB = mult(d2, d2);

        sum = dAtB / (dAtA + dBtB - dAtB);

        return sum;
     }

    private static double mult(double [] arr1, double [] arr2) {

        double sum = 0.0;
        for (int i = 0; i < arr1.length; i++) {
            sum += arr1[i] * arr2[i];
        }

        return sum;
     }

     public static double getTanimotoSimilarity(float[] d1, float[] d2) {

        double sum = 0;
        double dAtB = mult(d1, d2);
        double dAtA = mult(d1, d1);
        double dBtB = mult(d2, d2);

        sum = dAtB / (dAtA + dBtB - dAtB);

        return sum;
    }

    private static double mult(float [] arr1, float [] arr2) {

        double sum = 0.0;
        for (int i = 0; i < arr1.length; i++) {
            sum += arr1[i] * arr2[i];
        }

        return sum;
    }

     public static double getTanimotoSimilarity(int[] d1, int[] d2) {

        double sum = 0;
        double dAtB = mult(d1, d2);
        double dAtA = mult(d1, d1);
        double dBtB = mult(d2, d2);

        sum = dAtB / (dAtA + dBtB - dAtB);

        return sum;
    }

    public static double getCosine(int[] d1, int[] d2){
        DoubleVec dB = new DoubleVec(d1);
        DoubleVec dC = new DoubleVec(d2);
        dB.norm2One();
        dC.norm2One();
        return DoubleVec.getCosine(dB, dC);
    }

    private static double mult(int [] arr1, int [] arr2) {

        double sum = 0.0;
        for (int i = 0; i < arr1.length; i++) {
            sum += arr1[i] * arr2[i];
        }

        return sum;
    }

    public static double getMinMaxSimilarity(int[] d1, int[] d2){

         double sumMin = 0;
         double sumMax = 0;

        for (int i = 0; i < d1.length; i++) {
            int v1 = d1[i];
            int v2 = d2[i];
            sumMin += Math.min(v1, v2);
            sumMax += Math.max(v1, v2);
        }

        double score = sumMin / sumMax;

        return score;

    }


}
