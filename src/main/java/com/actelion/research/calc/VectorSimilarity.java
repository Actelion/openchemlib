package com.actelion.research.calc;

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

    private static double mult(int [] arr1, int [] arr2) {

        double sum = 0.0;
        for (int i = 0; i < arr1.length; i++) {
            sum += arr1[i] * arr2[i];
        }

        return sum;
    }


}
