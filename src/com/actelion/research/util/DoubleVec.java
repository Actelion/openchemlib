package com.actelion.research.util;

import java.text.DecimalFormat;
import java.util.Random;
import java.util.Vector;

/**
 * <p>Title: </p>
 * <p>Description: Vector with double values </p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: Axovan Ltd. </p>
 * @author Modest von Korff
 * @version 1.0
     * 18.09.2003 EuclideanDistanceFast runtime exception added cause of error in SOM
 */
// Please rename this class
public class DoubleVec implements Comparable<DoubleVec> {
	
	private static final DecimalFormat NF = new DecimalFormat("0.000");
	
    private static final Random RAND  = new Random();

    public static final int COSINE = 1;
    public static final int EUCLIDEAN = 2;
    public static final int EUCLIDEAN_FAST = 3;
    public static final int TANIMOTO = 4;
    public static final int TANIMOTO_INV = 5;

    private double data[];

    // Dot product for Tanimoto
    private double mDotProd;

    public DoubleVec(DoubleVec dVec) {
        init();
        data = new double[dVec.data.length];
        for (int ii = 0; ii < data.length; ii++) {
            data[ii] = dVec.data[ii];
        }
    }

    public DoubleVec(DoubleVec dVec, boolean bDotProd) {
        init();
        data = new double[dVec.data.length];
        for (int ii = 0; ii < data.length; ii++) {
            data[ii] = dVec.data[ii];
        }
        if(bDotProd)
            mDotProd = mult(data, data);
    }

    public DoubleVec(int iSize) {
        init();
        data = new double[iSize];
    }

    public DoubleVec(double[] dVec) {
        init();
        data = new double[dVec.length];
        for (int ii = 0; ii < data.length; ii++) {
            data[ii] = dVec[ii];
        }
    }

    public DoubleVec(int[] dVec) {
        init();
        data = new double[dVec.length];
        for (int ii = 0; ii < data.length; ii++) {
            data[ii] = dVec[ii];
        }
    }

    public DoubleVec(double[] dVec, boolean bDotProd) {
        init();
        data = new double[dVec.length];
        for (int ii = 0; ii < data.length; ii++) {
            data[ii] = dVec[ii];
        }
        if(bDotProd)
            mDotProd = mult(data, data);
    }

    public DoubleVec(Vector<Double> vec) {
        init();
        data = new double[vec.size()];
        for (int i = 0; i < data.length; i++) {
            data[i] = vec.get(i);
        }
    }

    public DoubleVec(Vector<Double> vec, boolean bDotProd) {
        init();
        
        data = new double[vec.size()];
        
        for (int i = 0; i < data.length; i++) {
            data[i] = vec.get(i);
        }
        
        if(bDotProd)
            mDotProd = mult(data, data);
    }

    /**
     * Adds aor subtracts a random value to the original value. The range is
     * specified by the percentage.
     *
     * @param dPercentage maximum range for the value change
     */
    public void addRNDvalue(double dPercentage) {
        Random rnd = new Random();

        for (int ii = 0; ii < data.length; ii++) {
            double dVal = data[ii] * (dPercentage / 100.0) * rnd.nextDouble();

            if (rnd.nextBoolean()) {
                data[ii] = data[ii] + dVal;
            }
            else {
                data[ii] = data[ii] - dVal;
            }
        }
    }

    public DoubleVec add(DoubleVec dvVec) {
        DoubleVec ret = new DoubleVec(data.length);

        if (data.length != dvVec.data.length) {
            throw new RuntimeException();
        }

        for (int ii = 0; ii < ret.data.length; ii++) {
            ret.data[ii] = data[ii] + dvVec.data[ii];
        }

        return ret;
    }

    public void addNoise(double fracNoise, double min, double max) {
        double d = max - min;
        for (int ii = 0; ii < data.length; ii++)
            if(RAND.nextDouble() < (fracNoise)) {
                double rnd = min + (RAND.nextDouble() * d);
                data[ii] = rnd;
            }
    }

    /**
     *
     * @param obj DoubleVec
     * @return -1 if the first different value is smaller than the corresponding
     * value in dv. 0 if bot vectors are equal. 1 if the first different value
     * is bigger than the corresponding value in dv.
     */
    public Object clone() {
        DoubleVec vec = new DoubleVec(this);
        return vec;
    }

    public int compareTo(DoubleVec dv) {
        
        int cmp = 0;
        for (int i = 0; i < data.length; i++) {
            if (data[i] > dv.data[i]) {
                cmp = 1;
                break;
            }
            else if (data[i] < dv.data[i]) {
                cmp = -1;
                break;
            }
        }
        return cmp;
    }

    public static double distance(DoubleVec dVec1, DoubleVec dVec2, int metric) throws
        Exception {
        double dDist = 0;
        // Vectors have to be normed.
        if (metric == DoubleVec.COSINE) {
            dDist = getCosine(dVec1, dVec2);
        }
        else if (metric == DoubleVec.EUCLIDEAN) {
            dDist = euclideanDistance(dVec1, dVec2);
        }
        else if (metric == DoubleVec.EUCLIDEAN_FAST) {
            dDist = getEuclideanDistanceFast(dVec1, dVec2);
        }
        else if (metric == DoubleVec.TANIMOTO) {
            dDist = getTanimotoDistance(dVec1, dVec2);
        }
        else if (metric == DoubleVec.TANIMOTO_INV) {
            dDist = getTanimotoDistInvDotProd(dVec1, dVec2);
        }
        else {
            throw new Exception("Unknown distance metric.");
        }

        return dDist;
    }

    public static DoubleVec devide(DoubleVec dVec1, DoubleVec dVec2) {

        DoubleVec dVecDev = new DoubleVec(dVec1.data.length);

        for (int ii = 0; ii < dVec1.data.length; ii++) {
            dVecDev.data[ii] = dVec1.data[ii] / dVec2.data[ii];
        }

        return dVecDev;
    }

    public boolean equal(DoubleVec dv) {
        boolean bEq = true;
        for (int ii = 0; ii < data.length; ii++) {
          if(data[ii] != dv.data[ii]) {
              bEq = false;
              break;
          }
        }
        return bEq;
    }

    public boolean equals(DoubleVec dv) {
        boolean bEq = true;
        for (int ii = 0; ii < data.length; ii++) {
          if(data[ii] != dv.data[ii]) {
              bEq = false;
              break;
          }
        }
        return bEq;
    }

    /**
     * Euclidean distance
     * @param dVec1
     * @param dVec2
     * @return
     */
    static public double euclideanDistance(DoubleVec dVec1, DoubleVec dVec2) {
        double dist = 0;

        double sum = 0;
        for (int i = 0; i < dVec1.data.length; i++) {
            sum += (dVec1.data[i] - dVec2.data[i]) * (dVec1.data[i] - dVec2.data[i]);
        }

        dist = Math.sqrt(sum);

        return dist;
    }
    
    static public double euclideanDistance(double [] arr1, double [] arr2) {
        double dist = 0;

        double sum = 0;
        for (int i = 0; i < arr1.length; i++) {
            sum += (arr1[i] - arr2[i]) * (arr1[i] - arr2[i]);
        }

        dist = Math.sqrt(sum);

        return dist;
    }

    static public double overlapDistance(DoubleVec dVec1, DoubleVec dVec2) {
        double dDist = 0;

        double occ1 = 0;
        double occ2 = 0;
        for (int ii = 0; ii < dVec1.data.length; ii++) {
        	if(dVec1.data[ii] != 0 )
        		occ1++;
        	if(dVec2.data[ii] != 0)
        		occ2++;
        }

        double dSum = 0;
        double occ = 0;
        for (int ii = 0; ii < dVec1.data.length; ii++) {
        	if(dVec1.data[ii] != 0 && dVec2.data[ii] != 0){
        		dSum += Math.abs(dVec1.data[ii] - dVec2.data[ii]);
        		occ++;
        	}
        	
        	
            
        }

        // dDist = Math.sqrt(dSum);
        dDist = 1.0 - (occ / Math.min(occ1,occ2));

        return dDist;
    }

    /**
     * Euclidean distance without sqrt
     * @param dVec1
     * @param dVec2
     * @return
     */
    static public double getEuclideanDistanceFast(DoubleVec dVec1, DoubleVec dVec2) throws
        ArrayIndexOutOfBoundsException {

        double dSum = 0;
        try {
			for (int ii = 0; ii < dVec1.data.length; ii+=4) {
			    dSum += (dVec1.data[ii] - dVec2.data[ii]) * (dVec1.data[ii] - dVec2.data[ii]);
			    dSum += (dVec1.data[ii+1] - dVec2.data[ii+1]) * (dVec1.data[ii+1] - dVec2.data[ii+1]);
			    dSum += (dVec1.data[ii+2] - dVec2.data[ii+2]) * (dVec1.data[ii+2] - dVec2.data[ii+2]);
			    dSum += (dVec1.data[ii+3] - dVec2.data[ii+3]) * (dVec1.data[ii+3] - dVec2.data[ii+3]);
			}
		} catch (RuntimeException e) {
			e.printStackTrace();
		}

        return dSum;
    }

    public double [] get() {
        return data;
    }

    public double get(int col) {
        return data[col];
    }

    public double getNorm() {
        double dNorm = 0;

        for (int ii = 0; ii < data.length; ii++) {
            dNorm += (data[ii] * data[ii]);
        }
        dNorm = Math.sqrt(dNorm);
        return dNorm;
    }

    /**
     * Vectors have to be normed!
     * @param dVec1 normed vector1
     * @param dVec2 normed vector2
     * @return Cosine
     */
    static public double getCosine(DoubleVec dVec1, DoubleVec dVec2) {
        double cosine = 0;
        for (int ii = 0; ii < dVec1.data.length; ii++) {
            cosine += dVec1.data[ii] * dVec2.data[ii];
        }
        return cosine;
    }

    static public double cubicDistance(DoubleVec dVec1, DoubleVec dVec2) {

        double dSum = 0;
        for (int ii = 0; ii < dVec1.data.length; ii++) {
            double dDist = Math.abs( (dVec1.data[ii] - dVec2.data[ii]));
            dSum += dDist * dDist * dDist;
        }

        return dSum;
    }

    public void initRND(double dMin, double dMax) {

        double dRange = dMax - dMin;

        for (int kk = 0; kk < data.length; kk++) {
            double dVal = dRange * Math.random() + dMin;
            data[kk] = dVal;
        }
    }

    private void init() {
        mDotProd = Double.NaN;
    }

    /**
         * The array contains the maximum and the minimum values for the initialisation
     * of each field in the double vector.
     * @param dArrMaxMin row 0: max val, row 1 min val
     */
    public void initRND(double[][] dArrMaxMin) {

        for (int kk = 0; kk < data.length; kk++) {
            // double dVal = (dRange * (new Random()).nextDouble()) + dMin;
            double dRange = dArrMaxMin[0][kk] - dArrMaxMin[1][kk];
            double dVal = dRange * Math.random() + dArrMaxMin[1][kk];
            data[kk] = dVal;
        }
    }

    static public double getManhattanBlockDistance(DoubleVec dVec1,
                                                DoubleVec dVec2) {
        double dDist = 0;

        double dSum = 0;
        for (int ii = 0; ii < dVec1.data.length; ii++) {
            dSum += Math.abs(dVec1.data[ii] - dVec2.data[ii]);
        }

        dDist = Math.sqrt(dSum);

        return dDist;
    }

    public DoubleVec mult(double dFactor) {

        DoubleVec ret = new DoubleVec(data.length);
        for (int ii = 0; ii < ret.data.length; ii++) {
            ret.data[ii] = data[ii] * dFactor;
        }

        return ret;
    }

    public static double mult(DoubleVec dVec1, DoubleVec dVec2) {

        double dSum = 0.0;

        for (int ii = 0; ii < dVec1.data.length; ii++) {
            dSum += dVec1.data[ii] * dVec2.data[ii];
        }

        return dSum;
    }

    private static double mult(double [] arr1, double [] arr2) {

        double dSum = 0.0;

        for (int ii = 0; ii < arr1.length; ii++) {
            dSum += arr1[ii] * arr2[ii];
        }

        return dSum;
    }

    /**
     * Elementwise multiplication
     * @param dVec1 input vector
     * @param dVec2 input vector
     * @return DoubleVec
     */
    public static DoubleVec multEl(DoubleVec dVec1, DoubleVec dVec2) {

        DoubleVec dVecMult = new DoubleVec(dVec1.data.length);

        for (int ii = 0; ii < dVec1.data.length; ii++) {
            dVecMult.data[ii] = dVec1.data[ii] * dVec2.data[ii];
        }

        return dVecMult;
    }

    public void norm2One() {
        double norm = getNorm();
        for (int ii = 0; ii < data.length; ii++) {
            data[ii] /= norm;
        }
    }

    public static DoubleVec minus(DoubleVec dVec1, DoubleVec dVec2) {

        DoubleVec dVecSub = new DoubleVec(dVec1.data.length);

        for (int i = 0; i < dVec1.data.length; i++) {
            dVecSub.data[i] = dVec1.data[i] - dVec2.data[i];
        }

        return dVecSub;
    }

    public void reduce(Vector<Integer> vecIndices) {
    	
        double [] arr = new double [vecIndices.size()];
        
        for (int i = 0; i < vecIndices.size(); i++) {
            int iIndex = vecIndices.get(i);
            arr[i] = data[iIndex];
        }
        
        data = arr;
    }

    public static DoubleVec plus(DoubleVec dVec1, DoubleVec dVec2) {

        DoubleVec dVecSum = new DoubleVec(dVec1.data.length);

        for (int ii = 0; ii < dVec1.data.length; ii++) {
            dVecSum.data[ii] = dVec1.data[ii] + dVec2.data[ii];
        }

        return dVecSum;
    }

    public void set(double dVal) {
        for (int ii = 0; ii < data.length; ii++) {
            data[ii] = dVal;
        }
    }

    public void set(double [] arr) {
        data = new double[arr.length];
        for (int ii = 0; ii < arr.length; ii++) {
            data[ii] = arr[ii];
        }
    }

    public void set(int [] arr) {
        data = new double[arr.length];
        for (int ii = 0; ii < arr.length; ii++) {
            data[ii] = arr[ii];
        }
    }

    public void set(int col, double val) {
        data[col] = val;
    }

    public void set(int start, int end, double val) {
        for (int ii = start; ii < end; ii++) {
            data[ii] = val;
        }
    }

/*
    public void set(int iSize, double dVal) {
        data = new double[iSize];
        set(dVal);
    }
*/
    public void setRNDvalue(double dCenter, double dRange) {
        double dMin = dCenter - (dRange / 2);
        for (int ii = 0; ii < data.length; ii++) {
            double dVal = dRange * Math.random();
            data[ii] = dMin + dVal;
        }
    }

    public int size() {
        return data.length;
    }

    public void setRNDvalue(double dRange) {
        for (int ii = 0; ii < data.length; ii++) {
            double dMin = data[ii] - (dRange / 2);
            double dVal = dRange * Math.random();
            data[ii] = dMin + dVal;
        }
    }

    /**
     * Substraction
     * @param dvSub
     * @return
     */
    public DoubleVec sub(DoubleVec dvSub) {

        DoubleVec ret = new DoubleVec(data.length);
        for (int ii = 0; ii < ret.data.length; ii++) {
            ret.data[ii] = data[ii] - dvSub.data[ii];
        }

        return ret;
    }

    /**
     * Calculates the Tanimoto coefficient according
     * http://www.pnylab.com/pny/papers/nmet/nmet/
     * @param dVec1 vector1
     * @param dVec2 vector2
     * @return Tanimoto: 1.0: maximum similarity, 0: maximum dissimilarity.
     *
     */
    static public double getTanimotoDistance(DoubleVec dVec1, DoubleVec dVec2) {

        double dSum = 0;
        double dAtB = mult(dVec1, dVec2);
        double dAtA = mult(dVec1, dVec1);
        double dBtB = mult(dVec2, dVec2);

        dSum = dAtB / (dAtA + dBtB - dAtB);

        return dSum;
    }

    static public double getTanimotoDistance(double[] d1, double[] d2) {

        double dSum = 0;
        double dAtB = mult(d1, d2);
        double dAtA = mult(d1, d1);
        double dBtB = mult(d2, d2);

        dSum = dAtB / (dAtA + dBtB - dAtB);

        return dSum;
    }
    
    /**
     * Calculates the Inverse Tanimoto coefficient
     * @param dVec1 vector1
     * @param dVec2 vector2
     * @return Tanimoto: 0.0: maximum similarity, 1.0: maximum dissimilarity.
     *
     */
    static public double getTanimotoDistInv(DoubleVec dVec1, DoubleVec dVec2) {

        double dSum = 0;
        double dAtB = mult(dVec1, dVec2);
        double dAtA = mult(dVec1, dVec1);
        double dBtB = mult(dVec2, dVec2);

        dSum = 1.0 - (dAtB / (dAtA + dBtB - dAtB));

        return dSum;
    }
    
    static public double getTanimotoDistInv(double[] d1, double[] d2) {
    	return 1.0 - getTanimotoDistance(d1,d2);
    }

    static public double getTanimotoDistInvDotProd(DoubleVec dVec1, DoubleVec dVec2) {

        double dSum = 0;
        double dAtB = mult(dVec1, dVec2);
        double dAtA = dVec1.mDotProd;
        double dBtB = dVec2.mDotProd;

        dSum = 1.0 - (dAtB / (dAtA + dBtB - dAtB));

        return dSum;
    }

    public String toString() {
    	StringBuilder str = new StringBuilder();
        for (int ii = 0; ii < data.length; ii++) {
            String sVal = NF.format(data[ii]);
            str.append(sVal + " ");
        }
        return str.toString();
    }

    /**
     *
     * @param iNumDigits number of decimal places
     * @return String with doubles in 0.0 notation
     * 16.10.2003 MvK
     */
    public String toString(int iNumDigits) {
        StringBuffer str = new StringBuffer();

        String sFormat = "0";
        if(iNumDigits > 0)
            sFormat += ".";

        for (int ii = 0; ii < iNumDigits; ii++) {
            sFormat = sFormat + "0";
        }

        DecimalFormat nf = new DecimalFormat(sFormat);

        for (int ii = 0; ii < data.length; ii++) {
            String sVal = nf.format(data[ii]);
            str.append(sVal + " ");
        }

        return str.toString();
    }

    public double[] toArray() {

        return data;
    }
    

}