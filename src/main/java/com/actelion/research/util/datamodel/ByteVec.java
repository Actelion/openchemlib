/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

package com.actelion.research.util.datamodel;


import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Vector;

import com.actelion.research.util.BitUtils;
import com.actelion.research.util.BurtleHasher;

/**
 * <p>ByteVec: </p>
 * <p>Description: Vector with byte values </p>
 * @author Modest von Korff
 * 05.04.2005 start implementation
 */
public class ByteVec implements Comparable<ByteVec>, Serializable {

	private static final long serialVersionUID = 27052009;
	
    public static final Random RANDOM = new Random();


    public static final int COSINE = 1;
    public static final int EUCLIDEAN = 2;
    public static final int EUCLIDEAN_FAST = 3;
    public static final int TANIMOTO = 4;
    public static final int TANIMOTO_INV = 5;

    

    private byte data[];
    
    private int hash;

    /**
     * Deep copy
     * @param v
     */
    public ByteVec(ByteVec v) {
        data = new byte[v.data.length];
        for (int ii = 0; ii < data.length; ii++) {
            data[ii] = v.data[ii];
        }
        calcHashCode();
    }

    /**
     * Hash code is calculated.
     * @param s
     */
    public ByteVec (String s){
    	data = s.getBytes();
    	calcHashCode();
    }
    
    public ByteVec(int iSize) {
        data = new byte[iSize];
        hash = -1;
    }

    public ByteVec(byte [] arr) {
        data = new byte[arr.length];
        for (int i = 0; i < data.length; i++) {
            data[i] = arr[i];
        }
    	calcHashCode();
    }

    public ByteVec(int [] arr) {
    	int fac = (Integer.SIZE / Byte.SIZE);
    	
    	data = new byte[arr.length * fac];
    	 
    	 for (int i = 0; i < arr.length; i++) {
			byte [] ab = IntVec.getByteVec(arr[i]);
			
			System.arraycopy(ab, 0, data, i*fac, ab.length);
			
		}
    	 
    }

    public ByteVec add(ByteVec dvVec) {
        ByteVec ret = new ByteVec(data.length);

        if (data.length != dvVec.data.length) {
            throw new RuntimeException();
        }

        for (int i = 0; i < ret.data.length; i++) {
            ret.data[i] = (byte)(data[i] + dvVec.data[i]);
        }

    	ret.calcHashCode();
        
       return ret;
    }

    public void append(String s) {

        int n = data.length;

        byte [] a = s.getBytes();

        int newSize = n + a.length;

        resize(newSize);

        System.arraycopy(a, 0, data, n, a.length);

    }

    private void resize(int newSize) {

        byte [] arrTmp = new byte[newSize];

        int n = Math.min(data.length, newSize);

        System.arraycopy(data, 0, arrTmp, 0, n);

        data = arrTmp;
    }

    /**
     *
     * @param o ByteVec
     * @return -1 if the first different value is smaller than the corresponding
     * value in dv. 0 if bot vectors are equal. 1 if the first different value
     * is bigger than the corresponding value in dv.
     */
    public int compareTo(ByteVec o) {
        
        int cmp = 0;
        
        if(data.length < o.data.length){
        	return -1;
        } else if(data.length > o.data.length){
        	return 1;
        } 
        
        for (int ii = 0; ii < data.length; ii++) {
            if (data[ii] > o.data[ii]) {
                cmp = 1;
                break;
            }
            else if (data[ii] < o.data[ii]) {
                cmp = -1;
                break;
            }
        }
        return cmp;
    }

    public static double distance(ByteVec dVec1, ByteVec dVec2, int metric) throws
        Exception {
        double dDist = 0;
        // Vectors have to be normed.
        if (metric == ByteVec.COSINE) {
            dDist = Cosine(dVec1, dVec2);
        }
        else if (metric == ByteVec.EUCLIDEAN) {
            dDist = euclideanDistance(dVec1, dVec2);
        }
        else if (metric == ByteVec.EUCLIDEAN_FAST) {
            dDist = EuclideanDistanceFast(dVec1, dVec2);
        }
        else if (metric == ByteVec.TANIMOTO) {
            dDist = getTanimotoDist(dVec1, dVec2);
        }
        else if (metric == ByteVec.TANIMOTO_INV) {
            dDist = getTanimotoDistInv(dVec1, dVec2);
        }
        else {
            throw new Exception("Unknown distance metric.");
        }

        return dDist;
    }

    public static ByteVec devide(ByteVec dVec1, ByteVec dVec2) {

        ByteVec dVecDev = new ByteVec(dVec1.data.length);

        for (int ii = 0; ii < dVec1.data.length; ii++) {
            dVecDev.data[ii] = (byte)(dVec1.data[ii] / dVec2.data[ii]);
        }

        dVecDev.calcHashCode();
        
        return dVecDev;
    }

    public boolean equals(Object o) {
    	
    	ByteVec bv = (ByteVec)o;
    	
    	if(size()!= bv.size())
    		return false;

    	if(hash != bv.hash){
    	    return false;
        }

        boolean equal = true;
        for (int i = 0; i < data.length; i++) {
          if(data[i] != bv.data[i]) {
              equal = false;
              break;
          }
        }
        return equal;
    }

    public boolean equal(ByteVec dv) {
        return equals(dv);
    }
    /**
     * Euclidean distance
     * @param dVec1
     * @param dVec2
     * @return
     */
    static public double euclideanDistance(ByteVec dVec1, ByteVec dVec2) {
        double dDist = 0;

        double dSum = 0;
        for (int ii = 0; ii < dVec1.data.length; ii++) {
            dSum += (dVec1.data[ii] - dVec2.data[ii]) *
                (dVec1.data[ii] - dVec2.data[ii]);
        }

        dDist = Math.sqrt(dSum);

        return dDist;
    }

    /**
     * Euclidean distance without sqrt
     * @param dVec1
     * @param dVec2
     * @return
     */
    static public double EuclideanDistanceFast(ByteVec dVec1, ByteVec dVec2) throws
        ArrayIndexOutOfBoundsException {

        if (dVec1.data.length != dVec2.data.length) {
            String sMessage = "Length double vector 1: " + dVec1.data.length +
                "Length double vector 2: " + dVec2.data.length + "\n";
            throw (new ArrayIndexOutOfBoundsException(sMessage));
        }

        double dSum = 0;
        for (int ii = 0; ii < dVec1.data.length; ii++) {
            dSum += (dVec1.data[ii] - dVec2.data[ii]) *
                (dVec1.data[ii] - dVec2.data[ii]);
        }

        return dSum;
    }

    public byte [] get() {
        return data;
    }

    public byte get(int col) {
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

    private void calcHashCode(){
    	hash = getHashCode(data);
    }
    
    public static int getHashCode(byte [] a){
    	return BurtleHasher.hashlittle(a, 13);
    }
    
    /**
     * The hash code has to be calculated before. 
     * Should be done automatically in any of the ByteVec functions changing a value, but please check.
     */
    public int hashCode(){
    	return hash;
    }
    
    /**
     * Vectors have to be normed!
     * @param dVec1 normed vector1
     * @param dVec2 normed vector2
     * @return Cosine
     */
    static public double Cosine(ByteVec dVec1, ByteVec dVec2) {
        double cosine = 0;
        for (int ii = 0; ii < dVec1.data.length; ii++) {
            cosine += dVec1.data[ii] * dVec2.data[ii];
        }
        return cosine;
    }

    static public double cubicDistance(ByteVec dVec1, ByteVec dVec2) {

        double dSum = 0;
        for (int ii = 0; ii < dVec1.data.length; ii++) {
            double dDist = Math.abs( (dVec1.data[ii] - dVec2.data[ii]));
            dSum += dDist * dDist * dDist;
        }

        return dSum;
    }

    public void initRND(byte min, byte max) {
        byte dRange = (byte)(max - min);
        for (int kk = 0; kk < data.length; kk++) {
            byte dVal = (byte)(dRange * RANDOM.nextInt(dRange) + min);
            data[kk] = dVal;
        }
        calcHashCode();
    }

    public static List<ByteVec> initRND(byte min, byte max, int lenByteVec, int numelements) {
        List<ByteVec> li = new ArrayList<ByteVec>(numelements);
        for (int ii = 0; ii < numelements; ii++) {
            ByteVec bv = new ByteVec(lenByteVec);
            bv.initRND(min,max);
            li.add(bv);
        }
        return li;
    }


    static public double manhattanBlockDistance(ByteVec dVec1,
                                                ByteVec dVec2) {
        double dDist = 0;

        double dSum = 0;
        for (int ii = 0; ii < dVec1.data.length; ii++) {
            dSum += Math.abs(dVec1.data[ii] - dVec2.data[ii]);
        }

        dDist = Math.sqrt(dSum);

        return dDist;
    }

    public ByteVec mult(double dFactor) {

        ByteVec ret = new ByteVec(data.length);
        for (int ii = 0; ii < ret.data.length; ii++) {
            ret.data[ii] = (byte)(data[ii] * dFactor);
        }

        ret.calcHashCode();
        
        return ret;
    }

    public static final double mult(ByteVec bv1, ByteVec bv2) {
        return mult(bv1.data, bv2.data);
    }
    
    public static final double mult(byte [] a1, byte [] a2) {

        double sum = 0.0;

        for (int i = 0; i < a1.length; i++) {
            sum += a1[i] * a2[i];
        }

        return sum;
    }

    /**
     * Elementwise multiplication
     * @param dVec1 input vector
     * @param dVec2 input vector
     * @return ByteVec
     */
    public static ByteVec multEl(ByteVec dVec1, ByteVec dVec2) {

        ByteVec bvMult = new ByteVec(dVec1.data.length);

        for (int ii = 0; ii < dVec1.data.length; ii++) {
            bvMult.data[ii] = (byte)(dVec1.data[ii] * dVec2.data[ii]);
        }

        bvMult.calcHashCode();
        
        return bvMult;
    }

    public void norm2One() {
        double norm = getNorm();
        for (int ii = 0; ii < data.length; ii++) {
            data[ii] /= norm;
        }
        calcHashCode();
    }

    public static ByteVec minus(ByteVec dVec1, ByteVec dVec2) {

        ByteVec bvSub = new ByteVec(dVec1.data.length);

        for (int ii = 0; ii < dVec1.data.length; ii++) {
            bvSub.data[ii] = (byte)(dVec1.data[ii] - dVec2.data[ii]);
        }

        bvSub.calcHashCode();
        
        return bvSub;
    }

    public void reduce(Vector<Integer> vecIndices) {
        byte [] arr = new byte [vecIndices.size()];
        for (int ii = 0; ii < vecIndices.size(); ii++) {
            int iIndex = ((Integer) vecIndices.get(ii)).intValue();
            arr[ii] = data[iIndex];
        }
        data = arr;
        
        calcHashCode();
    }

    public static ByteVec plus(ByteVec dVec1, ByteVec dVec2) {

        ByteVec bvSum = new ByteVec(dVec1.data.length);

        for (int ii = 0; ii < dVec1.data.length; ii++) {
            bvSum.data[ii] = (byte)(dVec1.data[ii] + dVec2.data[ii]);
        }

        bvSum.calcHashCode();
        
        return bvSum;
    }

    public void set(byte dVal) {
        for (int ii = 0; ii < data.length; ii++) {
            data[ii] = dVal;
        }
        calcHashCode();
    }

    /**
     * Hash value has to be recalculated!
     * @param col
     * @param val
     */
    public void set(int col, byte val) {
        data[col] = val;
        hash = -1;
    }

    public void setRNDvalue(double dCenter, double dRange) {
        double dMin = dCenter - (dRange / 2);
        for (int ii = 0; ii < data.length; ii++) {
            double dVal = dRange * Math.random();
            data[ii] = (byte)(dMin + dVal);
        }
        calcHashCode();
    }

    public int size() {
        return data.length;
    }

    public void setRNDvalue(double dRange) {
        for (int ii = 0; ii < data.length; ii++) {
            double dMin = data[ii] - (dRange / 2);
            double dVal = dRange * Math.random();
            data[ii] = (byte)(dMin + dVal);
        }
        calcHashCode();
    }

    /**
     * Substraction
     * @param dvSub
     * @return
     */
    public ByteVec sub(ByteVec dvSub) {

        ByteVec ret = new ByteVec(data.length);
        for (int ii = 0; ii < ret.data.length; ii++) {
            ret.data[ii] = (byte)(data[ii] - dvSub.data[ii]);
        }

        ret.calcHashCode();
        
        return ret;
    }

    /**
     * Calculates the Tanimoto coefficient
     * @param bv1 vector1
     * @param bv2 vector2
     * @return Tanimoto:  0.0: maximum distance, 1.0: maximum similarity.
     * Calculation according http://www.pnylab.com/pny/papers/nmet/nmet/
     * Congruent with DoubleVec
     */
    static public double getTanimotoDistBitWise(ByteVec bv1, ByteVec bv2) {

        int bitsOR = 0, bitsXOR = 0;
        for (int i = 0; i < bv1.data.length; i++) {
            bitsOR += Integer.bitCount(bv1.data[i] | bv2.data[i]);
            bitsXOR += Integer.bitCount(bv1.data[i] & bv2.data[i]);
        }

        if(bitsXOR == 0)
            return 0;

        double dSum = (double)(bitsXOR) / (double)(bitsOR);

        return dSum;
    }

    public int getBitsSet(){
    	
    	int n=0;
    	
    	for (int i = 0; i < data.length; i++) {
    		int t = data[i];
    		n += Integer.bitCount(t);
		}
    	
    	return n;
    }
    
    
    /**
     * Calculates the Tanimoto coefficient
     * @param bv1 vector1
     * @param bv2 vector2
     * @return Tanimoto: 0.0: maximum distance, 1.0: maximum similarity.
     *
     */
    static public double getTanimotoDist(ByteVec bv1, ByteVec bv2) {

        double sum = 0;
        double dAtB = mult(bv1, bv2);
        double dAtA = mult(bv1, bv1);
        double dBtB = mult(bv2, bv2);

        sum = dAtB / (dAtA + dBtB - dAtB);

        return sum;
    }
    
    static public double getTanimotoDist(byte [] a1, byte [] a2) {

        double sum = 0;
        double dAtB = mult(a1, a2);
        double dAtA = mult(a1, a1);
        double dBtB = mult(a2, a2);

        sum = dAtB / (dAtA + dBtB - dAtB);

        return sum;
    }
    
    /**
     * Byte wise inverse Tanimoto distance.
     * @param bv1
     * @param bv2
     * @return
     */
    static public double getTanimotoDistInv(ByteVec bv1, ByteVec bv2) {

        double dSum = 0;
        double dAtB = mult(bv1, bv2);
        double dAtA = mult(bv1, bv1);
        double dBtB = mult(bv2, bv2);

        dSum = 1.0 - (dAtB / (dAtA + dBtB - dAtB));

        return dSum;
    }

    /**
     * Converts a byte into its decimal.
     */
    public String toString() {
    	
        StringBuilder sb = new StringBuilder();
        
        DecimalFormat nf = new DecimalFormat("0");

        for (int i = 0; i < data.length; i++) {
            sb.append(nf.format(data[i]) + " ");
        }

        return sb.toString();
    }

    /**
     * Converts the bytes into chars.
     * The original string is reconstructed.
     * @return
     */
    public String toStringString() {
    	
    	StringBuilder sb = new StringBuilder();
        
        for (int i = 0; i < data.length; i++) {
        	
        	char c = (char)data[i];
        	
            sb.append(c);
        }

        return sb.toString();
    }
    
    public String toBinaryString() {
    	
    	StringBuilder sb = new StringBuilder();
        
        for (int i = 0; i < data.length; i++) {
        	
            String sVal = Integer.toBinaryString(data[i]);
            
            String sValFormat = "";
            
            for (int j = 0; j < sVal.length() - 1; j++) {
            	
              sValFormat += sVal.charAt(j) + " ";
            }
            
            sValFormat += sVal.charAt(sVal.length() - 1);

            int rest = Byte.SIZE - sVal.length();
            
            for (int j = 0; j < rest; j++) {
            	
              sValFormat = "0 " + sValFormat;
            }

            sb.append(sValFormat + " ");
        }

        return sb.toString();
    }
    
    public String toBinaryStringDense() {
    	
    	StringBuilder sb = new StringBuilder();
        
        for (int i = 0; i < data.length; i++) {
        	
            String sVal = Integer.toBinaryString(data[i]);
            
            String sValFormat = "";
            
            for (int j = 0; j < sVal.length() - 1; j++) {
            	
              sValFormat += sVal.charAt(j) + "";
            }
            
            sValFormat += sVal.charAt(sVal.length() - 1);

            int rest = Byte.SIZE - sVal.length();
            
            for (int j = 0; j < rest; j++) {
            	
              sValFormat = "0" + sValFormat;
            }

            sb.append(sValFormat + "");
        }

        return sb.toString();
    }

    /**
     *
     * @param iNumDigits number of decimal places
     * @return String with doubles in 0.0 notation
     * 16.10.2003 MvK
     */
    public String toString(int iNumDigits) {
    	StringBuilder str = new StringBuilder();

        String sFormat = "0.";
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

    public byte[] toArray() {

        return data;
    }

    
	public static byte [] toByteArray(long data) {

		return new byte[] {
			(byte) ((data >> 56) & 0xff),
			(byte) ((data >> 48) & 0xff),
			(byte) ((data >> 40) & 0xff),
			(byte) ((data >> 32) & 0xff),
			(byte) ((data >> 24) & 0xff),
			(byte) ((data >> 16) & 0xff),
			(byte) ((data >> 8) & 0xff),
			(byte) ((data >> 0) & 0xff),
		};
	}
	
	public static byte [] toByteArray(double [] data) {
		
		if (data == null) 
			return null;
		
		byte[] byts = new byte[data.length * 8];
		
		for (int i = 0; i < data.length; i++)
			System.arraycopy(toByteArray(Double.doubleToRawLongBits(data[i])), 0, byts, i * 8, 8);
		
		return byts;
	}
	
	public static byte [] toByteArray(int [] data) {
		
		if (data == null) 
			return null;
		
		byte[] byts = new byte[data.length * 4];
		
		for (int i = 0; i < data.length; i++)
			System.arraycopy(toByteArray(data[i]), 0, byts, i * 4, 4);
		
		return byts;
	}

	public static byte [] toByteArray(int data) {
		
		return new byte[] {
		
			(byte)((data >> 24) & 0xff),
			
			(byte)((data >> 16) & 0xff),
			
			(byte)((data >> 8) & 0xff),
			
			(byte)((data >> 0) & 0xff),
		
		};
	}

    public static byte [] reverse(byte [] a) {

        byte [] b = new byte[a.length];

        for (int i = 0; i < a.length; i++) {
            b[a.length-i-1] = a[i];
        }

        return b;
    }

    public static String toString(byte [] a, String seperator) {
    	    			
		StringBuilder sb = new StringBuilder();
		
		for (int i = 0; i < a.length; i++) {
			
			sb.append(a[i]);
			
			if(i < a.length-1){
				sb.append(seperator);
			}
			
		}
					
		return sb.toString();

    }

    public static String toString(byte [] a) {
		return toString(a, " ");
    }

    static public ByteVec testS1() {
       byte [] arr = {0,0,0,0,1,0,0,1,1,0,0,0,0,2,1,0,1,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,3,2,3,4,1,2,3,3,3,3,2,16,4,4,1,2,1,1,3,3,2,3,2,4,10,8,8,6,7,7,6,7,7,6,7,26,0,2,4,5,7,8,7,4,3,4,4,6,0,1,4,7,8,7,5,3,1,0,0,0,5,15,24,29,27,23,18,13,10,7,7,92,18,18,9,0,1,4,9,15,20,22,19,18};
       ByteVec dv = new ByteVec(arr);
       return dv;
    }
    static public ByteVec testTakeda1() {
       byte [] arr = {0,0,0,0,1,1,0,1,1,0,1,1,0,1,2,0,1,2,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,4,3,4,7,5,2,5,7,6,3,4,13,4,3,2,1,4,2,2,2,4,2,2,4,11,11,10,9,5,6,7,9,9,11,8,24,1,4,6,6,7,8,7,4,3,5,6,15,0,2,5,6,5,5,6,5,2,0,0,0,5,13,21,28,31,29,21,13,9,9,14,95,18,18,9,0,1,4,9,15,20,22,19,18};
       ByteVec dv = new ByteVec(arr);
       return dv;
    }

    static public ByteVec test01() {
       byte [] arr = {10,0,126};
       ByteVec dv = new ByteVec(arr);
       return dv;
    }

    static public ByteVec meanClust() {
       byte [] arr = {15,0,126};
       ByteVec dv = new ByteVec(arr);
        return dv;
    }

}