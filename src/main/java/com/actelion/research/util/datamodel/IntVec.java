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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import com.actelion.research.util.BitUtils;
import com.actelion.research.util.BurtleHasher;

public class IntVec implements Comparable<IntVec> {


    @SuppressWarnings("unused")
	private static final int UPPER_MASK = 0x80000000; // most significant bit
    
    public static final int LEN_INTEGER_BYTES = Integer.SIZE / Byte.SIZE;

    
    public static final int MASK_FIRST_BYTE  = 0x000000FF;
    public static final int MASK_SEC_BYTE    = 0x0000FF00;
    public static final int MASK_THIRD_BYTE  = 0x00FF0000;
    public static final int MASK_FOURTH_BYTE = 0xFF000000;

    public static final int MASK_INVERSE_FIRST_BYTE  = 0xFFFFFF00;
    public static final int MASK_INVERSE_SEC_BYTE    = 0xFFFF00FF;
    public static final int MASK_INVERSE_THIRD_BYTE  = 0xFF00FFFF;
    public static final int MASK_INVERSE_FOURTH_BYTE = 0x00FFFFFF;

    private int data[];
    
    private int hash;

    public IntVec() {
        
    }
    
    /**
     * Deep copy, hash code is calculated.
     * @param iv
     */
    public IntVec(IntVec iv) {
        this(iv.data);
    }

    /**
     * Don't forget to calculate the hash code after setting the bits!
     * @param size
     */
    public IntVec(int size) {
        init();
        data = new int[size];
    }

    /**
     * Deep copy, hash code is calculated.
     * @param arr
     */
    public IntVec(int [] arr) {
        init();
        
        data = new int[arr.length];
        
        System.arraycopy(arr, 0, data, 0, arr.length);
        
        calculateHashCode();
        
    }
	/**
	 * If true each int contains only a 0 or a 1. The binary information is summarized by factor 32 (Integer.size).
	 * @param arr
	 * @param binary
	 * @throws Exception
	 */
    public IntVec(int[] arr, boolean binary) {
        init();
        int len = (arr.length + Integer.SIZE-1)/ Integer.SIZE;
        
        data = new int[len];
        for (int i = 0; i < arr.length; i++) {
			if(arr[i]==1){
				setBit(i);
			} else if(arr[i]>1 || arr[i]<0){
				throw new RuntimeException("abs(value) in array larger than one, no binary data");
			}
		}
    }
    
    /**
     * Puts an boolean array into an int array.
     * @param arr
     * @throws Exception
     */
    public IntVec(boolean [] arr) {
        init();
        int len = (arr.length + Integer.SIZE-1)/ Integer.SIZE;
        
        data = new int[len];
        
        for (int i = 0; i < arr.length; i++) {
			if(arr[i]==true){
				setBit(i);
			} 
		}
    }
    public IntVec(List<Integer> vec) {
        init();
        data = new int[vec.size()];
        for (int i = 0; i < data.length; i++) {
            data[i] = vec.get(i).intValue();
        }
    }


    public IntVec add(IntVec dvVec) {
        IntVec ret = new IntVec(data.length);

        if (data.length != dvVec.data.length) {
            throw new RuntimeException();
        }

        for (int ii = 0; ii < ret.data.length; ii++) {
            ret.data[ii] = (int)(data[ii] + dvVec.data[ii]);
        }

        return ret;
    }

    public void copy(IntVec ivOrigin){
        System.arraycopy(ivOrigin.data, 0, data, 0, data.length);
    }

    public void clear() {
        for (int ii = 0; ii < data.length; ii++)
            data[ii] = 0;
        hash = -1;
    }

    /**
     *
     * @param iv IntVec
     * @return -1 if the first different value is smaller than the corresponding
     * value in dv. 0 if bot vectors are equal. 1 if the first different value
     * is bigger than the corresponding value in dv.
     */
    public int compareTo(IntVec iv) {
        
        int cmp = 0;
        for (int i = 0; i < data.length; i++) {
            if (data[i] > iv.data[i]) {
                cmp = 1;
                break;
            }
            else if (data[i] < iv.data[i]) {
                cmp = -1;
                break;
            }
        }
        return cmp;
    }

    
    public static IntVec devide(IntVec iv1, IntVec iv2) {

        IntVec dVecDev = new IntVec(iv1.data.length);

        for (int i = 0; i < iv1.data.length; i++) {
            dVecDev.data[i] = (int)(iv1.data[i] / iv2.data[i]);
        }

        return dVecDev;
    }

    public boolean equal(IntVec iv) {
        boolean eq = true;
        if(size() != iv.size()){
        	return false;
        }
        for (int i = 0; i < data.length; i++) {
          if(data[i] != iv.data[i]) {
              eq = false;
              break;
          }
        }
        return eq;
    }

    public boolean equals(Object o) {
        return equal((IntVec)o);
    }

    /**
     * Euclidean distance
     * @param iv1
     * @param iv2
     * @return
     */
    static public double getEuclidDist(IntVec iv1, IntVec iv2) {
        double dDist = 0;

        double dSum = 0;
        for (int i = 0; i < iv1.data.length; i++) {
            dSum += (iv1.data[i] - iv2.data[i]) *
                (iv1.data[i] - iv2.data[i]);
        }

        dDist = Math.sqrt(dSum);

        return dDist;
    }

    static public double getEuclidDistBitWise(IntVec iv1, IntVec iv2) {
        int bitsXOR = 0;
        for (int i = 0; i < iv1.data.length; i++) {
            bitsXOR += Integer.bitCount(iv1.data[i] ^ iv2.data[i]);
        }
        return Math.sqrt(bitsXOR);
    }

    public static int getOverlappingBitCount(int [] a, int [] b){
        int bits = 0;
        for (int i = 0; i < a.length; i++) {
            bits += Integer.bitCount(a[i] & b[i]);
        }
        return bits;
    }

    /**
     * Euclidean distance without sqrt
     * @param dVec1
     * @param dVec2
     * @return
     */
    static public double getEuclidDistFast(IntVec dVec1, IntVec dVec2) throws
        ArrayIndexOutOfBoundsException {

        if (dVec1.data.length != dVec2.data.length) {
            String sMessage = "Length double vector 1: " + dVec1.data.length +
                "Length double vector 2: " + dVec2.data.length + "\n";
            throw (new ArrayIndexOutOfBoundsException(sMessage));
        }

        double dSum = 0;
        for (int i = 0; i < dVec1.data.length; i++) {
            dSum += (dVec1.data[i] - dVec2.data[i]) *
                (dVec1.data[i] - dVec2.data[i]);
        }

        return dSum;
    }
    
    /**
     * Makes an logical OR and calculates the hash code.
     * @return
     */
    static public IntVec OR(IntVec iv1, IntVec iv2) {

    	IntVec iv = new IntVec(iv1.data.length);
    
	    for (int i = 0; i < iv1.data.length; i++) {
	    	iv.data[i] = iv1.data[i] | iv2.data[i]; 
		}
	    
	    iv.calculateHashCode();
	    
	    return iv;
    }
    
    
    /**
     * Makes an logical OR and calculates the hash code.
     * @param li
     * @return
     */
    static public IntVec OR(List<IntVec> li) {

    	IntVec iv = new IntVec(li.get(0));
    
    	int len = iv.data.length;
    	
	    for (int i = 1; i < li.size(); i++) {
	    	for (int j = 0; j < len; j++) {
	    		iv.data[j] = iv.data[j] | li.get(i).data[j];	
			}
		}
    
	    iv.calculateHashCode();
	    
	    return iv;
    }
    
    /**
     * Makes an logical AND and calculates the hash code.
     * @param iv1
     * @param iv2
     * @return
     */
    static public IntVec AND(IntVec iv1, IntVec iv2) {

    	IntVec iv = new IntVec(iv1.data.length);
    
	    for (int i = 0; i < iv1.data.length; i++) {
	    	iv.data[i] = iv1.data[i] & iv2.data[i]; 
		}
    
	    iv.calculateHashCode();
	    
	    return iv;
    }

    
    
    public static int [] getRND(int size) {
        int [] arr = new int [size];


        return arr;
    }

    public int [] get() {
        return data;
    }

    public int getByte(int indexBytes) {
    	
        return getByte(data, indexBytes);
    }
    
    public static int getByte(int [] data, int indexBytes) {
        int val = 0;
        int ind = indexBytes / 4;
        int indInInt = indexBytes % 4;

        switch(indInInt) {
            case 3:
                val = data[ind] & MASK_FIRST_BYTE;
                break;
            case 2:
                val = data[ind] & MASK_SEC_BYTE;
                val = val >>> 8;
                break;
            case 1:
                val = data[ind] & MASK_THIRD_BYTE;
                val = val >>>16;
                break;
            case 0:
                val = data[ind] & MASK_FOURTH_BYTE;
                val = val >>> 24;
                break;
        }

        return val;
    }
   
   
    public byte [] getByteVec() {
    	return getByteVec(data);
    }

    public static byte [] getByteVec(int intVal) {
    	byte [] arr = new byte[4];
        int val = 0;
       
    	arr[0] = (byte)(intVal & MASK_FIRST_BYTE);
    
        val = intVal & MASK_SEC_BYTE;
        arr[1] = (byte)(val >>> 8);
    
        val = intVal & MASK_THIRD_BYTE;
        arr[2] = (byte)(val >>>16);
    
        val = intVal & MASK_FOURTH_BYTE;
        arr[3] = (byte)(val >>> 24);

        return arr;
    }

    public static byte [] getByteVec(int [] a) {

        int fac = (Integer.SIZE / Byte.SIZE);

        byte [] arr = new byte[a.length * fac];

        for (int i = 0; i < a.length; i++) {

            final int intVal = a[i];

            int val = 0;

            final int indexByteVec = i*fac;

            arr[indexByteVec+0] = (byte)(intVal & MASK_FIRST_BYTE);

            val = intVal & MASK_SEC_BYTE;
            arr[indexByteVec+1] = (byte)(val >>> 8);


            val = intVal & MASK_THIRD_BYTE;
            arr[indexByteVec+2] = (byte)(val >>>16);


            val = intVal & MASK_FOURTH_BYTE;
            arr[indexByteVec+3] = (byte)(val >>> 24);
        }


        return arr;
    }

    public static int getSizeForBits(int bits){
    	
    	int sizeInteger = 0;
    	
    	sizeInteger = (bits + Integer.SIZE-1)/ Integer.SIZE;
    	
    	return sizeInteger;
    }
    
    public static int getNumberAbove(int [] a, int val){
    	
    	int n=0;
    	
    	for (int i = 0; i < a.length; i++) {
			if(a[i]>val){
				n++;
			}
		}
    	
    	return n;
    }
    
    public static int getInt(byte [] arr) {
    	
        int t1 = MASK_FIRST_BYTE & arr[0];
        int t2 = MASK_FIRST_BYTE & arr[1];
        int t3 = MASK_FIRST_BYTE & arr[2];
        int t4 = MASK_FIRST_BYTE & arr[3];
        
        
        t2 = (int)(t2 << 8);
        t3 = (int)(t3 << 16);
        t4 = (int)(t4 << 24);
        
        int v1 = t1;
        int v2 = (v1 | t2);
        int v3 = (v2 | t3);
        int v4 = (v3 | t4);
                
        return v4;
    }
    
    public int getBitsSet(){
    	int sum = 0;
	    for (int i = 0; i < data.length; i++) {
	        sum += Integer.bitCount(data[i]);
	    }
	    return sum;
    }

    public int [] getByteWise() {
    	int bytes = size()*4;
    	int [] arr = new int [bytes];
    	for (int i = 0; i < bytes; i++) {
			arr[i]=getByte(i);
		}
    	return arr;
    	
    }
    
    public int [] getBitWise() {
    	int bits = size()* Integer.SIZE;
    	int [] arr = new int [bits];
    	for (int i = 0; i < bits; i++) {
    		if(isBitSet(i))
    			arr[i]=1;
    		else
    			arr[i]=0;
		}
    	return arr;
    	
    }

    /**
     * Converts the IntVec into an array. Each field in the array contains a value that was derived with the given
     * resolution from IntVec.
     * @param iv
     * @param nValues so many values are encoded in IntVec.
     * @param bitsResolution So many bits were used to encode a single value.
     * @return
     */
    public static int [] extractForGivenResolution(IntVec iv, int nValues, int bitsResolution){

        int [] a = new int[nValues];

        int indexIntVec = 0;
        for (int i = 0; i < nValues; i++) {

            int v = 0;
            for (int j = 0; j < bitsResolution; j++) {
                if(iv.isBitSet(indexIntVec++)){
                    v |= 1 << j;
                }
            }
            a[i]=v;
        }
        return a;
    }
    
    public boolean allFieldsEquals(int v){
    	boolean b = true;

    	for (int i = 0; i < data.length; i++) {
	        if(data[i] != v) {
	        	b = false;
	        	break;
	        }
	    }
    	
    	
    	return b;
    }
    
    public int get(int col) {
        return data[col];
    }

    public double getNorm() {
        double dNorm = 0;

        for (int i = 0; i < data.length; i++) {
            dNorm += (data[i] * data[i]);
        }
        dNorm = Math.sqrt(dNorm);
        return dNorm;
    }
    
    public int hashCode(){
    	
    	if(hash==-1){
    		calculateHashCode();
    	}
    	
    	return hash;
    }
    
    public void calculateHashCode(){
    	if(data.length==1){
    		hash=data[0];
    	} else {
    		hash = BurtleHasher.hashlittle(data, 13);	
    	}
    }
    
    /**
     *
     * @param iv1 normed vector1
     * @param iv2 normed vector2
     * @return Cosine
     */
    static public double getCosine(IntVec iv1, IntVec iv2) {

        double ab = 0;

        double a=0;
        double b=0;

        for (int i = 0; i < iv1.data.length; i++) {
            ab += iv1.data[i] * iv2.data[i];

            a += iv1.data[i] * iv1.data[i];
            b += iv2.data[i] * iv2.data[i];
        }

        double c = Math.sqrt(ab)/(Math.sqrt(a)*Math.sqrt(b));

        return c;
    }

    static public double cubicDistance(IntVec dVec1, IntVec dVec2) {

        double dSum = 0;
        for (int i = 0; i < dVec1.data.length; i++) {
            double dDist = Math.abs( (dVec1.data[i] - dVec2.data[i]));
            dSum += dDist * dDist * dDist;
        }

        return dSum;
    }
    
    public static int calculateHashCode(IntVec iv){
    	int h = 0;
    	
    	int l = iv.size() * 4;
    	
    	byte [] a = new byte[l];
    	for (int i = 0; i < l; i++) {
			a[i] = (byte)iv.getByte(i);
		}
    	
    	h = BurtleHasher.hashlittle(a, 13);
    	
    	return h;
    }
    
    private static int[] convert(String sLine) {

        StringTokenizer st = new StringTokenizer(sLine);
        int[] dArray = new int[st.countTokens()];
        int i = 0;
        while (st.hasMoreTokens()) {
            String sNumber = st.nextToken();
            // The formatting sign "'" is not recognized by the Double.valueOf(...)
            // function.
            sNumber = sNumber.replaceAll("'", "");
            try {
                int val = (int) Double.parseDouble(sNumber);
                dArray[i] = val;
            }

            catch (NumberFormatException ex1) {
                System.err.println("No number: " + sNumber + ".");
                ex1.printStackTrace();
            }
            i++;
        }
        return dArray;
    }
    
    public static IntVec readBitStringDense(String s) {

    	int bits = s.length();
    	
    	if(bits%Integer.SIZE!=0){
    		throw new RuntimeException("Wrong size ("+ bits +") of string for coversion.");
    	}
    	
    	int size = bits/Integer.SIZE;
    	
        IntVec iv = new IntVec(size);

        for (int i = 0; i < bits; i++) {
        	
        	int index = bits-i-1;
        	
			if(s.charAt(index)=='1'){
				iv.setBit(i);
			} else if(s.charAt(index)!='0'){
				throw new RuntimeException("Illegal character.");
			} 
		}
       
        return iv;
    }

    private void init() {
    	hash = -1;
    }

    
    static public double manhattanBlockDistance(IntVec iv1,
                                                IntVec iv2) {
        double dDist = 0;

        double dSum = 0;
        for (int i = 0; i < iv1.data.length; i++) {
            dSum += Math.abs(iv1.data[i] - iv2.data[i]);
        }

        dDist = Math.sqrt(dSum);

        return dDist;
    }

    public IntVec mult(double factor) {

        IntVec ret = new IntVec(data.length);
        for (int iv = 0; iv < ret.data.length; iv++) {
            ret.data[iv] = (int)(data[iv] * factor);
        }

        return ret;
    }

    public static double mult(IntVec iv1, IntVec iv2) {
        double dSum = 0.0;
        for (int i = 0; i < iv1.data.length; i++) {
            dSum += iv1.data[i] * iv2.data[i];
        }
        return dSum;
    }

    public static double mult(int [] a, int [] b) {
        double dSum = 0.0;
        for (int i = 0; i < a.length; i++) {
            dSum += a[i] * b[i];
        }
        return dSum;
    }


    public static double multByteWise(IntVec iv1, IntVec iv2) {
        return multByteWise(iv1.data, iv2.data);
    }

    public static double multByteWise(int[] iv1, int[] iv2) {

        double dSum = 0.0;

        int b11,b12,b21,b22,b31,b32,b41,b42;
        for (int i = 0; i < iv1.length; i++) {
            b11 = iv1[i] & MASK_FIRST_BYTE;
            b12 = iv2[i] & MASK_FIRST_BYTE;

            b21 = (iv1[i] & MASK_SEC_BYTE) >> 8;
            b22 = (iv2[i] & MASK_SEC_BYTE) >> 8;

            b31 = (iv1[i] & MASK_THIRD_BYTE) >> 16;
            b32 = (iv2[i] & MASK_THIRD_BYTE) >> 16;

            b41 = (iv1[i] & MASK_FOURTH_BYTE) >> 24;
            b42 = (iv2[i] & MASK_FOURTH_BYTE) >> 24;

            dSum += b11 * b12 + b21 * b22 + b31 * b32 + b41 * b42;
        }

        return dSum;
    }
    
    public static double getSimilarityBytewiseOverlap(int[] a1, int[] a2) {

        int b11,b12,b21,b22,b31,b32,b41,b42;
        
		float similarity = 0;
		
		float ccOverlap=0;
		
		float ccTotal=0;
       
        for (int i = 0; i < a1.length; i++) {
        	
            b11 = a1[i] & MASK_FIRST_BYTE;
            b12 = a2[i] & MASK_FIRST_BYTE;

            b21 = (a1[i] & MASK_SEC_BYTE) >>> 8;
            b22 = (a2[i] & MASK_SEC_BYTE) >>> 8;

            b31 = (a1[i] & MASK_THIRD_BYTE) >>> 16;
            b32 = (a2[i] & MASK_THIRD_BYTE) >>> 16;

            b41 = (a1[i] & MASK_FOURTH_BYTE) >>> 24;
            b42 = (a2[i] & MASK_FOURTH_BYTE) >>> 24;

			ccOverlap += Math.min(b11, b12);
			ccTotal += Math.max(b11, b12);

			ccOverlap += Math.min(b21, b22);
			ccTotal += Math.max(b21, b22);

			ccOverlap += Math.min(b31, b32);
			ccTotal += Math.max(b31, b32);

			ccOverlap += Math.min(b41, b42);
			ccTotal += Math.max(b41, b42);

        }

        similarity = ccOverlap / ccTotal;
		
		return similarity;
    }
    
    public static IntVec subtractByteWise(IntVec iv1, IntVec iv2) {

    	IntVec ivSub = new IntVec(iv1.size());

        for (int i = 0; i < iv1.sizeBytes(); i++) {
        	int sub = iv1.getByte(i) - iv2.getByte(i);
        	ivSub.setByte(i,sub);
        }
        return ivSub;
    }

    public static IntVec maskByteWise(IntVec mask, IntVec query) {
    	IntVec ivSub = new IntVec(query);
        for (int i = 0; i < mask.sizeBytes(); i++) {
        	if(mask.getByte(i) > 0)
        		ivSub.setByte(i,0);
        }
        return ivSub;
    }

    
    /**
     * Elementwise multiplication
     * @param iv1 input vector
     * @param iv2 input vector
     * @return IntVec
     */
    public static IntVec multEl(IntVec iv1, IntVec iv2) {

        IntVec dVecMult = new IntVec(iv1.data.length);

        for (int i = 0; i < iv1.data.length; i++) {
            dVecMult.data[i] = (int)(iv1.data[i] * iv2.data[i]);
        }

        return dVecMult;
    }

    public void norm2One() {
        double norm = getNorm();
        for (int i = 0; i < data.length; i++) {
            data[i] /= norm;
        }
    	hash = -1;
    }

    public static IntVec minus(IntVec dVec1, IntVec dVec2) {

        IntVec dVecSub = new IntVec(dVec1.data.length);

        for (int i = 0; i < dVec1.data.length; i++) {
            dVecSub.data[i] = (int)(dVec1.data[i] - dVec2.data[i]);
        }

        return dVecSub;
    }

    public void read(String s) {
    	data = convert(s);
    	hash = -1;
    }

    public static IntVec [] read(File file) {
    	
    	List<IntVec> li = new ArrayList<IntVec>();
    	
    	try {
			BufferedReader buf = new BufferedReader(new InputStreamReader(new FileInputStream(file), StandardCharsets.UTF_8));
			
			while(buf.ready()) {
				String s = buf.readLine();
				int [] a = convert(s);
				IntVec v = new IntVec(a);
				li.add(v);
			}
			
			buf.close();
			
		} catch (FileNotFoundException e) {
			
			e.printStackTrace();
		} catch (IOException e) {
			
			e.printStackTrace();
		}
    	
    	IntVec [] arrIV = new IntVec [li.size()];
    	for (int i = 0; i < arrIV.length; i++) {
			arrIV[i]=(IntVec)li.get(i);
		}
    	
    	return arrIV;
    	
    }
    
    /**
     * Removes all values with no corresponding index in the list.
     * @param liIndices list with indices for values that will be kept.
     */
    public void reduce(List<Integer> liIndices) {
        int [] arr = new int [liIndices.size()];
        for (int i = 0; i < liIndices.size(); i++) {
            int iIndex = liIndices.get(i).intValue();
            arr[i] = data[iIndex];
        }
        data = arr;
    }

    public void resize(int newlen){

        if(data.length == newlen){
            return;
        }

        int intNewlen = 0;

        long max = Integer.MAX_VALUE;

        if(newlen >= max) {

            intNewlen = Integer.MAX_VALUE;

            new RuntimeException("Warning! Maximum length of integer array reached.").printStackTrace();

        } else {
            intNewlen = newlen;
        }

        int [] arr = new int [intNewlen];

        System.arraycopy(data, 0, arr, 0, Math.min(data.length, intNewlen));

        data = arr;

    }

    public static IntVec plus(IntVec dVec1, IntVec dVec2) {

        IntVec dVecSum = new IntVec(dVec1.data.length);

        for (int i = 0; i < dVec1.data.length; i++) {
            dVecSum.data[i] = (int)(dVec1.data[i] + dVec2.data[i]);
        }

        return dVecSum;
    }

    /**
     * Don't forget to set the hash code!
     * @param val
     */
    public void set(int val) {
        for (int i = 0; i < data.length; i++) {
            data[i] = val;
        }
    	hash = -1;
    }

    /**
     * Don't forget to set the hash code!
     * @param v
     */
    public void set(IntVec v) {
        for (int i = 0; i < data.length; i++) {
            data[i] = v.data[i];
        }
    	hash = -1;
    }

    /**
     * 
     * Don't forget to set the hash code!
     * @param col
     * @param val
     */
    public void set(int col, int val) {
        data[col] = val;
    	hash = -1;
    }
    /**
     * Counts from the left to the right
     * Don't forget to set the hash code!
     * @param i index bit
     */
    public void setBit(int i) {
        int ind = data.length - (i / Integer.SIZE) - 1;
        int indInInt = i % Integer.SIZE;
        int mask = 1;
        mask = mask << indInInt;
        data[ind] = data[ind] | mask;
    	hash = -1;
    }
    
    /**
     * Don't forget to set the hash code!
     * @param i
     */
    public void unsetBit(int i) {
        int ind = data.length - (i / Integer.SIZE) - 1;
        int indInInt = i % Integer.SIZE;
        int mask = 1;
        mask = mask << indInInt;
        data[ind] = data[ind] & ~mask;
    	hash = -1;
    }

    public boolean isValidBitIndex(int i) {
        boolean valid = false;

        int ind = (i / Integer.SIZE);

        if(ind < data.length){
            valid=true;
        }

        return valid;
    }

     
    
    /**
     * Counts from the left to the right
     * Don't forget to set the hash code!
     * @param i index byte
     * @param val value
     */

    public void setByte(int i, int val) {
        
        setByte(data, i, val);
        
    	hash = -1;
    }
    
    public static void setBytes(int data[], int val) {
    	
    	int n = data.length * LEN_INTEGER_BYTES;
    	
    	for (int i = 0; i < n; i++) {
    		setByte(data, i, val);
		}
    }
    
    public static void setByte(int data[], int i, int val) {
    	
        int ind = i / 4;
        int indInInt = i % 4;

        int mask = 0;
        switch(indInInt) {
            case 3:
                mask = MASK_INVERSE_FIRST_BYTE;
                break;
            case 2:
                mask = MASK_INVERSE_SEC_BYTE;
                val = val << 8;
                break;
            case 1:
                mask = MASK_INVERSE_THIRD_BYTE;
                val = val << 16;
                break;
            case 0:
                mask = MASK_INVERSE_FOURTH_BYTE;
                val = val << 24;
                break;
        }

        data[ind] = (data[ind] & mask) | val;
    	
    }
    

    public List<Integer> getIndicesBitsSet(){
    	List<Integer> li = new ArrayList<Integer>();
    	for (int i = 0; i < sizeBits(); i++) {
			if(isBitSet(i)){
				li.add(new Integer(i));
			}
		}
    	return li;
    }
    
    public boolean isBitSet(int i) {
    	return isBitSetNotStatic(data, i);
    }
    
    private boolean isBitSetNotStatic(int [] a, int i) {
        int ind = a.length - (i / Integer.SIZE) - 1;
        int indInInt = i % Integer.SIZE;
        int mask = 1;
        mask = mask << indInInt;
        if((a[ind] & mask) != 0)
            return true;
        else
            return false;
    }


    public void switchBit(int i) {
        int ind = data.length - (i / Integer.SIZE) - 1;
        int indInInt = i % Integer.SIZE;
        int mask = 1;
        mask = mask << indInInt;
        if((data[ind] & mask) != 0) {
            data[ind] = data[ind] & ~mask;
        } else {
            data[ind] = data[ind] | mask;
        }

    	hash = -1;
    }

    public void setBits(int iStart, int num) {
        for (int i = iStart; i < iStart + num; i++) {
            setBit(i);
        }
    }

    public void setBytes(int iStart, int num, int val) {
        for (int i = iStart; i < iStart + num; i++) {
            setByte(i, val);
        }
    	hash = -1;
    }


    public void setRNDvalue(double dCenter, double dRange) {
        double dMin = dCenter - (dRange / 2);
        for (int i = 0; i < data.length; i++) {
            double dVal = dRange * Math.random();
            data[i] = (int)(dMin + dVal);
        }
    	hash = -1;
    }

    public int size() {
        return data.length;
    }

    public int sizeBits() {
        return data.length * Integer.SIZE;
    }

    public int sizeBytes() {
        return data.length * LEN_INTEGER_BYTES;
    }


    public void setRNDvalue(double dRange) {
        for (int i = 0; i < data.length; i++) {
            double dMin = data[i] - (dRange / 2);
            double dVal = dRange * Math.random();
            data[i] = (int)(dMin + dVal);
        }
    	hash = -1;
    }

    /**
     * Substraction
     * @param dvSub
     * @return
     */
    public IntVec sub(IntVec dvSub) {

        IntVec ret = new IntVec(data.length);
        for (int i = 0; i < ret.data.length; i++) {
            ret.data[i] = (int)(data[i] - dvSub.data[i]);
        }

        return ret;
    }



    public String toString() {
    	
    	StringBuilder sb = new StringBuilder();
    	
        DecimalFormat nf = new DecimalFormat("0");

        for (int i = 0; i < data.length; i++) {
            sb.append(nf.format(data[i]));
            
            if(i<data.length-1){
                sb.append(" ");
            }
        }

        return sb.toString();
    }
    
    public String toStringHex() {
    	
    	StringBuilder sb = new StringBuilder();
    	

        for (int i = 0; i < data.length; i++) {
        	
            sb.append(Integer.toHexString(data[i]));
            
            if(i<data.length-1){
                sb.append(" ");
            }
        }

        return sb.toString();
    }

    public String toStringBinary() {
        StringBuilder sb = new StringBuilder();
        int si = size();
        for (int i = 0; i < si; i++) {
            String s = toStringBinary(get(i)) + " ";
            sb.append(s);
        }
        return sb.toString().trim();
    }
    
    public String toStringBinaryDense() {
        StringBuilder sb = new StringBuilder();
        int si = size();
        for (int i = 0; i < si; i++) {
            String s = toStringBinary(get(i), false);
            sb.append(s);
        }
        return sb.toString().trim();
    }
   

    public String toStringBytes() {
    	StringBuilder str = new StringBuilder();
        int si = sizeBytes();
        for (int ii = 0; ii < si; ii++) {
            String s = getByte(ii) + " ";
            str.append(s);
        }
        return str.toString().trim();
    }

    public double [] toDoubleBitWise() {

        double [] arr = new double [size() * Integer.SIZE];
        int cc = 0;
        for (int i = 0; i < data.length; i++) {
            int v = data[i];
            int mask = 1;
            for (int j = 0; j < Integer.SIZE; j++) {
                if((v & mask) != 0) {
                    arr[cc] = 1;
                } else {
                    arr[cc] = 0;
                }
                mask = mask << 1;
                cc++;
            }
        }

        return arr;
    }

    public int [] toIntByteWise() {
    	int [] arr = new int [size() * 4];
        
        for (int i = 0; i < sizeBytes(); i++) {
            arr[i] = getByte(i);
        }

        return arr;
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
        for (int i = 0; i < iNumDigits; i++) {
            sFormat = sFormat + "0";
        }

        DecimalFormat nf = new DecimalFormat(sFormat);

        for (int i = 0; i < data.length; i++) {
            String sVal = nf.format(data[i]);
            str.append(sVal + " ");
        }

        return str.toString();
    }

    public String write2String() throws IOException{

        StringBuilder sb = new StringBuilder();

        sb.append(data.length);
        sb.append(" ");
        sb.append(hash);
        sb.append(" ");
        sb.append(toString());

        return sb.toString();
    }


    public int[] toArray() {
        return data;
    }


    /**
     * Calculates the Tanimoto coefficient
     * @param iv1 vector1
     * @param iv2 vector2
     * @return Tanimoto: 1.0: identical bits are set. 0.0 no matching bits are set.
     * Calculation according http://www.pnylab.com/pny/papers/nmet/nmet/
     * Congruent with DoubleVec
     */
    static public double getTanimotoDistBitWise(IntVec iv1, IntVec iv2) {

        int bitsOR = 0, bitsAND = 0;

        for (int i = 0; i < iv1.data.length; i++) {

            bitsOR += Integer.bitCount(iv1.data[i] | iv2.data[i]);

            bitsAND += Integer.bitCount(iv1.data[i] & iv2.data[i]);
        }

        if(bitsAND == 0)
            return 0;

        double sum = (double)(bitsAND) / (double)(bitsOR);

        return sum;
    }

    /**
     *
     * @param arr1
     * @param arr2
     * @return Tanimoto: 1: identical bits are set. 0 no matching bits are set.
     */
    static public final double getTanimotoDistBitWise(int [] arr1, int [] arr2) {

        int bitsOR = 0, bitsAND = 0;
        for (int i = 0; i < arr1.length; i++) {
            bitsOR += Integer.bitCount(arr1[i] | arr2[i]);
            bitsAND += Integer.bitCount(arr1[i] & arr2[i]);
        }

        if(bitsAND == 0)
            return 0;

        double dSum = (double)(bitsAND) / (double)(bitsOR);

        return dSum;
    }
    /**
     * Calculates the Inverse Tanimoto coefficient
     * @param iv1 vector1
     * @param iv2 vector2
     * @return Inverse Tanimoto: 0: identical bits are set. 1 no matching bits are set.
     * Calculation according Duda, Hart, Stork; Pattern Classification;
     * Wiley 2001 p188.
     *
     */
    static public double getTanimotoDistInvBitWise(IntVec iv1, IntVec iv2) {
        return 1.0 - getTanimotoDistBitWise(iv1, iv2);
    }
    static public double getTanimotoDistInvBitWise(int [] arr1, int [] arr2) {
        return 1.0 - getTanimotoDistBitWise(arr1, arr2);
    }

    /**
     * Calculates the Inverse Tanimoto coefficient
     * @param iv1 vector1
     * @param iv2 vector2
     * @return Tanimoto: 1: identical bits are set. 0 no matching bits are set.
     *
     */
    static public double getTanimotoDist(IntVec iv1, IntVec iv2) {

        double sum = 0;
        double dAtB = mult(iv1, iv2);
        double dAtA = mult(iv1, iv1);
        double dBtB = mult(iv2, iv2);

        sum = dAtB / (dAtA + dBtB - dAtB);

        return sum;
    }

    static public double getTanimotoDist(int [] a, int [] b) {

        double sum = 0;
        double dAtB = mult(a, b);
        double dAtA = mult(a, a);
        double dBtB = mult(b, b);

        sum = dAtB / (dAtA + dBtB - dAtB);

        return sum;
    }

    /**
     *
     * @param iv1
     * @param iv2
     * @return Inverse Tanimoto: 0: identical bits are set. 1 no matching bits are set.
     */
    static final public double getTanimotoDistInv(IntVec iv1, IntVec iv2) {

        double sum = 0;
        double dAtB = mult(iv1, iv2);
        double dAtA = mult(iv1, iv1);
        double dBtB = mult(iv2, iv2);

        sum = 1.0 - (dAtB / (dAtA + dBtB - dAtB));

        return sum;
    }

    public static double getScoreQueryInBaseByteWise(IntVec query, IntVec base) {

        double dSumPosDiff = 0.0;
        double denominator = 0;

        int b11,b12,b21,b22,b31,b32,b41,b42;
        for (int i = 0; i < query.data.length; i++) {
            b11 = query.data[i] & MASK_FIRST_BYTE;
            b12 = base.data[i] & MASK_FIRST_BYTE;

            denominator += b11;
            int diff = b11 - b12;
            if(diff > 0) {
                dSumPosDiff += diff;
            }

            b21 = (query.data[i] & MASK_SEC_BYTE) >> 8;
            b22 = (base.data[i] & MASK_SEC_BYTE) >> 8;

            denominator += b21;
            diff = b21 - b22;
            if(diff > 0) {
                dSumPosDiff += diff;
            }

            b31 = (query.data[i] & MASK_THIRD_BYTE) >> 16;
            b32 = (base.data[i] & MASK_THIRD_BYTE) >> 16;

            denominator += b31;
            diff = b31 - b32;
            if(diff > 0) {
                dSumPosDiff += diff;
            }

            b41 = (query.data[i] & MASK_FOURTH_BYTE) >> 24;
            b42 = (base.data[i] & MASK_FOURTH_BYTE) >> 24;

            denominator += b41;
            diff = b41 - b42;
            if(diff > 0) {
                dSumPosDiff += diff;
            }

        }
        if(denominator > 0)
            dSumPosDiff /= denominator;

        return dSumPosDiff;
    }

    static public double getScoreQueryInBase(IntVec query, IntVec base) {

        double dSumPosDiff = 0;
        double denominator = 0;
        for (int i = 0; i < query.size(); i++) {
            int diff = query.get(i) - base.get(i);
            denominator += query.get(i);
            if(diff > 0) {
                dSumPosDiff += diff;

            }
        }
        if(denominator > 0)
            dSumPosDiff /= denominator;
        return dSumPosDiff;
    }
    public static double getScoreQueryInBaseBitWise(IntVec query, IntVec base) {
        return getScoreQueryInBaseBitWise(query.data, base.data);
    }

    /**
     * (sum bits set only in query) / (sum bits set in query)
     * @param query
     * @param base
     * @return
     */
    public static double getScoreQueryInBaseBitWise(int[] query, int[] base) {
        double sc = 0;

        for (int i = 0; i < query.length; i++) {

            int bitsCommon = 0;
            int bitsOnlyInQuery = 0;

            bitsCommon = query[i] | base[i];
            bitsOnlyInQuery = bitsCommon ^ base[i];

            sc += (double) Integer.bitCount(bitsOnlyInQuery) / (double) Integer.bitCount(query[i]);
        }
        sc /= (double)query.length;

        return sc;
    }

    /**
     * (sum bits set common) / (sum bits set in both)
     *
     *
     * @return
     */
    public static double getScoreFracBitsInCommonBitWise(IntVec v1, IntVec v2) {

        double sc = 0;

        double cc=0;
        for (int i = 0; i < v1.data.length; i++) {

            int bitsCommon = v1.data[i] & v2.data[i];
            int bitsTotalInV1 = v1.data[i] | v2.data[i];

            if(bitsTotalInV1 != 0) {
                sc += (double) Integer.bitCount(bitsCommon) / (double) Integer.bitCount(bitsTotalInV1);
            }
            cc++;
        }

        sc /= cc;

        return 1.0-sc;
    }




    /**
     *
     * @param iv1
     * @param iv2
     * @return Inverse Tanimoto: 0: identical bits are set. 1 no matching bits are set.
     */
    static public double getTanimotoDistInvByteWise(IntVec iv1, IntVec iv2) {

        double sum = 0;
        double dAtB = multByteWise(iv1, iv2);
        double dAtA = multByteWise(iv1, iv1);
        double dBtB = multByteWise(iv2, iv2);

        sum = 1.0 - (dAtB / (dAtA + dBtB - dAtB));

        return sum;
    }

    /**
     *
     * @param iv1
     * @param iv2
     * @return Inverse Tanimoto: 0: identical bits are set. 1 no matching bits are set.
     */
    static public double getTanimotoDistInvByteWise(int[] iv1, int[] iv2) {

        double sum = 0;
        double dAtB = multByteWise(iv1, iv2);
        double dAtA = multByteWise(iv1, iv1);
        double dBtB = multByteWise(iv2, iv2);

        sum = 1.0 - (dAtB / (dAtA + dBtB - dAtB));

        return sum;
    }

    /**
     * (sum bits set common) / (sum bits set in query)
     * @param query
     * @param base
     * @return
     */
    public static double getScoreFracBitsCommonQuery(IntVec query, IntVec base) {
        double sc = 0;

        for (int i = 0; i < query.data.length; i++) {

            int bitsCommon = query.data[i] & base.data[i];

            if(query.data[i] != 0) {
                sc += (double) Integer.bitCount(bitsCommon) / (double) Integer.bitCount(query.data[i]);
            }
        }

        sc /= (double)query.data.length;

        return 1.0-sc;
    }

    /**
     *
     * @param iv1
     * @param iv2
     * @return number of overlapping bits.
     */
    static public int getOverlap(IntVec iv1, IntVec iv2) {

        IntVec iv = IntVec.AND(iv1, iv2);

        return iv.getBitsSet();
    }

    public static void incrementByte(int data[], int i) {

        int ind = i / LEN_INTEGER_BYTES;

        int indInInt = i % LEN_INTEGER_BYTES;

        int valOld = 0;
        int increment = 1;

        int maskInverse = 0;

        switch(indInInt) {
            case 3:
                maskInverse = MASK_INVERSE_FIRST_BYTE;
                valOld = data[ind] & MASK_FIRST_BYTE;
                break;
            case 2:
                maskInverse = MASK_INVERSE_SEC_BYTE;
                valOld = data[ind] & MASK_SEC_BYTE;
                increment = increment << 8;
                break;
            case 1:
                maskInverse = MASK_INVERSE_THIRD_BYTE;
                valOld = data[ind] & MASK_THIRD_BYTE;
                increment = increment << 16;
                break;
            case 0:
                maskInverse = MASK_INVERSE_FOURTH_BYTE;
                valOld = data[ind] & MASK_FOURTH_BYTE;
                increment = increment << 24;
                break;
        }

        data[ind] = (data[ind]&maskInverse) | (valOld+increment);

    }

    public static void setAllBits(int [] a){
        for (int i = 0; i < a.length; i++) {
            a[i]=0xFFFFF;
        }
    }
    public static void writeBitStringDense(File fi, List<IntVec> li) throws IOException{
    	
    	BufferedWriter bw = new BufferedWriter(new FileWriter(fi));
    	
    	for (int i = 0; i < li.size(); i++) {
    		bw.append(li.get(i).toStringBinaryDense());
    		
    		if(i<li.size()-1){
    			bw.append("\n");
    		}
		}
    	
    	bw.close();
    }
    
    public static List<IntVec> readBitStringDense(File fi) throws IOException{
    	List<IntVec> li = new ArrayList<IntVec>();
    	
    	BufferedReader br = new BufferedReader(new FileReader(fi));
    	
    	String line="";
    	
    	int cc=0;
    	while((line=br.readLine())!=null){
    		
    		
    		try {
				IntVec iv = readBitStringDense(line);
				
				li.add(iv);
			} catch (Exception e) {
				System.err.println("Error in line " + cc + ".");
				e.printStackTrace();
			}
			cc++;
    	}
    	
    	br.close();
    	
    	return li;
    }

    public static IntVec read(InputStream s) throws IOException{
    	    	
    	int size = IntArray.parseInteger(s);
    	
    	int hash = IntArray.parseInteger(s);
    	
    	int [] a = new int [size];
    	
    	for (int i = 0; i < size; i++) {
    		a[i] = IntArray.parseInteger(s);
		}
    	
    	IntVec iv = new IntVec();
    	
    	iv.data = a;
    	
    	iv.hash = hash;
    	
    	return iv;
    	
    }
    
    public static boolean isBitSet(int [] a, int i) {
        int ind = a.length - (i / Integer.SIZE) - 1;
        int indInInt = i % Integer.SIZE;
        int mask = 1;
        mask = mask << indInInt;
        if((a[ind] & mask) != 0)
            return true;
        else
            return false;
    }

    /**
     *
     * @param val
     * @param index
     * @return
     */
    public static boolean isBitSet(int val, int index) {
        int indInInt = index % Integer.SIZE;
        int mask = 1;
        mask = mask << indInInt;
        if((val & mask) != 0)
            return true;
        else
            return false;
    }



    public static void setBit(int [] data, int i) {
        int ind = data.length - (i / Integer.SIZE) - 1;
        int indInInt = i % Integer.SIZE;
        int mask = 1;
        mask = mask << indInInt;
        data[ind] = data[ind] | mask;

    }

    public static String toStringBinary(int v) {

        return toStringBinary(v, true);
    }



    public static String toStringBinary(int v, boolean space) {
        StringBuilder sb = new StringBuilder();

        int len = Integer.SIZE;

        for (int i = 0; i < len; i++) {
            if ((v & 1) == 1) {
                sb.insert(0, "1");
            } else {
                sb.insert(0, "0");
            }
            if(space)
                sb.insert(0, " ");

            v = v >> 1;
        }

        return sb.toString().trim();
    }

}

