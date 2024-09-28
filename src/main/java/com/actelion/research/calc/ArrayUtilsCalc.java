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

package com.actelion.research.calc;


import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.StringTokenizer;
import java.util.TreeSet;

import com.actelion.research.calc.statistics.median.ModelMedianDouble;
import com.actelion.research.calc.statistics.median.ModelMedianInteger;
import com.actelion.research.util.datamodel.PointDouble;

/* TODO
*  PLEASE USE ANOTHER CLASSNAME FOR THIS!
   CHRISTIAN
*/

public class ArrayUtilsCalc {

	public final static int [] cat(int [] a, int [] b) {
		int [] c = new int [a.length + b.length];
		for (int i = 0; i < a.length; i++) {
			c[i]=a[i];
		}
		for (int i = 0; i < b.length; i++) {
			c[a.length+i]=b[i];
		}
		return c;
	}


	public final static boolean contains(int [] a, int b) {
		boolean bFound = false;
		for (int i = 0; i < a.length; i++) {
			if(a[i]==b){
				bFound=true;
				break;
			}
		}
		return bFound;
	}
	public final static boolean containsAll(int [] a, int [] b) {
		boolean bFound = true;
		for (int i = 0; i < b.length; i++) {
			if(!contains(a, b[i])){
				bFound=false;
				break;
			}
		}
		return bFound;
	}

	public final static int [] copy(int [] a) {
		int [] b = new int [a.length];
		for (int i = 0; i < b.length; i++) {
			b[i]=a[i];

		}
		return b;
	}

	public final static byte [] copy(byte [] a) {
		byte [] b = new byte [a.length];
		for (int i = 0; i < b.length; i++) {
			b[i]=a[i];

		}
		return b;
	}

	/**
	 *
	 * @param li
	 * @return deep copy.
	 */
	public final static List<int[]> copyIntArray(List<int[]> li) {

		List<int[]> liC = new ArrayList<int[]>(li.size());

		for (int[] a : li) {
			int [] c = new int [a.length];
			System.arraycopy(a, 0, c, 0, a.length);
			liC.add(c);
		}

		return liC;
	}

	public final static Object [] copy(Object [] a) {
		Object b [] = new Object [a.length];
		for (int ii = 0; ii < a.length; ii++) {
			b[ii] = a[ii];
		}
		return b;
	}

	public static boolean equals(int [] a, int [] b){
		if(a.length!=b.length){
			return false;
		}
		for (int i = 0; i < b.length; i++) {
			if(a[i]!=b[i]){
				return false;
			}
		}
		return true;
	}

	public final static double [] extractCol(double [][] a, int col) {
		double b [] = new double [a.length];
		for (int ii = 0; ii < a.length; ii++) {
			b[ii] = a[ii][col];
		}
		return b;
	}

	public final static double [] filter(int [] arrData, double [] arrFilter){
		double [] arr = new double [arrData.length];

		for (int i = 0; i < arr.length; i++) {

			double val = 0;
			for (int j = 0; j < arrFilter.length; j++) {
				int indexFilter = (-arrFilter.length / 2) + j;
				int indexData = indexFilter + i;
				if(indexData >= 0 && indexData < arr.length){
					val += arrData[indexData] * arrFilter[j];
				}
			}
			arr[i]=val;
		}
		return arr;
	}

	public final static double [] filter(byte [] arrData, double [] arrFilter){
		double [] arr = new double [arrData.length];

		for (int i = 0; i < arr.length; i++) {

			double val = 0;
			for (int j = 0; j < arrFilter.length; j++) {
				int indexFilter = (-arrFilter.length / 2) + j;
				int indexData = indexFilter + i;
				if(indexData >= 0 && indexData < arr.length){
					val += arrData[indexData] * arrFilter[j];
				}
			}
			arr[i]=val;
		}
		return arr;
	}

	public final static boolean findIdentical(int [] a, int [] b) {
		boolean bFound = false;

		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b.length; j++) {
				if(a[i]==b[j]){
					bFound=true;
					break;
				}
			}
		}

		return bFound;
	}


    public static final double getCorrPearson(List<PointDouble> li) {

    	final double [] a = new double [li.size()];
    	final double [] b = new double [li.size()];
    	for (int i = 0; i < li.size(); i++) {
			a[i]=li.get(i).x;
			b[i]=li.get(i).y;
		}
        
        return getCorrPearson(a, b);
    }

    public static final double getCorrPearson(Matrix A, Matrix B) {
        
        
        final double [] a = A.toArray();
        
        final double [] b = B.toArray();
        
        
        return getCorrPearson(a, b);
    }
    
    public static final double getCorrPearson(double [] a, double [] b) {
        
        final double [] aCent = ArrayUtilsCalc.getCentered(a);
        
        final double [] aCentNorm = ArrayUtilsCalc.getNormalized(aCent);
        
        final double [] bCent = ArrayUtilsCalc.getCentered(b);
        
        final double [] bCentNorm = ArrayUtilsCalc.getNormalized(bCent);
        
        final double val = ArrayUtilsCalc.getCorrPearsonStandardized(aCentNorm,bCentNorm);
        
        return val;
    }

	public static final double getCorrPearsonStandardized(double [] a, double [] b) {

		final double covXY = getCovarianceCentered(a,b);
		
		final double varA = getVariance(a);
        
		final double varB = getVariance(b);
        
		final double val = covXY / (varA * varB);
        
        return val;
	}
	public static double [] getCentered(double [] arr) {
		double [] arrCent = new double [arr.length];
		final double mean = getMean(arr);
		for (int i = 0; i < arr.length; i++) {
			arrCent[i] = arr[i]-mean;
		}
		return arrCent;
	}

	public static double getCovarianceCentered(double [] a, double [] b) {
		
		double sum = 0;
		for (int i = 0; i < a.length; i++) {
			sum += a[i]*b[i];
		}
		
		final double covXY = sum / (a.length - 1);
		
		return covXY;
	}

	public static double getGiniCoefficient(double [] a){
		double sum=0, sumDiff=0;
		for (int i = 0; i < a.length; i++) {
			sum += a[i];
			for (int j = 0; j < a.length; j++) {
				sumDiff = Math.abs(a[i]-a[j]);
			}
		}
		double gini = sumDiff/(2* a.length * sum);
		return gini;
	}

	public static final double [] getNormalized(double [] arr) {
		
		final double [] arrNorm = new double [arr.length];
		
		double sdv = getStandardDeviation(arr);

		for (int i = 0; i < arr.length; i++) {
			arrNorm[i] = arr[i]/sdv;
		}
		
		return arrNorm;

	}

	
	public static final double getMean(double [] arr) {
		double sum = 0;
		for (int i = 0; i < arr.length; i++) {
			sum += arr[i];
		}
		return sum/arr.length;
	}

	public static final double getMean(int [] arr) {
		double sum = 0;
		for (int i = 0; i < arr.length; i++) {
			sum += arr[i];
		}
		return sum/arr.length;
	}

	public static ModelMedianDouble getMedian(double [] arr) {
		
		Arrays.sort(arr);
		
		ModelMedianDouble m = new ModelMedianDouble();
		
		m.lowerQuartile = getPercentileFromSorted(arr, 0.25);
		
		m.median = getPercentileFromSorted(arr, 0.5);
		
		m.upperQuartile = getPercentileFromSorted(arr, 0.75);
		
		m.size = arr.length;
		
		return m;
		
	}

	public static ModelMedianInteger getMedian(int [] arr) {

		Arrays.sort(arr);

		ModelMedianInteger m = new ModelMedianInteger();

		m.lowerQuartile = getPercentileFromSorted(arr, 0.25);

		m.median = getPercentileFromSorted(arr, 0.5);

		m.upperQuartile = getPercentileFromSorted(arr, 0.75);

		m.size = arr.length;

		return m;

	}



	/**
	 * 
	 * @param arr list has to be sorted in ascending order.
	 * @param fraction 0.25 lower quartile, 0,5 median and 0.75 upper quartile.
	 * @return
	 */
	public static double getPercentileFromSorted(double [] arr, double fraction) {
		
		if(arr.length==1){
			return arr[0];
		}

		double percentile=0;
		
		int len = arr.length;
		
		if(((int)(len*fraction))==(len*fraction)) {
			int index1 = (int)(len*fraction)-1;
			int index2 = index1+1;
			
			if(index1<0){
				throw new RuntimeException("Fraction to small.");
			}
			
			percentile = (arr[index1] +  arr[index2])/2.0;
			
		} else {
			int index1 = (int)(len*fraction);
			
			percentile = arr[index1];
		}
		
		return percentile;
	}

	public static int getPercentileFromSorted(int [] arr, double fraction) {

		if(arr.length==1){
			return arr[0];
		}

		int percentile=0;

		int len = arr.length;

		if(((int)(len*fraction))==(len*fraction)) {
			int index1 = (int)(len*fraction)-1;
			int index2 = index1+1;

			if(index1<0){
				throw new RuntimeException("Fraction to small.");
			}

			percentile = (int)((arr[index1] +  arr[index2])/2.0 + 0.5);

		} else {
			int index1 = (int)(len*fraction);

			percentile = arr[index1];
		}

		return percentile;
	}

	public static final double getStandardDeviation(double [] arr) {
		double sum=0;
		double mean = getMean(arr);
		for (int i = 0; i < arr.length; i++) {
			sum += (arr[i]-mean)*(arr[i]-mean);
		}
		double sdv = Math.sqrt(sum / (arr.length - 1));
		return sdv;
	}


	public static final double getVariance(double [] arr) {

		double sum=0;

		final double mean = getMean(arr);

		for (int i = 0; i < arr.length; i++) {
			sum += (arr[i]-mean)*(arr[i]-mean);
		}

		final double var = sum / (arr.length - 1);

		return var;
	}
	public static final double getVariance(int [] arr) {
		double sum=0;
		final double mean = getMean(arr);
		for (int i = 0; i < arr.length; i++) {
			sum += (arr[i]-mean)*(arr[i]-mean);
		}
		final double var = sum / (arr.length - 1);
		return var;
	}


	public final static int sum(int [] a) {
        int b = 0;
        for (int ii = 0; ii < a.length; ii++) {
          b += a[ii];
        }
        return b;
    }

	public final static int sum(long [] a) {
        int b = 0;
        for (int ii = 0; ii < a.length; ii++) {
          b += a[ii];
        }
        return b;
    }

	public final static int sum(byte [] a) {
        int b = 0;
        for (int ii = 0; ii < a.length; ii++) {
          b += a[ii];
        }
        return b;
    }

	public final static double sum(double [] a) {
        double b = 0;
        for (int ii = 0; ii < a.length; ii++) {
          b += a[ii];
        }
        return b;
    }

	/**
	 * Resize an array of Object
	 */
//	public final static Object resize(Object a, int newSize) {
//		Class cl = a.getClass();
//		if (!cl.isArray()) return null;
//		int size = Array.getLength(a);
//		Class componentType = a.getClass().getComponentType();
//		Object newArray = Array.newInstance(componentType, newSize);
//		System.arraycopy(a, 0, newArray, 0, Math.min(size, newSize));
//		return newArray;
//	}

	public final static String [] resize(String [] arr, int newSize) {
		
		String [] tmp = new String[newSize];
		int size = Math.min(arr.length, newSize);
		System.arraycopy(arr, 0, tmp, 0, size);
		return tmp;
	}

	public final static int [] resize(int [] arr, int newSize) {
		
		int [] tmp = new int[newSize];
		int size = Math.min(arr.length, newSize);
		System.arraycopy(arr, 0, tmp, 0, size);
		return tmp;
	}

	public final static byte [] resize(byte [] arr, int newSize) {
		
		byte [] tmp = new byte[newSize];
		
		int size = Math.min(arr.length, newSize);
		
		System.arraycopy(arr, 0, tmp, 0, size);
		
		return tmp;
	}
	
	public final static boolean [] resize(boolean [] arr, int newSize) {
		
		boolean [] tmp = new boolean[newSize];
		
		int size = Math.min(arr.length, newSize);
		
		System.arraycopy(arr, 0, tmp, 0, size);
		
		return tmp;
	}


	public final static double [] resize(double [] arr, int newSize) {
		
		double [] tmp = new double[newSize];
		int size = Math.min(arr.length, newSize);
		
		System.arraycopy(arr, 0, tmp, 0, size);
		
		return tmp;
	}
	
	public final static void removeDoubletsInt(List<int[]> li) {
		
		for(int i = 0; i<li.size();i++ ){
			for (int j = li.size() - 1; j > i; j--) {
				int [] a1 = li.get(i);
				int [] a2 = li.get(j);
				boolean bEq = true;
				if(a2.length != a1.length) {
					bEq = false;
					break;
				}
				for (int k = 0; k < a2.length; k++) {
					if(a1[k]!=a2[k]){
						bEq = false;
						break;
					}
				}
				if(bEq)
					li.remove(j);
			}
		}
	}
	/**
	 * Removes arrays which contains identical integer. The integer comparison 
	 * is independend from the order of the integer in the array. 
	 * @param li list with int [] as elements.
	 */
	public final static void removeDoubletsIntOrderIndepend(List<int []> li) {
		
		for(int i = 0; i<li.size();i++ ){
			for (int j = li.size() - 1; j > i; j--) {
				int [] a1 = li.get(i);
				int [] a2 = li.get(j);
				boolean bEq = true;
				
				for (int k = 0; k < a1.length; k++) {
					boolean bFound = false;
					for (int l = 0; l < a2.length; l++) {
						if(a1[k]==a2[l]){
							bFound = true;
							break;
						}
					}
					if(!bFound){
						bEq = false;
						break;
					}
				}
				if(bEq)
					li.remove(j);
			}
		}
	}

	public final static double [][] resize(double mData [][], int iNumberRowsNew, int iNumberColsNew) {
        double [][] dTmp = new double [iNumberRowsNew][iNumberColsNew];

        int iRows = iNumberRowsNew;
        int iCols = iNumberColsNew;

        if(iNumberRowsNew > mData.length)
            iRows = mData.length;
        if(iNumberColsNew > mData[0].length)
            iCols = mData[0].length;
        for (int ii = 0; ii < iRows; ii++) {
            for (int jj = 0; jj < iCols; jj++) {
                dTmp[ii][jj] = mData[ii][jj];
            }
        }
        return dTmp;
    }
	
	public final static boolean [][] resize(boolean mData [][], int rows) {
        return resize(mData, rows, mData[0].length);
    }
	
	public final static boolean [][] resize(boolean mData [][], int iNumberRowsNew, int iNumberColsNew) {
		boolean [][] dTmp = new boolean [iNumberRowsNew][iNumberColsNew];

        int iRows = iNumberRowsNew;
        int iCols = iNumberColsNew;

        if(iNumberRowsNew > mData.length)
            iRows = mData.length;
        if(iNumberColsNew > mData[0].length)
            iCols = mData[0].length;
        for (int ii = 0; ii < iRows; ii++) {
            for (int jj = 0; jj < iCols; jj++) {
                dTmp[ii][jj] = mData[ii][jj];
            }
        }
        return dTmp;
    }
	
	public final static int[] reverse(int [] arr) {
		int[] res = new int[arr.length];

		for (int i = 0; i < res.length; i++) {
			res[res.length - i - 1] = arr[i];
		}
		return res;
	}

	public final static void reverse(Object [] mArrResult) {
		for (int i = 0; i < mArrResult.length/2; i++) {
			Object res = mArrResult[mArrResult.length - i - 1];
			mArrResult[mArrResult.length - i - 1] = mArrResult[i];
			mArrResult[i] = res;
		}
	}
	
	public static List<Integer> getOverlap(int [] a1, int [] a2){
		
		TreeSet<Integer> ts = new TreeSet<Integer>();
		
		for (int i = 0; i < a1.length; i++) {
			ts.add(a1[i]);
		}
		List<Integer> li = new ArrayList<Integer>();
		for (int i = 0; i < a2.length; i++) {
			if(!ts.add(a2[i])){
				li.add(a2[i]);
			}
		}
		
		return li;
	}
	
	public static List<Integer> getUnique(int [] a1, int [] a2){
		
		TreeSet<Integer> ts = new TreeSet<Integer>();
		
		for (int i = 0; i < a1.length; i++) {
			ts.add(a1[i]);
		}
		
		for (int i = 0; i < a2.length; i++) {
			ts.add(a2[i]);
		}
		
		List<Integer> li = new ArrayList<Integer>(ts);
		
		return li;
	}

	
	public final static int [] getUnique(int [] arr) {
		
		TreeSet<Integer> ts = new TreeSet<Integer>();
		for (int i = 0; i < arr.length; i++) {
			ts.add(arr[i]);
		}
		
		int[] res = ArrayUtilsCalc.toIntArray(ts);
		
		return res;
	}

	/**
	 * Converts a List of Integer to an int[]
	 * @param list
	 * @return an array of int
	 */
	public final static int[] toIntArray(Collection<Integer> list) {
		int[] res = new int[list.size()];
		int index = 0;
		Iterator<Integer> iter = list.iterator();
		while(iter.hasNext()) {
			Integer i = (Integer) iter.next();
			res[index++] = i.intValue();
		}
		return res;
	}


	public final static String [] toArray(List<String> list) {
		String [] res = new String[list.size()];
		for (int i = 0; i < list.size(); i++) {
			res[i]=list.get(i);
		}
		return res;
	}
	
	public final static String [][] toArrayStrStr(List<List<String>> list) {
		String [][] res = new String[list.size()][list.get(0).size()];
		
		for (int i = 0; i < list.size(); i++) {
			for (int j = 0; j < list.get(i).size(); j++) {
				res[i][j] = list.get(i).get(j);
			}
		}
		return res;
	}

	public final static double[] toDoubleArray(List<Double> list) {
		double[] res = new double[list.size()];
		int index = 0;
		for (double d : list) {
			res[index++] = d;
		}
		return res;
	}
	
	
	public final static double[] toDoubleArray(int [] a) {
		double[] res = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			res[i] = a[i];
		}
		return res;
	}
	
	public final static int[] toIntArray(double [] a) {
		int[] res = new int[a.length];
		for (int i = 0; i < a.length; i++) {
			res[i] = (int)a[i];
		}
		return res;
	}


	public final static List<Integer> toList(int [] a) {
		if(a==null)
			return null;
		
		List<Integer> li = new ArrayList<Integer>(a.length);
		for (int i = 0; i < a.length; i++) {
			li.add(a[i]);
		}
		
		return li;
	}
	
	public final static List<String> toList(String [] a) {
		if(a==null)
			return null;
		
		List<String> li = new ArrayList<String>(a.length);
		for (int i = 0; i < a.length; i++) {
			li.add(a[i]);
		}
		
		return li;
	}

	public final static int indexOf(Object[] array, Object obj) {
		for (int i = 0; i < array.length; i++) {
			if(array[i].equals(obj)) return i;
		}
		return -1;
	}

	public final static int indexOf(int[] array, int obj) {
		for (int i = 0; i < array.length; i++) {
			if(array[i] == obj) return i;
		}
		return -1;
	}

	public final static int lastIndexOf(int[] array, int obj) {
		for (int i = array.length - 1; i >= 0; i--) {
			if(array[i] == obj) return i;
		}
		return -1;
	}

	public final static int lastIndexOfNot(byte [] array, int obj) {
		for (int i = array.length - 1; i >= 0; i--) {
			if(array[i] != obj) return i;
		}
		return -1;
	}
	
	public final static double min(double[] array) {
		if(array.length==0) return 0;
		double res = array[0];
		for(int i=1; i<array.length; i++) {
			res = Math.min(res, array[i]);
		}
		return res;
	}
	
	public final static int min(int[] array) {
		if(array.length==0) return 0;
		int res = array[0];
		for(int i=1; i<array.length; i++) {
			res = Math.min(res, array[i]);
		}
		return res;
	}


	public final static double min(double[][] array, int col) {
		if(array.length==0) return 0;
		double res = array[0][col];
		for(int i=1; i<array.length; i++) {
			res = Math.min(res, array[i][col]);
		}
		return res;
	}

	public final static double min(float[][] array, int col) {
		if(array.length==0) return 0;
		float res = array[0][col];
		for(int i=1; i<array.length; i++) {
			res = Math.min(res, array[i][col]);
		}
		return res;
	}

	public final static byte max(byte[] array) {
		if(array.length==0) return 0;
		byte res = array[0];
		for(int i=1; i<array.length; i++) {
			res = (byte)Math.max(res, array[i]);
		}
		return res;
	}
	public final static double max(double[] array) {
		if(array.length==0) return 0;
		double res = array[0];
		for(int i=1; i<array.length; i++) {
			res = Math.max(res, array[i]);
		}
		return res;
	}
	
	public final static double maxDouble(List<Double> array) {
		if(array.size()==0) return 0;
		double res = array.get(0);
		for(int i=1; i<array.size(); i++) {
			res = Math.max(res, array.get(i));
		}
		return res;
	}
	
	public final static int maxInt(List<Integer> array) {
		if(array.size()==0) return 0;
		int res = array.get(0);
		for(int i=1; i<array.size(); i++) {
			res = Math.max(res, array.get(i));
		}
		return res;
	}
	
	public final static double max(double[][] array, int col) {
		if(array.length==0) return 0;
		double res = array[0][col];
		for(int i=1; i<array.length; i++) {
			res = Math.max(res, array[i][col]);
		}
		return res;
	}
	public final static double max(float[][] array, int col) {
		if(array.length==0) return 0;
		float res = array[0][col];
		for(int i=1; i<array.length; i++) {
			res = Math.max(res, array[i][col]);
		}
		return res;
	}

	public final static int max(int [] array) {
		if(array.length==0) return 0;
		int res = array[0];
		for(int i=1; i<array.length; i++) {
			res = Math.max(res, array[i]);
		}
		return res;
	}

	/**
	 * Separator is ','. Vector can start with '[' and end with ']'.
	 * @param s
	 * @return
	 */
	public static int [] readIntArray(String s){
		return readIntArray(s, ",");
	}
	
	public static int [] readIntArray(String s, String seperator){
		
		s = s.replace('[', ' ');
		s = s.replace(']', ' ');
		
		s = s.trim();
		
		StringTokenizer st = new StringTokenizer(s, seperator);
		
		List<Integer> li = new ArrayList<Integer>();
		
		while(st.hasMoreTokens()){
			li.add(Integer.parseInt(st.nextToken().trim()));
		}
		
		return toIntArray(li);
	}
	
	public static double [] readDoubleArray(String s){
		return readDoubleArray(s, ",");
	}

	public static double [] readDoubleArray(String s, String seperator){
		
		s = s.replace('[', ' ');
		s = s.replace(']', ' ');
		
		s = s.trim();
		
		StringTokenizer st = new StringTokenizer(s, seperator);
		
		List<Double> li = new ArrayList<Double>();
		
		while(st.hasMoreTokens()){
			li.add(Double.parseDouble(st.nextToken().trim()));
		}
		
		return toDoubleArray(li);
	}
	
	public final static void set(int [] array, int val) {
		for(int i=0; i < array.length; i++) {
			array[i] = val;
		}
	}
	public final static void set(float [] array, float val) {
		for(int i=0; i < array.length; i++) {
			array[i] = val;
		}
	}
	public final static void set(double [] array, double val) {
		for(int i=0; i < array.length; i++) {
			array[i] = val;
		}
	}

	public final static void set(int [][] array, int val) {
		for(int i=0; i < array.length; i++)
            for(int j=0; j < array[0].length; j++)
                array[i][j] = val;
	}

	public final static void set(short [][] array, short val) {
		for(int i=0; i < array.length; i++)
            for(int j=0; j < array[0].length; j++)
                array[i][j] = val;
	}

	public final static void set(double [][] array, double val) {
		for(int i=0; i < array.length; i++)
            for(int j=0; j < array[0].length; j++)
                array[i][j] = val;
	}

	public final static void set(float [][] array, float val) {
		for(int i=0; i < array.length; i++)
            for(int j=0; j < array[0].length; j++)
                array[i][j] = val;
	}

	public final static String toStringBinary(int[] v) {
		StringBuilder sb = new StringBuilder();
		
		sb.append("[");
		for(int i=0; i<v.length; i++) {
			sb.append(toStringBinary(v[i]));
			if(i<v.length-1){
				sb.append(" ");
			}
		}
		sb.append("]");
		
		return sb.toString();
	}
	
	public static String toStringBinary(int v) {
		String str = "";

		int len = Integer.SIZE;

		for (int ii = 0; ii < len; ii++) {
			if ((v & 1) == 1) {
				str = "1 " + str;
			} else {
				str = "0 " + str;
			}
			v = v >> 1;
		}

		return str.trim();
	}

	public final static String toString(int[] v) {
		
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for(int i=0; i<v.length; i++) {
			sb.append((i>0?",":"") + v[i] );
		}
		sb.append("]");
		return sb.toString();
	}
	
	public final static String toString(Collection<Integer> li) {
		
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		int cc=0;
		for (int i : li) {
			sb.append((cc>0?",":"") +i );
			cc++;
		}
		sb.append("]");
		return sb.toString();
	}
	
	public final static String toStringNoBrackets(Collection<Integer> li, String sep) {
		
		StringBuilder sb = new StringBuilder();
		int cc=0;
		for (int i : li) {
			sb.append((cc>0?sep:"") +i );
			cc++;
		}
		return sb.toString();
	}
	
	/**
	 * Writes a list into a string with line terminators.
	 * @param li
	 * @param step number of numbers in one line
	 * @return
	 */
	public final static String toStringIntegerList(List<Integer> li, int step) {
		
		StringBuilder sb = new StringBuilder();
		sb = new StringBuilder();
		for (int i = 0; i < li.size(); i+=step) {
			if(i+step>li.size())
				step = li.size() - i;
			for (int j = i; j < i+step; j++) {
				sb.append(li.get(j));
				if(j<i+step-1){
					sb.append(" ");
				}
			}
			sb.append("\n");
		}
		return sb.toString();
	}
	
	public final static String toStringLongList(List<Long> li, int step) {
		
		StringBuilder sb = new StringBuilder();
		sb = new StringBuilder();
		for (int i = 0; i < li.size(); i+=step) {
			if(i+step>li.size())
				step = li.size() - i;
			for (int j = i; j < i+step; j++) {
				sb.append(li.get(j));
				if(j<i+step-1){
					sb.append(" ");
				}
			}
			sb.append("\n");
		}
		return sb.toString();
	}
	
	public final static String [] toStringArray(List<Integer> li) {
		String [] a = new String [li.size()];
		for (int i = 0; i < li.size(); i++) {
			a[i]=Integer.toString(li.get(i));
		}
		return a;
	}
	
	
	public final static String toString(byte[] v) {
		
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for(int i=0; i<v.length; i++) {
			
			int val = v[i] & 0xFF;
			
			sb.append((i>0?", ":"") + val );
		}
		sb.append("]");
		return sb.toString();
	}

	public final static String toStringPure(int[] v) {
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<v.length; i++) 
			sb.append((i>0?" ":"") + v[i] );
		return sb.toString();
	}

	public final static String toStringIntArrays(List<int []> li) {
		StringBuilder buff = new StringBuilder();
		for (int[] element : li) {
			buff.append(toString(element) + "\n");
		}
		return buff.toString();
	}

	public final static String toString(int [][] v) {
		StringBuilder sb = new StringBuilder();
		for(int i=0; i< v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                sb.append( (j > 0 ? "," : "") + v[i][j]);
            }
            sb.append("\n");
        }
		return sb.toString();
	}
	
	public final static String toStringFormatted(int [] arrTop, int [] arrBottom) {
		int [][] v = new int [2][];
		
		v[0]=arrTop;
		v[1]=arrBottom;
		
		return toStringFormatted(v);
	}
	
	public final static String toStringFormatted(int [][] v) {
		
		
		int maxAbs = 0;
		for(int i=0; i< v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
            	if(Math.abs(v[i][j])>maxAbs){
            		maxAbs = Math.abs(v[i][j]);
            	}
            }
        }
		
		int len = Integer.toString(maxAbs).length()+1;

		StringBuilder sb = new StringBuilder();
		
		for(int i=0; i< v.length; i++) {
			
            for (int j = 0; j < v[0].length; j++) {
            	
            	StringBuilder sbVal = new StringBuilder(Integer.toString(v[i][j]));
            	
            	while(sbVal.length()<len){
            		sbVal.insert(0, " ");
            	}
            	
            	sb.append(sbVal);
            	
            	if(j<v[0].length-1){
            		sb.append(" ");
            	}
            }
            
            sb.append("\n");
            
        }
		
		return sb.toString();
	}

	public final static String toString(double[] v) {
		StringBuilder sb = new StringBuilder("[");
		for(int i=0; i<v.length; i++) {
			sb.append((i>0?",":"") + v[i]);
		}
        sb.append("]");
		return sb.toString();
	}

	public final static String toStringPure(double[] v) {
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<v.length; i++) {
			sb.append((i>0?"\t":"") + v[i]);
		}
		return sb.toString();
	}

	public final static String toString(double[] v, NumberFormat nf) {
		StringBuilder sb = new StringBuilder("[");
		for(int i=0; i<v.length; i++) {
			sb.append((i>0?",":"") + nf.format(v[i]));
		}
        sb.append("]");
		return sb.toString();
	}

	public final static String toString(double [][] v) {
		StringBuilder res = new StringBuilder();
		for(int i=0; i< v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                res.append( (j > 0 ? "," : "") + v[i][j]);
            }
            res.append("\n");
        }
		return res.toString();
	}
	
    public final static String toString(double[] v, int iDigits) {
        
        String sFormat = "";

        sFormat += "0";
        int iCounter = 0;
        if(iDigits > 0)
            sFormat += ".";

        while(iCounter < iDigits) {
          sFormat += "0";
          iCounter++;
        }
        NumberFormat nf = new DecimalFormat(sFormat);
        StringBuilder sb = new StringBuilder();
		for(int i=0; i<v.length; i++) {
			sb.append((i>0?",":"") + nf.format(v[i])) ;
		}
		
		return sb.toString();

    }
    public final static String toString(float[] v, int iDigits) {
        
        String sFormat = "";

        sFormat += "0";
        int iCounter = 0;
        if(iDigits > 0)
            sFormat += ".";

        while(iCounter < iDigits) {
          sFormat += "0";
          iCounter++;
        }
        NumberFormat nf = new DecimalFormat(sFormat);
        StringBuilder sb = new StringBuilder();
		for(int i=0; i<v.length; i++) {
			sb.append((i>0?",":"") + nf.format(v[i])) ;
		}
		
		return sb.toString();

    }
        
    public final static String toString(double[][] v, int iDigits) {
        
        String sFormat = "";

        sFormat += "0";
        int iCounter = 0;
        if(iDigits > 0)
            sFormat += ".";

        while(iCounter < iDigits) {
          sFormat += "0";
          iCounter++;
        }
        NumberFormat nf = new DecimalFormat(sFormat);
        StringBuilder res = new StringBuilder();
		for(int i=0; i<v.length; i++) {
			for(int j=0; j<v[i].length; j++) {
				res.append((j>0?" ":"") + nf.format(v[i][j])) ;
			}
			res.append((i<v.length-1?"\n":""));
		}
		
		return res.toString();

    }

    public final static String toStringNoDigits(double[] v) {

        NumberFormat nf = new DecimalFormat("0");
		String res = "";
		for(int i=0; i<v.length; i++) {
			res += (i>0?",":"") + nf.format(v[i]) ;
		}
		return res + "";
	}

	public final static String toString(Object[] v) {
		String res = "[";
		for(int i=0; i<v.length; i++) {
			res += (i>0?",":"") + v[i] ;
		}
		return res + "]";
	}

	public final static void shift(int[] v, int n) {
		int[] copy = new int[v.length];
		for(int i=0; i<v.length; i++) copy[i] = v[(i+n+v.length)%v.length];
		System.arraycopy(copy, 0, v, 0, v.length);
	}
	
	
	public final static List<Integer> parseInteger(String s, String sep) {
		
    	StringTokenizer st = new StringTokenizer(s, sep);
    	    	
    	List<Integer> li = new ArrayList<Integer>();
    	
    	while(st.hasMoreTokens()){
    		int v = Integer.parseInt(st.nextToken());
    		li.add(v);
    	}
    	
    	return li;
	}
	
}
