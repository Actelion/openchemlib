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

package com.actelion.research.util;

import java.lang.reflect.Array;
import java.text.DecimalFormat;
import java.util.*;

public class ArrayUtils {


	public static boolean equals(int [] a, int [] b){

		boolean eq = true;

		if(a==null && b==null){
			return true;
		} else if(a==null){
			return false;
		} else if(b==null){
			return false;
		}

		if(a.length!=b.length){
			return false;
		}

		for (int i = 0; i < a.length; i++) {
			if(a[i]!=b[i]){
				eq=false;
				break;
			}
		}

		return eq;
	}

	/**
	 * Resize an array 
	 */
	public final static Object resize(Object a, int newSize) {
		Class cl = a.getClass();
		if (!cl.isArray()) return null;
		int size = Array.getLength(a);
		Class componentType = a.getClass().getComponentType();
		Object newArray = Array.newInstance(componentType, newSize);
		System.arraycopy(a, 0, newArray, 0, Math.min(size, newSize));
		return newArray;
	}
	
	/**
	 * Resize an array of Object
	 */
	public final static double[] cut(double a[], int off, int len) {
		double[] res = new double[a.length-len];
		for(int i=0; i<off; i++) {
			res[i] = a[i];
		}
		for(int i=off; i<res.length; i++) {
			res[i] = a[i+len];			
		}
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
		Iterator iter = list.iterator();
		while(iter.hasNext()) {
			Integer i = (Integer) iter.next();
			res[index++] = i.intValue();
		}
		return res;
	}

	public final static<T> List<T> toList(T [] arr) {
		List<T> li = new ArrayList<T>(arr.length);

		for (T t : arr) {
			li.add(t);
		}

		return li;
	}

	public final static List<String> toList(String [] arr) {
		List<String> li = new ArrayList<>(arr.length);
		for (String t : arr) {
			li.add(t);
		}
		return li;
	}

	public final static List<Integer> toList(int [] arr) {
		List<Integer> li = new ArrayList<>(arr.length);
		for (int t : arr) {
			li.add(t);
		}
		return li;
	}

	public final static void toList(int [] arr, List<Integer> li) {
		for (int t : arr) {
			li.add(t);
		}
	}


	public final static String[] toStringArray(List<String> list) {
		String[] res = new String[list.size()];
		int index = 0;
		Iterator<String> iter = list.iterator();
		while(iter.hasNext()) {
			res[index++] = iter.next();
		}
		return res;
	}

	public final static byte [] toByteArray(List<Byte> li) {
		byte[] res = new byte[li.size()];
		int index = 0;
		Iterator<Byte> iter = li.iterator();
		while(iter.hasNext()) {
			res[index++] = iter.next();
		}
		return res;
	}

	public final static int indexOf(Object[] array, Object obj) {
		for (int i = 0; i < array.length; i++) {
			if(array[i].equals(obj)) return i;
		}
		return -1;		
	}

	public final static int indexOf(int[] array, int obj) {
		return indexOf(array, obj, 0);
	}

	public final static int indexOf(int[] array, int obj, int startIndex) {
		for (int i = startIndex; i < array.length; i++) {
			if(array[i] == obj) return i;
		}
		return -1;		
	}

	public final static int indexOf(byte[] array, byte obj) {
		return indexOf(array, obj, 0);
	}

	public final static int indexOf(byte[] array, byte obj, int startIndex) {
		for (int i = startIndex; i < array.length; i++) {
			if(array[i] == obj) return i;
		}
		return -1;
	}

	public final static int sum(int[] array) {
		int res = 0;
		for(int i=0; i<array.length; i++) {
			res += array[i];  
		}
		return res;
	}
	
	public final static double sum(double[] array) {
		double res = 0;
		for(int i=0; i<array.length; i++) {
			res += array[i];  
		}
		return res;
	}
	
	public final static double min(double[] array) {
		if(array.length==0) return 0;
		double res = array[0];
		for(int i=1; i<array.length; i++) {
			res = Math.min(res, array[i]);  
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

	public final static int max(int[] array) {
		if(array.length==0) return 0;
		int res = array[0];
		for(int i=1; i<array.length; i++) {
			res = Math.max(res, array[i]);
		}
		return res;
	}
	public final static int max(int[] [] array) {
		if(array.length==0) return 0;

		int r = array.length;
		int c = array[0].length;

		int res = array[0][0];
		for(int i=0; i<r; i++) {
			for(int j=0; j<c; j++) {
				res = Math.max(res, array[i][j]);
			}
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

	public final static float max(float [] array) {

		if(array.length==0)
			return 0;

		float res = array[0];

		for(int i=1; i<array.length; i++) {
			res = Math.max(res, array[i]);
		}

		return res;
	}

	public final static String toString(int[] v) {
		String res = "[";
		for(int i=0; i<v.length; i++) {
			res += (i>0?", ":"") + v[i] ;			 
		}
		return res + "]";
	}
	
	public final static String toString(byte[] v) {
		String res = "[";
		for(int i=0; i<v.length; i++) {
			res += (i>0?", ":"") + v[i] ;			 
		}
		return res + "]";
	}
	
	public final static String toString(double[] v) {
		String res = "[";
		for(int i=0; i<v.length; i++) {
			res += (i>0?", ":"") + v[i] ;			 
		}
		return res + "]";
	}

	public final static String toString(float[] v) {
		String res = "[";
		for(int i=0; i<v.length; i++) {
			res += (i>0?", ":"") + v[i] ;
		}
		return res + "]";
	}

	public final static String toString(double[] v, DecimalFormat df) {
		StringBuilder sb =  new StringBuilder();
		sb.append("[");

		for(int i=0; i<v.length; i++) {
			sb.append(df.format(v[i]));
			if(i<v.length-1){
				sb.append(", ");
			}
		}
		sb.append("]");
		return sb.toString();
	}

	public final static String toString(Object[] v) {
		String res = "[";
		for(int i=0; i<v.length; i++) {
			res += (i>0?", ":"") + v[i] ;			 
		}
		return res + "]";
	}
	
	public final static void shift(int[] v, int n) {
		int[] copy = new int[v.length];
		for(int i=0; i<v.length; i++) copy[i] = v[(i+n+v.length)%v.length];
		System.arraycopy(copy, 0, v, 0, v.length);
	}

	public static void shuffle(int [] arr, Random random){

		int cycles = 7;

		int size = arr.length;

		for (int i = 0; i < cycles; i++) {

			for (int j = 0; j < size; j++) {
				int dest = random.nextInt(size);

				if(dest==j){
					continue;
				}

				int v = arr[j];

				arr[j] = arr[dest];

				arr[dest] = v;
			}
		}
	}

	public static void shuffle(int [] arr){

		Random rnd = new Random();

		shuffle(arr, rnd);
	}

	public static <T> void shuffle(T [] arr){

		Random rnd = new Random();

		int cycles = 7;

		int size = arr.length;

		for (int i = 0; i < cycles; i++) {

			for (int j = 0; j < size; j++) {
				int dest = rnd.nextInt(size);

				if(dest==j){
					continue;
				}

				T v = arr[j];

				arr[j] = arr[dest];

				arr[dest] = v;
			}
		}
	}

	
	/**
	 * Copy an array 
	 */
	public final static Object copy(Object a) {
		Class cl = a.getClass();
		if (!cl.isArray()) return null;
		int size = Array.getLength(a);
		Class componentType = a.getClass().getComponentType();
		Object newArray = Array.newInstance(componentType, size);
		System.arraycopy(a, 0, newArray, 0, size);
		return newArray;
	}
	
	public final static boolean contains(List<int[]> list, int[] arr) {
		for (int[] arr2: list ) {
			if(arr.length!=arr2.length) continue;
			for (int i = 0; i < arr2.length; i++) if(arr2[i]!=arr[i]) continue;			
			return true;
		}		
		return false;	
	}

	public static void set(boolean [] a, boolean v){
		for (int i = 0; i < a.length; i++) {
			a[i]=v;
		}
	}

	public final static boolean isOverlap(int[] a1, int[] a2) {
		boolean ov = false;
		all:
		for (int v1 : a1) {
			for (int v2 : a2) {
				if(v1==v2){
					ov=true;
					break all;
				}
			}
		}
		return ov;
	}

	/**
	 *
	 * @param s string with single digits '123456789'.
	 * @return
	 */
	public static int [] parseSingleDigitString(String s) {
		if (s == null)
			return null;
		int[] arr = new int[s.length()];
		for (int i = 0; i < s.length(); i++) {
			int c = Integer.parseInt(Character.toString(s.charAt(i)));
			arr[i] = c;
		}
		return arr;
	}
}
