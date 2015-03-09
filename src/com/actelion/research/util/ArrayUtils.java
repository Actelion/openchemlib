/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Joel Freyss
 */

package com.actelion.research.util;

import java.lang.reflect.Array;
import java.util.*;

public class ArrayUtils {

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
	public final static int[] toIntArray(List list) {
		int[] res = new int[list.size()];
		int index = 0;
		Iterator iter = list.iterator();
		while(iter.hasNext()) {
			Integer i = (Integer) iter.next();
			res[index++] = i.intValue();
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
		for (int i = 0; i < array.length; i++) {
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
}
