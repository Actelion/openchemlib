/*
* Copyright (c) 1997 - 2015
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.List;
import java.util.StringTokenizer;




@SuppressWarnings({"unchecked", "rawtypes"})
public class CompareUtils {


	public static <T> boolean contains(T[] array, T object) {
		for (T t : array) {
			if(t==object) return true;
		}
		return false;
	}

	
//	private static final Pattern pattern;
//	static {
//		pattern = Pattern.compile("[\\_\\/\\,\\-\\.\\:\n\\s]");
//	}
	
//	private static List<String> split(String s){
//		List<String> res = new ArrayList<String>();
//		StringTokenizer st = new StringTokenizer(s, "_/,-.:\n\t ", true);
//		while(st.hasMoreTokens()) {
//			res.add(st.nextToken());
//		}
//		return res;
//	}
	/**
	 * Compare 2 strings by splitting the chains into blocks and comparing each blocks individually. Each Block is then compared as integer or string
	 * Example
	 * Dose 10 ->  [Dose] [10]
	 * Dose 100 -> [Dose] [100]
	 * 
	 * Caveat:
	 * The sorting does not consider "-" as negative but as a separator, so the sorting goes:
	 * _1, -1, -2, -10, 1, 2, 10
	 * 
	 * @param o1
	 * @param o2
	 * @return
	 */
	public static final int compare(String o1, String o2) {
		if(o1==null && o2==null) return 0; 
		if(o1==null) return 1; 
		if(o2==null) return -1;
		
		StringTokenizer st1 = new StringTokenizer(o1, "_/,-.:\n\t ", true);
		StringTokenizer st2 = new StringTokenizer(o2, "_/,-.:\n\t ", true);

		String s1, s2;
		boolean allDigits1, allDigits2;
		int c;
		while(st1.hasMoreTokens() && st2.hasMoreTokens()) {
			s1 = st1.nextToken();
			s2 = st2.nextToken();
			
			if(s1.length()==0) {
				if(s2.length()==0) continue;
				return -1;
			}
			if(s2.length()==0) return 1;
			
			//Compare first the numeric value if possible
			allDigits1 = true;
			allDigits2 = true;
			for(int j=0; allDigits1 && j<s1.length();j++) {
				if(!Character.isDigit(s1.charAt(j))) allDigits1 = false;
			}
			for(int j=0; allDigits2 && j<s2.length();j++) {
				if(!Character.isDigit(s2.charAt(j))) allDigits2 = false;
			}
			if(allDigits1 && allDigits2) {
				c = fastIntValueOf(s1) - fastIntValueOf(s2);
				if(c!=0) return c;
			} else if(allDigits1 && !allDigits2) {
				return -1;
			} else if(!allDigits1 && allDigits2) {
				return 1;
			} else {
				//Compare by string first case insensitive
				c = s1.compareToIgnoreCase(s2);
				if(c!=0) return c;
			}

		}

		if(!st1.hasMoreTokens() && st2.hasMoreTokens()) return -1;
		if(st1.hasMoreTokens() && !st2.hasMoreTokens()) return 1;
		
		return o1.compareTo(o2);
	}

	public static int fastIntValueOf( String str ) {
	    int ival = 0;
	    for(int i=0; i<str.length(); i++) {
	    	ival = ival*10 + (str.charAt(i)-'0');
	    }
	    return ival;
	}
	
	public static int compare(Object[] a1, Object[] a2) {
		for (int i = 0; i < a1.length || i < a2.length; i++) {
			Object o1 = i < a1.length? a1[i]: null;
			Object o2 = i < a2.length? a2[i]: null;
			int c = CompareUtils.compare(o1, o2);
			if(c!=0) return c;
		}
		return 0;
	}
	public static int compare(Object o1, Object o2) {
		if(o1==null && o2==null) return 0; 
		if(o1==null) return 1; //Null at the end
		if(o2==null) return -1;

		if((o1 instanceof String) && (o2 instanceof String)) {
			return compare((String) o1, (String) o2);
		} else if((o1 instanceof Object[]) && (o2 instanceof Object[])) {
			return compare((Object[]) o1, (Object[]) o2);
		} else if((o1 instanceof Comparable) && (o2 instanceof Comparable)) {			
			return ((Comparable) o1).compareTo((Comparable) o2);
		} else {
			return compare(o1.toString(), o2.toString());			
		}
	}
	
	public static final int compareAsDate(Object o1, Object o2) {
		if(o1==null && o2==null) return 0; 
		if(o1==null) return 1; //Null at the end
		if(o2==null) return -1;
		
		Date d1 = Formatter.parseDateTime(o1.toString());
		Date d2 = Formatter.parseDateTime(o2.toString());
		if(d1!=null && d2!=null) {
			return d1.compareTo(d2);
		} else if(d1==null && d2!=null) {
			return 1;
		} else if(d1!=null && d2==null) {
			return -1;
		} else {
			return compare(o1, o2);
		}
	}

	
	/**
	 * Comparator to allow null values
	 */
	public static final Comparator<Object> OBJECT_COMPARATOR = new Comparator<Object>() {
		@Override
		public int compare(Object o1, Object o2) {
			if(o1==null && o2==null) return 0; 
			if(o1==null) return 1; //Null at the end
			if(o2==null) return -1;
			return CompareUtils.compare(o1.toString(), o2.toString());
		}
	};
	
	
	
	/**
	 * Comparator for strings, comparing blocks separately so that the order becomes SLIDE-1, SLIDE-2, SLIDE-10, SLIDE-10-1, ...
	 */
	public static final Comparator<Object> STRING_COMPARATOR = new Comparator<Object>() {
		@Override
		public int compare(Object o1, Object o2) {
			return CompareUtils.compare(o1, o2);
		}
	};
	
	/**
	 * Comparator for dates, allowing formats such as yyyy, mm.yyyy, ...
	 */
	public static final Comparator<Object> DATE_COMPARATOR = new Comparator<Object>() {		
		@Override
		public int compare(Object o1, Object o2) {
			return CompareUtils.compareAsDate(o1, o2);
		}
	};

	public static boolean equals(Object obj1, Object obj2) {
		return obj1==null? obj2==null: obj1.equals(obj2);
			
	}

	/**
	 * Test Speed
	 * @param args
	 */
	public static void main(String[] args) {
//		List<String> initial = Arrays.asList(new String[] {"heart", "", "lung", "lung/left", "1", "10", "2", "3", "11", "Box1","Box10", "Box2", "Box 2", "Box 1", "Box 10", "10.9.2012", "11.9.12", "2012", "Genomics", "_2", "_3", "_10", " 1", "_1", "-10", "-1. 0", "-2", "-1.00", "-1. 00", "-1.  00", "1-1", "1-10", "Proteomics", "Clinical Analysis", "Lung/Right", "Lung", "Lung/Left", "1", "-1", "-1.1", "-1.10", "2.A", "10.B-1", "10.B-10", "10.B-2", "1.C", "21.D", "3.B-3", "Heart","11.C", "11. C", "10. B-1", "2.C", "1.B","d030; Heart","d030; Heart/Left ventricle + Septum","d030; Heart/Right ventricle"});
		List<String> initial = Arrays.asList(new String[] {"1", "4", "2B", "2A", "3A", "3B", "a", "A", "b", "B", "c", "C", "1", "11", "12", "2", "21", "22"});
		List<String> l = new ArrayList<String>();
		l.addAll(initial);
//		l.addAll(initial);
//		l.addAll(initial);
//		l.addAll(initial);
//		l.addAll(initial);
//		l.addAll(initial);
//		l.addAll(initial);
//		l.addAll(initial);
		
		
		List<String> l2 = new ArrayList<String>(l);
		long st = System.currentTimeMillis();
		Collections.sort(l2);
		long t1 = System.currentTimeMillis()-st;
		
		List<String> l3 = new ArrayList<String>(l);
		st = System.currentTimeMillis();
		Collections.sort(l3, DATE_COMPARATOR);
		long t2 = System.currentTimeMillis()-st;
		
		st = System.currentTimeMillis();
		Collections.sort(l, STRING_COMPARATOR);
		long t3 = System.currentTimeMillis()-st;
		
		for (String string : l) {
			System.out.println(string);
		}
		
		System.out.println();
		System.out.println("CompareUtils.normal: "+t1);
		System.out.println("CompareUtils.date: "+t2);
		System.out.println("CompareUtils.string: "+t3);
	}


	
}
