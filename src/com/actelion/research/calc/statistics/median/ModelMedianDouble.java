/*
 * Project: DD_core
 * @(#)ModelMedianDouble.java
 *
 * Copyright (c) 1997- 2015
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * Author: MvK
 */

package com.actelion.research.calc.statistics.median;

import java.util.List;

import com.actelion.research.util.Formatter;
import com.actelion.research.util.StringFunctions;


public class ModelMedianDouble {
	
	private static final int LENGTH = 7;
	
	public double lowerQuartile;
	
	public double median;
	
	public double upperQuartile;
	
	public int id;
	
	public int size;
	
	
	/**
	 * 
	 */
	public ModelMedianDouble() {
		// TODO Auto-generated constructor stub
	}
	
	/**
	 * @param lowerQuartile
	 * @param median
	 * @param upperQuartile
	 * @param id
	 * @param size
	 */
	public ModelMedianDouble(double lowerQuartile, double median, double upperQuartile, int id, int size) {
		super();
		this.lowerQuartile = lowerQuartile;
		this.median = median;
		this.upperQuartile = upperQuartile;
		this.id = id;
		this.size = size;
	}
	
	public ModelMedianDouble(ModelMedianDouble m) {
		super();
		this.lowerQuartile = m.lowerQuartile;
		this.median = m.median;
		this.upperQuartile = m.upperQuartile;
		this.id = m.id;
		this.size = m.size;
	}



	public double range(){
		return upperQuartile-lowerQuartile;
	}
	
	
	public String toString() {
		StringBuilder sb = new StringBuilder(); 
		
		sb.append(Formatter.format4(lowerQuartile));
		sb.append("\t");
		sb.append(Formatter.format4(median));
		sb.append("\t");
		sb.append(Formatter.format4(upperQuartile));
		
		return sb.toString();
	}
	
	public static String toString(List<ModelMedianDouble> liModelMedian) {
		
		StringBuilder sb = new StringBuilder(); 
		
		int length = LENGTH;

		int lengthText = 15;
		
		sb.append(StringFunctions.format2DefinedLengthTrailing("Id", lengthText));
		for (int i = 0; i < liModelMedian.size(); i++) {
			ModelMedianDouble m = liModelMedian.get(i);
			
			String s = StringFunctions.format2DefinedLengthLeading(Integer.toString(m.id), length);
			
			sb.append(s);
		}
		
		sb.append("\n");
		
		sb.append(StringFunctions.format2DefinedLengthTrailing("Upper quartile", lengthText));
		for (int i = 0; i < liModelMedian.size(); i++) {
			ModelMedianDouble m = liModelMedian.get(i);
			
			String s = StringFunctions.format2DefinedLengthLeading(Formatter.format2(m.upperQuartile), length);
			
			sb.append(s);
		}
		
		sb.append("\n");
		
		sb.append(StringFunctions.format2DefinedLengthTrailing("Median", lengthText));
		for (int i = 0; i < liModelMedian.size(); i++) {
			ModelMedianDouble m = liModelMedian.get(i);
			
			String s = StringFunctions.format2DefinedLengthLeading(Formatter.format2(m.median), length);
			
			sb.append(s);
		}
		
		sb.append("\n");
		
		sb.append(StringFunctions.format2DefinedLengthTrailing("Lower quartile", lengthText));
		for (int i = 0; i < liModelMedian.size(); i++) {
			ModelMedianDouble m = liModelMedian.get(i);
			
			String s = StringFunctions.format2DefinedLengthLeading(Formatter.format2(m.lowerQuartile), length);
			
			sb.append(s);
		}
		
		sb.append("\n");
		
		sb.append(StringFunctions.format2DefinedLengthTrailing("Size", lengthText));
		for (int i = 0; i < liModelMedian.size(); i++) {
			ModelMedianDouble m = liModelMedian.get(i);
			
			String s = StringFunctions.format2DefinedLengthLeading(Integer.toString(m.size), length);
			
			sb.append(s);
		}
		
		return sb.toString();
	}

}
