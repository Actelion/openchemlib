/*
 * Copyright (c) 2018.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
 */

package com.actelion.research.calc.statistics.median;

import com.actelion.research.util.Formatter;
import com.actelion.research.util.StringFunctions;

import java.util.List;

public class ModelMedianFloat {

	private static final int LENGTH = 7;

	public float lowerQuartile;

	public float median;

	public float upperQuartile;

	public int id;

	public int size;


	/**
	 *
	 */
	public ModelMedianFloat() {

	}

	/**
	 * @param lowerQuartile
	 * @param median
	 * @param upperQuartile
	 * @param id
	 * @param size
	 */
	public ModelMedianFloat(float lowerQuartile, float median, float upperQuartile, int id, int size) {
		super();
		this.lowerQuartile = lowerQuartile;
		this.median = median;
		this.upperQuartile = upperQuartile;
		this.id = id;
		this.size = size;
	}

	public ModelMedianFloat(ModelMedianFloat m) {
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
		
		sb.append(Formatter.format4((double)lowerQuartile));
		sb.append("\t");
		sb.append(Formatter.format4((double)median));
		sb.append("\t");
		sb.append(Formatter.format4((double)upperQuartile));
		
		return sb.toString();
	}
	
	public static String toString(List<ModelMedianFloat> liModelMedian) {
		
		StringBuilder sb = new StringBuilder(); 
		
		int length = LENGTH;

		int lengthText = 15;
		
		sb.append(StringFunctions.format2DefinedLengthTrailing("Id", lengthText));
		for (int i = 0; i < liModelMedian.size(); i++) {
			ModelMedianFloat m = liModelMedian.get(i);
			
			String s = StringFunctions.format2DefinedLengthLeading(Integer.toString(m.id), length);
			
			sb.append(s);
		}
		
		sb.append("\n");
		
		sb.append(StringFunctions.format2DefinedLengthTrailing("Upper quartile", lengthText));
		for (int i = 0; i < liModelMedian.size(); i++) {
			ModelMedianFloat m = liModelMedian.get(i);
			
			String s = StringFunctions.format2DefinedLengthLeading(Formatter.format2((double)m.upperQuartile), length);
			
			sb.append(s);
		}
		
		sb.append("\n");
		
		sb.append(StringFunctions.format2DefinedLengthTrailing("Median", lengthText));
		for (int i = 0; i < liModelMedian.size(); i++) {
			ModelMedianFloat m = liModelMedian.get(i);
			
			String s = StringFunctions.format2DefinedLengthLeading(Formatter.format2((double)m.median), length);
			
			sb.append(s);
		}
		
		sb.append("\n");
		
		sb.append(StringFunctions.format2DefinedLengthTrailing("Lower quartile", lengthText));
		for (int i = 0; i < liModelMedian.size(); i++) {
			ModelMedianFloat m = liModelMedian.get(i);
			
			String s = StringFunctions.format2DefinedLengthLeading(Formatter.format2((double)m.lowerQuartile), length);
			
			sb.append(s);
		}
		
		sb.append("\n");
		
		sb.append(StringFunctions.format2DefinedLengthTrailing("Size", lengthText));
		for (int i = 0; i < liModelMedian.size(); i++) {
			ModelMedianFloat m = liModelMedian.get(i);
			
			String s = StringFunctions.format2DefinedLengthLeading(Integer.toString(m.size), length);
			
			sb.append(s);
		}
		
		return sb.toString();
	}

}
