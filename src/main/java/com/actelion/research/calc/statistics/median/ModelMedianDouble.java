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
