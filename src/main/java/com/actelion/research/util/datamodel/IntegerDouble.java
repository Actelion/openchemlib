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
 * @author Modest v. Korff
 */

/**
 * 
 */
package com.actelion.research.util.datamodel;

import com.actelion.research.util.ConstantsDWAR;
import com.actelion.research.util.Formatter;

import java.text.NumberFormat;
import java.util.Comparator;
import java.util.List;

/**
 * IntegerDouble
 * Oct 24, 2011 MvK: Start implementation
 */
public class IntegerDouble {
	
	private int iv;
	
	private double dv;
	
	public IntegerDouble() {
		
	}
	
	/**
	 * @param iv
	 * @param dv
	 */
	public IntegerDouble(int iv, double dv) {
		super();
		this.iv = iv;
		this.dv = dv;
	}
	
	public IntegerDouble(IntegerDouble id) {
		this.iv = id.iv;
		this.dv = id.dv;
	}
	
	
	public int hashCode() {
		return iv;
	}
	
	public int getInt() {
		return iv;
	}
	
	public double getDouble() {
		return dv;
	}
		
	public void setInteger(int iv) {
		this.iv = iv;
	}
	
	public void setDouble(double dv) {
		this.dv = dv;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(iv);
		sb.append("\t");
		sb.append(Formatter.format3(dv));
		
		return sb.toString();
	}

	public static String toStringValues(List<IntegerDouble> li, NumberFormat df) {

		StringBuilder sb = new StringBuilder();

		for (IntegerDouble id : li) {
			if(sb.length()>0){
				sb.append(ConstantsDWAR.SEP_VALUE);
			}
			sb.append(df.format(id.getDouble()));
		}

		return sb.toString();
	}

	public static String toStringValues(IntegerDouble [] arr, NumberFormat df) {

		StringBuilder sb = new StringBuilder();

		for (IntegerDouble id : arr) {
			if(sb.length()>0){
				sb.append(ConstantsDWAR.SEP_VALUE);
			}
			sb.append(df.format(id.getDouble()));
		}

		return sb.toString();
	}


	
	public static Comparator<IntegerDouble> getComparatorDouble(){
		
		return new Comparator<IntegerDouble>() {
			
			public int compare(IntegerDouble id1, IntegerDouble id2) {
				
				if(id1.dv>id2.dv){
					return 1;
				}else if(id1.dv<id2.dv){
					return -1;
				}
				
				return 0;
			}
		};
	}
	
	public static Comparator<IntegerDouble> getComparatorInt(){

		return new Comparator<IntegerDouble>() {

			public int compare(IntegerDouble id1, IntegerDouble id2) {

				if(id1.iv>id2.iv){
					return 1;
				}else if(id1.iv<id2.iv){
					return -1;
				}

				return 0;
			}
		};
	}




}