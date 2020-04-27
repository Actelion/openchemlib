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

/**
 * 
 */
package com.actelion.research.util.datamodel;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Comparator;
import java.util.List;

import com.actelion.research.util.ConstantsDWAR;
import com.actelion.research.util.Formatter;

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