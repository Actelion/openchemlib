/*
 * Copyright (c) 2020.
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

package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.util.Formatter;

import java.util.Comparator;

public class PPNodeVizTriangle implements Comparable<PPNodeVizTriangle> {
	
	private PPNodeViz n1;
	
	private PPNodeViz n2;
	
	private PPNodeViz n3;
	
	public double area;
	
	/**
	 * @param n1
	 * @param n2
	 * @param n3
	 */
	public PPNodeVizTriangle(PPNodeViz n1, PPNodeViz n2, PPNodeViz n3) {
		super();
		this.n1 = n1;
		this.n2 = n2;
		this.n3 = n3;

		this.area = -1;
	}
	
	/**
	 * @return the area
	 */
	public double getArea() {
		return area;
	}
	
	/**
	 * @param area the area to set
	 */
	public void setArea(double area) {
		this.area = area;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	
	public int compareTo(PPNodeVizTriangle ta) {
		
		if(area > ta.area){
			return 1;
		} 
		
		if(area < ta.area){
			return -1;
		}
		
		if(area == ta.area){
			return 0;
		}
				
		return 0;
	};

	
	
	public PPNodeViz getNodeA(){
		return n1; 
	}
	
	public PPNodeViz getNodeB(){
		return n2; 
	}
	
	public PPNodeViz getNodeC(){
		return n3; 
	}
	
	public Coordinates getCoordinatesA(){
		return n1.getCoordinates();
	}
	
	public Coordinates getCoordinatesB(){
		return n2.getCoordinates();
	}
	
	public Coordinates getCoordinatesC(){
		return n3.getCoordinates();
	}
	
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("area=");
		builder.append(Formatter.format3(area));
		return builder.toString();
	}

	public static Comparator<PPNodeVizTriangle> getComparatorArea(){
		
		Comparator<PPNodeVizTriangle> comp = new Comparator<PPNodeVizTriangle>() {
			
		
			public int compare(PPNodeVizTriangle ta1, PPNodeVizTriangle ta2) {
				
				if(ta1.area > ta2.area){
					return 1;
				} 
				
				if(ta1.area < ta2.area){
					return -1;
				}
				
				return 0;

			}
		};
		
		return comp;
		
		
	}
	
	
}