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