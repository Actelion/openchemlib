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

package com.actelion.research.calc.regression.linear.simple;

import java.awt.Point;
import java.util.List;
import java.util.Vector;

import com.actelion.research.util.datamodel.DoubleArray;
import com.actelion.research.util.datamodel.PointDouble;

/**
 * LinearRegression
 * 2009 MvK: Start implementation
 */
public class LinearRegression {
	
	private double intercept; // Intercept
	
	private double slope; // Slope
	
	//vector of points
	private Vector<PointDouble> v;
	
	//vector of residuals
	private Vector<PointDouble> residuals;

	private double xMean;
	
	private double yMean;
	
	/**
	 * 
	 */
	public LinearRegression() {
		v = new Vector<PointDouble>();  	
		residuals = new Vector<PointDouble>();   
		xMean = 0;
		yMean = 0;
	}
	
	/**  
	*  add a point to the vector of points
	* @param p the point to be added
	*/

	public void addPoint(Point p){
		v.addElement(new PointDouble(p));
	}
	
	public void addPoint(PointDouble p){
		v.addElement(p);
	}
	
	public void addPoint(double x, double y){
		v.addElement(new PointDouble(x, y));
	}


	public void clear(){
		v = new Vector<PointDouble>();
	}

	public List<PointDouble> getValues(){
		return v;
	}
	
	public double [][] getValuesAsArray(){
		
		double [][] arr = new double [2][v.size()];
				
		for (int i = 0; i < v.size(); i++) {
			arr[0][i]=v.get(i).x;
			arr[1][i]=v.get(i).y;
		}
		
		return arr;
	}
	
	public DoubleArray getValuesAsArrayX(){
		
		DoubleArray arr = new DoubleArray(v.size());
				
		for (int i = 0; i < v.size(); i++) {
			arr.add(v.get(i).x);
		}
		
		return arr;
	}
	
	public DoubleArray getValuesAsArrayY(){
		
		DoubleArray arr = new DoubleArray(v.size());
				
		for (int i = 0; i < v.size(); i++) {
			arr.add(v.get(i).y);
		}
		
		return arr;
	}
	
	
	public Vector<PointDouble> regress() {
		
		xMean = 0;
		yMean = 0;
		
		for (PointDouble p : v) {
			xMean += p.x;
			yMean += p.y;
		}
		
		xMean /= v.size();
		yMean /= v.size();

		double sxy2=0;
		double sxx2=0;
		for (PointDouble p : v) {
			sxy2 += (p.x-xMean)*(p.y-yMean);
			sxx2 += (p.x-xMean)*(p.x-xMean);
		}
		
		slope = sxy2/sxx2;
		
		intercept = yMean - slope * xMean;
		
		Vector<PointDouble> resid = new Vector<PointDouble>();
		PointDouble q;
		
		for (PointDouble p : v) {
			double currentResidual = p.y-(intercept+slope*p.x);
			q = new PointDouble(p.x, (int)currentResidual);
			resid.addElement(q);	
		}
		
		return resid;

	}
		
	public void calculate() {
		if (v.size() >1) {
			residuals=regress();
		}
	}
	
	public Vector<PointDouble> getResiduals(){
		return residuals;
	}
	
	public PointDouble getMean(){
		return new PointDouble(xMean+0.5, yMean+0.5);
	}


	public double getIntercept() {
		return intercept;
	}


	public double getSlope() {
		return slope;
	}


	public double getXMean() {
		return xMean;
	}


	public double getYMean() {
		return yMean;
	}
	
	public double getY(double x) {
		return intercept + slope*x;
	}


	public void setSlope(double slope) {
		this.slope = slope;
	}
	
}
