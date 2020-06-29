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
package com.actelion.research.chem.interactionstatistics;

import com.actelion.research.util.ArrayUtils;
import com.actelion.research.util.FastSpline;

/**
 * Class used to represent a Protein Ligand Function
 * 
 */
public class SplineFunction {
	
	private int[] occurencesArray;
	private FastSpline spline, derivate;	
	private double[] discreteFunction;


	public int[] getOccurencesArray() {
		return occurencesArray;
	}
	
	public void setOccurencesArray(int[] occurencesArray) {
		this.occurencesArray = occurencesArray;
	}
	
	public void setDiscreteFunction(double[] discreteFunction) {
		this.discreteFunction = discreteFunction;
	}
	
	public void setSplineFunction(FastSpline spline) {
		this.spline = spline;
		this.derivate = spline.derivative();
	}
	
	public int getTotalOccurences() {
		return ArrayUtils.sum(occurencesArray);		
	}
	
	public double[] getFGValue(double v) {
		if(spline==null) return new double[]{0.0, 0.0};
		try{		
			double value = spline.value(v);
			double dev = derivate.value(v);
			return new double[]{value, dev};
		} catch (ArrayIndexOutOfBoundsException e) {
			return new double[]{0.0, 0.0};
		}
	}
	
	public double getDiscreteValue(double d) {
		int index = (int) (0.5+d/InteractionDistanceStatistics.BIN_SIZE);
		return discreteFunction[index];
	}
	
	public FastSpline getSpline() {
		return spline;
	}


}
