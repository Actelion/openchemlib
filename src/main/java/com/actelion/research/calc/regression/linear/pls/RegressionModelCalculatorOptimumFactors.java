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
 */

package com.actelion.research.calc.regression.linear.pls;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.regression.ModelError;
import com.actelion.research.util.datamodel.ModelXYIndex;

/**
 * RegressionModelCalculatorOptimumFactors
 * @author Modest von Korff
 * @version 1.0
 * Aug 14, 2015 MvK Start implementation
 */
public class RegressionModelCalculatorOptimumFactors {

	private static final double FRACTION_LEAVE_OUT = 0.25;
	
	private static final int LMO_REPETITIONS = 7;
	
	
	private boolean centerData;
	
	private Matrix B;
		
	private Matrix YHat;
	
	private int factorsMin;
	
	/**
	 * 
	 */
	public RegressionModelCalculatorOptimumFactors() {
		centerData = false;
	}
	
	
	
	/**
	 * @param centerData the centerData to set
	 */
	public void setCenterData(boolean centerData) {
		this.centerData = centerData;
	}


	/**
	 * Calculates the PLS regression model for the given data set. 
	 * A Leave Multiple Out estimator is used to assess the optimum number of factors.
	 * The calculation starts with factorsStart factor up to factorsEnd.
	 * @param dataXYTrain
	 * @param factorsStart
	 * @param factorsEnd
	 * @return
	 */
	public ModelError calculateModel(ModelXYIndex dataXYTrain, int factorsStart, int factorsEnd){
		
		
		SimPLSLMOValidation simPLSLMOValidation = new SimPLSLMOValidation(dataXYTrain.X, dataXYTrain.Y);
				
		simPLSLMOValidation.setFractionLeaveOut(FRACTION_LEAVE_OUT);
		
		simPLSLMOValidation.setNumRepetitions(LMO_REPETITIONS);
		
		double errTestMin = Integer.MAX_VALUE;
		
		factorsMin = 1;
		
		for (int i = factorsStart; i < factorsEnd+1; i++) {
			
			simPLSLMOValidation.setNumFactors(i);
			
			double errTest = simPLSLMOValidation.calculateMedianTestError();
			
			if(errTest < errTestMin){
				errTestMin = errTest;
				factorsMin = i;
			}
			
		}
		
		System.out.println("Optimum number of factors " + factorsMin);
		
		PLSRegressionModelCalculator rmc = new PLSRegressionModelCalculator();
		
		rmc.setCenterData(centerData);

		rmc.setFactors(factorsMin);

		Matrix yHat = rmc.createModel(dataXYTrain);

		ModelError modelErrorTrain = ModelError.calculateError(dataXYTrain.Y, yHat);


		B = rmc.getB();
		
		YHat = rmc.getYHat();
		       
        return modelErrorTrain;
	}

	/**
	 * @return the b
	 */
	public Matrix getB() {
		return B;
	}


	/**
	 * @return the yHat
	 */
	public Matrix getYHat() {
		return YHat;
	}



	/**
	 * @return the factorsMin
	 */
	public int getFactorsMin() {
		return factorsMin;
	}

	
	
}
