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
import com.actelion.research.calc.statistics.median.MedianStatisticFunctions;
import com.actelion.research.calc.statistics.median.ModelMedianDouble;
import com.actelion.research.util.datamodel.ModelXY;
import com.actelion.research.util.datamodel.ModelXYCrossValidation;
import com.actelion.research.util.datamodel.ModelXYIndex;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * 
 * 
 * SimPLSLMOValidation
 * @author Modest von Korff
 * @version 1.0
 * Jul 26, 2011 MvK: Start implementation
 */
@Deprecated // Replaced by LMOCV for multiple regression techniques
public class SimPLSLMOValidation {
	
    public static final DecimalFormat NF = new DecimalFormat("0.000");
		
	private int nRepetitions;
	
	private int nFactors;
	
	private List<ModelError> liModelErrorTest;
	
	private List<ModelError> liModelErrorTrain;
	
	private boolean centerData;

	private ModelXYCrossValidation modelXYCrossValidation;
	
	public SimPLSLMOValidation(Matrix X, Matrix Y) {
		
		modelXYCrossValidation = new ModelXYCrossValidation(new ModelXY(X, Y));
		
	}

	public void setFractionLeaveOut(double fracOut) {
		modelXYCrossValidation.setFractionLeaveOut(fracOut);
	}
	
	public void setNumRepetitions(int nRepetitions) {
		this.nRepetitions = nRepetitions;
	}
	
	public void setNumFactors(int nFactors) {
		this.nFactors = nFactors;
	}
	
	public void setCenterData(boolean centerData) {
		this.centerData = centerData;
	}

	/**
	 * 
	 * @return median error of all repetitions.
	 */
	public double calculateMedianTestError() {
		
		liModelErrorTest = new ArrayList<ModelError>();
		
		liModelErrorTrain = new ArrayList<ModelError>();
		
		PLSRegressionModelCalculator rmc = new PLSRegressionModelCalculator();
		
		rmc.setCenterData(centerData);
		
		for (int i = 0; i < nRepetitions; i++) {
			
			modelXYCrossValidation.next();

			ModelXYIndex modelXYIndex = new ModelXYIndex();

			modelXYIndex.X = modelXYCrossValidation.getXtrain();
			modelXYIndex.Y = modelXYCrossValidation.getYtrain();

			Matrix yHat = rmc.createModel(modelXYIndex);

			ModelError modelErrorTrain = ModelError.calculateError(modelXYIndex.Y, yHat);

			ModelError modelErrorTest = rmc.calculateModelErrorTest(modelXYCrossValidation.getXtest(), modelXYCrossValidation.getYtest());
			
			liModelErrorTrain.add(modelErrorTrain);
			
			liModelErrorTest.add(modelErrorTest);
			
		}
		
		List<Double> liErrorTest = ModelError.getError(liModelErrorTest);
		
		ModelMedianDouble modelMedian = MedianStatisticFunctions.getMedianForDouble(liErrorTest);
		
		return modelMedian.median;
		
	}

	
	
}
