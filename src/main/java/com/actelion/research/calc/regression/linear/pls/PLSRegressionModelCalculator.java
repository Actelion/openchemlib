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
import com.actelion.research.calc.regression.ARegressionMethod;
import com.actelion.research.calc.regression.ModelError;
import com.actelion.research.util.datamodel.ModelXYIndex;

/**
 * PLSRegressionModelCalculator
 * @author Modest von Korff
 * @version 1.0
 * Aug 14, 2015 MvK Start implementation
 */
public class PLSRegressionModelCalculator extends ARegressionMethod<ParameterPLS> {

	public static final int FACTORS = 15;

	private SimPLS simPLS;

	private Matrix B;
	
	private Matrix Xvar;
	
	private Matrix YHat;

	private Matrix X, Y;

	private Matrix XtrainPreprocessed;
	
	private Matrix YtrainPreprocessed;

	/**
	 * 
	 */
	public PLSRegressionModelCalculator() {
		setParameterRegressionMethod(new ParameterPLS(FACTORS));
	}

	public PLSRegressionModelCalculator(ParameterPLS parameterPLS) {
		setParameterRegressionMethod(parameterPLS);
	}

	/**
	 * @param centerData the centerData to set
	 */
	public void setCenterData(boolean centerData) {
		getParameter().setCenterData(centerData);
	}


	public void setFactors(int factors){
		getParameter().setFactors(factors);
	}

	/**
	 *
	 * @param dataXYTrain
	 * @return Yhat from the train data
	 */
	public Matrix createModel(ModelXYIndex dataXYTrain){

		X = dataXYTrain.X;
		Y = dataXYTrain.Y;

		XtrainPreprocessed = dataXYTrain.X;
		
		YtrainPreprocessed = dataXYTrain.Y;

		
		if(getParameter().isCenterData()){
			
			XtrainPreprocessed = dataXYTrain.X.getCenteredMatrix();
			
			YtrainPreprocessed = dataXYTrain.Y.getCenteredMatrix();
			
			// System.out.println("Calculate PLS with centered data.");
			
		} else {
			// System.out.println("Calculate PLS with raw data.");
		}
		
		simPLS = new SimPLS();
		
		simPLS.simPlsSave(XtrainPreprocessed, YtrainPreprocessed, getParameter().getFactors());
		
		Matrix R = simPLS.getR();

		if(R.cols() == 1 && R.rows() == 1 && R.get(0,0)==0){
			System.out.println("RegressionModelCalculator R = 0.");
		}

		Matrix Q = simPLS.getQ();

		B = R.multiply(false, true, Q);
		
		Xvar = XtrainPreprocessed.getVarianceCols();
								
        YHat = SimPLS.invLinReg_Yhat(B, X, dataXYTrain.X, Y);
        
        return YHat;
	}

	/**
	 * With centering of Xtest with Xtrain.
	 * @param Xtest
	 * @return
	 */
	public Matrix calculateYHat(Matrix Xtest){
		
		Matrix YHatTest = SimPLS.invLinReg_Yhat(B, X, Xtest, Y);
		
		return YHatTest;
	}

	@Override
	public double calculateYHat(double[] arrRow) {
		Matrix YHatTest = SimPLS.invLinReg_Yhat(B, X, new Matrix(true, arrRow), Y);
		return YHatTest.get(0,0);
	}

	public double calculateYHat(byte[] arrRow) {
		Matrix YHatTest = SimPLS.invLinReg_Yhat(B, X, new Matrix(true, arrRow), Y);
		return YHatTest.get(0,0);
	}

	public double calculateYHat(int[] arrRow) {
		Matrix YHatTest = SimPLS.invLinReg_Yhat(B, X, new Matrix(true, arrRow), Y);
		return YHatTest.get(0,0);
	}

	public Matrix calculateYHatWithoutDeCentering(Matrix Xtest){

		Matrix YHatTest = SimPLS.invLinReg_Yhat(B, Xtest);

		return YHatTest;
	}

	public ModelError calculateModelErrorTest(Matrix Xtest, Matrix Ytest){
		
		Matrix YHatTest = SimPLS.invLinReg_Yhat(B, XtrainPreprocessed, Xtest, YtrainPreprocessed);
		
		return ModelError.calculateError(Ytest, YHatTest);
	}
	
	
	/**
	 * @return the b
	 */
	public Matrix getB() {
		return B;
	}



	/**
	 * @return the xvar
	 */
	public Matrix getXvar() {
		return Xvar;
	}



	/**
	 * @return the yHat
	 */
	public Matrix getYHat() {
		return YHat;
	}


	public Matrix getT(Matrix XPreprocessed) {
		return simPLS.getT(XPreprocessed);
	}

	
	
}
