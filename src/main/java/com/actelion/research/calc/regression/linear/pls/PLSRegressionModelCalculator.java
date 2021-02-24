package com.actelion.research.calc.regression.linear.pls;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.regression.ARegressionMethod;
import com.actelion.research.calc.regression.ModelError;
import com.actelion.research.util.datamodel.ModelXYIndex;

/**
 * PLSRegressionModelCalculator
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
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
