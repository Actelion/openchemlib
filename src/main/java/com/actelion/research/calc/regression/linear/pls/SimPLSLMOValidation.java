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
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
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
