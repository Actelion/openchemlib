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

package com.actelion.research.chem.properties.complexity;

import com.actelion.research.calc.CorrelationCalculator;
import com.actelion.research.calc.Logarithm;
import com.actelion.research.calc.regression.linear.simple.LinearRegression;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.datamodel.DoubleArray;

import java.util.List;


public class ObjectiveExhaustiveStatistics {


	public static boolean VERBOSE = false;
	
	public static boolean SCALE = false;

//	private static final int MIN_NUM_REGRESSION_POINTS = 5;
//	
//	private static final int NUM_BONDS_REGRESSION_STARTS = 3; 
//	
//	private static final int DELTA_BONDS_UP = 2; 
	
	
	public static final int NUM_BONDS_START_REGRESSION_SMALL_MOLECULE = 1; 
	
	public static final int DELTA_BONDS = 1; 
	
	public static final double BASE_LOG = Math.E;
	
	public static final double ALMOST_ZERO = 1E-6; 
	
	
	/**
	 * 
	 */
	public ObjectiveExhaustiveStatistics() {
		
	}
	
//	/**
//	 * 03.09.2013
//	 * Test function 
//	 * @param resultFragmentsStatistic
//	 * @return
//	 */
//	public ResultObjective calculateScore(ResultFragmentsStatistic resultFragmentsStatistic){
//		
//		ResultObjective resultObjective = new ResultObjective();
//		
//		List<ModelExhaustiveStatistics> liModelExhaustiveStatistics = resultFragmentsStatistic.getExhaustiveStatistics();
//				
//		int nQuarterBonds = (int)(resultFragmentsStatistic.getBonds() / 4.0);
//				
//		double slopeUnique = getMaxSlopeAroundPointOfInflection(liModelExhaustiveStatistics, nQuarterBonds);
//						
//		double ratio = SymmetryCalculator.getRatioNonSymmetricAtoms(resultFragmentsStatistic.getMol());
//
//		double score = slopeUnique * (1.0-ratio);
//		
//		resultObjective.setScore(score);
//		
//		
//		// System.out.println("ObjectiveExhaustiveStatistics slopeUnique " +  Formatter.format3(slopeUnique) + " slopeCorrectionFactor " +  Formatter.format3(slopeCorrectionFactor));					
//
//		
//		return resultObjective;
//	}
//	
//	private static double getMaxSlopeAroundPointOfInflection(List<ModelExhaustiveStatistics> liModelExhaustiveStatistics, int startBondSearch){
//		
//		int indexStart = 0;
//		for (int i=0; i < liModelExhaustiveStatistics.size(); i++) {
//			
//			ModelExhaustiveStatistics model = liModelExhaustiveStatistics.get(i);
//			
//			if(model.getNumBondsInFragment()==startBondSearch){
//				indexStart = i;
//				break;
//			}
//			
//		}
//		
//		
//		// int [] arrOffsetIndex = {0, -1, 1, -2, 2, -3, 3, -4, 4};
//		int [] arrOffsetIndex = {-2, -1, 0, 1, 2};
//
//		double slopeMax = 0;
//		
//		for (int i = 0; i < arrOffsetIndex.length; i++) {
//			
//			LinearRegression linearRegressionUnique = new LinearRegression();
//			
//			int index = indexStart+arrOffsetIndex[i];
//			
//			if(index-1<0){
//				continue;
//			}
//			if(index+1>liModelExhaustiveStatistics.size()-1){
//				break;
//			}
//			
//			
//			ModelExhaustiveStatistics modelLow = liModelExhaustiveStatistics.get(index-1);
//			double vl = Math.log(modelLow.getUnique());
//			
//			ModelExhaustiveStatistics model = liModelExhaustiveStatistics.get(index);
//			double v = Math.log(model.getUnique());
//			
//			ModelExhaustiveStatistics modelHigh = liModelExhaustiveStatistics.get(index+1);
//			
//			double vh = Math.log(modelHigh.getUnique());
//
//			
//			linearRegressionUnique.addPoint(modelLow.getNumBondsInFragment(), vl);
//			linearRegressionUnique.addPoint(model.getNumBondsInFragment(), v);
//			linearRegressionUnique.addPoint(modelHigh.getNumBondsInFragment(), vh);
//			
//			linearRegressionUnique.calculate();
//			
//			double slopeUnique = linearRegressionUnique.getSlope();
//			
//			System.out.println("bonds " + model.getNumBondsInFragment() + " slope " + Formatter.format3(slopeUnique));
//			
//			if(slopeUnique>slopeMax){
//				slopeMax=slopeUnique;
//			}
//
//		}
//		
//		return slopeMax;
//		
//	}
	
	
	
	
	
	/**
	 * 28.08.2013
	 * @param resultFragmentsStatistic
	 * @return
	 */
	public ResultObjective calculateScore(ResultFragmentsStatistic resultFragmentsStatistic){
		
		ResultObjective resultObjective = new ResultObjective();
		
		List<ModelExhaustiveStatistics> liModelExhaustiveStatistics = resultFragmentsStatistic.getExhaustiveStatistics();
					
		
		
		int nQuarterBonds = (int)(resultFragmentsStatistic.getBonds() / 4.0);
		
		int numBondsRegressionStarts = nQuarterBonds - DELTA_BONDS;
		
		int numBondsRegressionStops = nQuarterBonds + DELTA_BONDS;
		
			
		LinearRegression linearRegressionUnique = getLinearRegressionUnique(liModelExhaustiveStatistics, numBondsRegressionStarts, numBondsRegressionStops);
		
		
		double slopeUnique = linearRegressionUnique.getSlope();
				
		DoubleArray arrX =  linearRegressionUnique.getValuesAsArrayX();
		
		DoubleArray arrY =  linearRegressionUnique.getValuesAsArrayY();
		
		double r = new CorrelationCalculator().calculateCorrelation(arrX, arrY, CorrelationCalculator.TYPE_BRAVAIS_PEARSON);
		
		double r2 = r*r;
		
		resultObjective.setSlopeR2(r2);
		
		resultObjective.setSlope(slopeUnique);
		
		resultObjective.setNumRegressionPoints(linearRegressionUnique.getValues().size());
		
		
		double ratio = SymmetryCalculator.getRatioSymmetricAtoms(resultFragmentsStatistic.getMol());

		double score = slopeUnique * (1.0-ratio);
		
		resultObjective.setScore(score);
		
		
		// System.out.println("ObjectiveExhaustiveStatistics slopeUnique " +  Formatter.format3(slopeUnique) + " slopeCorrectionFactor " +  Formatter.format3(slopeCorrectionFactor));					

		
		return resultObjective;
	}
	
//	/**
//	 * 21.11.2014 
//	 * Just the highest difference in fragments.
//	 * @param resultFragmentsStatistic
//	 * @return
//	 */
//	public ResultObjective calculateScore(ResultFragmentsStatistic resultFragmentsStatistic){
//		
//		ResultObjective resultObjective = new ResultObjective();
//		
//		List<ModelExhaustiveStatistics> liModelExhaustiveStatistics = resultFragmentsStatistic.getExhaustiveStatistics();
//					
//		double maxSlope = 0;
//		for (int i = 1; i < liModelExhaustiveStatistics.size(); i++) {
//			int frags = liModelExhaustiveStatistics.get(i-1).getFragments();
//			int frags1 = liModelExhaustiveStatistics.get(i).getFragments();
//			
//			double slope = ((double)frags1)/frags;
//			
//			if(slope>maxSlope){
//				maxSlope=slope;
//			}
//		}
//		
//		resultObjective.setScore(maxSlope);
//		
//		return resultObjective;
//	}
	
	
	private static LinearRegression getLinearRegressionUnique(List<ModelExhaustiveStatistics> liModelExhaustiveStatistics, int numBondsRegressionStarts, int numBondsRegressionStops){
		
		LinearRegression linearRegressionUnique = new LinearRegression();
				
		for (int i=0; i < liModelExhaustiveStatistics.size(); i++) {
						
			ModelExhaustiveStatistics model = liModelExhaustiveStatistics.get(i);

			if(model.getNumBondsInFragment() < numBondsRegressionStarts)
				continue;
			
			if(model.getNumBondsInFragment() > numBondsRegressionStops)
				break;
			
			double yUnique = Logarithm.get(model.getUnique(), BASE_LOG);
			
			linearRegressionUnique.addPoint(model.getNumBondsInFragment(), yUnique);
			
		}
		
		linearRegressionUnique.calculate();
		
		return linearRegressionUnique;
	}
	
	public static int getNeededNumberOfBondsInFragment(StereoMolecule mol){
		
		int b = mol.getBonds();
		
		int n = b / 4 + DELTA_BONDS;
		
		return n;
	}
	
	
	

}
