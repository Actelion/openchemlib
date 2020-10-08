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

import com.actelion.research.calc.geometry.Triangle;
import com.actelion.research.calc.regression.linear.simple.LinearRegression;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.ConstantsDWAR;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


public class SummaryFragments {

	public static final String TAG_SUM_UNIQUE_FRAGMENTS_UNTIL_MAX = "SumUniqueFragmentsUntilMaxFrags";

	public static final String TAG_NUM_BONDS_AT_MAXIMUM_FREQ = "NumBondsAtMaxFreq";

	public static final String TAG_NUM_BONDS_AT_MAXIMUM_SLOPE = "NumBondsAtMaxSlope";

	public static final String TAG_MAXIMUM_SLOPE = "MaxSlope";

	public static final String TAG_RATIO_USED_BONDS_UP_TO_MAX_FREQ = "RatioUsedBondsUpToMaxFreq";

	public static final String TAG_RATIO_NON_SYMMETRIC_ATOMS = "RatioNonSymmetricAtoms";

    public static final String TAG_COMPLEXITY = "Complexity";



	private static final int DELTA_SLOPE = 1;

	private int nAtomsMolecule;
	
	private int nBondsMolecule;

	private StereoMolecule mol;

	private int sumFragsUntilMaximumFrequency;

	private int numBondsAtMaximumFrequency;

	private int numBondsAtMaxSlope;

	private double maximumSlope;

	private double ratioUsedBondsUpToMaxFreq;

	private double ratioNonSymmetricAtoms;

	private double score;

	/**
	 * 
	 */
	public SummaryFragments(ResultFragmentsStatistic resultFragmentsStatistic) {
				
		this.nAtomsMolecule=resultFragmentsStatistic.getAtoms();
		
		this.nBondsMolecule=resultFragmentsStatistic.getBonds();
		
		mol = resultFragmentsStatistic.getMol();

		summary(resultFragmentsStatistic.getExhaustiveStatistics());

	}
	
	/**
	 * Working version 16.04.2012
	 * Good results for the five sample datasets.
	 * @param liModelExhaustiveStatistics
	 */
	private void summary(List<ModelExhaustiveStatistics> liModelExhaustiveStatistics){

		int delta = DELTA_SLOPE;

		Collections.sort(liModelExhaustiveStatistics, ModelExhaustiveStatistics.getComparatorNumBonds());

		double [] arrUniqueFragsLn = new double[liModelExhaustiveStatistics.size()];

		int [] arrBondsInFrag = new int[liModelExhaustiveStatistics.size()];

		double maxValue = 0;
		numBondsAtMaximumFrequency = -1;
		int indexEndInclusively = -1;

		sumFragsUntilMaximumFrequency = 0;

		for (int i=0; i<liModelExhaustiveStatistics.size();i++) {
			
			ModelExhaustiveStatistics model = liModelExhaustiveStatistics.get(i);

			double uniqueFragsLn = Math.log(model.getUnique());

			arrUniqueFragsLn[i] = uniqueFragsLn;

			arrBondsInFrag[i] = model.getNumBondsInFragment();

			if(arrUniqueFragsLn[i]>maxValue){

				maxValue=arrUniqueFragsLn[i];

				numBondsAtMaximumFrequency=arrBondsInFrag[i];

				indexEndInclusively = i;

				sumFragsUntilMaximumFrequency += model.getUnique();
			}
		}

		double areaUntilMaxValue = getIntegratedAreaWeighted(arrBondsInFrag, arrUniqueFragsLn, indexEndInclusively);

		int start = delta;
		int end = liModelExhaustiveStatistics.size() - delta;

		double [] arrSlope = new double[liModelExhaustiveStatistics.size()];

		maximumSlope = -10000;
		numBondsAtMaxSlope = 0;
		for (int i = start; i < end; i++) {

			int startSlide = i - delta;
			int endSlide = i + delta+1;

			LinearRegression linearRegression = new LinearRegression();

			for (int j = startSlide; j < endSlide; j++) {

				linearRegression.addPoint(j, arrUniqueFragsLn[j]);
			}

			linearRegression.regress();

			double slope = linearRegression.getSlope();

			arrSlope[i] = slope;

			if(slope>maximumSlope){
				numBondsAtMaxSlope=arrBondsInFrag[i];
				maximumSlope=slope;
			}
		}


		ratioUsedBondsUpToMaxFreq = (numBondsAtMaximumFrequency)/(double)nBondsMolecule;

		ratioNonSymmetricAtoms = 1.0 - SymmetryCalculator.getRatioSymmetricAtoms(mol);

		score = (areaUntilMaxValue) * (ratioUsedBondsUpToMaxFreq) * (maximumSlope) * (ratioNonSymmetricAtoms);


	}

	public int getAtomsMolecule() {
		return nAtomsMolecule;
	}

	public int getBondsMolecule() {
		return nBondsMolecule;
	}

	public StereoMolecule getMol() {
		return mol;
	}

	public int getNumBondsAtMaximumFrequency() {
		return numBondsAtMaximumFrequency;
	}

	public double getMaximumSlope() {
		return maximumSlope;
	}

	public int getNumBondsAtMaxSlope() {
		return numBondsAtMaxSlope;
	}

	public double getRatioUsedBondsUpToMaxFreq() {
		return ratioUsedBondsUpToMaxFreq;
	}

	public double getRatioNonSymmetricAtoms() {
		return ratioNonSymmetricAtoms;
	}

    public int getSumFragsUntilMaximumFrequency() {
        return sumFragsUntilMaximumFrequency;
    }

    public double getComplexityScore() {
		return score;
	}

	/**
	 * Takes the natural logarithm for each area fraction.
	 * @param arrXBondsInFrag
	 * @param arrY
	 * @param indexEndInclusively
	 * @return
	 */
	public static double getIntegratedAreaWeighted(int [] arrXBondsInFrag, double [] arrY, int indexEndInclusively) {

		double areaWeighted = 0;

		int end = indexEndInclusively+1;

		for (int i = 1; i < end; i++) {

			double a = arrXBondsInFrag[i]-arrXBondsInFrag[i-1];
			double b = arrY[i]-arrY[i-1];

			double areaTriangle = Triangle.getAreaRightTriangle(a,b);

			double areaRectangle = arrY[i-1] * a;

			double area = areaTriangle+areaRectangle;

			if(area<1){
				area=1;
			}

			areaWeighted += Math.log(area);

		}

		return areaWeighted;
	}

	public static List<String> getHeaderTags(){

		List<String> li = new ArrayList<String>();

        li.add(ConstantsDWAR.TAG_ATOMS);

        li.add(ConstantsDWAR.TAG_BONDS);

        li.add(TAG_SUM_UNIQUE_FRAGMENTS_UNTIL_MAX);
        li.add(TAG_NUM_BONDS_AT_MAXIMUM_FREQ);
        li.add(TAG_NUM_BONDS_AT_MAXIMUM_SLOPE);
        li.add(TAG_MAXIMUM_SLOPE);
        li.add(TAG_RATIO_USED_BONDS_UP_TO_MAX_FREQ);
        li.add(TAG_RATIO_NON_SYMMETRIC_ATOMS);
        li.add(TAG_COMPLEXITY);

		return li;
	}
}