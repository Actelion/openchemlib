package com.actelion.research.chem.descriptor.flexophore.completegraphmatcher;

import com.actelion.research.calc.filter.SlidingWindow;
import com.actelion.research.chem.descriptor.flexophore.IMolDistHist;
import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;

/**
 * 
 * 
 * HistogramMatchCalculator
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Oct 2, 2012 MvK: Start implementation
 * May 15 2013 MvK: Heavy bug detected. Wrong similarity results. reset() added.
 * Mar 01 2016 MvK sliding filter added.
 */
public class HistogramMatchCalculator {
	
	private static final boolean DEBUG = false;

	// private static final double [] FILTER = {1};

	// private static final double [] FILTER = {0.125,0.75,0.125};

	// Filter 07.04.2020
	private static final double [] FILTER = {0.25,0.5,0.25};

	// private static final double [] FILTER = {0.05, 0.1, 0.2, 0.3, 0.2, 0.1, 0.5};


	
	private byte [] arrTmpHistMol;

	private int [] arrTmpHistMolBlurred;

	private byte [] arrTmpHistFrag;

	private int [] arrTmpHistFragBlurred;

	private SlidingWindow slidingWindow;

	public HistogramMatchCalculator() {
		
		arrTmpHistMol =  new byte[ConstantsFlexophoreGenerator.BINS_HISTOGRAM];
		
		arrTmpHistMolBlurred =  new int[ConstantsFlexophoreGenerator.BINS_HISTOGRAM];
		
		arrTmpHistFrag =  new byte[ConstantsFlexophoreGenerator.BINS_HISTOGRAM];
		
		arrTmpHistFragBlurred =  new int[ConstantsFlexophoreGenerator.BINS_HISTOGRAM];

		slidingWindow = new SlidingWindow(FILTER);
	}
	
	private void reset(){
		for (int i = 0; i < ConstantsFlexophoreGenerator.BINS_HISTOGRAM; i++) {
			arrTmpHistMol[i]=0;
			arrTmpHistMolBlurred[i]=0;
			arrTmpHistFrag[i]=0;
			arrTmpHistFragBlurred[i]=0;
		}
	}


	/**
	 *
	 * @param query
	 * @param indexQueryPPPoint1 index for pharmacophore point 1 in the query
	 * @param indexQueryPPPoint2 index for pharmacophore point 2 in the query
	 * @param base
	 * @param indexBasePPPoint1 index for pharmacophore point 1 in the base
	 * @param indexBasePPPoint2 index for pharmacophore point 2 in the base
	 * @return
	 */
	public double getSimilarity(IMolDistHist query, int indexQueryPPPoint1, int indexQueryPPPoint2, IMolDistHist base, int indexBasePPPoint1, int indexBasePPPoint2) {
		
		double sc = 0;
		
		byte [] arr1 = query.getDistHist(indexQueryPPPoint1, indexQueryPPPoint2, arrTmpHistMol);
		
		byte [] arr2 = base.getDistHist(indexBasePPPoint1, indexBasePPPoint2, arrTmpHistFrag);

		byte [] arr1Filt = slidingWindow.filter(arr1);

		byte [] arr2Filt = slidingWindow.filter(arr2);

		sc = getSimilarity(arr1Filt, arr2Filt);

		reset();
		
		return sc;

	}

	private void blurrHistogram(byte [] arr, int [] arrBlurred){

		slidingWindow.filter(arr);

	}


	
	/**
	 * 
	 * @param arr1
	 * @param arr2
	 * @return Similarity value between 0 and 1. 1: histograms are identical.
	 */
	private double getSimilarity(int [] arr1, int [] arr2) {

		if(DEBUG){
			System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			System.out.println("!!!HistogramMatchCalculator DEBUG mode!!!");
			System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			
			return 1.0;
		}

		double sumOverlap = 0;

		double sum = 0;

		for (int i = 0; i < arr1.length; i++) {
			
			sumOverlap += Math.min(arr1[i], arr2[i]); 
			
		}

		double score = sumOverlap / (sum/2.0);

		if(score>1)
			score = 1.0;
		
		return score;
	}

	private double getSimilarity(byte [] arr1, byte [] arr2) {

		if(DEBUG){
			System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			System.out.println("!!!HistogramMatchCalculator DEBUG mode!!!");
			System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");

			return 1.0;
		}

		double sumOverlap = 0;

		double sum = 0;

		for (int i = 0; i < arr1.length; i++) {

			sumOverlap += Math.min(arr1[i], arr2[i]);

			sum += arr1[i] + arr2[i];

		}

		double score = sumOverlap / (sum/2.0);

		if(score>1)
			score = 1.0;

		return score;
	}


}
