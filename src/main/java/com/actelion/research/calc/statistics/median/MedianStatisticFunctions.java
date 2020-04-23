/*
 * Copyright (c) 2018.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
 */

package com.actelion.research.calc.statistics.median;

import java.util.Collections;
import java.util.List;

/**
 * MedianStatisticFunctions
 * Aug 5, 2011 MvK: Start implementation
 */
public class MedianStatisticFunctions {
	
	/**
	 * 
	 * @param liScore list has to be sorted in ascending order.
	 * @param fraction 0.25 lower quartile, 0,5 median and 0.75 upper quartile.
	 * @return
	 */
	public static double getPercentileFromSorted(List<Double> liScore, double fraction) {
		
		if(liScore.size()==1){
			return liScore.get(0);
		}

		double percentile=0;
		
		int len = liScore.size();
		
		if(((int)(len*fraction))==(len*fraction)) {
			int index1 = (int)(len*fraction)-1;
			int index2 = index1+1;
			
			if(index1<0){
				throw new RuntimeException("Fraction to small.");
			}
			
			percentile = (liScore.get(index1)+liScore.get(index2))/2.0;
			
		} else {
			int index1 = (int)(len*fraction);
			
			percentile = liScore.get(index1);
		}
		
		return percentile;
	}
	
	public static double getPercentileFromSorted(double [] arr, double fraction) {
		
		return getPercentileFromSorted(arr, fraction, 0, arr.length);
	}
	
	public static double getPercentileFromSorted(float [] arr, double fraction) {

		return getPercentileFromSorted(arr, fraction, 0, arr.length);
	}

	public static double getPercentileFromSorted(double [] arr, double fraction, int indexStart, int length) {
		
		if(arr.length==1){
			return arr[0];
		}

		double percentile=0;
				
		if(((int)(length*fraction))==(length*fraction)) {
			
			int index1 = (int)(length*fraction)-1 + indexStart;
			
			int index2 = index1+1 + indexStart;
			
			if(index1<0){

				throw new RuntimeException("Fraction to small.");
			}
			
			percentile = (arr[index1]+arr[index2])/2.0;
			
		} else {
			int index1 = (int)(length*fraction) + indexStart;
			
			percentile = arr[index1];
		}
		
		return percentile;
	}

	public static double getPercentileFromSorted(float [] arr, double fraction, int indexStart, int length) {

		if(arr.length==1){
			return arr[0];
		}

		double percentile=0;

		if(((int)(length*fraction))==(length*fraction)) {

			int index1 = (int)(length*fraction)-1 + indexStart;

			int index2 = index1+1 + indexStart;

			if(index1<0){

				throw new RuntimeException("Fraction to small.");
			}

			percentile = (arr[index1]+arr[index2])/2.0;

		} else {
			int index1 = (int)(length*fraction) + indexStart;

			percentile = arr[index1];
		}

		return percentile;
	}

	public static double getPercentileFromSortedInt(List<Integer> liScore, double fraction) {
		
		if(liScore.size()==1){
			return liScore.get(0);
		}
		
		double percentile=0;
		
		int len = liScore.size();
		
		if(((int)(len*fraction))==(len*fraction)) {
			int index1 = (int)((len*fraction)+0.5)-1;
			int index2 = index1+1;
			
			if(index1<0){
				index1=0;
			}
			
			percentile = (liScore.get(index1)+liScore.get(index2))/2.0;
			
		} else {
			int index1 = (int)(len*fraction);
			
			percentile = liScore.get(index1);
		}
		
		return percentile;
	}

	public static double getPercentileFromSortedLong(List<Long> liScore, double fraction) {

		if(liScore.size()==1){
			return liScore.get(0);
		}

		long percentile=0;

		int len = liScore.size();

		if(((int)(len*fraction))==(len*fraction)) {
			int index1 = (int)((len*fraction)+0.5)-1;
			int index2 = index1+1;

			if(index1<0){
				index1=0;
			}

			percentile = (long)((liScore.get(index1)+liScore.get(index2))/2.0);

		} else {
			int index1 = (int)(len*fraction);

			percentile = liScore.get(index1);
		}

		return percentile;
	}

	/**
	 * 
	 * @param liScore the list is sorted in the method.
	 * @return
	 */
	public static ModelMedianInteger getMedianForInteger(List<Integer> liScore) {
		
		Collections.sort(liScore);
		
		ModelMedianInteger modelMedian = new ModelMedianInteger();
		
		modelMedian.lowerQuartile = (int)(MedianStatisticFunctions.getPercentileFromSortedInt(liScore, 0.25) + 0.5);
		
		modelMedian.median = (int)(MedianStatisticFunctions.getPercentileFromSortedInt(liScore, 0.5) + 0.5);
		
		modelMedian.upperQuartile = (int)(MedianStatisticFunctions.getPercentileFromSortedInt(liScore, 0.75) + 0.5);
		
		modelMedian.size = liScore.size();
		
		return modelMedian;

	}
	
	public static ModelMedianLong getMedianForLong(List<Long> liScore) {

		Collections.sort(liScore);

		ModelMedianLong modelMedian = new ModelMedianLong();

		modelMedian.lowerQuartile = (long)(MedianStatisticFunctions.getPercentileFromSortedLong(liScore, 0.25) + 0.5);

		modelMedian.median = (long)(MedianStatisticFunctions.getPercentileFromSortedLong(liScore, 0.5) + 0.5);

		modelMedian.upperQuartile = (long)(MedianStatisticFunctions.getPercentileFromSortedLong(liScore, 0.75) + 0.5);

		modelMedian.size = liScore.size();

		return modelMedian;

	}

	/**
	 * 
	 * @param liScore the list is sorted in the method.
	 * @return
	 */
	public static ModelMedianDouble getMedianForDouble(List<Double> liScore) {
		
		Collections.sort(liScore);
		
		ModelMedianDouble modelMedian = new ModelMedianDouble();
		
		modelMedian.lowerQuartile = MedianStatisticFunctions.getPercentileFromSorted(liScore, 0.25);
		
		modelMedian.median = MedianStatisticFunctions.getPercentileFromSorted(liScore, 0.5);
		
		modelMedian.upperQuartile = MedianStatisticFunctions.getPercentileFromSorted(liScore, 0.75);
		
		modelMedian.size = liScore.size();
		
		return modelMedian;

	}

	

}
