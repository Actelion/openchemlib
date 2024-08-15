/*
 * Copyright (c) 1997 - 2018
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * @author Modest v. Korff
 */

package com.actelion.research.calc.statistics.median;

import com.actelion.research.util.ArrayUtils;

import java.util.ArrayList;
import java.util.Collection;
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
	public static ModelMedianInteger getMedianForInteger(int [] a) {
		return getMedianForInteger(ArrayUtils.toList(a));
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
	 * @param
	 * @return
	 */
	public static ModelMedianDouble getMedianForDouble(Collection<Double> liScoreRaw) {

		if(liScoreRaw.size()<1)
			return null;

		List<Double> liScore = new ArrayList<>(liScoreRaw);

		Collections.sort(liScore);
		
		ModelMedianDouble modelMedian = new ModelMedianDouble();
		
		modelMedian.lowerQuartile = MedianStatisticFunctions.getPercentileFromSorted(liScore, 0.25);
		
		modelMedian.median = MedianStatisticFunctions.getPercentileFromSorted(liScore, 0.5);
		
		modelMedian.upperQuartile = MedianStatisticFunctions.getPercentileFromSorted(liScore, 0.75);
		
		modelMedian.size = liScore.size();
		
		return modelMedian;

	}

	

}
