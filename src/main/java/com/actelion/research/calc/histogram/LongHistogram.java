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

package com.actelion.research.calc.histogram;

import java.util.List;
import java.util.Random;

/**
 * LongHistogram
 * @author Modest von Korff
  * Mar 08, 2023 MvK Start implementation
 */
public class LongHistogram {


	public static final long [][] ARR_BINS_EXAMPLE = {{1,2}, {2,3}, {3,4}, {4,5}, {5,6}, {6,8}, {8,10}, {10,20}, {20,50}, {50,100}, {100,1000}};


	private long [][] arrBins;

	private long [] arrCounts;

	public LongHistogram(long [][] arrBins) {
		this(arrBins, true);
	}

	public LongHistogram(long [][] arrBins, boolean consecutive) {
		
		this.arrBins = arrBins;
		
		arrCounts = new long [arrBins.length];
		
		if(consecutive){
			for (int i = 1; i < arrBins.length; i++) {
				if(arrBins[i-1][1] != arrBins[i][0]){
					throw new RuntimeException("Non consecutive bins!");
				}
			}
		}
	}
	
	/**
	 * Added to the bin where fullfilling the criteria v >= lower bound and v < higher bound. 
	 * @param v
	 */
	public void add(long v){
		for (int i = 0; i < arrBins.length; i++) {
			if((v >= arrBins[i][0]) && (v < arrBins[i][1])){
				arrCounts[i]++;
			}
		}
	}
	
	public void add(long [] a){
		for (int i = 0; i < a.length; i++) {
			add(a[i]);
		}
	}

	public void add(List<Long> li){
		for (int i = 0; i < li.size(); i++) {
			add(li.get(i));
		}
	}

	public int getTotalCounts(){
		int c = 0;
		for (int i = 0; i < arrBins.length; i++) {
			c += arrCounts[i];
		}
		return c;
	}

	public long [] getBinWithNPercentOfAllCounts(int percent){

		long nTotal = getTotalCounts();

		long nPercent = (long) (nTotal * (percent/100.0));

		long c = 0;
		
		int index = -1;
		for (int i = 0; i < arrBins.length; i++) {
			c += arrCounts[i];
			
			if(c >= nPercent){
				index = i;
				break;
			}
		}

		long [] bin = new long [2];
		
		bin[0] = arrBins[index][0]; 
		bin[1] = arrBins[index][1]; 
		
		return bin;
	}
	

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		
		String sep = "\t";
		int [] arrLenOut = new int [arrCounts.length];
		for (int i = 0; i < arrBins.length; i++) {
			long max = arrCounts[i];
			if(arrBins[i][0]>max){
				max = arrBins[i][0];
			}
			if(arrBins[i][1]>max){
				max = arrBins[i][1];
			}
			String sVal = Long.toString(max);
			arrLenOut[i]=sVal.length()+1;
		}
		
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < arrBins.length; i++) {
			String s = Long.toString(arrBins[i][0]);
			StringBuilder sb0 = new StringBuilder(s);
			while(sb0.length()<arrLenOut[i]){
				sb0.append(" ");
			}
			sb.append(sb0);
		}

		sb.append("\n");
		for (int i = 0; i < arrBins.length; i++) {
			String s = Long.toString(arrBins[i][1]);
			StringBuilder sb0 = new StringBuilder(s);
			while(sb0.length()<arrLenOut[i]){
				sb0.append(" ");
			}
			sb.append(sb0);
		}
		
		sb.append("\n");
		
		for (int i = 0; i < arrCounts.length; i++) {
			String s = Long.toString(arrCounts[i]);
			StringBuilder sb0 = new StringBuilder(s);
			while(sb0.length()<arrLenOut[i]){
				sb0.append(" ");
			}
			sb.append(sb0);
		}
		return sb.toString();
	}
	
	public static void main(String[] args) {
		
		int n = 10000;
		
		int max = 1000;
		
		// int [][] arrBins = ARR_BINS_EXAMPLE;
		
		int bins = 20;
		
		long [][] arrBins = getBinsEquallyDistributed(bins, max);
		LongHistogram ih = new LongHistogram(arrBins);

		Random rnd = new Random();
		
		for (int i = 0; i < n; i++) {
			ih.add(rnd.nextInt(max));
		}
		
		System.out.println(ih.toString());
		
		
		
	}
	
	public static long [][] getBinsEquallyDistributed(int bins, long maxValue){
		return getBinsEquallyDistributed(bins, 0, maxValue);
	}
	
	public static long [][] getBinsEquallyDistributed(int bins, long minValue, long maxValue){

		long rangeBin = (maxValue-minValue) / bins;
		
		if(bins > maxValue-minValue){
			rangeBin = 1;
			bins = (int)(maxValue-minValue + 1);
		}

		long [][] arr = new long [bins][];

		long binStart = minValue;
		
		for (int i = 0; i < bins; i++) {

			long binEnd = binStart + rangeBin;

			long [] bin = new long [2];
			
			bin[0]=binStart;
			bin[1]=binEnd;
			
			binStart = binEnd;

			arr[i]=bin;
		}
		
		if(arr[arr.length-1][1] <= maxValue){
			arr[arr.length-1][1] = maxValue+1;
		}
		
		return arr;
		
	}
	
	
	

}
