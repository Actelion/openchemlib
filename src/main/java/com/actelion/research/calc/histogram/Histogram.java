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

import com.actelion.research.calc.Matrix;
import com.actelion.research.util.Formatter;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.Locale;
import java.util.Random;


/**
 * 
 * Histogram
 * @author Modest von Korff
 * 29 Feb 2008 MvK: Start implementation
 * 06 Jun 2016 MvK: added cumulative histogram to output.
 * 03.12.2018 improvements toString()
 */
public class Histogram {



	private double minimumX;

	private double maximumX;

	private float [] arrBins;

	private int [] arrFrequencies;

	private double maximumY;

	private double sumY;

	private double binwidth;

	private double [] arrRaw;

	private int digits;

	/**
	 *
	 * @param arrRaw data
	 * @param min minimum value considered in the histogram
	 * @param max maximum value considered in the histogram
	 * @param bins
	 */
	public Histogram(double [] arrRaw, double min, double max, int bins) {
		initialize(arrRaw,min,max,bins);
	}

	public Histogram(Collection<Double> values, double min, double max, int bins) {
		double [] v = new double[values.size()];
		int c=0;
		for (Double value : values) {
			v[c++]=value;
		}
		initialize(v,min,max,bins);
	}

	public Histogram(double [] arrRaw, int bins) {

		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;

		for (double v : arrRaw) {
			if(v>max){
				max=v;
			}
			if(v<min){
				min=v;
			}
		}
		initialize(arrRaw,min,max,bins);
	}



	public void initialize(double [] arrRaw, double min, double max, int bins) {

		if(bins==0)
			throw new RuntimeException("Number of bins is 0.");

		this.arrRaw = arrRaw;

		minimumX = min - min * ConstantsHistogram.TINY_FACTOR;

		maximumX = max + max * ConstantsHistogram.TINY_FACTOR;

		arrFrequencies = new int [bins];

		arrBins = new float [bins];

		digits = -1;

		calcBins();

		calcHistogram();
	}

    private void calcHistogram(){
    	int bins = arrFrequencies.length;
    	
    	binwidth =  bins / (maximumX-minimumX);
    	
    	maximumY = 0;
    	for (int i = 0; i < arrRaw.length; i++) {
			if(arrRaw[i]>= minimumX && arrRaw[i]< maximumX){
				int pos = (int)((arrRaw[i]-minimumX)*binwidth);
				arrFrequencies[pos]++;
				sumY++;
				if(arrFrequencies[pos]>maximumY)
					maximumY=arrFrequencies[pos];
			}
		}
    }

    public int getBinIndex(double v) {
		int pos = (int)((v-minimumX)*binwidth);
		return pos;
	}

    private void calcBins(){
		double incr = (maximumX-minimumX)/arrBins.length;
		arrBins[0]=(float)minimumX;
		for (int i = 1; i < arrBins.length; i++) {
			arrBins[i]=(float)(arrBins[i-1]+incr);
		}
    }

    public int [] getFrequencies(){
    	return arrFrequencies;
    }

	public int getFrequency(int index){
		return arrFrequencies[index];
	}

    public float [] getBins(){
    	return arrBins;
    }

    public float getBin(int index){
    	return arrBins[index];
    }

    public int getNumBins(){
    	return arrBins.length;
    }
    public double getMaximumY(){
    	return maximumY;
    }
    public double getSumY(){
    	return sumY;
    }
    
    /**
     * Parts of the range outside the bins are not considered. 
     * @param low
     * @param up
     * @return
     */
    public double getSumFromRange(double low, double up){
    	int bins = arrFrequencies.length;
    	
    	double devider =  bins / (maximumX-minimumX);
    	
    	int posStart = (int)Math.max(0, ((low-minimumX)*devider));
    	int posEnd = (int)Math.min(bins-1, ((up-minimumX)*devider));
    	double sum=0;
    	
    	for (int i = posStart; i < posEnd+1; i++) {
			sum += arrFrequencies[i];
		}
    	return sum;
    	
    }
    
	public Matrix getHistogram(){

		Matrix maBins = MatrixBasedHistogram.getHistogramBins(minimumX, maximumX, arrFrequencies.length);

		Matrix maHist = MatrixBasedHistogram.getHistogram(arrRaw, maBins);

		return maHist;

	}

	public double getMinimumX() {
		return minimumX;
	}

	public double getMaximumX() {
		return maximumX;
	}

	public void setDigits(int digits) {
		this.digits = digits;
	}

	/**
     * First row: bins, upper border.
     * Sec row: frequencies
	 * Third: row cumulative fraction
     */
    public String toString() {
    	StringBuilder sb = new StringBuilder();
    	double range = maximumX-minimumX;
    	int logRange = (int)Math.log10(range);
    	int bins = arrFrequencies.length;
    	int logBins = (int)Math.log10(bins);

		if(digits<0)
    		digits = logBins-logRange+1;
    	
    	String sFormatDigits="";
    	for (int i = 0; i < digits; i++) {
    		sFormatDigits += "0";
		}
    	
    	String sFormat="0";
    	
    	if(sFormatDigits.length()>0)
    		sFormat += "." + sFormatDigits;

    	NumberFormat nfX = new DecimalFormat(sFormat, new DecimalFormatSymbols(Locale.US));
    	NumberFormat nfY = new DecimalFormat("0", new DecimalFormatSymbols(Locale.US));
    	NumberFormat nfFractionCumulative = new DecimalFormat("0.00", new DecimalFormatSymbols(Locale.US));

    	String [] arrStrX = new String [bins];
    	String [] arrStrY = new String [bins];
    	String [] arrStrValCumulative = new String [bins];

		int sumFreq = 0;

    	for (int i = 0; i < arrStrY.length; i++) {
    		String sX = nfX.format(arrBins[i]);
    		String sY = nfY.format(arrFrequencies[i]);
			sumFreq += arrFrequencies[i];
			double fractionCumulative = (double)sumFreq / arrRaw.length;
    		String sFractionCumulative = nfFractionCumulative.format(fractionCumulative);
			int maxLen = sX.length();
			if(sY.length() > maxLen) {
				maxLen = sY.length();
			}

			if(sFractionCumulative.length() > maxLen) {
				maxLen = sFractionCumulative.length();
			}

    		while(sX.length() < maxLen){
    			sX = " " + sX ;
    		}

    		while(sY.length() < maxLen){
    			sY = " " + sY ;
    		}
    		
    		while(sFractionCumulative.length() < maxLen){
				sFractionCumulative = " " + sFractionCumulative ;
    		}

			arrStrX[i] = sX;
    		
    		arrStrY[i] = sY;

			arrStrValCumulative[i] = sFractionCumulative;
 		}
    	
    	for (int i = 0; i < arrStrX.length; i++) {
			sb.append(arrStrX[i]);
			if(i<arrStrX.length-1){
				sb.append(" ");
			}
		}
    	sb.append("\n");
    	for (int i = 0; i < arrStrY.length; i++) {
			sb.append(arrStrY[i]);
			if(i<arrStrY.length-1){
				sb.append(" ");
			}
		}

    	sb.append("\n");
    	for (int i = 0; i < arrStrValCumulative.length; i++) {
			sb.append(arrStrValCumulative[i]);
			if(i<arrStrValCumulative.length-1){
				sb.append(" ");
			}
		}

    	return sb.toString();
    }

	/**
	 * Replaces the spaces with tab. So, the output can be used in Excel by copy/paste.
	 * @return
	 */
	public String toStringWithTabs() {
		String s = toString();
		String sTab = s.replaceAll("[ ]{1,}", "\t");
		return sTab;
	}


		public double getBinWidth() {
		return binwidth;
	}

	public static void main(String[] args) {
		int n=1000;

		double min=0;
		double max=0.1;
		double minBin=0;
		double maxBin=1;
		int bins=20;

		double [] arr = new double [n];
		Random rnd = new Random();
		for (int i = 0; i < n; i++) {
			// arr[i]= min + rnd.nextInt((int)(max-min));
			arr[i]= min + rnd.nextDouble();
		}



		Histogram h = new Histogram(arr,minBin,maxBin,bins);
		System.out.println(h);
		System.out.println("sum " + h.getSumY());

	}

}
