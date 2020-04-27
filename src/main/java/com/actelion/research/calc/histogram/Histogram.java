package com.actelion.research.calc.histogram;

import com.actelion.research.calc.Matrix;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Random;


/**
 * 
 * Histogram
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
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

	public Histogram(double [] arrRaw, double min, double max, int bins) {
		initialize(arrRaw,min,max,bins);
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
    	
    	int digits = logBins-logRange+1;
    	
    	String sFormatDigits="";
    	for (int i = 0; i < digits; i++) {
    		sFormatDigits += "0";
		}
    	
    	String sFormat="0";
    	
    	if(sFormatDigits.length()>0)
    		sFormat += "." + sFormatDigits;
    	
    	
    	NumberFormat nfX = new DecimalFormat(sFormat);
    	NumberFormat nfY = new DecimalFormat("0");
    	NumberFormat nfFractionCumulative = new DecimalFormat("0.00");

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

	public double getBinWidth() {
		return binwidth;
	}

	public static void main(String[] args) {
		int n=1000;

		double min=0;
		double max=0.1;
		double minBin=0;
		double maxBin=0.1;
		int bins=20;

		double [] arr = new double [n];
		Random rnd = new Random();
		for (int i = 0; i < n; i++) {
			// arr[i]= min + rnd.nextInt((int)(max-min));
			arr[i]= min + rnd.nextDouble()*0.1;
		}



		Histogram h = new Histogram(arr,minBin,maxBin,bins);
		System.out.println(h);
		System.out.println("sum " + h.getSumY());

	}

}
