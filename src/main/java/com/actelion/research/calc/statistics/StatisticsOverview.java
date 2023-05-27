package com.actelion.research.calc.statistics;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Locale;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.histogram.ConstantsHistogram;
import com.actelion.research.calc.histogram.IntegerHistogram;
import com.actelion.research.calc.histogram.MatrixBasedHistogram;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.datamodel.DoubleArray;
import com.actelion.research.util.datamodel.IntArray;

/**
 * 
 * StatisticsOverview
 * Some basic statistics about the given data
 * 28 Jun 2010 MvK: Start implementation
 */
public class StatisticsOverview {

	public static final String TAG_MEAN = "Avr";
	public static final String TAG_SDV = "SDV";
	public static final String TAG_MEDIAN = "Median";
	public static final String TAG_PERCENTILE05 = "Percentile.05";
	public static final String TAG_PERCENTILE95 = "Percentile.95";

	public static final NumberFormat DF1 = new DecimalFormat("0.0", new DecimalFormatSymbols(Locale.US));
	public static final NumberFormat DF3 = new DecimalFormat("0.000", new DecimalFormatSymbols(Locale.US));
	public static final NumberFormat DF4 = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.US));
	private static final NumberFormat DF3Plus = new DecimalFormat("0.000##", new DecimalFormatSymbols(Locale.US));

	private static final int BINS = 20;
	
	private static final int WIDTH = 8;
	
	private static final int DIGITS = 2;

	private String name;
	
	private double min;
	private double mean;
	private double max;

	private double sdv;
	
	private double percentile05;
	private double leftQuartile;

	private double median;

	private double rightQuartile;
	private double percentile95;

	
	private Matrix histogram;

	private int valsBelowHistMin;
	
	private int valsAboveHistMax;
	
	private int bins;
	
	private DoubleArray data;
	
	private boolean evaluated;
	
	public StatisticsOverview() {
		data = new DoubleArray();
		bins = BINS;
		evaluated = false;
	}
	
	public StatisticsOverview(DoubleArray da) {
		data = da;
		bins = BINS;
		evaluated = false;
	}
	
	public void add(double value){
		data.add(value);
		evaluated = false;
		
	}
	
	public void add(double [] arr){
		
		for (int i = 0; i < arr.length; i++) {
			data.add(arr[i]);	
		}
		
		evaluated = false;
		
	}
	
	public void add(int [] arr){
		
		for (int i = 0; i < arr.length; i++) {
			data.add(arr[i]);	
		}
		
		evaluated = false;
		
	}

	public DoubleArray getData() {
		return data;
	}

	public double getMean() {
		return mean;
	}

	public double getSdv() {
		return sdv;
	}

	public double getMedian() {
		return median;
	}

	public ModelStatisticsOverview evaluate(){

		double histMin = data.min() - data.min() * ConstantsHistogram.TINY_FACTOR;
		double histMax = data.max() + data.max() * ConstantsHistogram.TINY_FACTOR;

		ModelStatisticsOverview modelStatisticsOverview = new ModelStatisticsOverview();

		if(Math.abs(histMin-histMax) > Matrix.TINY08) {

			evaluate(histMin, histMax);
			modelStatisticsOverview.min = data.min();
			modelStatisticsOverview.avr = mean;
			modelStatisticsOverview.max = data.max();
			modelStatisticsOverview.sdv = sdv;

		} else {
			modelStatisticsOverview.min = 0;
			modelStatisticsOverview.avr = 0;
			modelStatisticsOverview.max = 0;
			modelStatisticsOverview.sdv = 0;
		}

		return modelStatisticsOverview;
	}
	
	public void evaluate(double histMin, double histMax){

		if(Math.abs(histMin-histMax) < Matrix.TINY08) {
			throw new RuntimeException("Equal histogram boundaries! histMin " + histMin + " histMax" + histMax + ".");
		}

		calculateMedianStatistics();

		double [] arr = data.get();

		min = ArrayUtilsCalc.min(arr);

		mean = ArrayUtilsCalc.getMean(arr);

		max = ArrayUtilsCalc.max(arr);

		sdv = ArrayUtilsCalc.getStandardDeviation(arr);
		
		Matrix maBins = MatrixBasedHistogram.getHistogramBins(histMin, histMax, bins);

		Matrix ma = new Matrix(true, arr);

		histogram = MatrixBasedHistogram.getHistogram(ma, maBins);

		for (int i = 0; i < ma.cols(); i++) {
			if(ma.get(0, i)<histMin)
				valsBelowHistMin++;
			else if(ma.get(0, i)>histMax)
				valsAboveHistMax++;
		}
		
		evaluated = true;
	}
	
	public void evaluateIntegerBins(double histMin, double histMax){
		
		calculateMedianStatistics();
		
		Matrix ma = new Matrix(true, data.get());
		
		mean = ma.getMean();
		
		sdv = ma.getStandardDeviation();
		
		int [][] arrBins = IntegerHistogram.getBinsEquallyDistributed(bins, (int)histMin, (int)histMax);
			
		Matrix maBins = new Matrix(arrBins);
		
		maBins = maBins.getTranspose();

		histogram = MatrixBasedHistogram.getHistogram(ma, maBins);
		
		for (int i = 0; i < ma.cols(); i++) {
			if(ma.get(0, i)<histMin)
				valsBelowHistMin++;
			else if(ma.get(0, i)>histMax)
				valsAboveHistMax++;
		}
		
		evaluated = true;
	}
	
    private double calculateMedianStatistics() {
    	
    	double [] arr = data.get();
    	
    	Arrays.sort(arr);
    	
		percentile05 = getQuartile(arr, 0.05);
		leftQuartile = getQuartile(arr, 0.25);

		median = getQuartile(arr, 0.5);
		
		rightQuartile = getQuartile(arr, 0.75);
		percentile95 = getQuartile(arr, 0.95);

		return median;
    }

	public double getPercentile05() {
		return percentile05;
	}

	public double getPercentile95() {
		return percentile95;
	}

	private static double getQuartile(double [] arr, double q) {
    	double v = 0;

		if(q<0){
			throw new RuntimeException("Negative values are not allowed!");
		}

		if(arr.length % 2==0) {
			if(((int)(q * arr.length))==0){
				v = arr[0];
			} else {
				int p1 = (int) ((q * arr.length) - 1);
				int p2 = (int) (q * arr.length);
				v = (arr[p1] + arr[p2]) / 2.0;
			}
		} else {
			int p = (int)(arr.length * q);
			v = arr[p];
		}
    	
    	return v;
    }

    public double getQuartile(double q) {
    	double [] arr = data.get();
    	Arrays.sort(arr);
    	return getQuartile(arr, q);
    }

    public int getNumData(){
    	return data.size();
    }
	
	@Override
	public String toString() {

		if(data.size()==0){
			return "";
		}

		if(!evaluated)
			evaluate();
		
		StringBuilder sb = new StringBuilder();
		
		sb.append("Name\t" + name);
		sb.append("\n");
		sb.append("Values\t" + data.size());
		sb.append("\n");
		sb.append("min\t" + min);
		sb.append("\n");
		sb.append("Mean\t" + DF3Plus.format(mean));
		sb.append("\n");
		sb.append("max\t" + max);
		sb.append("\n");
		sb.append("SDV\t" + DF3Plus.format(sdv));
		sb.append("\n");
		sb.append("Quartile 0.25\t" + DF3Plus.format(leftQuartile));
		sb.append("\n");
		sb.append("Median\t" + DF3Plus.format(median));
		sb.append("\n");
		sb.append("Quartile 0.75\t" + DF3Plus.format(rightQuartile));
		sb.append("\n");
		sb.append("Histogram values below hist min " + valsBelowHistMin + ", values above hist max " + valsAboveHistMax);
		sb.append("\n");
		sb.append(MatrixBasedHistogram.histogram2String(histogram, DIGITS, WIDTH));
		
		Matrix histTrans = histogram.getTranspose();
		
		DecimalFormat dfBins = Matrix.format(DIGITS);
		
		sb.append("\n");
		for (int i = 0; i < histTrans.rows(); i++) {
			sb.append(Matrix.format(histTrans.get(i,0), dfBins, WIDTH));
			sb.append("\t");
			sb.append(Matrix.format(histTrans.get(i,1), dfBins, WIDTH));
			sb.append("\t");
			sb.append(Matrix.format(histTrans.get(i,2), dfBins, WIDTH));
			sb.append("\n");
		}

		return sb.toString();
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public void setBins(int bins) {
		this.bins = bins;
	}

	public int getValsBelowHistMin() {
		return valsBelowHistMin;
	}

	public int getValsAboveHistMax() {
		return valsAboveHistMax;
	}

	public Matrix getHistogram() {
		return histogram;
	}


	public static ModelStatisticsOverview get(DoubleArray da){

		StatisticsOverview statisticsOverview = new StatisticsOverview(da);

		return statisticsOverview.evaluate();

	}

	public static ModelStatisticsOverviewMedian getMedianOverview(DoubleArray da){
		StatisticsOverview statisticsOverview = new StatisticsOverview(da);
		statisticsOverview.evaluate();
		ModelStatisticsOverviewMedian model =
				new ModelStatisticsOverviewMedian(
						statisticsOverview.percentile05,
						statisticsOverview.leftQuartile,
						statisticsOverview.median,
						statisticsOverview.rightQuartile,
						statisticsOverview.percentile95);

		return model;

	}

	public static String toString(DoubleArray da, String text1){

		StringBuilder sb = new StringBuilder();

		StatisticsOverview so = new StatisticsOverview(da);
		so.evaluate();
		sb.append(text1);
		sb.append("\t");
		sb.append(DF3.format(so.getMean()));
		sb.append("\t");
		sb.append(DF3.format(so.getSdv()));
		sb.append("\t");
		sb.append(DF3.format(so.getMedian()));
		sb.append("\t");
		sb.append(DF3.format(so.getPercentile05()));
		sb.append("\t");
		sb.append(DF3.format(so.getPercentile95()));

		return sb.toString();
	}

	public static String toString(IntArray ia, String text){

		StringBuilder sb = new StringBuilder();

		DoubleArray da = new DoubleArray(ia);

		StatisticsOverview so = new StatisticsOverview(da);
		so.evaluate();
		sb.append(text);
		sb.append("\t");
		sb.append(DF1.format(so.getMean()));
		sb.append("\t");
		sb.append(DF1.format(so.getSdv()));
		sb.append("\t");
		sb.append(DF1.format(so.getMedian()));
		sb.append("\t");
		sb.append(DF1.format(so.getPercentile05()));
		sb.append("\t");
		sb.append(DF1.format(so.getPercentile95()));

		return sb.toString();
	}

	public static String toStringHeader(){

		StringBuilder sb = new StringBuilder();
		sb.append("text");
		sb.append("\t");
		sb.append("mean");
		sb.append("\t");
		sb.append("sdv");
		sb.append("\t");
		sb.append("median");
		sb.append("\t");
		sb.append("perc05");
		sb.append("\t");
		sb.append("perc95");

		return sb.toString();
	}


	public static String toString(DoubleArray da, String text1, String text2){

		StringBuilder sb = new StringBuilder();

		StatisticsOverview so = new StatisticsOverview(da);
		so.evaluate();
		sb.append(text1);
		sb.append("\t");
		sb.append(text2);
		sb.append("\t");
		sb.append(DF3.format(so.getMean()));
		sb.append("\t");
		sb.append(DF3.format(so.getSdv()));
		sb.append("\t");
		sb.append(DF3.format(so.getMedian()));
		sb.append("\t");
		sb.append(DF3.format(so.getPercentile05()));
		sb.append("\t");
		sb.append(DF3.format(so.getPercentile95()));

		return sb.toString();
	}
}
