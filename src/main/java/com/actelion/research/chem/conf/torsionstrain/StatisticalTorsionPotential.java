package com.actelion.research.chem.conf.torsionstrain;

import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.interactionstatistics.SplineFunction;
import com.actelion.research.util.FastSpline;
import com.actelion.research.util.SmoothingSplineInterpolator;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class StatisticalTorsionPotential {
	
	private volatile Map<String,SplineFunction> torsionPotentials;
	private final Map<String,int[]> torsionStatistics;
	
	private static volatile StatisticalTorsionPotential instance;
	private static final String DATABASE_COD = "cod/";
	private static final String DATABASE_CSD = "csd/";
	private static final String TORSION_BINS_FILE = "torsionBins.txt";
	private static final String TORSION_IDS_FILE = "torsionID.txt";
	private static final String BASE_PATH = "/resources/";
	private static String database;
	
	private List<String> torsionIDs;
	
	public static final double OCCURENCE_CUTOFF = 500;
	public static final double BIN_SIZE = 5.0;
	public static final double CHI = 1e-04;
	public static final double MAX = 10;

	
	public static StatisticalTorsionPotential getInstance() {

		if (instance==null) { 
			synchronized(StatisticalTorsionPotential.class) {
				if (instance==null) {
					instance = new StatisticalTorsionPotential();
				}
			}
		}
		return instance;
	}
	
	private StatisticalTorsionPotential() {
		torsionIDs = new ArrayList<String>();
		torsionStatistics = new ConcurrentHashMap<String,int[]>();
		TorsionDB.initialize(TorsionDB.MODE_BINS);
		if (database == null) {
			InputStream is = TorsionDB.class.getResourceAsStream(BASE_PATH+DATABASE_CSD+TORSION_BINS_FILE);
			if (is != null) {
				database = DATABASE_CSD;
				}
			else {
				database = DATABASE_COD;
			}
		}
		
		
	initialize();
		


	}
	
	private void initialize() {
		BufferedReader torsionIDReader = new BufferedReader(new InputStreamReader(TorsionDB.class.getResourceAsStream(
				BASE_PATH+database+TORSION_IDS_FILE), StandardCharsets.UTF_8));
		
		try {
			readTorsionIDs(torsionIDReader);
			initBins();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException();
		}
	calculatePotentials();
	}
	
	private synchronized void calculatePotentials() {
		splineCalculation();
	}
	
	private void readTorsionIDs(BufferedReader reader) throws IOException {
		try {

			String line;
			while((line = reader.readLine())!=null && line.length()!=0) {
				String s = line.trim();
				torsionIDs.add(s);
			}
			reader.close();
		}
			
	 	catch(Exception e) {
	 		throw new RuntimeException(e);
	}
		
	}
	
	private void initBins()  {
		for(String torsionID : torsionIDs) {
			byte [] counts = TorsionDB.getTorsionBinCounts(torsionID);
			int[] occurences = new int[counts.length];
			IntStream.range(0, occurences.length).forEach(i -> occurences[i] = (int) counts[i]);
			torsionStatistics.putIfAbsent(torsionID, occurences);
		}
	}
	
	public  SplineFunction getFunction(String torsionID) {
		return torsionPotentials.get(torsionID);
	}
	

	
	private void splineCalculation() {
		torsionPotentials = new HashMap<String,SplineFunction>();

		double[] referenceSum = new double[(int)(360.0/BIN_SIZE)];


		
		
		torsionStatistics.entrySet().stream().forEach(e -> {
			SplineFunction potential = new SplineFunction();
			potential.setOccurencesArray(e.getValue());
			torsionPotentials.putIfAbsent(e.getKey(), potential);
		});

		
		Map<String, double[]> discreteFunctions = torsionStatistics.entrySet().stream()
			    .collect(Collectors.toMap(e -> e.getKey(), e -> normalization(e.getValue()))); //check this line
		
		AtomicInteger runCount = new AtomicInteger(0);
		discreteFunctions.entrySet().stream().
		forEach(statistics -> { 
			runCount.getAndIncrement();
			IntStream.range(0,statistics.getValue().length).forEach(i -> {
				referenceSum[i]+=statistics.getValue()[i];
			}
			);
			
		});
		for(int i=0;i<referenceSum.length;i++) {
			referenceSum[i]/=runCount.get();
		}
		//double[] refSum = distanceNormalization(referenceSum);

		discreteFunctions.replaceAll((k,v)-> normalize(v,referenceSum));
		double[] X = new double[referenceSum.length];
		IntStream.range(0, X.length).forEach(i-> X[i]= (i+0.5)*BIN_SIZE);
		//System.out.println(Arrays.toString(X));
		for(String l : discreteFunctions.keySet()) {
			
			double[] sigma = new double[X.length];
			Arrays.fill(sigma, 1);

			//
			//  Smoothing Spline
			//
			SmoothingSplineInterpolator interpolator = new SmoothingSplineInterpolator();
			interpolator.setLambda(0.005);
			interpolator.setSigma(sigma);

			FastSpline ff = interpolator.interpolate(X, discreteFunctions.get(l));
			
			double[] Y = new double[X.length];
			for(int i=0;i<X.length;i++) {
				try {
				Y[i] = ff.value(X[i]);
				}
				catch(Exception e) {
					e.printStackTrace();
					Y[i] = 0;
				}
			}
			

			SplineFunction potential = torsionPotentials.get(l);
			potential.setSplineFunction(ff);
			potential.setDiscreteFunction(discreteFunctions.get(l));

		}
	}
				


	
	private double[] normalize(double[] arr, double[] reference) {
		IntStream.range(0,arr.length).forEach(i -> arr[i]=-Math.log((arr[i] + CHI) / (reference[i] + CHI)));
		return arr;
	}
	
	private double[] normalization(int[] Y) {
		double[] YNorm = new double[Y.length];
		double normalizedSum = 0;
		for(int index=0; index<Y.length; index++) {
		
			normalizedSum += Y[index];
		}				
		if(normalizedSum==0) return YNorm;
		for(int index=0; index<Y.length; index++) {
			YNorm[index] = (Y[index]/normalizedSum);
		}
		return YNorm;
	}
	
	public SplineFunction getFunction(long l) {
		return torsionPotentials.get(l);
	}
	


}
