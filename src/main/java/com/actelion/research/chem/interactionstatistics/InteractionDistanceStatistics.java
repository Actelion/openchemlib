package com.actelion.research.chem.interactionstatistics;

import com.actelion.research.util.FastSpline;
import com.actelion.research.util.SmoothingSplineInterpolator;

import java.io.*;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class InteractionDistanceStatistics {
	
	private static volatile InteractionDistanceStatistics instance = new InteractionDistanceStatistics(); //eager initialization
	private static final String BASE_PATH = "/resources/interactionstatistics/";
	private volatile Map<Long,SplineFunction> pairPotentials;
	public static final double OCCURENCE_CUTOFF = 500;
	public static final double CUTOFF_RADIUS = 6.0;
	public static final double BIN_SIZE = 0.2;
	public static final double CHI = 1e-04;
	public static final double MAX = 10;
	private final Map<Long,int[]> interactionStatistics;

	
	public static InteractionDistanceStatistics getInstance() {

		if (instance==null) { 
			synchronized(InteractionDistanceStatistics.class) {
				if (instance==null) {
					instance = new InteractionDistanceStatistics();
				}
			}
		}
		return instance;
	}
	
	public int getInteractionClasses() {
		return pairPotentials.keySet().size();
	}

	public Map<Long, SplineFunction> getPairPotentials() {
		return pairPotentials;
	}

	public List<Integer> getAtomTypes(){
		HashSet<Integer> hs = new HashSet<>();

		for (long pair : pairPotentials.keySet()) {
			int [] a = splitLongToInts(pair);
			hs.add(a[0]);
			hs.add(a[1]);
		}
		return new ArrayList<>(hs);
	}

	public Set<Integer> getAtomKeySet() {
		Set<Integer> atomKeySet = new HashSet<Integer>();
		for(long l : pairPotentials.keySet()) {
			int[] pair = splitLongToInts(l);
			int a = getKey(pair[0]);
			int b = getKey(pair[1]);
			atomKeySet.add(a);
			atomKeySet.add(b);
		}
		return atomKeySet;
	}


	
	private synchronized void calculatePotentials() {
		splineCalculation();
	}
	
	private InteractionDistanceStatistics() {
		interactionStatistics = new ConcurrentHashMap<Long,int[]>();
		initialize();

	}
	
	private void initialize() {
	try {
		readFromFile();
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
	calculatePotentials();
	}
	
	
	
	
	public long combineIntsToLong(int a, int b) {
		int a1;
		int a2;
		if(a<b) {
			a1 = a;
			a2 = b;
		}
		
		else {
			a1 = b;
			a2 = a;
		}
		long l = (((long) a1<<32) | a2 );
		
		return l;
	}
	
	private int[] splitLongToInts(long l) {
		int a = (int)(l >> 32);
		int b = (int)l;
		
		return new int[] {a,b};
	}
	
	private boolean isGenericAtomPair(int a, int b) {
		return (isGenericAtomType(a) && isGenericAtomType(b));
	}
	
	private boolean isGenericAtomType(int a) {
		return Integer.toBinaryString(a).length()<=InteractionAtomTypeCalculator.AtomFlagCount.BASIC_ATOM_FLAG_COUNT.getCount();
	}
	
	//private boolean isExtendedAtomType(int a) {
	//	return Integer.toBinaryString(a).length()<=InteractionAtomTypeCalculator.AtomFlagCount.EXTENDED_ATOM_FLAG_COUNT.getCount();
	//}
	
	private boolean isSpecificAtomType(int a) {
		return Integer.toBinaryString(a).length()<=InteractionAtomTypeCalculator.AtomFlagCount.FUNC_GROUP_FLAG_COUNT.getCount();
	}
	
	public int getKey(int atomType)  {
		if (isSpecificAtomType(atomType)) {
			return (atomType & InteractionAtomTypeCalculator.AtomPropertyMask.SPECIFIC.getMask());
		}
		else return atomType;
	}
	
	private boolean isGenericAtomPair(long l) {
		int[] pair = splitLongToInts(l);
		return isGenericAtomPair(pair[0], pair[1]);
	}

 	
	public void addInteraction(int atom1, int atom2, double dist) {
		if (dist>=CUTOFF_RADIUS) return;
		//InteractionAtomPair ap =new InteractionAtomPair(atom1,atom2);
		long l = combineIntsToLong(atom1,atom2);
		interactionStatistics.putIfAbsent(l, new int[(int)(CUTOFF_RADIUS/BIN_SIZE)]);
		int index = (int) (0.5+dist/BIN_SIZE);
		int[] occurences = interactionStatistics.get(l);
		if(index<occurences.length)
		occurences[index]++;
	}
	

	

	
	public void write(String file) throws IOException {
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), StandardCharsets.UTF_8));
		//Write the distance pair interactions
		for(long l: interactionStatistics.keySet()) {
			int[] occurences = interactionStatistics.get(l);
			writer.write(Long.toString(l));
			for(int index=0; index<occurences.length; index++) {
				writer.write(" " + occurences[index]);
			}
			writer.write(System.getProperty("line.separator"));
		}

		writer.write(System.getProperty("line.separator"));
		writer.close();

	}
	
	public void readFromFile() throws IOException {
		try {
			//String file = "/home/joel/PL_stat_no_neighbours.txt";
			String file = BASE_PATH + "InteractionStatistics.txt";
			URL url =  InteractionDistanceStatistics.class.getResource(file);
			if(url==null) {
				throw new RuntimeException("Could not find the interactions parameter file in the classpath: " + file);
			}
			InputStream is = url.openStream();
			//InputStream is = new FileInputStream(file);
			BufferedReader reader = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8));
			String line;
			while((line = reader.readLine())!=null && !line.isEmpty()) {
				String s[] = line.split(" ");

				long l = Long.parseLong(s[0]);
				
				int[] occurences = new int[(int)(CUTOFF_RADIUS/BIN_SIZE)];
				for(int i=1;i<s.length;i++) {
					occurences[i-1] = Integer.parseInt(s[i]); 
				}

				interactionStatistics.putIfAbsent(l, occurences);

			
			}
			reader.close();
		}
			
	 	catch(Exception e) {
	 		throw new RuntimeException(e);
	}
	}

		

		
		private void splineCalculation() {
			pairPotentials = new HashMap<Long,SplineFunction>();

			double[] referenceSum = new double[(int)(CUTOFF_RADIUS/BIN_SIZE)];


			
			
			interactionStatistics.entrySet().stream().forEach(e -> {
				SplineFunction potential = new SplineFunction();
				potential.setOccurencesArray(e.getValue());
				pairPotentials.putIfAbsent(e.getKey(), potential);
			});

			
			Map<Long, double[]> discreteFunctions = interactionStatistics.entrySet().stream()
				    .collect(Collectors.toMap(e -> e.getKey(), e -> distanceNormalization(e.getValue()))); //check this line
			
			AtomicInteger runCount = new AtomicInteger(0);
			discreteFunctions.entrySet().stream().filter(e -> isGenericAtomPair(e.getKey())).
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
			for(long l : discreteFunctions.keySet()) {
				
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
						Y[i] = 0;
					}
				}
				

				/*
				 * add repulsive part for short distances
				 */
				
				int globalMinIndex = -1;
				int firstMaxIndex = -1;
				for(int i=1;i<Y.length-1;i++) {
					if(Y[i-1]<Y[i] && Y[i+1]<Y[i]) {//maximum
						if(firstMaxIndex==-1) {
							firstMaxIndex = i;
						}
					}
					else if(Y[i-1]>Y[i] && Y[i+1]>Y[i]) { //minimum
						if(globalMinIndex==-1) globalMinIndex = i;
						else if(Y[i]<Y[globalMinIndex]) globalMinIndex = i;
				}
				}

				if(firstMaxIndex>=0 && firstMaxIndex<globalMinIndex) {
					double derivative = MAX/(firstMaxIndex*BIN_SIZE);
					//double derivative = (Y[firstMaxIndex]-Y[firstMaxIndex+1])/BIN_SIZE;
					int i = firstMaxIndex-1;
					while(i>=0) {
						Y[i] = Y[i+1]+derivative*BIN_SIZE;
						i--;
					}
					ff = interpolator.interpolate(X, Y);
				}
				
				SplineFunction potential = pairPotentials.get(l);
				potential.setSplineFunction(ff);
				potential.setDiscreteFunction(discreteFunctions.get(l));

			}
		}
					


		
		private double[] normalize(double[] arr, double[] reference) {
			IntStream.range(0,arr.length).forEach(i -> arr[i]=-Math.log((arr[i] + CHI) / (reference[i] + CHI)));
			return arr;
		}
		
		private double[] distanceNormalization(int[] Y) {
			double[] YNorm = new double[Y.length];
			double normalizedSum = 0;
			for(int index=0; index<Y.length; index++) {
				double v = Math.PI * 4.0/3 * (Math.pow(InteractionDistanceStatistics.BIN_SIZE*(index+0.5),3)-Math.pow(InteractionDistanceStatistics.BIN_SIZE*(index-0.5),3));
				normalizedSum += Y[index] / v;
			}				
			if(normalizedSum==0) return YNorm;
			for(int index=0; index<Y.length; index++) {
				double v = Math.PI * 4.0/3 * (Math.pow(InteractionDistanceStatistics.BIN_SIZE*(index+0.5),3)-Math.pow(InteractionDistanceStatistics.BIN_SIZE*(index-0.5),3));
				YNorm[index] = (Y[index]) / (v * normalizedSum);
			}
			return YNorm;
		}
		
		public SplineFunction getFunction(long l) {
			return pairPotentials.get(l);
		}
		
		public SplineFunction getFunction(int a, int b) {

			
			SplineFunction spline;

			int a1 = a & InteractionAtomTypeCalculator.AtomPropertyMask.SPECIFIC.getMask();
			int b1 = b & InteractionAtomTypeCalculator.AtomPropertyMask.SPECIFIC.getMask();
			long l = combineIntsToLong(a1,b1);

			spline = pairPotentials.get(l);
			if(spline!=null && spline.getTotalOccurences()>OCCURENCE_CUTOFF) return spline;
			int a2 =  a & InteractionAtomTypeCalculator.AtomPropertyMask.EXTENDED.getMask();
			int b2 =  b & InteractionAtomTypeCalculator.AtomPropertyMask.EXTENDED.getMask();
			l = combineIntsToLong(a1,b2);

			spline = pairPotentials.get(l);
			if(spline!=null && spline.getTotalOccurences()>OCCURENCE_CUTOFF ) return spline;
			l = combineIntsToLong(a2,b1);

			spline = pairPotentials.get(l);
			if(spline!=null && spline.getTotalOccurences()>OCCURENCE_CUTOFF) return spline;
			l = combineIntsToLong(a2,b2);
			spline = pairPotentials.get(l);

			if(spline!=null && spline.getTotalOccurences()>OCCURENCE_CUTOFF) return spline;
			
			int a3 =  a & InteractionAtomTypeCalculator.AtomPropertyMask.BASIC.getMask();
			int b3 =  b & InteractionAtomTypeCalculator.AtomPropertyMask.BASIC.getMask();
			
			l = combineIntsToLong(a3,b1);
			spline = pairPotentials.get(l);

			if(spline!=null && spline.getTotalOccurences()>OCCURENCE_CUTOFF) return spline;
			l = combineIntsToLong(a1,b3);
			spline = pairPotentials.get(l);

			if(spline!=null && spline.getTotalOccurences()>OCCURENCE_CUTOFF) return spline;
			l = combineIntsToLong(a3,b2);
			spline = pairPotentials.get(l);

			if(spline!=null && spline.getTotalOccurences()>OCCURENCE_CUTOFF) return spline;
			l = combineIntsToLong(a2,b3);
			spline = pairPotentials.get(l);

			if(spline!=null && spline.getTotalOccurences()>OCCURENCE_CUTOFF) return spline;
			
			l = combineIntsToLong(a3,b3);

			spline = pairPotentials.get(l);
			return spline;
			

		}
	
			
			
	
	
	
	
	

	
	
	

}
