package com.actelion.research.calc;

import com.actelion.research.chem.descriptor.ISimilarityCalculator;
import com.actelion.research.util.Pipeline;
import com.actelion.research.util.datamodel.IIdentifiedObject;
import com.actelion.research.util.datamodel.IdentifiedObject;
import com.actelion.research.util.datamodel.ScorePoint;

import java.awt.*;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicLong;


/**
 * 
 * SimilarityMulticore
 * T is the descriptor object class
 * @author Modest von Korff
 * @version 1.0
 * 10 Dec 2010 MvK: Start implementation
 * Nov 2011 MvK: Generalization via interface definitions.
 * 24 Apr 2013 MvK: Some improvements on the thread handling.
 * 04 Dec 2014 MvK: Some improvements on the thread handling.
 * 09.05.2016 MvK: Calculates now the similarity matrix.
 * 26.04.2017 MvK: Using ExecutorServices now.
 * 26.11.2018 code changed. Similarity for identically labeled descriptors will now be calculated.
 */
public class SimilarityMulticore<T> {
	
	private static final int MAX_KERNELS = 80;
	
	// private static boolean VERBOSE = false;
	
	// private static final double DEFAULT_SIMILARITY = 1.0;
	
	private static final double DEFAULT_MINIMUM_SIMILARITY = 0.01;
	
	private static final long SLEEP_SHORT = 10;
	
	private static final long SLEEP_ULTRA_SHORT = 1;
	
	private ISimilarityCalculator<T> similarityCalculator;
	private List<IdentifiedObject<T>> liDescriptor1;
	private List<IdentifiedObject<T>> liDescriptor2;
	private int kernels;
	private AtomicLong sleep;
	private Pipeline<Point> queueIndices;
	private ConcurrentLinkedQueue<ScorePoint> queueScore;
	private List<RunSimilarityCalc> liRun;
	private int similarities2Calculate;
	private AtomicLong calculationsPerSecond;
	private Matrix maSimilarity;
	private boolean verbose;

	/**
	 *
	 * @param similarityCalculator
	 */
	public SimilarityMulticore(ISimilarityCalculator<T> similarityCalculator) {
		this(similarityCalculator, Math.min(Runtime.getRuntime().availableProcessors()-1, MAX_KERNELS));
	}
	
	public SimilarityMulticore(ISimilarityCalculator<T> similarityCalculator, int kernels) {
		this.similarityCalculator = similarityCalculator;
		this.kernels = kernels;
		queueIndices = new Pipeline<Point>();
		queueScore = new ConcurrentLinkedQueue<ScorePoint>();
		sleep = new AtomicLong();
		calculationsPerSecond = new AtomicLong();
		verbose = false;
	}

	public void setVerbose() {
		this.verbose = true;
	}

	public void run(IdentifiedObject<T> descriptor, List<IdentifiedObject<T>> liDescriptor2) {
		List<IdentifiedObject<T>> liOneSample = new ArrayList<IdentifiedObject<T>>();
		liOneSample.add(descriptor);
		run(liOneSample, liDescriptor2);
	}

	public void run(List<IdentifiedObject<T>> liDescriptor) {
		run(liDescriptor, liDescriptor, true);
	}

	public void run(List<IdentifiedObject<T>> liDescriptor1, List<IdentifiedObject<T>> liDescriptor2) {
		run(liDescriptor1, liDescriptor2, false);
	}

	/**
	 *
	 * @param liDescriptor1 list with descriptors,
	 * @param liDescriptor2 list with descriptors,
	 * liDescriptor1 will be compared with liDescriptor2 via ISimilarityCalculator given in constructor.
	 * RFesulting is a similarity matrix with rows = liDescriptor1.size() and cols = liDescriptor2.size()
	 */
	private void run(List<IdentifiedObject<T>> liDescriptor1, List<IdentifiedObject<T>> liDescriptor2, boolean singleList) {

		calculationsPerSecond.set(-1);

		long t1 = new Date().getTime();
		
		sleep.set(SLEEP_ULTRA_SHORT);
		
		if(verbose) {
			System.out.println("SimilarityMulticore start.");
			System.out.println("SimilarityMulticore kernels\t" + kernels);
		}

		this.liDescriptor1 = liDescriptor1;
		
		this.liDescriptor2 = liDescriptor2;

		if(verbose){
			System.out.println("liDescriptor1 " + liDescriptor1.size() + " liDescriptor2 " + liDescriptor2.size() + ".");
		}
		
		queueScore.clear();

		maSimilarity = new Matrix(liDescriptor1.size(), liDescriptor2.size());

		if(singleList){
			fillCalculationIndexQueueSingleList();
		} else {
			fillCalculationIndexQueueTwoLists();
		}

		liRun = new ArrayList();

		ExecutorService executorService = Executors.newFixedThreadPool(kernels);

		for (int i = 0; i < kernels; i++) {
			RunSimilarityCalc rsc = new RunSimilarityCalc(i, similarityCalculator, queueIndices, liDescriptor1, liDescriptor2, maSimilarity, singleList, queueScore);
			liRun.add(rsc);
			executorService.execute(rsc);
		}
		executorService.shutdown();
		while(!executorService.isTerminated()){
			try {Thread.sleep(1);} catch (InterruptedException e) {}
		}
		
		long t2 = new Date().getTime();
		long sec = (t2-t1) / 1000;

		if(sec!=0){
			calculationsPerSecond.set(getCalculatedSimilarityValues() / sec);
		}
		
		if(verbose){
			System.out.println("Similarity calculations " + getCalculatedSimilarityValues());
			System.out.println("Similarity calculations per second " + calculationsPerSecond.get());
									
			int sumCalc = 0;
			for (int i = 0; i < liRun.size(); i++) {
				RunSimilarityCalc rsc = liRun.get(i);
				sumCalc += rsc.getNSimilarityCalculations();
				System.out.println("Thread " + rsc.getIndexThread() + " calcs " + rsc.getNSimilarityCalculations());
			}
			
			System.out.println("Sum calcs " + sumCalc + ".");
		}
		
		sleep.set(SLEEP_SHORT);

	}
	
	public long getCalculationsPerSecond(){
		return calculationsPerSecond.get();
	}
	
	public int getSimilarities2Calculate(){
		return similarities2Calculate;
	}
	

	public long getCalculatedSimilarityValues(){
		
		long ccCalc = 0;
		
		for (RunSimilarityCalc rsc : liRun) {
			ccCalc += rsc.getNSimilarityCalculations();
		}
		
		return ccCalc;
	}
	
	private boolean isFinished() {
		
		if(!queueIndices.isAllDataIn()){
			return false;
		}
		
		if(!queueIndices.isEmpty()){
			return false;
		}

		boolean finished = true;
		
		if(queueScore.size() != similarities2Calculate){
			finished=false;
		}
				
		return finished;
	}

	public boolean hasMoreResults() {
		return !queueScore.isEmpty();
	}
	
	/**
	 * 
	 * @return similarity score with the id numbers of the compared input objects.
	 * The x value is the identifier from the object from liDescriptor1 and the y value from liDescriptor2.
	 */
	public ScorePoint getNextResult() {
		return queueScore.poll();
	}

	private void fillCalculationIndexQueueTwoLists(){
		
		queueIndices.setAllDataIn(false);
		
		similarities2Calculate = liDescriptor1.size() * liDescriptor2.size();
					
		for (int i = 0; i < liDescriptor1.size(); i++) {
			for (int j = 0; j < liDescriptor2.size(); j++) {
				Point p = new Point(i,j);
				queueIndices.addData(p);
			}
		}

		queueIndices.setAllDataIn(true);
				
		if(verbose){
			System.out.println("SimilarityMulticore sim to calc " + similarities2Calculate + ".");
		}

	}

	private void fillCalculationIndexQueueSingleList(){

		queueIndices.setAllDataIn(false);

		similarities2Calculate = ((liDescriptor1.size() * liDescriptor1.size()) - liDescriptor1.size()) / 2;

		for (int i = 0; i < liDescriptor1.size(); i++) {
			for (int j = i; j < liDescriptor1.size(); j++) {
				Point p = new Point(i,j);
				queueIndices.addData(p);
			}
		}

		queueIndices.setAllDataIn(true);

		if(verbose){
			System.out.println("SimilarityMulticore sim to calc " + similarities2Calculate + ".");
		}
	}


	public Matrix getSimilarityMatrix() {
		return maSimilarity;
	}

	private static class RunSimilarityCalc<T> implements Runnable {
		
		private ISimilarityCalculator iSimilarityCalculator;
		private Pipeline<Point> queueIndices;
		private List<IIdentifiedObject<T>> liDescriptor1;
		private List<IIdentifiedObject<T>> liDescriptor2;
		private Matrix maSimilarity;
		private boolean singleList;
		private ConcurrentLinkedQueue<ScorePoint> queueScore;
		private AtomicLong calculatedSimilarities;
		private int indexThread;
		
		public RunSimilarityCalc(int indexThread,
								 ISimilarityCalculator similarityCalculator,
								 Pipeline<Point> queueIndices,
								 List<IIdentifiedObject<T>> liDescriptor1,
								 List<IIdentifiedObject<T>> liDescriptor2,
								 Matrix maSimilarity,
								 boolean singleList,
								 ConcurrentLinkedQueue<ScorePoint> queueScore) {
			
			this.indexThread = indexThread;
			this.iSimilarityCalculator = similarityCalculator.getThreadSafeCopy();
			this.queueIndices = queueIndices;
			this.liDescriptor1 = liDescriptor1;
			this.liDescriptor2 = liDescriptor2;
			this.maSimilarity = maSimilarity;
			this.singleList = singleList;
			this.queueScore = queueScore;
			calculatedSimilarities = new AtomicLong();
		}
		
		public void run() {
			
			while(!queueIndices.wereAllDataFetched()) {

				Point p = queueIndices.pollData();

				if(p == null) {
					try {Thread.sleep(SLEEP_SHORT);} catch (InterruptedException e) {}
					continue;
				}

				int indexX = p.x;
				int indexY = p.y;

				IIdentifiedObject<T> idObj1 = null;
				idObj1 = liDescriptor1.get(indexX);
				IIdentifiedObject<T> idObj2 = null;
				idObj2 = liDescriptor2.get(indexY);
				
				ScorePoint sp = new ScorePoint((int)idObj1.getId(), (int)idObj2.getId());
					
				try {
					double sc = iSimilarityCalculator.getSimilarity(idObj1.getData(), idObj2.getData());
					if(sc<DEFAULT_MINIMUM_SIMILARITY){
						sc=DEFAULT_MINIMUM_SIMILARITY;
					}
					
					calculatedSimilarities.incrementAndGet();
					
					sp.setScore(sc);
					maSimilarity.set(indexX, indexY, sc);
					if(singleList) {
						maSimilarity.set(indexY, indexX, sc);
					}

					queueScore.add(sp);

				} catch (Exception e) {
					sp.setScore(Double.NaN);
					queueScore.add(sp);
					e.printStackTrace();
				}
			}
		}

		/**
		 * @return the indexThread
		 */
		protected int getIndexThread() {
			return indexThread;
		}

		public long getNSimilarityCalculations() {
			return calculatedSimilarities.get();
		}
		
	}
}

