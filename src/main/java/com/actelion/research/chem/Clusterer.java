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
 *
 * @author Thomas Sander
 */

package com.actelion.research.chem;

import com.actelion.research.calc.DataProcessor;
import com.actelion.research.chem.descriptor.DescriptorHandler;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class Clusterer<T> extends DataProcessor {
	private volatile int[]			mClusterNo,mNoOfMembers;
	private volatile int			mNoOfCompounds;
	private volatile float[][]		mSimilarityMatrix;
	private volatile T[]			mDescriptor;
	private volatile DescriptorHandler<T,?> mDescriptorHandler;
	private volatile AtomicInteger	mSMPCompoundIndex;

	private boolean[]				mIsRepresentative;
	private int						mNoOfClusters,mThreadCount;
	private ExecutorService			mExecutor;
	private ClusterWorker<T>[]		mClusterWorker;

    @SuppressWarnings("unchecked")
    public Clusterer(DescriptorHandler<T,?> descriptorHandler, T[] descriptor) {
		mDescriptorHandler = descriptorHandler;
		mDescriptor = descriptor;
		mNoOfCompounds = mDescriptor.length;

		mSimilarityMatrix = new float[mNoOfCompounds][];
		for (int i=1; i<mNoOfCompounds; i++)
			mSimilarityMatrix[i] = new float[i];

		mThreadCount = Runtime.getRuntime().availableProcessors();
		if (mThreadCount != 1) {
			mExecutor = Executors.newFixedThreadPool(mThreadCount);
			mClusterWorker = new ClusterWorker[mThreadCount];
			for (int t=0; t<mThreadCount; t++)
				mClusterWorker[t] = new ClusterWorker<T>();
			}
		}

	/**
	 * Defines the criteria for stopping the clustering.
	 * At least one of the two limits must be in the applicable valid range.
	 * @param similarityLimit >0...<=1.0 or 0.0 if not applied
	 * @param clusterCountLimit >=2...objectCount or -1 if not applied
	 */
	public void cluster(double similarityLimit, int clusterCountLimit) {
		calculateSimilarityMatrix(false);
		if (threadMustDie()) {
		    stopProgress("clustering cancelled");
			return;
			}

		mNoOfMembers = new int[mNoOfCompounds];	// initialize no of cluster members
		mClusterNo = new int[mNoOfCompounds];	// initialize compound's cluster numbers
		for (int i=0; i<mNoOfCompounds; i++) {
			mNoOfMembers[i] = 1;
			mClusterNo[i] = i;
			}

		if (clusterCountLimit < 1)
			clusterCountLimit = 1;

		if (similarityLimit != 0.0)
			startProgress("Clustering Compounds...", 0, (int)(5000.0*(1.0-similarityLimit)));
		else
			startProgress("Clustering Compounds...", 0, mNoOfCompounds - clusterCountLimit);

		mNoOfClusters = mNoOfCompounds;
		while (mNoOfClusters > clusterCountLimit) {
			float maxSimValue = 0;		// find highest similarity level
			int maxCluster1 = -1;
			int maxCluster2 = -1;
			if (mThreadCount == 1) {
				for (int cluster2=1; cluster2<mNoOfCompounds; cluster2++) {
					if (mNoOfMembers[cluster2] > 0) {
						for (int cluster1=0; cluster1<cluster2; cluster1++) {
							if (mNoOfMembers[cluster1] != 0) {
								if (maxSimValue < mSimilarityMatrix[cluster2][cluster1]) {
									maxSimValue = mSimilarityMatrix[cluster2][cluster1];
									maxCluster1 = cluster1;
									maxCluster2 = cluster2;
									}
								}
							}
						}
					}
				}
			else {
				runInParallel(ClusterWorker.FIND_MAXIMUM_SIMILARITY);
				for (ClusterWorker<T> worker:mClusterWorker) {
					if (maxSimValue < worker.getMaxSimilarity()) {
						maxSimValue = worker.getMaxSimilarity();
						maxCluster1 = worker.getCluster1();
						maxCluster2 = worker.getCluster2();
						}
					}
				}

			if (maxSimValue < similarityLimit)
				break;

			for (int i=0; i<maxCluster1; i++)	// calculate new weighted similarity values
				if (mNoOfMembers[i] != 0)
					mSimilarityMatrix[maxCluster1][i]
							= (mNoOfMembers[maxCluster1] * mSimilarityMatrix[maxCluster1][i]
							 + mNoOfMembers[maxCluster2] * mSimilarityMatrix[maxCluster2][i])
							/ (mNoOfMembers[maxCluster1] + mNoOfMembers[maxCluster2]);
			for (int i=maxCluster1+1; i<maxCluster2; i++)
				if (mNoOfMembers[i] != 0)
					mSimilarityMatrix[i][maxCluster1]
							= (mNoOfMembers[maxCluster1] * mSimilarityMatrix[i][maxCluster1]
							 + mNoOfMembers[maxCluster2] * mSimilarityMatrix[maxCluster2][i])
							/ (mNoOfMembers[maxCluster1] + mNoOfMembers[maxCluster2]);
			for (int i=maxCluster2+1; i<mNoOfCompounds; i++)
				if (mNoOfMembers[i] != 0)
					mSimilarityMatrix[i][maxCluster1]
							= (mNoOfMembers[maxCluster1] * mSimilarityMatrix[i][maxCluster1]
							 + mNoOfMembers[maxCluster2] * mSimilarityMatrix[i][maxCluster2])
							/ (mNoOfMembers[maxCluster1] + mNoOfMembers[maxCluster2]);

			mNoOfMembers[maxCluster1] += mNoOfMembers[maxCluster2];
			mNoOfMembers[maxCluster2] = 0;

			for (int i=0; i<mNoOfCompounds; i++)
				if (mClusterNo[i] == maxCluster2)
					mClusterNo[i] = maxCluster1;

			mNoOfClusters--;

			if (threadMustDie()) {
			    stopProgress("clustering cancelled");
				return;
				}

//System.out.println("Clusters: "+mNoOfClusters+", similarity value: "+maxSimValue);
		    if (similarityLimit != 0.0)
			    updateProgress((int)(5000.0*(1.0 - maxSimValue)));
			else
			    updateProgress(mNoOfCompounds - mNoOfClusters);
			}

		findRepresentatives();

		mExecutor.shutdown();

		while(!mExecutor.isTerminated()){
			try {Thread.sleep(500);} catch (InterruptedException e) {e.printStackTrace();}
		}

	    stopProgress("clustering finished");
		}


	public boolean isRepresentative(int compound) {
		return mIsRepresentative[compound];
		}


	public int getClusterNo(int compound) {
		return mClusterNo[compound];
		}


	public int getClusterCount() {
		return mNoOfClusters;
		}


    /**
     *  Renumber cluster numbers starting at 1 to eliminate unused numbers.
     */
	public void regenerateClusterNos() {
		int[] newClusterNo = new int[mNoOfCompounds];
		int clusterNo = 1;
		for (int i=0; i<mNoOfCompounds; i++) {
			if (newClusterNo[mClusterNo[i]] == 0)
				newClusterNo[i] = clusterNo++;
			mClusterNo[i] = newClusterNo[mClusterNo[i]];
			}
		}


	private void calculateSimilarityMatrix(boolean withinClustersOnly) {
		startProgress("Calculating Similaries...", 0, 1000);
		if (mThreadCount == 1) {
			for (int compound2=1; compound2<mNoOfCompounds && !threadMustDie(); compound2++) {
				for (int compound1=0; compound1<compound2; compound1++)
					if (!withinClustersOnly || mClusterNo[compound1] == mClusterNo[compound2])
						mSimilarityMatrix[compound2][compound1] = (float)mDescriptorHandler.getSimilarity(mDescriptor[compound1], mDescriptor[compound2]);

				updateProgress((int)(1000.0*compound2*compound2/mNoOfCompounds/mNoOfCompounds));
				}
			}
		else {
			if (withinClustersOnly)
				runInParallel(ClusterWorker.CALC_CLUSTER_SIMILARITIES);
			else
				runInParallel(ClusterWorker.CALC_ALL_SIMILARITIES);
			}
		}

	public float getSimilarity(int index1, int index2){
		if(index1==index2) return 1.0f;
		return mSimilarityMatrix[Math.max(index1,index2)][Math.min(index1,index2)];
	}

	private void findRepresentatives() {
		calculateSimilarityMatrix(true);
		if (threadMustDie())
			return;

		float[] lowSim = new float[mNoOfCompounds];
		float[] simSum = new float[mNoOfCompounds];
		for (int i=0; i<mNoOfCompounds; i++) {
			lowSim[i] = 1;
			simSum[i] = 0;
			}

		startProgress("Locating Representatives...", 0, mNoOfCompounds);
		for (int cluster2=1; cluster2<mNoOfCompounds; cluster2++) {
			if (threadMustDie())
				return;

			updateProgress(cluster2);

			for (int cluster1=0; cluster1<cluster2; cluster1++) {
				if (mClusterNo[cluster1] == mClusterNo[cluster2]) {
					if (lowSim[mClusterNo[cluster1]] > mSimilarityMatrix[cluster2][cluster1])
						lowSim[mClusterNo[cluster1]] = mSimilarityMatrix[cluster2][cluster1];
					simSum[cluster1] += mSimilarityMatrix[cluster2][cluster1];
					simSum[cluster2] += mSimilarityMatrix[cluster2][cluster1];
					}
				}
			}

		int[] representative = new int[mNoOfCompounds];
		for (int i=0; i<mNoOfCompounds; i++)
			representative[i] = -1;

		for (int index=0; index<mNoOfCompounds; index++)	// locate cluster representative compounds
			if (representative[mClusterNo[index]] == -1
			 || simSum[representative[mClusterNo[index]]] < simSum[index])
				representative[mClusterNo[index]] = index;

		mIsRepresentative = new boolean[mNoOfCompounds];
		for (int index=0; index<mNoOfCompounds; index++)
			if (representative[mClusterNo[index]] == index)
				mIsRepresentative[index] = true;
		}


	private void runInParallel(int whatToDo) {
		CountDownLatch doneSignal = new CountDownLatch(mThreadCount);
		for (ClusterWorker<T> worker:mClusterWorker) {
			worker.initJob(whatToDo, doneSignal);
			mExecutor.execute(worker);
			}
		try {
			doneSignal.await();
			}
		catch (InterruptedException e) {}
		}


	private class ClusterWorker<U> implements Runnable {
		private static final int CALC_ALL_SIMILARITIES = 1;
		private static final int CALC_CLUSTER_SIMILARITIES = 2;
		private static final int FIND_MAXIMUM_SIMILARITY = 3;

		private CountDownLatch mDoneSignal;
		private int mWhatToDo,mCluster1,mCluster2;
		private float mMaxSimilarity;
		private DescriptorHandler<T,?> mThreadSafeDH;

		public void initJob(int whatToDo, CountDownLatch doneSignal) {
			mWhatToDo = whatToDo;
			mDoneSignal = doneSignal;
	    	mSMPCompoundIndex = new AtomicInteger(mNoOfCompounds);
			mThreadSafeDH = mDescriptorHandler.getThreadSafeCopy();
			}

		public void run() {
			switch (mWhatToDo) {
			case CALC_ALL_SIMILARITIES:
				int compound2 = mSMPCompoundIndex.decrementAndGet();
				while (compound2 >= 1 && !threadMustDie()) {
					for (int compound1=0; compound1<compound2; compound1++)
						mSimilarityMatrix[compound2][compound1] = (float)mThreadSafeDH.getSimilarity(mDescriptor[compound1], mDescriptor[compound2]);

    				compound2 = mSMPCompoundIndex.decrementAndGet();
    				updateProgress(1000-(int)(1000.0*compound2*compound2/mNoOfCompounds/mNoOfCompounds));
					}
				break;
			case CALC_CLUSTER_SIMILARITIES:
				compound2 = mSMPCompoundIndex.decrementAndGet();
				while (compound2 >= 1 && !threadMustDie()) {
					for (int compound1=0; compound1<compound2; compound1++)
						if (mClusterNo[compound1] == mClusterNo[compound2])
							mSimilarityMatrix[compound2][compound1] = (float)mThreadSafeDH.getSimilarity(mDescriptor[compound1], mDescriptor[compound2]);

    				compound2 = mSMPCompoundIndex.decrementAndGet();
    				updateProgress(1000-(int)(1000.0*compound2*compound2/mNoOfCompounds/mNoOfCompounds));
					}
				break;
			case FIND_MAXIMUM_SIMILARITY:
				mMaxSimilarity = 0;		// find highest similarity level
				mCluster1 = -1;
				mCluster2 = -1;
				int cluster2 = mSMPCompoundIndex.decrementAndGet();
				while (cluster2 >= 1 && !threadMustDie()) {
					if (mNoOfMembers[cluster2] > 0) {
						for (int cluster1=0; cluster1<cluster2; cluster1++) {
							if (mNoOfMembers[cluster1] != 0) {
								if (mMaxSimilarity < mSimilarityMatrix[cluster2][cluster1]) {
									mMaxSimilarity = mSimilarityMatrix[cluster2][cluster1];
									mCluster1 = cluster1;
									mCluster2 = cluster2;
									}
								}
							}
						}
					cluster2 = mSMPCompoundIndex.decrementAndGet();
					}
				break;
				}
			mDoneSignal.countDown();
			}

		public float getMaxSimilarity() {
			return mMaxSimilarity;
			}

		public int getCluster1() {
			return mCluster1;
			}

		public int getCluster2() {
			return mCluster2;
			}
		}
	}
