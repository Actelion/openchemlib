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

package com.actelion.research.calc;

import java.awt.Point;
import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public abstract class SelfOrganizedMap extends DataProcessor {
	private static final int cModeNeighbourhoodMask = 7;
	public static final int cModeNeighbourhoodGaussean = 0;
	public static final int cModeNeighbourhoodMexicanHat = 1;
	public static final int cModeNeighbourhoodLinear = 2;
	public static final int cModeTopologyUnlimited = 8;
	public static final int cModeGrowDuringOptimization = 16;
	public static final int cModeFastBestMatchFinding = 32;

	protected SOMController	mController;
	protected Object[][]	mReferenceVector;

	protected int			mNX,mNY,mMode;
	private int				mCycle,mCyclesPerNode,mConstantInfluenceCycles,mThreadCount,
							mInputVectorCount,mInputVectorIndex,mCorrectQuickBestMatches;
	private boolean			mFindBestMatchQuickly;
	private double			mMaxRange,mDiagonal;
	protected double[][]	mInfluence;
	private Point[]			mLastBestMatch;
	private int[][]			mSMPSOMIndex,mSMPInfluenceIndex;
	private Rectangle		mSMPInfluenceRect;
	private ExecutorService	mExecutor;
	private SOMWorker[]		mSOMWorker;

	/**
	 * Constructor to be used if SOM interna are read from a SOM file with read()
	 */
	public SelfOrganizedMap() {
		initializeSMP();
		}

	public SelfOrganizedMap(int nx, int ny, int mode) {
		initializeSMP();
		initializeReferenceVectors(nx, ny, mode);
		}

	public void initializeReferenceVectors(int nx, int ny, int mode) {
		mNX = nx;
		mNY = ny;
		mMode = mode;
		mReferenceVector = new Object[nx][ny];

		if (mThreadCount != 1)
			mSMPSOMIndex = getSMPArraySplitting(nx, ny);
		}

	private void initializeSMP() {
		mThreadCount = Runtime.getRuntime().availableProcessors();
		if (mThreadCount != 1) {
			mExecutor = Executors.newFixedThreadPool(mThreadCount);
			mSOMWorker = new SOMWorker[mThreadCount];
			for (int t=0; t<mThreadCount; t++)
				mSOMWorker[t] = new SOMWorker(t);
			}
		}

	public void setController(SOMController sc) {
		mController = sc;
		}

	public int getWidth() {
		return mNX;
		}

	public int getHeight() {
		return mNY;
		}

	public int getCreationMode() {
		return mMode;
		}

	public void organize() {
		if (mController.getInputVectorCount() == 0)
			return;

		initializeNormalization();

		mDiagonal = Math.sqrt((double)((mNX-1)*(mNX-1) + (mNY-1)*(mNY-1)));
		mCyclesPerNode = 16;
        mConstantInfluenceCycles = Math.max(1, mNX * mNY * mCyclesPerNode / 2560);

		if ((mMode & cModeGrowDuringOptimization) != 0) {
			mNX /= 8;
			mNY /= 8;
			}

		mReferenceVector = new Object[mNX][mNY];
		for (int x=0; x<mNX; x++)
			for (int y=0; y<mNY; y++)
				mReferenceVector[x][y] = getRandomVector();

		mInputVectorCount = mController.getInputVectorCount();

			// don't use fast best match finding if input vectors are randomly generated
		if (mInputVectorCount == -1)
			mMode &= ~cModeFastBestMatchFinding;

		if ((mMode & cModeFastBestMatchFinding) != 0) {
			mLastBestMatch = new Point[mInputVectorCount];
			mFindBestMatchQuickly = false;	// start out with slow full search
			}

		if ((mMode & cModeGrowDuringOptimization) != 0) {
			int cyclesPhaseOne = mCyclesPerNode * mNX * mNY;
			int overallCycles = cyclesPhaseOne
							  + cyclesPhaseOne / 3 * 4
							  + cyclesPhaseOne / 3 * 16
							  + cyclesPhaseOne / 3 * 64;
			mCycle = 0;
			startProgress("Self-Organizing map...", 0, overallCycles);

			optimize(0, cyclesPhaseOne);
			int cyclesInPhase = cyclesPhaseOne;
			for (int phase=2; phase<=4; phase++) {
				grow();
				cyclesInPhase *= 4;
				optimize(cyclesInPhase - cyclesInPhase / 3, cyclesInPhase);
				}
			stopProgress("Map completed.");
			}
		else {
			int overallCycles = mCyclesPerNode * mNX * mNY;
			mCycle = 0;

			startProgress("Self-Organizing map...", 0, overallCycles);
			optimize(0, overallCycles);
			stopProgress("Map completed.");
			}
		}

	public Object getReferenceVector(int x, int y) {
		return mReferenceVector[x][y];
		}

	public double[][] getInfluence() {
		return mInfluence;
		}

	public double getChaos() {
		double sum = 0.0;
		for (int x=0; x<mNX; x++) {
			for (int y=0; y<mNY; y++) {
				sum += Math.sqrt((x == 0) ? getDissimilarity(mReferenceVector[0][y], mReferenceVector[mNX-1][y])
										  : getDissimilarity(mReferenceVector[x][y], mReferenceVector[x-1][y]));
				sum += Math.sqrt((y == 0) ? getDissimilarity(mReferenceVector[x][0], mReferenceVector[x][mNY-1])
										  : getDissimilarity(mReferenceVector[x][y], mReferenceVector[x][y-1]));
				}
			}
		return sum / (double)(mNX * mNY);
		}

	public double getMatchScore() {
		final int randomCount = 1;
		double sum = 0.0;
		for (int i=0; i<randomCount; i++) {
			Object randomVector = getRandomVector();
			Object bestMatch = null;
			double minDissimilarity = Double.POSITIVE_INFINITY;
			for (int x=0; x<mNX; x++) {
				for (int y=0; y<mNY; y++) {
					double dissimilarity = getDissimilarity(mReferenceVector[x][y], randomVector);
					if (minDissimilarity > dissimilarity) {
						minDissimilarity = dissimilarity;
						bestMatch = mReferenceVector[x][y];
						}
					}
				}
			sum += Math.sqrt(getDissimilarity(randomVector, bestMatch));
			}
		return sum / (double)randomCount;
		}

	protected double getTimeInfluence(double time) {
		return 1.0 - time;
		}

	protected double getNeighbourInfluence(int dx, int dy, double time) {
											// distance is normalized by map diagonal
		double distance = Math.sqrt(dx*dx+dy*dy)/mDiagonal;
		double f = 0.0;

		switch (mMode & cModeNeighbourhoodMask) {
		case cModeNeighbourhoodGaussean:
			f = Math.exp(distance * distance * Math.log(0.001) / (mMaxRange * mMaxRange));
			return (f < 0.001) ? 0.0 : f;
		case cModeNeighbourhoodMexicanHat:
			f = 1.0 - distance * distance / (mMaxRange * mMaxRange);
			return (f < 0.0) ? 0.0 : f;
		case cModeNeighbourhoodLinear:
			f = 1.0 - distance / mMaxRange;
			return (f < 0.0) ? 0.0 : f;
			}

		return f;
		}

	/**
	 * Calculates start and end array indices for individual threads.
	 * @param xCount
	 * @param yCount
	 */
	private int[][] getSMPArraySplitting(int xCount, int yCount) {
		int nodeCount = xCount * yCount;
		int nodesPerThread = nodeCount / mThreadCount;
		int remainingNodes = nodeCount % mThreadCount;
		int[][] index = new int[mThreadCount+1][2];
		for (int i=0; i<mThreadCount; i++) {
			int nodes = (i < remainingNodes) ? nodesPerThread + 1 : nodesPerThread;
			int dx = nodes % xCount;
			int dy = nodes / xCount;
			if (index[i][0] + dx >= xCount) {
				dx -= xCount;
				dy++;
				}
			index[i+1][0] = index[i][0] + dx;
			index[i+1][1] = index[i][1] + dy;
			}
		return index;
		}

	private void optimize(int startCycle, int cycles) {
		mInfluence = null;

		for (int cycle=startCycle; cycle<cycles; cycle++) {
			if (((cycle-startCycle) % mConstantInfluenceCycles) == 0)
				calculateInfluences((double)cycle/(double)cycles);

			updateProgress(mCycle++);
			if (threadMustDie())
				break;

			Object inputVector = normalizeVector(mController.getInputVector(mInputVectorIndex));

			if ((mMode & cModeFastBestMatchFinding) != 0) {
				if (mFindBestMatchQuickly)
					mLastBestMatch[mInputVectorIndex] = findBestMatchLocationQuickly(inputVector);
				else {
					Point quickBestMatch = (mLastBestMatch[mInputVectorIndex] == null) ?
							null : findBestMatchLocationQuickly(inputVector);

					mLastBestMatch[mInputVectorIndex] = findBestMatchLocation(inputVector);

						// switch to quick best match search after 2*mInputVectorCount correct quick best match runs
					if (quickBestMatch != null
					 && quickBestMatch.equals(mLastBestMatch[mInputVectorIndex]))
						mCorrectQuickBestMatches++;
					else
						mCorrectQuickBestMatches = 0;
					if (mCorrectQuickBestMatches == 2*mInputVectorCount)
						mFindBestMatchQuickly = true;
//System.out.println("mCycle:"+mCycle+", mCorrectQuickBestMatches:"+mCorrectQuickBestMatches+"["+mLastBestMatch[mInputVectorIndex].x+","+mLastBestMatch[mInputVectorIndex].y+"]");
					}

				applyInfluences(inputVector, mLastBestMatch[mInputVectorIndex]);
				}
			else {
				applyInfluences(inputVector, findBestMatchLocation(inputVector));
				}

				// mInputVectorCount may be -1 if input vectors are generated on the fly
			if (++mInputVectorIndex >= mInputVectorCount)
				mInputVectorIndex = 0;
			}
		}

	protected void applyInfluences(Object inputVector, Point location) {
			// apply to every node shifts caused by assigned input vector
		int maxRange = (int)(mDiagonal * mMaxRange);
		int x1 = location.x-maxRange;
		int x2 = location.x+maxRange;
		int y1 = location.y-maxRange;
		int y2 = location.y+maxRange;
		if ((mMode & cModeTopologyUnlimited) != 0) {
			if (x2 - x1 >= mNX) {
				x1 = 0;
				x2 = mNX;
				}
			if (y2 - y1 >= mNY) {
				y1 = 0;
				y2 = mNY;
				}

			if (mThreadCount != 1 && maxRange != 0) {
				if (mSMPInfluenceRect == null || mSMPInfluenceRect.width != x2-x1 || mSMPInfluenceRect.height != y2-y1)
					mSMPInfluenceIndex = getSMPArraySplitting(x2-x1, y2-y1);
				mSMPInfluenceRect = new Rectangle(x1, y1, x2-x1, y2-y1);
				applyInfluencesSMP(inputVector, location);
				}
			else {
				for (int x=x1; x<x2; x++) {
					for (int y=y1; y<y2; y++) {
						int dx = Math.abs(location.x - x);
						int dy = Math.abs(location.y - y);
						if (dx > mNX / 2)
							dx = mNX - dx;
						if (dy > mNY / 2)
							dy = mNY - dy;
						if (mInfluence[dx][dy] > 0.0)
							updateReference(inputVector, mReferenceVector[x<0?x+mNX:x<mNX?x:x-mNX][y<0?y+mNY:y<mNY?y:y-mNY], mInfluence[dx][dy]);
						}
					}
				}
			}
		else {
			if (x1 < 0)
				x1 = 0;
			if (x2 > mNX)
				x2 = mNX;
			if (y1 < 0)
				y1 = 0;
			if (y2 > mNY)
				y2 = mNY;

			if (mThreadCount != 1 && maxRange != 0) {
				mSMPInfluenceRect = new Rectangle(x1, y1, x2-x1, y2-y1);
				mSMPInfluenceIndex = getSMPArraySplitting(x2-x1, y2-y1);
				applyInfluencesSMP(inputVector, location);
				}
			else {
				for (int x=x1; x<x2; x++) {
					for (int y=y1; y<y2; y++) {
						int dx = Math.abs(location.x - x);
						int dy = Math.abs(location.y - y);
						if (mInfluence[dx][dy] > 0.0)
							updateReference(inputVector, mReferenceVector[x][y], mInfluence[dx][dy]);
						}
					}
				}
			}
		}

	protected void applyInfluencesSMP(Object inputVector, Point location) {
		CountDownLatch doneSignal = new CountDownLatch(mThreadCount);
		for (SOMWorker worker:mSOMWorker) {
			worker.initApplyInfluences(inputVector, location, doneSignal);
			mExecutor.execute(worker);
			}
		try {
			doneSignal.await();
			}
		catch (InterruptedException e) {}
		}

	public Point findBestMatchLocation(Object inputVector) {
		if (inputVector == null)
			return null;

		// used internally during map optimization
			//  and externally to assign any vector to location on completed map
		if (mThreadCount != 1)
			return findBestMatchLocationSMP(inputVector);

		Point minLocation = new Point(-1, -1);
		double minDissimilarity = Double.POSITIVE_INFINITY;
		for (int x=0; x<mNX; x++) {
			for (int y=0; y<mNY; y++) {
				double dissimilarity = getDissimilarity(mReferenceVector[x][y], inputVector);
				if (minDissimilarity > dissimilarity) {
					minDissimilarity = dissimilarity;
					minLocation.x = x;
					minLocation.y = y;
					}
				}
			}
		return minLocation;
		}

	public Point findBestMatchLocationSMP(Object inputVector) {
		if (inputVector == null)
			return null;

		CountDownLatch doneSignal = new CountDownLatch(mThreadCount);
		for (SOMWorker worker:mSOMWorker) {
			worker.initFindBestMatch(inputVector, doneSignal);
			mExecutor.execute(worker);
			}
		try {
			doneSignal.await();
			}
		catch (InterruptedException e) {}

	    Point minLocation = new Point(-1, -1);
		double minDissimilarity = Double.POSITIVE_INFINITY;
		for (SOMWorker worker:mSOMWorker) {
			double dissimilarity = worker.getBestMatchDissimilarity();
			if (minDissimilarity > dissimilarity) {
				minDissimilarity = dissimilarity;
				minLocation = worker.getBestMatchLocation();
				}
			}

		return minLocation;
		}

	public Point findBestMatchLocationQuickly(Object inputVector) {
		if (inputVector == null)
			return null;

		// used internally during map optimization
			//  and externally to assign any vector to location on completed map
		Point minLocation = mLastBestMatch[mInputVectorIndex];
		double minDissimilarity = getDissimilarity(mReferenceVector[minLocation.x][minLocation.y], inputVector);
		boolean[][] locationChecked = new boolean[mNX][mNY];
		locationChecked[minLocation.x][minLocation.y] = true;
		boolean found;
		do {
//System.out.print("["+minLocation.x+","+minLocation.y+"];");
			Point p = minLocation;
			found = false;
			for (int xdif=-1; xdif<2; xdif++) {
				int x = p.x + xdif;
				if ((mMode & cModeTopologyUnlimited) != 0) {
					if (x < 0)
						x = mNX - 1;
					else if (x >= mNX)
						x = 0;
					}
				else if (x < 0 || x >= mNX)
					continue;
				for (int ydif=-1; ydif<2; ydif++) {
					int y = p.y + ydif;
					if ((mMode & cModeTopologyUnlimited) != 0) {
						if (y < 0)
							y = mNY - 1;
						else if (y >= mNY)
							y = 0;
						}
					else if (y < 0 || y >= mNY)
						continue;

					if (!locationChecked[x][y]) {
						locationChecked[x][y] = true;
						double dissimilarity = getDissimilarity(mReferenceVector[x][y], inputVector);
						if (minDissimilarity > dissimilarity) {
							minDissimilarity = dissimilarity;
							minLocation = new Point(x, y);
							found = true;
							}
						}
					}
				}
			} while (found);

//System.out.println();
		return minLocation;
		}

	public double[] findExactMatchLocation(Object inputVector) {
		if (inputVector == null)
			return null;

		Point p = findBestMatchLocation(inputVector);
		double[] location = new double[3];
		location[0] = p.x;
		location[1] = p.y;

		int x1 = p.x - 1;
		int x2 = p.x + 1;
		int y1 = p.y - 1;
		int y2 = p.y + 1;
		if ((mMode & cModeTopologyUnlimited) != 0) {
			if (x1 == -1)
				x1 += mNX;
			if (x2 == mNX)
				x2 = 0;
			if (y1 == -1)
				y1 += mNY;
			if (y2 == mNY)
				y2 = 0;
			}

		double dis0 = Math.sqrt(getDissimilarity(mReferenceVector[p.x][p.y], inputVector));
		if (dis0 > 0.0) {
			if (x1 == -1)
				location[0] += 0.5 * dis0 / Math.sqrt(getDissimilarity(mReferenceVector[x2][p.y], inputVector));
			else if (x2 == mNX)
				location[0] -= 0.5 * dis0 / Math.sqrt(getDissimilarity(mReferenceVector[x1][p.y], inputVector));
			else {
				double dis1 = Math.sqrt(getDissimilarity(mReferenceVector[x1][p.y], inputVector)) - dis0;
				double dis2 = Math.sqrt(getDissimilarity(mReferenceVector[x2][p.y], inputVector)) - dis0;
				if (dis1 + dis2 != 0)
					location[0] += dis1 / (dis1 + dis2) - 0.5;
				}

			if (y1 == -1)
				location[1] += 0.5 * dis0 / Math.sqrt(getDissimilarity(mReferenceVector[p.x][y2], inputVector));
			else if (y2 == mNY)
				location[1] -= 0.5 * dis0 / Math.sqrt(getDissimilarity(mReferenceVector[p.x][y1], inputVector));
			else {
				double dis1 = Math.sqrt(getDissimilarity(mReferenceVector[p.x][y1], inputVector)) - dis0;
				double dis2 = Math.sqrt(getDissimilarity(mReferenceVector[p.x][y2], inputVector)) - dis0;
				if (dis1 + dis2 != 0)
					location[1] += dis1 / (dis1 + dis2) - 0.5;
				}
			}

		location[2] = dis0;	// dissimilarity of match
		return location;
		}

	private void grow() {
			// works currently only with unlimited SOMs
		Object[][] oldReferenceVector = mReferenceVector;
		mReferenceVector = new Object[mNX*2][mNY*2];
		for (int x=0; x<mNX; x++) {
			for (int y=0; y<mNY; y++) {
				int nextX = (x == mNX-1) ? 0 : x+1;
				int nextY = (y == mNY-1) ? 0 : y+1;
				mReferenceVector[x*2][y*2] = oldReferenceVector[x][y];
				mReferenceVector[x*2+1][y*2] = getMeanVector(
								oldReferenceVector[x][y], oldReferenceVector[nextX][y]);
				mReferenceVector[x*2][y*2+1] = getMeanVector(
								oldReferenceVector[x][y], oldReferenceVector[x][nextY]);
				mReferenceVector[x*2+1][y*2+1] = getMeanVector(
								oldReferenceVector[x][y], oldReferenceVector[nextX][nextY]);
				}
			}
		mNX *= 2;
		mNY *= 2;

		if ((mMode & cModeFastBestMatchFinding) != 0) {
			for (int i=0; i<mLastBestMatch.length; i++) {
				if (mLastBestMatch[i] != null) {
					mLastBestMatch[i].x *= 2;
					mLastBestMatch[i].y *= 2;
					}
				}
			}
		}

	protected void calculateInfluences(double time) {
		final double cStartRange = 1.0;		// 1.0 := diagonal of entire map
		final double cFinalRange = 0.05;	// 0.05
		mMaxRange = cStartRange * Math.exp(Math.log(cFinalRange/cStartRange)*Math.pow(time, 0.3));
		double timeInfluence = getTimeInfluence(time);
		int dxmax = Math.min(1+(int)(mDiagonal * mMaxRange), mNX);
		int dymax = Math.min(1+(int)(mDiagonal * mMaxRange), mNY);
		mInfluence = new double[dxmax][dymax];
		for (int dx=0; dx<dxmax; dx++)
			for (int dy=0; dy<dymax; dy++)
				mInfluence[dx][dy] = 0.5 * timeInfluence * getNeighbourInfluence(dx, dy, time);
		}

	public void write(BufferedWriter writer) throws IOException {
		writer.write("<width=\""+mNX+"\">");
		writer.newLine();
		writer.write("<height=\""+mNY+"\">");
		writer.newLine();
		writer.write("<creationMode=\""+mMode+"\">");
		writer.newLine();

		startProgress("Writing SOM Vectors...", 0, mNY);
		for (int y=0; y<mNY; y++) {
			updateProgress(y);
			for (int x=0; x<mNX; x++) {
				writer.write("<reference["+x+"]["+y+"]=\""+referenceVectorToString(x, y)+"\">");
				writer.newLine();
				}
			}
		stopProgress("SOM Vectors Written");
		}

	public void read(BufferedReader reader) throws Exception {
		String theLine = reader.readLine();
		boolean error = !theLine.startsWith("<width=");
		if (!error) {
			mNX = Integer.parseInt(extractValue(theLine));
			theLine = reader.readLine();
			error = !theLine.startsWith("<height=");
			}

		if (!error) {
			mNY = Integer.parseInt(extractValue(theLine));
			theLine = reader.readLine();
			error = !theLine.startsWith("<creationMode=");
			}

		if (!error) {
			mMode = Integer.parseInt(extractValue(theLine));
			}

		if (!error && mThreadCount != 1)
			mSMPSOMIndex = getSMPArraySplitting(mNX, mNY);

		mReferenceVector = new Object[mNX][mNY];
		startProgress("Reading SOM Vectors...", 0, mNY);
		for (int y=0; y<mNY && !error; y++) {
			updateProgress(y);
			for (int x=0; x<mNX && !error; x++) {
				theLine = reader.readLine();
				error = !theLine.startsWith("<reference["+x+"]["+y+"]=");
				if (!error)
					setReferenceVector(x, y, extractValue(theLine));
				}
			}

		stopProgress("SOM Vectors Reading Done");
		if (error)
			throw new IOException("Invalid SOM file format");
		}

	public static String extractValue(String theLine) {
		int index1 = theLine.indexOf("=\"") + 2;
		int index2 = theLine.indexOf("\"", index1);
		return theLine.substring(index1, index2);
		}


	protected abstract String referenceVectorToString(int x, int y);
	protected abstract void setReferenceVector(int x, int y, String ref) throws Exception;
	protected abstract void initializeNormalization();
	protected abstract void updateReference(Object inputVector, Object referenceVector, double influence);
	protected abstract Object getMeanVector(Object vector1, Object vector2);
	protected abstract Object getRandomVector();
	public abstract Object normalizeVector(Object vector);
	public abstract double getDissimilarity(Object vector1, Object vector2);

	private class SOMWorker implements Runnable {
		private static final int FIND_BEST_MATCH = 1;
		private static final int APPLY_INFLUENCES = 2;

		private CountDownLatch mDoneSignal;
		private int mThreadIndex;
		private Object mInputVector;
		private int mWhatToDo;
	    private Point mLocation;
		private double mMinDissimilarity;

		private SOMWorker(int threadIndex) {
			mThreadIndex = threadIndex;
			}

		public void initFindBestMatch(Object inputVector, CountDownLatch doneSignal) {
			mWhatToDo = FIND_BEST_MATCH;
			mInputVector = inputVector;
			mDoneSignal = doneSignal;
			}

		public void initApplyInfluences(Object inputVector, Point location, CountDownLatch doneSignal) {
			mWhatToDo = APPLY_INFLUENCES;
			mInputVector = inputVector;
			mLocation = location;
			mDoneSignal = doneSignal;
			}

		public void run() {
			switch (mWhatToDo) {
			case FIND_BEST_MATCH:
			    mLocation = new Point(-1, -1);
				mMinDissimilarity = Double.POSITIVE_INFINITY;
				int y1 = mSMPSOMIndex[mThreadIndex][1];
				int y2 = mSMPSOMIndex[mThreadIndex+1][1];
				for (int y=y1; y<=y2; y++) {
					int x1 = (y == y1) ? mSMPSOMIndex[mThreadIndex][0] : 0;
					int x2 = (y == y2) ? mSMPSOMIndex[mThreadIndex+1][0] : mNX;
					for (int x=x1; x<x2; x++) {
						double dissimilarity = getDissimilarity(mReferenceVector[x][y], mInputVector);
						if (mMinDissimilarity > dissimilarity) {
							mMinDissimilarity = dissimilarity;
							mLocation.x = x;
							mLocation.y = y;
							}
						}
					}
				break;
			case APPLY_INFLUENCES:
				if ((mMode & cModeTopologyUnlimited) != 0) {
					y1 = mSMPInfluenceIndex[mThreadIndex][1];
					y2 = mSMPInfluenceIndex[mThreadIndex+1][1];
					for (int yy=y1; yy<=y2; yy++) {
						int x1 = (yy == y1) ? mSMPInfluenceIndex[mThreadIndex][0] : 0;
						int x2 = (yy == y2) ? mSMPInfluenceIndex[mThreadIndex+1][0] : mSMPInfluenceRect.width;
						int y = yy + mSMPInfluenceRect.y;
						for (int xx=x1; xx<x2; xx++) {
							int x = xx + mSMPInfluenceRect.x;
							int dx = Math.abs(mLocation.x - x);
							int dy = Math.abs(mLocation.y - y);
							if (dx > mNX / 2)
								dx = mNX - dx;
							if (dy > mNY / 2)
								dy = mNY - dy;
							if (mInfluence[dx][dy] > 0.0)
								updateReference(mInputVector, mReferenceVector[x<0?x+mNX:x<mNX?x:x-mNX][y<0?y+mNY:y<mNY?y:y-mNY], mInfluence[dx][dy]);
							}
						}
					}
				else {
					y1 = mSMPInfluenceIndex[mThreadIndex][1];
					y2 = mSMPInfluenceIndex[mThreadIndex+1][1];
					for (int yy=y1; yy<=y2; yy++) {
						int x1 = (yy == y1) ? mSMPInfluenceIndex[mThreadIndex][0] : 0;
						int x2 = (yy == y2) ? mSMPInfluenceIndex[mThreadIndex+1][0] : mSMPInfluenceRect.width;
						int y = yy + mSMPInfluenceRect.y;
						for (int xx=x1; xx<x2; xx++) {
							int x = xx + mSMPInfluenceRect.x;
							int dx = Math.abs(mLocation.x - x);
							int dy = Math.abs(mLocation.y - y);
							if (mInfluence[dx][dy] > 0.0)
								updateReference(mInputVector, mReferenceVector[x][y], mInfluence[dx][dy]);
							}
						}
					}
				break;
				}
			mDoneSignal.countDown();
			}

		public Point getBestMatchLocation() {
			return mLocation;
			}

		public double getBestMatchDissimilarity() {
			return mMinDissimilarity;
			}
		}
	}
