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

import com.actelion.research.calc.ProgressController;
import com.actelion.research.chem.descriptor.*;
import com.actelion.research.util.ByteArrayComparator;

import java.nio.charset.StandardCharsets;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

public class StructureSearch {
	public static final int SEARCH_RUNNING = -1;
	public static final int SEARCH_PENDING = 0;
	public static final int SEARCH_STOPPED = 1;
	public static final int QUERY_MISSING = 2;
	public static final int SEARCH_TYPE_NOT_SUPPORTED = 3;
	public static final int SUCCESSFUL_COMPLETION = 4;
	public static final int COUNT_LIMIT_EXCEEDED = 5;
	public static final int TIME_LIMIT_EXCEEDED = 6;
	public static final String[] COMPLETION_TEXT = { "not started", "stopped", "query missing", "unsupported search type", "successful", "count limit hit", "time limit hit" };

	private final StructureSearchSpecification mSpecification;
	private final StructureSearchDataSource mDataSource;
	private final StructureSearchController mSearchController;
	private final ProgressController mProgressController;
	private volatile StereoMolecule[] mQueryFragment,mDoubleQueryFragment;
	private volatile ByteArrayComparator mIDCodeComparator;
	private volatile DescriptorHandler mDescriptorHandler;
	private volatile Object[] mQueryDescriptor;
	private volatile long[] mQueryHashCode;
	private volatile byte[][] mQueryIDCode;
	private volatile int mDescriptorColumn;
	private volatile int mMaxSSSMatches,mMaxNonSSSMatches, mStatus;
	private volatile long mStopTime,mMaxMillis;
	private ConcurrentLinkedQueue<Integer> mResultQueue;
	private AtomicInteger mSMPIndex,mMatchCount;

	/**
	 * This contructs a new structure search, which upon calling start()
	 * runs a multithreaded structure search on the structure rows provided by dataSource.
	 * If a searchController is given, this is asked for every row, whether the row
	 * meets all preconditions and qualifies for the search. 
	 * @param specification
	 * @param dataSource
	 * @param searchController may be null, if all rows need to be searched
	 * @param progressController may be null
	 * @param dhFactory if null then the default DescriptorHandlerStandard2DFactory is used
	 */
	public StructureSearch(StructureSearchSpecification specification,
						   StructureSearchDataSource dataSource,
						   StructureSearchController searchController,
						   ProgressController progressController,
						   DescriptorHandlerFactory dhFactory) {
		mSpecification = specification;
		mDataSource = dataSource;
		mSearchController = searchController;
		mProgressController = progressController;
		mStatus = SEARCH_PENDING;

		if (mSpecification != null) {
			// define needed descriptor handlers
			if (mSpecification.isSimilaritySearch()) {
				DescriptorHandlerFactory factory = (dhFactory != null) ? dhFactory : DescriptorHandlerStandard2DFactory.getFactory();
				mDescriptorHandler = factory.getDefaultDescriptorHandler(specification.getDescriptorShortName());
				}
			else if (mSpecification.isSubstructureSearch()) {
				mDescriptorHandler = DescriptorHandlerLongFFP512.getDefaultInstance();
				}
			}
		}

	/**
	 * If the search shall be aborted once it exceeds a given number of matches,
	 * then define the maximum number of matches with this method before starting the search.
	 * In case a search would return more than the defined maximum of allowed matches,
	 * then the search would stop at the allowed maximum and return those matches.
	 * @param maxSSSMatches maximum number of allowed sub-reaction/retron search matches (0: no limit)
	 * @param maxNonSSSMatches maximum number of allowed matches for other search types (0: no limit)
	 */
	public void setMatchLimit(int maxSSSMatches, int maxNonSSSMatches) {
		mMaxSSSMatches = maxSSSMatches;
		mMaxNonSSSMatches = maxNonSSSMatches;
		}

	/**
	 * If the search shall be aborted once it exceeds a given elapsed time limit,
	 * then define the maximum allowed search time in milliseconds.
	 * If a search time limit is reached, then the search would return all matches found.
	 * @param maxMillis maximum allowed elapsed search milliseconds (0: no limit)
	 */
	public void setTimeLimit(long maxMillis) {
		mMaxMillis = maxMillis;
		}

	public String getCompletionStatus() {
		return COMPLETION_TEXT[mStatus];
		}

	public int[] start() {
		if (!mDataSource.isSupportedSearchType(mSpecification)) {
			mStatus = SEARCH_TYPE_NOT_SUPPORTED;
			return null;
			}

		mMatchCount = new AtomicInteger(0);

		if (!mSpecification.isNoStructureSearch()) {
			final int queryStructureCount = mSpecification.getStructureCount();
			if (queryStructureCount == 0) {
				mStatus = QUERY_MISSING;
				return null;
				}

			mDescriptorColumn = -1;
	        boolean largestFragmentOnly = mSpecification.isLargestFragmentOnly();

			if (mSpecification.isSubstructureSearch() || mSpecification.isSimilaritySearch()) {
				if (mSpecification.isSubstructureSearch()) {
					mDescriptorColumn = mDataSource.getDescriptorColumn(DescriptorConstants.DESCRIPTOR_FFP512.shortName);
					mQueryFragment = new StereoMolecule[queryStructureCount];
					for (int i=0; i<queryStructureCount; i++) {
						mQueryFragment[i] = new IDCodeParser(false).getCompactMolecule(mSpecification.getIDCode(i));
						mQueryFragment[i].ensureHelperArrays(Molecule.cHelperParities);
						}
					if (mSpecification.isSingleMatchOnly()) {
						mDoubleQueryFragment = new StereoMolecule[queryStructureCount];
						for (int i=0; i<queryStructureCount; i++) {
							mDoubleQueryFragment[i] = new StereoMolecule(mQueryFragment[i].getAtoms(), mQueryFragment[i].getBonds());
							mDoubleQueryFragment[i].addMolecule(mQueryFragment[i]);
							mDoubleQueryFragment[i].addMolecule(mQueryFragment[i]);
							mDoubleQueryFragment[i].ensureHelperArrays(Molecule.cHelperParities);
							}
						}
					}
				else {
					final String descriptorShortName = mSpecification.getDescriptorShortName();
					mDescriptorColumn = (descriptorShortName == null) ? -1 : mDataSource.getDescriptorColumn(descriptorShortName);
					}

				mQueryDescriptor = new Object[queryStructureCount];
				boolean missingDescriptorFound = false;
				for (int i=0; i<queryStructureCount; i++) {
					mQueryDescriptor[i] = mSpecification.getDescriptor(i);
					if (mQueryDescriptor[i] == null)
						missingDescriptorFound = true;
					}

				if (missingDescriptorFound)
					calculateQueryDescriptorsAndWait();
				}
			else if (mSpecification.isExactSearch()) {
				mIDCodeComparator = new ByteArrayComparator();
				mQueryIDCode = new byte[queryStructureCount][];
				for (int i=0; i<queryStructureCount; i++) {
					if (largestFragmentOnly) {
						StereoMolecule query = new IDCodeParser(true).getCompactMolecule(mSpecification.getIDCode(i));
						mQueryIDCode[i] = CanonizerUtil.getIDCode(query, CanonizerUtil.IDCODE_TYPE.NORMAL, largestFragmentOnly).getBytes(StandardCharsets.UTF_8);
						}
					else {
						mQueryIDCode[i] = mSpecification.getIDCode(i);
						}
					}
				}
			else if (mSpecification.isNoStereoSearch()) {
				mQueryHashCode = new long[queryStructureCount];
				for (int i=0; i<queryStructureCount; i++)
					mQueryHashCode[i] = CanonizerUtil.getNoStereoHash(
							new IDCodeParser(false).getCompactMolecule(mSpecification.getIDCode(i)), largestFragmentOnly);
				}
			else if (mSpecification.isTautomerSearch()) {
				mQueryHashCode = new long[queryStructureCount];
				for (int i=0; i<queryStructureCount; i++)
					mQueryHashCode[i] = CanonizerUtil.getTautomerHash(
							new IDCodeParser(true).getCompactMolecule(mSpecification.getIDCode(i)), largestFragmentOnly);
				}
			else if (mSpecification.isNoStereoTautomerSearch()) {
				mQueryHashCode = new long[queryStructureCount];
				for (int i=0; i<queryStructureCount; i++)
					mQueryHashCode[i] = CanonizerUtil.getNoStereoTautomerHash(
							new IDCodeParser(false).getCompactMolecule(mSpecification.getIDCode(i)), largestFragmentOnly);
				}
			else if (mSpecification.isBackboneSearch()) {
				mQueryHashCode = new long[queryStructureCount];
				for (int i=0; i<queryStructureCount; i++)
					mQueryHashCode[i] = CanonizerUtil.getBackboneHash(
							new IDCodeParser(false).getCompactMolecule(mSpecification.getIDCode(i)), largestFragmentOnly);
				}
			}

    	mSMPIndex = new AtomicInteger(mDataSource.getRowCount());

    	mResultQueue = new ConcurrentLinkedQueue<>();

		if (mProgressController != null && mSpecification.getStructureCount() > 1023)
			mProgressController.startProgress("Searching structures", 0, mSpecification.getStructureCount());

		mStopTime = (mMaxMillis == 0) ? Long.MAX_VALUE : System.currentTimeMillis() + mMaxMillis;
		mStatus = SEARCH_RUNNING;

		int threadCount = Runtime.getRuntime().availableProcessors();
    	SearchThread[] t = new SearchThread[threadCount];
    	for (int i=0; i<threadCount; i++) {
    		t[i] = new SearchThread("Structure Search "+(i+1));
    		t[i].setPriority(Thread.MIN_PRIORITY);
    		t[i].start();
    		}

    	// the controller thread must wait until all others are finished
    	// before the next task can begin or the dialog is closed
    	for (int i=0; i<threadCount; i++)
    		try { t[i].join(); } catch (InterruptedException e) {}

		if (mStatus == SEARCH_RUNNING)
			mStatus = SUCCESSFUL_COMPLETION;

		int[] result = new int[mResultQueue.size()];
    	int i=0;
    	for (Integer integer:mResultQueue)
    		result[i++] = integer;

    	return result;
		}

	private void calculateQueryDescriptorsAndWait() {
    	mSMPIndex = new AtomicInteger(mQueryDescriptor.length);
		int threadCount = Math.min(mQueryDescriptor.length, Runtime.getRuntime().availableProcessors());
    	Thread[] t = new Thread[threadCount];
    	for (int i=0; i<threadCount; i++) {
    		t[i] = new Thread("Query Descriptor Calculation "+(i+1)) {
    			public void run() {
    				while (true) {
        				int index = mSMPIndex.decrementAndGet();
        				if (index < 0)
        					break;

        				StereoMolecule mol = new IDCodeParser(false).getCompactMolecule(mSpecification.getIDCode(index));
        				mQueryDescriptor[index] = mDescriptorHandler.createDescriptor(mol);
    					}
    				}
    			};
    		t[i].setPriority(Thread.MIN_PRIORITY);
    		t[i].start();
    		}

    	for (int i=0; i<threadCount; i++)
    		try { t[i].join(); } catch (InterruptedException e) {}
		}

	private class SearchThread extends Thread {
		private SSSearcherWithIndex mSSSearcher;

		public SearchThread(String name) {
			super(name);
			if (mSpecification.isSubstructureSearch())
				mSSSearcher = new SSSearcherWithIndex();
			}

		public void run() {
			int row = mSMPIndex.decrementAndGet();
			while (row >= 0) {
				if ((mProgressController != null && mProgressController.threadMustDie())) {
					mStatus = SEARCH_STOPPED;
					break;
					}

				if (System.currentTimeMillis() > mStopTime) {
					mStatus = TIME_LIMIT_EXCEEDED;
					break;
					}

				if (mProgressController != null && row%1024==1023)
					mProgressController.updateProgress(mSpecification.getStructureCount()-row);

				if (mSearchController == null || mSearchController.rowQualifies(row)) {
					boolean isMatch = false;

					if (mSpecification.isSubstructureSearch()) {
						if (mMaxSSSMatches != 0 && mMatchCount.get() > mMaxSSSMatches) {
							mStatus = COUNT_LIMIT_EXCEEDED;
							break;
							}

						for (int s=0; !isMatch && s<mDataSource.getStructureCount(row); s++) {
							mSSSearcher.setMolecule(mDataSource.getIDCode(row, s, false), (long[])mDataSource.getDescriptor(mDescriptorColumn, row, s, false));
							for (int i=0; i<mQueryFragment.length; i++) {
								mSSSearcher.setFragment(mQueryFragment[i], (long[])mQueryDescriptor[i]);
								if (mSSSearcher.isFragmentInMolecule()) {
									if (mSpecification.isSingleMatchOnly()) {
										mSSSearcher.setFragment(mDoubleQueryFragment[i], (long[])mQueryDescriptor[i]);
										if (!mSSSearcher.isFragmentInMolecule()) {
											isMatch = true;
											break;
											}
										}
									else {
										isMatch = true;
										break;
										}
									}
								}
							}
						}
					else {
						if (mMaxNonSSSMatches != 0 && mMatchCount.get() > mMaxNonSSSMatches) {
							mStatus = COUNT_LIMIT_EXCEEDED;
							break;
							}

						if (mSpecification.isNoStructureSearch()) {
							isMatch = true;
							}
						else if (mSpecification.isSimilaritySearch()) {
							for (int s=0; !isMatch && s<mDataSource.getStructureCount(row); s++) {
								for (Object o : mQueryDescriptor) {
									if (mDescriptorHandler.getSimilarity(o, mDataSource.getDescriptor(mDescriptorColumn, row, s, mSpecification.isLargestFragmentOnly()))
											>=mSpecification.getSimilarityThreshold()) {
										isMatch = true;
										break;
									}
								}
								}
							}
						else if (mSpecification.isExactSearch()) {
							for (int s=0; !isMatch && s<mDataSource.getStructureCount(row); s++) {
								for (byte[] bytes : mQueryIDCode) {
									if (mIDCodeComparator.compare(bytes, mDataSource.getIDCode(row, s, mSpecification.isLargestFragmentOnly())) == 0) {
										isMatch = true;
										break;
									}
								}
								}
							}
						else if (mSpecification.isNoStereoSearch()) {
							for (int s=0; !isMatch && s<mDataSource.getStructureCount(row); s++) {
								for (int i=0; i<mQueryHashCode.length; i++) {
									if (mQueryHashCode[i] != 0 && mQueryHashCode[i] == mDataSource.getNoStereoCode(row, s, mSpecification.isLargestFragmentOnly())) {
										isMatch = true;
										break;
										}
									}
								}
							}
						else if (mSpecification.isTautomerSearch()) {
							for (int s=0; !isMatch && s<mDataSource.getStructureCount(row); s++) {
								for (long l : mQueryHashCode) {
									if (l != 0 && l == mDataSource.getTautomerCode(row, s, mSpecification.isLargestFragmentOnly())) {
										isMatch = true;
										break;
									}
								}
								}
							}
						else if (mSpecification.isNoStereoTautomerSearch()) {
							for (int s=0; !isMatch && s<mDataSource.getStructureCount(row); s++) {
								for (long l : mQueryHashCode) {
									if (l != 0 && l == mDataSource.getNoStereoTautomerCode(row, s, mSpecification.isLargestFragmentOnly())) {
										isMatch = true;
										break;
									}
								}
								}
							}
						else if (mSpecification.isBackboneSearch()) {
							for (int s=0; !isMatch && s<mDataSource.getStructureCount(row); s++) {
								for (long l : mQueryHashCode) {
									if (l != 0 && l == mDataSource.getBackboneCode(row, s, mSpecification.isLargestFragmentOnly())) {
										isMatch = true;
										break;
									}
								}
								}
							}
						}

					if (isMatch) {
						mResultQueue.add(row);
						mMatchCount.incrementAndGet();
						}
					}

				row = mSMPIndex.decrementAndGet();
				}

/*    				if (mSMPWorkingThreads.decrementAndGet() == 0) {
				if (!progressController.threadMustDie())
					// do cleanup some stuff

				// do cleanup some stuff
				}*/
			}
		}
	}
