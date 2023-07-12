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

package com.actelion.research.chem.reaction;

import com.actelion.research.calc.ProgressController;
import com.actelion.research.chem.*;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerReactionFP;
import com.actelion.research.util.IntArrayComparator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

public class ReactionSearch {
	private static final boolean MULTITHREADED_SEARCH = true;

	private volatile ReactionSearchSpecification mSpecification;
	private volatile ReactionSearchDataSource mDataSource;
	private volatile StructureSearchController mSearchController;
	private volatile ProgressController mProgressController;
	private volatile Reaction[] mQueryReaction;
	private volatile StereoMolecule[] mQueryReactant,mQueryProduct,mQueryRetron;
	private volatile long[] mQueryHash;
	private volatile DescriptorHandlerLongFFP512 mDescriptorHandlerFFP512;
	private volatile DescriptorHandlerReactionFP mDescriptorHandlerRxnFP;
	private volatile long[][] mQueryReactionDescriptor,mQueryReactantDescriptor,mQueryProductDescriptor,mQueryRetronDescriptor;
	private volatile int mMaxSSSMatches,mMaxNonSSSMatches,mStatus;
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
	 */
	public ReactionSearch(ReactionSearchSpecification specification,
	                      ReactionSearchDataSource dataSource,
						  StructureSearchController searchController,
						  ProgressController progressController) {
		mSpecification = specification;
		mDataSource = dataSource;
		mSearchController = searchController;
		mProgressController = progressController;
		mStatus = StructureSearch.SEARCH_PENDING;

		if (mSpecification != null) {
			// define needed descriptor handlers
			if (mSpecification.isSimilaritySearch()) {
				mDescriptorHandlerRxnFP = DescriptorHandlerReactionFP.getDefaultInstance();
				}
			else if (mSpecification.isSubreactionSearch()
				  || mSpecification.isRetronSearch()) {
				mDescriptorHandlerFFP512 = DescriptorHandlerLongFFP512.getDefaultInstance();
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
		return StructureSearch.COMPLETION_TEXT[mStatus];
		}

	public int[] start() {
		if (!mDataSource.isSupportedSearchType(mSpecification)) {
			mStatus = StructureSearch.SEARCH_TYPE_NOT_SUPPORTED;
			return null;
			}

		mMatchCount = new AtomicInteger(0);

		if (!mSpecification.isNoReactionSearch()) {
			final int queryReactionCount = mSpecification.getReactionCount();
			if (queryReactionCount == 0) {
				mStatus = StructureSearch.QUERY_MISSING;
				return null;
				}

			if (mSpecification.isSubreactionSearch()
			 || mSpecification.isSimilaritySearch()) {
				mQueryReaction = new Reaction[queryReactionCount];
				mQueryReactant = new StereoMolecule[queryReactionCount];
				mQueryProduct = new StereoMolecule[queryReactionCount];
				for (int i=0; i<queryReactionCount; i++) {
					mQueryReactant[i] = mergeMolecules(ReactionEncoder.decodeMolecules(mSpecification.getEncodedQuery(i), true, true, true, false));
					mQueryReactant[i].ensureHelperArrays(Molecule.cHelperParities);
					mQueryProduct[i] = mergeMolecules(ReactionEncoder.decodeMolecules(mSpecification.getEncodedQuery(i), true, true, false, true));
					mQueryProduct[i].ensureHelperArrays(Molecule.cHelperParities);
					mQueryReaction[i] = new Reaction();
					mQueryReaction[i].addReactant(mQueryReactant[i]);
					mQueryReaction[i].addProduct(mQueryProduct[i]);
					}
				if (mSpecification.isSubreactionSearch())
					ensureMoleculeDescriptors();
				else
					ensureReactionDescriptors();
				}
			else if (mSpecification.isRetronSearch()) {
				mQueryRetron = new StereoMolecule[queryReactionCount];
				for (int i=0; i<queryReactionCount; i++) {
					mQueryRetron[i] = new IDCodeParser(false).getCompactMolecule(mSpecification.getEncodedQuery(i));
					mQueryRetron[i].ensureHelperArrays(Molecule.cHelperParities);
					}
				ensureRetronDescriptors();
				}
			else if (mSpecification.isExactSearch()) {
				mQueryHash = new long[queryReactionCount];
				for (int i=0; i<queryReactionCount; i++) {
					int index = mSpecification.getEncodedQuery(i).indexOf(ReactionEncoder.OBJECT_DELIMITER);
					if (index > 0) {
						String[] idcodes = mSpecification.getEncodedQuery(i).substring(0, index).split(ReactionEncoder.MOLECULE_DELIMITER_STRING);
						for (String idcode:idcodes)
							mQueryHash[i] += CanonizerUtil.StrongHasher.hash(idcode);
						}
					}
				}
			else if (mSpecification.isNoStereoSearch()) {
				mQueryHash = new long[queryReactionCount];
				for (int i=0; i<queryReactionCount; i++) {
					int index = mSpecification.getEncodedQuery(i).indexOf(ReactionEncoder.OBJECT_DELIMITER);
					if (index > 0) {
						String[] idcodes = mSpecification.getEncodedQuery(i).substring(0, index).split(ReactionEncoder.MOLECULE_DELIMITER_STRING);
						for (String idcode:idcodes)
							mQueryHash[i] += CanonizerUtil.getNoStereoHash(
								new IDCodeParser(false).getCompactMolecule(idcode), false);
						}
					}
				}
			}

    	mSMPIndex = new AtomicInteger(mDataSource.getRowCount());

    	mResultQueue = new ConcurrentLinkedQueue<>();

		if (mProgressController != null && mSpecification.getReactionCount() > 1023)
			mProgressController.startProgress("Searching reactions", 0, mSpecification.getReactionCount());

		mStopTime = (mMaxMillis == 0) ? Long.MAX_VALUE : System.currentTimeMillis() + mMaxMillis;
		mStatus = StructureSearch.SEARCH_RUNNING;

		if (MULTITHREADED_SEARCH) {
			int threadCount = Runtime.getRuntime().availableProcessors();
			SearchThread[] t = new SearchThread[threadCount];
			for (int i = 0; i<threadCount; i++) {
				t[i] = new SearchThread("Reaction Search " + (i + 1));
				t[i].setPriority(Thread.MIN_PRIORITY);
				t[i].start();
				}

			// the controller thread must wait until all others are finished
			// before the next task can begin or the dialog is closed
			for (int i = 0; i<threadCount; i++)
				try {
					t[i].join();
					}
				catch (InterruptedException e) {}
			}
		else {
			new SearchThread("Reaction Search").run();
			}

		if (mStatus == StructureSearch.SEARCH_RUNNING)
			mStatus = StructureSearch.SUCCESSFUL_COMPLETION;

		int[] result = new int[mResultQueue.size()];
    	int i=0;
    	for (Integer integer:mResultQueue)
    		result[i++] = integer.intValue();

    	return result;
		}

	private void ensureMoleculeDescriptors() {
		final int queryReactionCount = mSpecification.getReactionCount();
		mQueryReactantDescriptor = new long[queryReactionCount][];
		mQueryProductDescriptor = new long[queryReactionCount][];
		boolean missingDescriptorFound = false;
		for (int i=0; i<queryReactionCount; i++) {
			mQueryReactantDescriptor[i] = mSpecification.getReactantDescriptor(i);
			if (mQueryReactantDescriptor[i] == null)
				missingDescriptorFound = true;
			mQueryProductDescriptor[i] = mSpecification.getProductDescriptor(i);
			if (mQueryProductDescriptor[i] == null)
				missingDescriptorFound = true;
		}

		if (!missingDescriptorFound)
			return;

		mSMPIndex = new AtomicInteger(mQueryReactantDescriptor.length);
		int threadCount = Math.min(queryReactionCount, Runtime.getRuntime().availableProcessors());
		Thread[] t = new Thread[threadCount];
		for (int i=0; i<threadCount; i++) {
			t[i] = new Thread("Query Molecule Descriptor Calculation "+(i+1)) {
				public void run() {
					while (true) {
						int index = mSMPIndex.decrementAndGet();
						if (index < 0)
							break;

						if (mQueryReactantDescriptor[index] == null)
							mQueryReactantDescriptor[index] = mDescriptorHandlerFFP512.createDescriptor(mQueryReactant[index]);
						if (mQueryProductDescriptor[index] == null)
							mQueryProductDescriptor[index] = mDescriptorHandlerFFP512.createDescriptor(mQueryProduct[index]);
						}
					}
				};
			t[i].setPriority(Thread.MIN_PRIORITY);
			t[i].start();
			}

		for (int i=0; i<threadCount; i++)
			try { t[i].join(); } catch (InterruptedException e) {}
		}

	private void ensureRetronDescriptors() {
		final int queryReactionCount = mSpecification.getReactionCount();
		mQueryRetronDescriptor = new long[queryReactionCount][];
		boolean missingDescriptorFound = false;
		for (int i=0; i<queryReactionCount; i++) {
			mQueryRetronDescriptor[i] = mSpecification.getRetronDescriptor(i);
			if (mQueryRetronDescriptor[i] == null)
				missingDescriptorFound = true;
			}

		if (!missingDescriptorFound)
			return;

		mSMPIndex = new AtomicInteger(queryReactionCount);
		int threadCount = Math.min(queryReactionCount, Runtime.getRuntime().availableProcessors());
		Thread[] t = new Thread[threadCount];
		for (int i=0; i<threadCount; i++) {
			t[i] = new Thread("Query Retron Descriptor Calculation "+(i+1)) {
				public void run() {
					while (true) {
						int index = mSMPIndex.decrementAndGet();
						if (index < 0)
							break;

						if (mQueryRetronDescriptor[index] == null)
							mQueryRetronDescriptor[index] = mDescriptorHandlerFFP512.createDescriptor(mQueryRetron[index]);
						}
					}
				};
			t[i].setPriority(Thread.MIN_PRIORITY);
			t[i].start();
			}

		for (int i=0; i<threadCount; i++)
			try { t[i].join(); } catch (InterruptedException e) {}
		}

	private void ensureReactionDescriptors() {
		final int queryReactionCount = mSpecification.getReactionCount();
		mQueryReactionDescriptor = new long[queryReactionCount][];
		boolean missingDescriptorFound = false;
		for (int i=0; i<queryReactionCount; i++) {
			mQueryReactionDescriptor[i] = mSpecification.getReactionDescriptor(i);
			if (mQueryReactionDescriptor[i] == null)
				missingDescriptorFound = true;
		}

		if (!missingDescriptorFound)
			return;

		if (MULTITHREADED_SEARCH && queryReactionCount > 1) {
	        mSMPIndex = new AtomicInteger(mQueryReactionDescriptor.length);
			int threadCount = Math.min(queryReactionCount, Runtime.getRuntime().availableProcessors());
	        Thread[] t = new Thread[threadCount];
	        for (int i=0; i<threadCount; i++) {
	            t[i] = new Thread("Query Reaction Descriptor Calculation "+(i+1)) {
	                public void run() {
	                    while (true) {
	                        int index = mSMPIndex.decrementAndGet();
	                        if (index < 0)
	                            break;

	                        if (mQueryReactionDescriptor[index] == null)
								mQueryReactionDescriptor[index] = mDescriptorHandlerRxnFP.createDescriptor(mQueryReaction[index]);
	                        }
	                    }
	                };
	            t[i].setPriority(Thread.MIN_PRIORITY);
	            t[i].start();
	            }

	        for (int i=0; i<threadCount; i++)
	            try { t[i].join(); } catch (InterruptedException e) {}
			}
		else {
			for (int i=0; i<mQueryReactionDescriptor.length; i++)
				mQueryReactionDescriptor[i] = mDescriptorHandlerRxnFP.createDescriptor(mQueryReaction[i]);
			}
		}

	private StereoMolecule mergeMolecules(StereoMolecule[] mol) {
		if (mol == null)
			return null;

		if (mol.length > 1) {	// gather all products within product[0]
			mol[0] = new StereoMolecule(mol[0]);
			for (int i=1; i<mol.length; i++)
				mol[0].addMolecule(mol[i]);
			}

		return mol[0];
		}

	private class SearchThread extends Thread {
		private SRSearcher mSRSearcher;
		private SSSearcherWithIndex mReactantSearcher,mProductSearcher;

		public SearchThread(String name) {
			super(name);
			if (mSpecification.isSubreactionSearch()) {
				mSRSearcher = new SRSearcher();
				}
			else if (mSpecification.isRetronSearch()) {
				mReactantSearcher = new SSSearcherWithIndex();
				mProductSearcher = new SSSearcherWithIndex();
				}
			}

		public void run() {
			int row = mSMPIndex.decrementAndGet();
			while (row >= 0) {
				if ((mProgressController != null && mProgressController.threadMustDie())) {
					mStatus = StructureSearch.SEARCH_STOPPED;
					break;
					}

				if (System.currentTimeMillis() > mStopTime) {
					mStatus = StructureSearch.TIME_LIMIT_EXCEEDED;
					break;
					}

				if (mProgressController != null && row%1024==1023)
					mProgressController.updateProgress(mSpecification.getReactionCount()-row);

				if (mSearchController == null || mSearchController.rowQualifies(row)) {
					boolean isMatch = false;

					if (mSpecification.isSubreactionSearch()
					 || mSpecification.isRetronSearch()) {
						if (mMaxSSSMatches != 0 && mMatchCount.get() > mMaxSSSMatches) {
							mStatus = StructureSearch.COUNT_LIMIT_EXCEEDED;
							break;
							}

						if (mSpecification.isSubreactionSearch()) {
							long[] reactantFFP = mDataSource.getMergedReactantDescriptor(row);
							long[] productFFP = mDataSource.getMergedProductDescriptor(row);
							mSRSearcher.setReaction(mDataSource.getReactionCode(row), mDataSource.getMapping(row), mDataSource.getCoordinates(row), reactantFFP, productFFP);

							for (int i = 0; i<mQueryReaction.length; i++) {
								mSRSearcher.setQuery(mQueryReaction[i], mQueryReactantDescriptor[i], mQueryProductDescriptor[i]);
								if (mSRSearcher.isQueryInReaction()) {
									isMatch = true;
									break;
									}
								}
							}
						else {  // retron search
							for (int i = 0; i<mQueryRetron.length; i++) {
								if (!mDescriptorHandlerFFP512.calculationFailed(mQueryRetronDescriptor[i])) {
									long[] productFFP = mDataSource.getMergedProductDescriptor(row);
									mProductSearcher.setFragment(mQueryRetron[i], mQueryRetronDescriptor[i]);
									mProductSearcher.setMolecule((StereoMolecule)null, productFFP);
									if (mProductSearcher.isFragmentIndexInMoleculeIndex()) {
										StereoMolecule product = mergeMolecules(ReactionEncoder.decodeMolecules(
												mDataSource.getReactionCode(row),
												mDataSource.getCoordinates(row),
												mDataSource.getMapping(row), false, true));
										mProductSearcher.setMolecule(product, productFFP);
										int inProductCount = mProductSearcher.findFragmentInMoleculeWithoutIndex(SSSearcher.cCountModeSeparated);
										if (inProductCount != 0) {
											long[] reactantFFP = mDataSource.getMergedReactantDescriptor(row);
											mReactantSearcher.setFragment(mQueryRetron[i], mQueryRetronDescriptor[i]);
											mReactantSearcher.setMolecule((StereoMolecule)null, reactantFFP);
											int inReactantCount = 0;
											if (mReactantSearcher.isFragmentIndexInMoleculeIndex()) {
												StereoMolecule reactant = mergeMolecules(ReactionEncoder.decodeMolecules(
														mDataSource.getReactionCode(row),
														mDataSource.getCoordinates(row),
														mDataSource.getMapping(row), true, false));
												mReactantSearcher.setMolecule(reactant, reactantFFP);
												inReactantCount = mReactantSearcher.findFragmentInMoleculeWithoutIndex(SSSearcher.cCountModeSeparated);
												if (inReactantCount != 0 && mDataSource.getMapping(row) != null) {
													inProductCount -= countEquivalentMatches(product, mProductSearcher.getGraphMatcher().getMatchList());
													if (inProductCount <= inReactantCount)
														inReactantCount -= countEquivalentMatches(reactant, mReactantSearcher.getGraphMatcher().getMatchList());
													}
												}
											// TODO check, whether we also have to take into account in catalyst occurrences
											if (inProductCount > inReactantCount) {
												isMatch = true;
												break;
												}
											}
										}
									}
								}
							}
						}
					else {
						if (mMaxNonSSSMatches != 0 && mMatchCount.get() > mMaxNonSSSMatches) {
							mStatus = StructureSearch.COUNT_LIMIT_EXCEEDED;
							break;
							}

						if (mSpecification.isNoReactionSearch()) {
							isMatch = true;
							}
						else if (mSpecification.isSimilaritySearch()) {
							for (int i=0; i<mQueryReactionDescriptor.length; i++) {
								if (mDescriptorHandlerRxnFP.getReactionCenterSimilarity(mQueryReactionDescriptor[i], mDataSource.getReactionDescriptor(row))
									>= mSpecification.getReactionCenterSimilarity()
								 && mDescriptorHandlerRxnFP.getPeripherySimilarity(mQueryReactionDescriptor[i], mDataSource.getReactionDescriptor(row))
									>= mSpecification.getPeripherySimilarity()) {
									isMatch = true;
									break;
									}
								}
							}
						else if (mSpecification.isExactSearch()) {
							for (int i=0; i<mQueryHash.length; i++) {
								if (mQueryHash[i] == mDataSource.getExactHash(row)) {
									isMatch = true;
									break;
									}
								}
							}
						else if (mSpecification.isNoStereoSearch()) {
							for (int i=0; i<mQueryHash.length; i++) {
								if (mQueryHash[i] == mDataSource.getNoStereoHash(row)) {
									isMatch = true;
									break;
									}
								}
							}
						}

					if (isMatch) {
						mResultQueue.add(new Integer(row));
						mMatchCount.incrementAndGet();
						}
					}

				row = mSMPIndex.decrementAndGet();
				}
			}
		}

	private int countEquivalentMatches(StereoMolecule mol, ArrayList<int[]> matchList) {
		int equivalentCount = 0;
		TreeSet<int[]> uniqueSet = new TreeSet<>(new IntArrayComparator());
		for (int[] match:matchList) {
			int[] mappingNo = new int[match.length];
			boolean foundZero = false;
			for (int i=0; i<match.length; i++) {
				mappingNo[i] = (match[i] == -1) ? -1 : mol.getAtomMapNo(match[i]);
				if (mappingNo[i] == 0) {
					foundZero = true;
					break;
				}
			}
			if (!foundZero) {
				Arrays.sort(mappingNo);
				if (uniqueSet.contains(mappingNo))
					equivalentCount++;
				else
					uniqueSet.add(mappingNo);
				}
			}
		return equivalentCount;
		}
	}
