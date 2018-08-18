package com.actelion.research.chem;

import com.actelion.research.calc.ProgressController;
import com.actelion.research.chem.descriptor.*;
import com.actelion.research.util.ByteArrayComparator;

import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

public class StructureSearch {
	private volatile StructureSearchSpecification mSpecification;
	private volatile StructureSearchDataSource mDataSource;
	private volatile StructureSearchController mSearchController;
	private volatile ProgressController mProgressController;
	private volatile StereoMolecule[] mQueryFragment;
	private volatile ByteArrayComparator mIDCodeComparator;
	private volatile DescriptorHandler mDescriptorHandler;
	private volatile Object[] mQueryDescriptor;
	private volatile long[] mQueryHashCode;
	private volatile byte[][] mQueryIDCode;
	private volatile int mDescriptorColumn;
	private volatile int mMaxSSSMatches,mMaxNonSSSMatches;
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
	 * If the search shall be aborted once it exceed a given number of matches,
	 * then define the maximum number of matches with this method before starting the search.
	 * Calling start with then return the first maximum count valid matches.
	 * @param maxSSSMatches maximum number of allowed sub-structure search matches (0: no limit)
	 * @param maxNonSSSMatches maximum number of allowed matches for other search types (0: no limit)
	 */
	public void setMatchLimit(int maxSSSMatches, int maxNonSSSMatches) {
		mMaxSSSMatches = maxSSSMatches;
		mMaxNonSSSMatches = maxNonSSSMatches;
		}

	public int[] start() {
		if (!mDataSource.isSupportedSearchType(mSpecification))
			return null;

		mMatchCount = new AtomicInteger(0);

		if (!mSpecification.isNoStructureSearch()) {
			final int queryStructureCount = mSpecification.getStructureCount();
			if (queryStructureCount == 0)
				return null;

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
						mQueryIDCode[i] = CanonizerUtil.getIDCode(query, CanonizerUtil.IDCODE_TYPE.NORMAL, largestFragmentOnly).getBytes();
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

    	mResultQueue = new ConcurrentLinkedQueue<Integer>();

		if (mProgressController != null && mSpecification.getStructureCount() > 1023)
			mProgressController.startProgress("Searching structures", 0, mSpecification.getStructureCount());

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

    	int[] result = new int[mResultQueue.size()];
    	int i=0;
    	for (Integer integer:mResultQueue)
    		result[i++] = integer.intValue();

    	return result;
		}

	private void calculateQueryDescriptorsAndWait() {
    	mSMPIndex = new AtomicInteger(mQueryDescriptor.length);
		int threadCount = Runtime.getRuntime().availableProcessors();
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
			while (row >= 0 && (mProgressController == null || !mProgressController.threadMustDie())) {
				if (mProgressController != null && row%1024==1023)
					mProgressController.updateProgress(mSpecification.getStructureCount()-row);

				if (mSearchController == null || mSearchController.rowQualifies(row)) {
					boolean isMatch = false;

					if (mSpecification.isSubstructureSearch()) {
						if (mMaxSSSMatches != 0 && mMatchCount.get() > mMaxSSSMatches)
							break;

						mSSSearcher.setMolecule(mDataSource.getIDCode(row, false), (long[])mDataSource.getDescriptor(mDescriptorColumn, row, false));
						for (int i=0; i<mQueryFragment.length; i++) {
							mSSSearcher.setFragment(mQueryFragment[i], (long[])mQueryDescriptor[i]);
							if (mSSSearcher.isFragmentInMolecule()) {
								isMatch = true;
								break;
								}
							}
						}
					else {
						if (mMaxNonSSSMatches != 0 && mMatchCount.get() > mMaxNonSSSMatches)
							break;

						if (mSpecification.isNoStructureSearch()) {
							isMatch = true;
							}
						else if (mSpecification.isSimilaritySearch()) {
							for (int i=0; i<mQueryDescriptor.length; i++) {
								if (mDescriptorHandler.getSimilarity(mQueryDescriptor[i], mDataSource.getDescriptor(mDescriptorColumn, row, mSpecification.isLargestFragmentOnly()))
									 >= mSpecification.getSimilarityThreshold()) {
									isMatch = true;
									break;
									}
								}
							}
						else if (mSpecification.isExactSearch()) {
							for (int i=0; i<mQueryIDCode.length; i++) {
								if (mIDCodeComparator.compare(mQueryIDCode[i], mDataSource.getIDCode(row, mSpecification.isLargestFragmentOnly())) == 0) {
									isMatch = true;
									break;
									}
								}
							}
						else if (mSpecification.isNoStereoSearch()) {
							for (int i=0; i<mQueryHashCode.length; i++) {
								if (mQueryHashCode[i] == mDataSource.getNoStereoCode(row, mSpecification.isLargestFragmentOnly())) {
									isMatch = true;
									break;
									}
								}
							}
						else if (mSpecification.isTautomerSearch()) {
							for (int i=0; i<mQueryHashCode.length; i++) {
								if (mQueryHashCode[i] == mDataSource.getTautomerCode(row, mSpecification.isLargestFragmentOnly())) {
									isMatch = true;
									break;
									}
								}
							}
						else if (mSpecification.isNoStereoTautomerSearch()) {
							for (int i=0; i<mQueryHashCode.length; i++) {
								if (mQueryHashCode[i] == mDataSource.getNoStereoTautomerCode(row, mSpecification.isLargestFragmentOnly())) {
									isMatch = true;
									break;
									}
								}
							}
						else if (mSpecification.isBackboneSearch()) {
							for (int i=0; i<mQueryHashCode.length; i++) {
								if (mQueryHashCode[i] == mDataSource.getBackboneCode(row, mSpecification.isLargestFragmentOnly())) {
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

/*    				if (mSMPWorkingThreads.decrementAndGet() == 0) {
				if (!progressController.threadMustDie())
					// do cleanup some stuff

				// do cleanup some stuff
				}*/
			}
		}
	}
