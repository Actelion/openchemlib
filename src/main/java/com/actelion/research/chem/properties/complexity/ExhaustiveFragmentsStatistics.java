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
*/

package com.actelion.research.chem.properties.complexity;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.mcs.ExhaustiveFragmentGeneratorBonds;
import com.actelion.research.chem.mcs.RunBondVector2IdCode;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.Pipeline;
import com.actelion.research.util.datamodel.ByteVec;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class ExhaustiveFragmentsStatistics {
	
	private static boolean ELUSIVE = false; 

	public static final int TOTAL_CAPACITY = (int)(10*Math.pow(10,6));

	public static final int MINLEN_FRAG = 1;
		
	private static final long SLEEP = 10;
	
	private static final int CAPACITY = 200 * 1000 * 1000;

	private static final int LIMIT_BONDS_WHERE_IDCODE_IS_STORED = 20;

	private int minNumBondsFragment;
	
	private ExhaustiveFragmentGeneratorBonds efg;
	
	private HashSet<ByteVec> hsIdCode;
	
	private List<List<ByteVec>> liliIdCode;
	
	private Pipeline<IBitArray> pipeInputFragIndexListsFromEFG;

	private Pipeline<FragmentDefinedByBondsIdCode> pipeOutputFragmentDefinedByBondsIdCode;

	private List<RunBondVector2IdCode> liRunBondVector2IdCode;

	private AtomicInteger processedFragments;

	private boolean collectFragmentIdCodes;
		
	
	/**
	 * 
	 * @param bits maximum number of bonds that can be stored in the bit arrays
	 * @param threads
	 * @param totalCapacity
	 */
	
	public ExhaustiveFragmentsStatistics(int bits, int threads, int totalCapacity) {
		init(bits, threads, totalCapacity);
		
	}
	
	public ExhaustiveFragmentsStatistics(int bits, int totalCapacity) {
		
		int threads = Runtime.getRuntime().availableProcessors()-1;
		
		threads = Math.max(1, threads);
		
		init(bits, threads, totalCapacity);
		
	}
	
	public ExhaustiveFragmentsStatistics(int bits) {

		int threads = Runtime.getRuntime().availableProcessors()-1;
		
		threads = Math.max(1, threads);

		init(bits, threads, TOTAL_CAPACITY);
	}
	
	private void init(int bits, int threadsBondVector2IdCode, int totalCapacity){

		efg  = new ExhaustiveFragmentGeneratorBonds(bits, totalCapacity);

		// Formatter

		System.out.println("ExhaustiveFragmentsStatistics init(...) totalCapacity " + Formatter.group(totalCapacity));

		minNumBondsFragment = MINLEN_FRAG;
				
		hsIdCode = new HashSet<ByteVec>(CAPACITY);
		
		processedFragments = new AtomicInteger();
		
		pipeInputFragIndexListsFromEFG = new Pipeline<IBitArray>();

		pipeOutputFragmentDefinedByBondsIdCode = new Pipeline<>();

		liRunBondVector2IdCode = new ArrayList<>();
		
		for (int i = 0; i < threadsBondVector2IdCode; i++) {
			RunBondVector2IdCode runBondVector2IdCode = new RunBondVector2IdCode(i, pipeInputFragIndexListsFromEFG, pipeOutputFragmentDefinedByBondsIdCode);
			liRunBondVector2IdCode.add(runBondVector2IdCode);
        	new Thread(runBondVector2IdCode).start();
		}
	}
	
	/**
	 * Creates a list with ModelExhaustiveStatistics. Each object contains the 
	 * number of fragments and the number of unique fragments for a certain number of atoms (length).
	 * @param mol
	 * @return
	 */
	public ResultFragmentsStatistic create(StereoMolecule mol, int maxNumBondsFragmentDesired){
		
		if(!pipeInputFragIndexListsFromEFG.isEmpty()){
			throw new RuntimeException("Error in algorithm!");
		}
		
		int maxNumBondsFragment = Math.min(efg.getMaximumCapacityBondsInFragment(), maxNumBondsFragmentDesired);

		List<ModelExhaustiveStatistics> liModelExhaustiveStatistics = new ArrayList<>();
		
		if(collectFragmentIdCodes) {
			for (int i = 0; i < maxNumBondsFragment+1; i++) {
				liliIdCode.get(i).clear();
			}
		}
		
		efg.set(mol, maxNumBondsFragment);
		
		efg.generateFragmentsAllBonds();
		
		for (RunBondVector2IdCode runBondVector2IdCode : liRunBondVector2IdCode) {
			runBondVector2IdCode.init(mol);
		}
		
		double bondsMol = mol.getBonds();
				
		for (int nBondsFrag = minNumBondsFragment; nBondsFrag < maxNumBondsFragment+1; nBondsFrag++) {

			if(isELUSIVE()){
				System.out.println("ExhaustiveFragmentsStatistics create(...) process for bond count " + nBondsFrag + ".");
			}

			hsIdCode.clear();

			pipeInputFragIndexListsFromEFG.reset();

			pipeOutputFragmentDefinedByBondsIdCode.reset();


			List<IBitArray> liIntegerListResultsEFG = efg.getFragments(nBondsFrag);

			if(isELUSIVE()){
				System.out.println("ExhaustiveFragmentsStatistics create(...) processing for bond count " + nBondsFrag + ". Retrieved list with fragments " +  Formatter.group(liIntegerListResultsEFG.size()) + ".");
			}

			processedFragments.set(0);

			if(!pipeInputFragIndexListsFromEFG.isEmpty()){
				throw new RuntimeException("Error in algorithm!");
			}

			// Feeds bond vector to idcode
			pipeInputFragIndexListsFromEFG.addData(liIntegerListResultsEFG);

			// Do not set all data in for pipeInputFragIndexListsFromEFG here. This is set in finalize().

			Runnable runAdd2HashSet = new Runnable() {
				@Override
				public void run() {

					try {
						while (!pipeOutputFragmentDefinedByBondsIdCode.wereAllDataFetched()) {
							FragmentDefinedByBondsIdCode fragmentIndexIdCode = pipeOutputFragmentDefinedByBondsIdCode.pollData();
							if (fragmentIndexIdCode == null) {
								try {
									Thread.sleep(SLEEP);
								} catch (InterruptedException e) {
									e.printStackTrace();
								}
								continue;
							}
							hsIdCode.add(new ByteVec(fragmentIndexIdCode.getIdCode()));
						}
					} catch (Throwable e){
						e.printStackTrace();
					} finally {
						if(isELUSIVE()){
							System.out.println("ExhaustiveFragmentsStatistics Runnable runAdd2HashSet finally reached. hsIdCode " + hsIdCode.size());
						}
					}
				}
			};

			ExecutorService executorServiceAdd2HashSet = Executors.newSingleThreadExecutor();

			executorServiceAdd2HashSet.submit(runAdd2HashSet);

			executorServiceAdd2HashSet.shutdown();

			while(!executorServiceAdd2HashSet.isTerminated()){
				try {Thread.sleep(SLEEP);} catch (InterruptedException e) {e.printStackTrace();}
				if(pipeOutputFragmentDefinedByBondsIdCode.getAdded() == liIntegerListResultsEFG.size()){
					pipeOutputFragmentDefinedByBondsIdCode.setAllDataIn();
				}
			}

			if(isELUSIVE()){
				System.out.println("ExhaustiveFragmentsStatistics create(...) add idcode to HashSet finished for bond count " + nBondsFrag + ".");
			}
			
			double ratioCoveredBonds =  (double)nBondsFrag / bondsMol;
						
			ModelExhaustiveStatistics modelLengthFragsUnique = 
				new ModelExhaustiveStatistics(
						nBondsFrag, 
						liIntegerListResultsEFG.size(), 
						hsIdCode.size(), 
						ratioCoveredBonds);
			
			
			
			liModelExhaustiveStatistics.add(modelLengthFragsUnique);
			
			if(collectFragmentIdCodes) {
				liliIdCode.get(nBondsFrag).addAll(hsIdCode);
			}

			if(isELUSIVE()){
				System.out.println("ExhaustiveFragmentsStatistics create(...) finished for bond count " + nBondsFrag + ".");
			}

		}
		
		ResultFragmentsStatistic fragmentsStatistic = new ResultFragmentsStatistic(mol, liModelExhaustiveStatistics);
		
		return fragmentsStatistic;
	}
	
	
		
	public boolean isCapacityLimitBreakes() {
		return efg.isCapacityLimitBreakes();
	}
	

	/**
	 * 
	 * @return List of lists with unique IdCodes.
	 * The index in the list equals the number of heavy atoms in the fragments.
	 */
	public List<List<ByteVec>> getLiLiIdCode() {
		return liliIdCode;
	}
	
	private boolean areAllReachedEndOfRun(){
		
		boolean endOfRunReached = true;
		
		for(RunBondVector2IdCode runBondVector2IdCode : liRunBondVector2IdCode){
			if(!runBondVector2IdCode.isEndOfRunReached()){
				endOfRunReached = false;
				break;
			}
		}
		
		return endOfRunReached;
	}
	
	public void roundUp() throws Throwable {
		
		pipeInputFragIndexListsFromEFG.setAllDataIn(true);
		
		while(!areAllReachedEndOfRun()){
			try {Thread.sleep(SLEEP);} catch (InterruptedException e) {e.printStackTrace();}
			System.out.println("ExhaustiveFragmentsStatistics roundUp() waiting for end of run.");
		}

	}

	public int getMaximumNumberBondsInMolecule(){
		return efg.getMaximumNumberBondsInMolecule();
	}


	public void setCollectFragmentIdCodes(boolean collectFragmentIdCodes) {
		
		this.collectFragmentIdCodes = collectFragmentIdCodes;
		
		if(collectFragmentIdCodes) {

			liliIdCode = new ArrayList<List<ByteVec>>();
			
			int [] arrCapacity = new int [LIMIT_BONDS_WHERE_IDCODE_IS_STORED];

			arrCapacity[0] = ContainerFragBondsSolutions.START_CAPACITY;
			
			for (int i = 1; i < arrCapacity.length; i++) {
				int capacity = (int)(arrCapacity[i-1] * ContainerFragBondsSolutions.FACTOR_CAPACITY);
				arrCapacity[i] = capacity;
			}

			for (int i = 0; i < MINLEN_FRAG; i++) {
				liliIdCode.add(new ArrayList<ByteVec>());
			}
					
			for (int i = 0; i < arrCapacity.length; i++) {
				liliIdCode.add(new ArrayList<ByteVec>(arrCapacity[i]));
			}
		}

	}

	
	public static Comparator<MultipleNonOverlapSolution> getComparatorCoverage() {
		return new Comparator<MultipleNonOverlapSolution>() {

			@Override
			public int compare(MultipleNonOverlapSolution o1, MultipleNonOverlapSolution o2) {
			
				if(o1.getContainer().getBitsSet() > o2.getContainer().getBitsSet()){
					return 1;
				} else if(o1.getContainer().getBitsSet() < o2.getContainer().getBitsSet()){
					return -1;
				}
				
				return 0;
			}
		};
	}

	/**
	 * @return the eLUSIVE
	 */
	public static boolean isELUSIVE() {
		return ELUSIVE;
	}

	/**
	 * @param elusive the ELUSIVE to set
	 */
	public static void setELUSIVE(boolean elusive) {
		ELUSIVE = elusive;
		ExhaustiveFragmentGeneratorBonds.setELUSIVE(elusive);
	}
	
	

	

}

