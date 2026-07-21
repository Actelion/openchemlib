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

import com.actelion.research.chem.Molecule;
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
	





	private AtomicInteger processedFragments;

	private boolean collectFragmentIdCodes;
	private int threads;
	
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
		if(isELUSIVE())
			System.out.println("ExhaustiveFragmentsStatistics init(...) totalCapacity " + Formatter.group(totalCapacity));
		minNumBondsFragment = MINLEN_FRAG;
		hsIdCode = new HashSet<ByteVec>(CAPACITY);
		processedFragments = new AtomicInteger();
		this.threads = threadsBondVector2IdCode;


	}
	
	/**
	 * Creates a list with ModelExhaustiveStatistics. Each object contains the 
	 * number of fragments and the number of unique fragments for a certain number of atoms (length).
	 * @param molIn
	 * @return
	 */
	public ResultFragmentsStatistic create(StereoMolecule molIn, int maxNumBondsFragmentDesired){

		StereoMolecule mol = new StereoMolecule(molIn);
		mol.ensureHelperArrays(Molecule.cHelperRings);

		int maxNumBondsFragment = Math.min(efg.getMaximumCapacityBondsInFragment(), maxNumBondsFragmentDesired);

		if(collectFragmentIdCodes) {
			for (int i = 0; i < maxNumBondsFragment+1; i++) {
				liliIdCode.get(i).clear();
			}
		}
		efg.set(mol, maxNumBondsFragment);
		efg.generateFragmentsAllBonds();
		ExecutorService exServe = Executors.newFixedThreadPool(threads);
		List<RunBondVector2IdCode> liRunBondVector2IdCode = new ArrayList<>();
		for (int i = 0; i < threads; i++) {
			RunBondVector2IdCode runBondVector2IdCode = new RunBondVector2IdCode(mol);
			liRunBondVector2IdCode.add(runBondVector2IdCode);
		}

		double bondsMol = mol.getBonds();
		processedFragments.set(0);
		List<ModelExhaustiveStatistics> liModelExhaustiveStatistics = new ArrayList<>();
		for (int nBondsFrag = minNumBondsFragment; nBondsFrag < maxNumBondsFragment+1; nBondsFrag++) {
			if(isELUSIVE()){
				System.out.println("ExhaustiveFragmentsStatistics create(...) process for bond count " + nBondsFrag + ".");
			}
			hsIdCode.clear();
			List<IBitArray> liIntegerListResultsEFG = efg.getFragments(nBondsFrag);
			// Feeds bond vector to idcode

			Pipeline<IBitArray> pipeInputFragIndexListsFromEFG = new Pipeline<>(liIntegerListResultsEFG);
			Pipeline<FragmentDefinedByBondsIdCode> pipeOutputFragmentDefinedByBondsIdCode = new Pipeline<>();
			for (int i = 0; i < threads; i++) {
				RunBondVector2IdCode run = liRunBondVector2IdCode.get(i);;
				run.init(pipeInputFragIndexListsFromEFG, pipeOutputFragmentDefinedByBondsIdCode);
				exServe.submit(run);
			}

			if(isELUSIVE()){
				System.out.println("ExhaustiveFragmentsStatistics create(...) processing for bond count " + nBondsFrag + ". Retrieved list with fragments " +  Formatter.group(liIntegerListResultsEFG.size()) + ".");
			}
			boolean allDone = false;
			while (!allDone){
				try {Thread.sleep(SLEEP);} catch (InterruptedException e) {e.printStackTrace();}
				allDone = true;
				for (RunBondVector2IdCode run : liRunBondVector2IdCode) {
					if(!run.isEndOfRunReached()){
						allDone = false;
						break;
					}
				}
			}
			pipeOutputFragmentDefinedByBondsIdCode.setAllDataIn();
			List<FragmentDefinedByBondsIdCode> fragmentDefinedByBondsIdCodeList = pipeOutputFragmentDefinedByBondsIdCode.pollAll();

			for (FragmentDefinedByBondsIdCode fragmentDefinedByBondsIdCode : fragmentDefinedByBondsIdCodeList) {
				hsIdCode.add(new ByteVec(fragmentDefinedByBondsIdCode.getIdCode()));
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
		exServe.shutdown();
		while (!exServe.isTerminated()){
			try {Thread.sleep(SLEEP);} catch (InterruptedException e) {e.printStackTrace();}
		}

		List<ModelExhaustiveStatistics> liModelExhaustiveStatisticsCpy = new ArrayList<>(liModelExhaustiveStatistics);

		ResultFragmentsStatistic fragmentsStatistic = new ResultFragmentsStatistic(mol, liModelExhaustiveStatisticsCpy);
		
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

