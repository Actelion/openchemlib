package com.actelion.research.chem.mcs;

import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.properties.complexity.IBitArray;
import com.actelion.research.chem.shredder.Fragment;
import com.actelion.research.util.Pipeline;

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
public class ExhaustiveFragmentGeneratorParallizer {

	private static final int CAPACITY = 10000;

	public static final int TOTAL_CAPACITY = (int)(10*Math.pow(10,6));

	/**
	 * Needed if wildcards are added. Each wildcard adds one bond to the fragment size.
	 */
	private static final int OFFSET_FRAG_SIZE = 24;
	
	private Pipeline<StereoMolecule> queuePipe;
	
	private List<ConcurrentHashMap<String, Fragment>> liHashMapIdCode_Fragment;
	
	private AtomicInteger ccProcessedRecords;

	private AtomicInteger ccMoleculeToLarge;

	
	private int nCores;
	
	private int bits;
	
	/**
	 * 
	 * @param bits maximum number of bonds in the fragment. 
	 * Actually it is the maximum number of bits to store the bond information.
	 */
	public ExhaustiveFragmentGeneratorParallizer(int bits) {
		
		this.bits = bits;
		
		liHashMapIdCode_Fragment = new ArrayList<ConcurrentHashMap<String,Fragment>>();

		ccProcessedRecords = new AtomicInteger();

		ccMoleculeToLarge = new AtomicInteger();

		nCores = Runtime.getRuntime().availableProcessors();
		
		if(nCores>1){
			nCores--;
		}
		// nCores = 1;
		
	}
	
	
	public void process(Pipeline<StereoMolecule> queuePipe, int maxSizeFrag, boolean cleaveRingBonds, boolean addWildcards) {
		
		
		for (ConcurrentHashMap<String, Fragment> hm : liHashMapIdCode_Fragment) {
			hm.clear();
		}
		
		int maxBondsPlusBondsWildcards = maxSizeFrag+OFFSET_FRAG_SIZE+1;
		
		for (int i = liHashMapIdCode_Fragment.size(); i < maxBondsPlusBondsWildcards; i++) {
			
			if(i < 4) {
				liHashMapIdCode_Fragment.add(new ConcurrentHashMap<String, Fragment>());	
			} else {
				liHashMapIdCode_Fragment.add(new ConcurrentHashMap<String, Fragment>(CAPACITY));
			}
			
		}
		
		
		this.queuePipe = queuePipe;
		
		ccProcessedRecords.set(0);


		int cores = Runtime.getRuntime().availableProcessors();

		int poolSize = cores-1;

		if(poolSize==0){
			poolSize=1;
		}

		ExecutorService executorService = Executors.newFixedThreadPool(poolSize);


		List<RunEFG> liRunEFG = new ArrayList<RunEFG>();


		for (int i = 0; i < poolSize; i++) {

			RunEFG extractor = new RunEFG(bits, maxSizeFrag, cleaveRingBonds, addWildcards);

			executorService.execute(extractor);

			liRunEFG.add(extractor);

		}

		executorService.shutdown();

		while (!executorService.isTerminated()){

			try {
				Thread.sleep(10000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}



			System.out.println(" Processed records " + ccProcessedRecords.get() + ".");

			System.out.println(" Molecules still in pipeline " + queuePipe.sizePipe() + ".");

			int ccFrags = 0;
			for (ConcurrentHashMap<String, Fragment> hm : liHashMapIdCode_Fragment) {
				ccFrags += hm.size();
			}

			System.out.println(" Unique fragments " + ccFrags + ".");

		}

		System.out.println(" Processed records " + ccProcessedRecords.get() + ".");

		System.out.println(" Molecules still in pipeline " + queuePipe.sizePipe() + ".");

		int ccFrags = 0;
		for (ConcurrentHashMap<String, Fragment> hm : liHashMapIdCode_Fragment) {
			ccFrags += hm.size();
		}

		System.out.println(" Unique fragments " + ccFrags + ".");

		System.out.println("Finished all " + new Date().toString());


	}

	
	
	/**
	 * 
	 * @return list of unique fragments.
	 */
	public List<Fragment> getFragmentList(){
		
		List<Fragment> li = new ArrayList<Fragment>();
				
		for (ConcurrentHashMap<String, Fragment> hm : liHashMapIdCode_Fragment) {
			li.addAll(hm.values());
		}

		return li;

	}
		
	class RunEFG implements Runnable {

		 private int maxSizeFrag;
		 
		 private boolean addWildcards;
		 
		 private boolean cleaveRingBonds;
		
		 private HashSet<String> hsIdCode;
		 
		 private ExhaustiveFragmentGeneratorBonds efg;
		
		/**
		 * @param maxSizeFrag
		 * @param cleaveRingBonds
		 * @param addWildcards
		 */
		public RunEFG(int bits, int maxSizeFrag, boolean cleaveRingBonds, boolean addWildcards) {
			
			super();
			
			this.maxSizeFrag = maxSizeFrag;
			
			this.addWildcards = addWildcards;
			
			this.cleaveRingBonds = cleaveRingBonds;
			
			hsIdCode = new HashSet<String>();
			
			efg = new ExhaustiveFragmentGeneratorBonds(bits, TOTAL_CAPACITY);
		}


		public void run() {
			
			while(!queuePipe.wereAllDataFetched()){
				
				StereoMolecule mol = queuePipe.pollData();

				if(mol == null){
					try {Thread.sleep(500);} catch (InterruptedException e) {e.printStackTrace();}
					continue;
				}

				if(mol.getBonds() <= bits) {

					try {
						process(mol);
					} catch (Exception e) {
						e.printStackTrace();
					}

				} else {

					ccMoleculeToLarge.incrementAndGet();
				}

				ccProcessedRecords.incrementAndGet();
			}
		}
		
		
		public void process(StereoMolecule mol) {
								
			hsIdCode.clear();
			
			BondVector2IdCode bondVector2IdCode = new BondVector2IdCode(mol);

			efg.set(mol, maxSizeFrag);
			
			efg.generateFragmentsAllBonds();
			
			int bonds = mol.getBonds();
			
			int n = Math.min(bonds, maxSizeFrag+1);
		
			for (int i = 0; i < n; i++) {
				List<IBitArray> liFragDefByBnds = efg.getFragments(i);
				
				for (IBitArray fragDefByBnds : liFragDefByBnds) {
					
					if(!cleaveRingBonds && bondVector2IdCode.containsFragmentOpenRing(fragDefByBnds)){
						continue;
					}


					int bnds = fragDefByBnds.getBitsSet();

					Fragment fragmentNew = bondVector2IdCode.getFragment(fragDefByBnds, addWildcards);

					if(fragmentNew.getMol().getBonds() >= liHashMapIdCode_Fragment.size()){

						Canonizer can = new Canonizer(fragmentNew.getMol());

						System.out.println("ExhaustiveFragmentGeneratorParallizer RunEFG bonds fragment " + bnds);

						System.out.println("ExhaustiveFragmentGeneratorParallizer RunEFG bonds fragment with wildcards " + fragmentNew.getMol().getBonds());

						System.out.println(can.getIDCode());

						continue;

					}

					
					ConcurrentHashMap<String, Fragment> hmIdCode_Fragment = liHashMapIdCode_Fragment.get(fragmentNew.getMol().getBonds());
					
					if(!hmIdCode_Fragment.containsKey(fragmentNew.getIdcode())){
						fragmentNew.setFrequencySumAll(0);
												
						hmIdCode_Fragment.put(fragmentNew.getIdcode(), fragmentNew);
					}
					
					Fragment fragment = hmIdCode_Fragment.get(fragmentNew.getIdcode());
					
					fragment.incrementFrequencySumAll();
					
					if(hsIdCode.add(fragment.getIdcode())){
						
						fragment.incrementFrequencyOnePerMol();
						
					}
				}
			}
		}
		
	}


}
