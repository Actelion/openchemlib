package com.actelion.research.chem.mcs;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.ErrorHashMap;
import com.actelion.research.util.Pipeline;
import com.actelion.research.util.datamodel.ByteVec;


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
public class ExhaustiveMCSGeneratorParallel {
	
	private static final int MIN_ATOMS_MCS = 4;
	
	private static final int MAX_THREADS = 20;
	
	private static final long SLEEP = 100;
	
	private Pipeline<ByteVec> pipeFragByteVec;
	
	private ConcurrentHashMap<ByteVec, ByteVec> hmMCS_MCS;
	
	private ArrayList<ByteVec> liMolByteVec;
	
	private AtomicLong ccMCSCalculations;
	
	private int maxThreads;
	
	private int ringStatus;
	
	public ExhaustiveMCSGeneratorParallel(int ringStatus) {
		
		this.ringStatus = ringStatus;
		
		maxThreads = MAX_THREADS;
				
		pipeFragByteVec = new Pipeline<ByteVec>();
		
		hmMCS_MCS = new ConcurrentHashMap<ByteVec, ByteVec>();
		
        liMolByteVec = new ArrayList<ByteVec>();
        
        ccMCSCalculations = new AtomicLong();
	}
	
	
	
	public List<StereoMolecule> process(List<StereoMolecule> liMol, File workdir) throws Exception {
				
		System.out.println("MCS process.");

		// LogHandler log = LogHandler.getLog("maximumCommonSubstructure.log");
		
		int nProcessors = Runtime.getRuntime().availableProcessors();
		
		int threads = Math.min(nProcessors, maxThreads);
		
		System.out.println("MCS calculations on " + threads + " processors.");

		
		
		liMolByteVec.clear();
			
		int cc=0; 
			
		for (StereoMolecule mol : liMol) {

			try {
				Canonizer can = new Canonizer(mol);
				
				String idcode = can.getIDCode();
				
				ByteVec bv = new ByteVec(idcode);
				
				liMolByteVec.add(bv);
				
				
				pipeFragByteVec.addData(bv);
				
				
				cc++;
			} catch (Exception e) {
				System.err.println("Molecule list size " + liMol.size() + ".");
				System.err.println("Index molecule " + cc + ".");
				e.printStackTrace();
			}
		}
		
		List<MCSThread> liMCSRunnable = new ArrayList<MCSThread>();
		
		List<Thread> liMCSThread = new ArrayList<Thread>();
		
        for (int i = 0; i < nProcessors; i++) {
        	
        	MCSThread mcsThread = new MCSThread(i, ringStatus);
        	
        	liMCSRunnable.add(mcsThread);
        	
        	Thread th = new Thread(mcsThread);
        	
        	liMCSThread.add(th);
        	
        	th.start();
		}
		
		System.out.println("Calculate MCS from " + liMol.size() + " molecules.");
		
        int ccCycles=0;
        while(!hasMCSCalculationFinished(liMCSThread,liMCSRunnable)) {
        	
        	try {Thread.sleep(SLEEP);} catch (InterruptedException e) {e.printStackTrace();};
        	
        	ccCycles++;
        	if(ccCycles%10==0)
        		System.out.println("Unique mcs " + hmMCS_MCS.size() + ", queue " + pipeFragByteVec.sizePipe() + ".");

        }
		
        List<StereoMolecule> liMolMCS = new ArrayList<StereoMolecule>();
        
        IDCodeParser parser = new IDCodeParser(false);
        
        List<ByteVec> liIdCodeMCS = Collections.list(hmMCS_MCS.keys());
        
        for (ByteVec bvIdCodeMCS : liIdCodeMCS) {
        	liMolMCS.add(parser.getCompactMolecule(bvIdCodeMCS.toStringString()));
		} 
        
        return liMolMCS;
	}
	
	private boolean hasMCSCalculationFinished(List<Thread> liThread, List<MCSThread> liMCSRunnable) {
		
		if(!pipeFragByteVec.isEmpty()){
			return false;
		}
		
		boolean calculating = false;
		for (MCSThread mcsThread : liMCSRunnable) {
			if(mcsThread.isCalculating()){
				calculating=true;
				break;
			}
		}
		
		if(calculating){
			return false;
		}
		
		
		long calculated1 = ccMCSCalculations.get();
		
		try {Thread.sleep(10000);} catch (InterruptedException e) {e.printStackTrace();}
		
		long calculated2 = ccMCSCalculations.get();
		
		
		if(calculated1 != calculated2) {
			return false;
		}
		
		pipeFragByteVec.setAllDataIn(true);
		
		boolean finished = true;
				
		for (Thread th : liThread) {
			
			if(th.isAlive()){
				finished = false;
				break;
			}
		}
				
		return finished;
	}

	
	static class MCSTask {
		
		ByteVec bcIdCode1;
		
		ByteVec bcIdCode2;
		
		public MCSTask(ByteVec bcIdCode1, ByteVec bcIdCode2) {
			
			this.bcIdCode1 = bcIdCode1;
			
			this.bcIdCode2 = bcIdCode2;
		}
		
		public MCSTask(MCSTask mcsTask) {
			
			this.bcIdCode1 = mcsTask.bcIdCode1;
			
			this.bcIdCode2 = mcsTask.bcIdCode2;
		}


	}
	
	static class MCSResult extends MCSTask {
		
		StereoMolecule molMCS;
		
		public MCSResult(MCSTask mcsTask, StereoMolecule molMCS) {
			super(mcsTask);
			this.molMCS = molMCS;
		}


	}
	
	class MCSThread implements Runnable {

		int indexThread;
		
		private MCS mcs;
		
		private IDCodeParser idCodeParser;
		
		private AtomicBoolean calculating;
		
		private ErrorHashMap ehm;
		
		private int nFailedSimilarityCalculations;
		
		private int added;
		
		
		public MCSThread(int indexThread, int ringStatus) {
			
			
			this.indexThread = indexThread;
			
			ehm = new ErrorHashMap();
			
			idCodeParser = new IDCodeParser(false);
			
			mcs = new MCS(ringStatus);
						
			calculating = new AtomicBoolean();
		}
		
		public void run() {
			
			
			while(!pipeFragByteVec.wereAllDataFetched()){
				
				ByteVec bvFrag = pipeFragByteVec.pollData();
				
				if(bvFrag==null){
					calculating.set(false);
					try {Thread.sleep(SLEEP);} catch (InterruptedException e) {e.printStackTrace();}
					continue;
				}
								
				calculating.set(true);
				for (int i = 0; i < liMolByteVec.size(); i++) {
					MCSTask mcsTask = new MCSTask(liMolByteVec.get(i), bvFrag);
					
					List<MCSResult> liResult = processMCS(mcsTask);
					
					for (MCSResult mcsResult : liResult) {
						
						Canonizer can = new Canonizer(mcsResult.molMCS);
						
						ByteVec bvMCS = new ByteVec(can.getIDCode());
						
						if(!hmMCS_MCS.containsKey(bvMCS)){
							
							hmMCS_MCS.put(bvMCS, bvMCS);
							
							pipeFragByteVec.addData(bvMCS);
							
						}
					}
				}
			}
		}

		private List<MCSResult> processMCS(MCSTask  mcsTask){
			
			List<MCSResult> liResult=null;
			try {
				
				StereoMolecule mol1 = idCodeParser.getCompactMolecule(mcsTask.bcIdCode1.toStringString());
				
				StereoMolecule mol2 = idCodeParser.getCompactMolecule(mcsTask.bcIdCode2.toStringString());
				
				
				if(mol1.getAtoms()>mol2.getAtoms()) {
					mcs.set(mol1, mol2);
				} else {
					mcs.set(mol2, mol1);
				}
				
				List<StereoMolecule> liMolMCS = mcs.getAllCommonSubstructures();
				
				ccMCSCalculations.incrementAndGet();
				
				if(liMolMCS==null){
					liResult=new ArrayList<MCSResult>();
				} else if(!liMolMCS.isEmpty()){
					
					liResult=new ArrayList<MCSResult>();
					
					for (StereoMolecule molMCS : liMolMCS) {
						liResult.add(new MCSResult(mcsTask, molMCS));	
					}
					
				}
				
			} catch (Exception e) {
				e.printStackTrace();
				nFailedSimilarityCalculations++;
			}
			
			return liResult;
		}
		
		public boolean isCalculating(){
			return calculating.get();
		}
		
		

		public ErrorHashMap getEhm() {
			return ehm;
		}

		public long getAdded() {
			return added;
		}
		
		public long getComparisons() {
			return ccMCSCalculations.get();
		}

		public int getArrFailedSimilarityCalculations() {
			return nFailedSimilarityCalculations;
		}
	}
}
