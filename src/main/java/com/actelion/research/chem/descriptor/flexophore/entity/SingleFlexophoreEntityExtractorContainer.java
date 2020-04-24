package com.actelion.research.chem.descriptor.flexophore.entity;

import com.actelion.research.chem.descriptor.flexophore.UnparametrizedAtomTypeException;
import com.actelion.research.util.Pipeline;
import com.actelion.research.util.datamodel.IDCodeCoord;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;


/**
 * SingleFlexophoreEntityExtractor
 * Extract pharmacophore points and linker.
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Jan 14, 2013 MvK Start implementation
 */
public class SingleFlexophoreEntityExtractorContainer {

	
	private static final int CAPACITY = 10000;
	
	private Pipeline<IDCodeCoord> queuePipe;
	
	private ConcurrentHashMap<FlexophorePoint, String> hmFlexophorePoint;
	
	private ConcurrentHashMap<Linker, String> hmLinker;
	
	private int nCores;
	
	public SingleFlexophoreEntityExtractorContainer() {
		
		hmFlexophorePoint = new ConcurrentHashMap<FlexophorePoint, String>(CAPACITY);
		
		hmLinker = new ConcurrentHashMap<Linker, String>(CAPACITY);
		
		nCores = Runtime.getRuntime().availableProcessors();
		
		if(nCores>1){
			nCores--;
		}
		// nCores = 1;
		
	}
	
	
	public void process(Pipeline<IDCodeCoord> queuePipe) throws UnparametrizedAtomTypeException{
		this.queuePipe = queuePipe;
		
		List<Thread> liThread = new ArrayList<Thread>();
		
		for (int i = 0; i < nCores; i++) {
			SingleFlexophoreEntityExtractorThread extractor = new SingleFlexophoreEntityExtractorThread();
			
			Thread th = new Thread(extractor);
			
			th.start();
			
			liThread.add(th);
		}
		
		boolean allThreadsFinished = false;
		
		while(!allThreadsFinished){
			
			allThreadsFinished = true;
			
			for (Thread th : liThread) {
				if(th.isAlive()){
					allThreadsFinished = false;
				}
			}
			
			try {Thread.sleep(500);} catch (InterruptedException e) {e.printStackTrace();}
		}
		
	}

	
	
		
	public List<FlexophorePoint> getPharmacophorePointList(){
		
		
		
		List<FlexophorePoint> li = Collections.list(hmFlexophorePoint.keys());
		
		return li;

	}
	
	public List<Linker> getLinkerList(){
			
		List<Linker> li = Collections.list(hmLinker.keys());
		
		return li;

	}
	
	class SingleFlexophoreEntityExtractorThread implements Runnable {

		public void run() {
			
			while(!queuePipe.wereAllDataFetched()){
				
				IDCodeCoord idcodeCoord = queuePipe.pollData();
				
				if(idcodeCoord == null){
					try {Thread.sleep(500);} catch (InterruptedException e) {e.printStackTrace();}
					continue;
				}
				
				try {
					process(idcodeCoord.idcode, idcodeCoord.coordinates);
				} catch (UnparametrizedAtomTypeException e) {
					e.printStackTrace();
				} 
				
			}
		}
		
		
		public void process(String idcode, String coordinates) throws UnparametrizedAtomTypeException{
			
			SingleFlexophoreEntityExtractor singleFlexophoreEntityExtractor = new SingleFlexophoreEntityExtractor(idcode, coordinates);
					
			singleFlexophoreEntityExtractor.process();
					
			List<FlexophorePoint> liFlexophorePoint = singleFlexophoreEntityExtractor.getPharmacophorePointList();
			
			for (FlexophorePoint fp : liFlexophorePoint) {
				fp.calculateHashCode();
				hmFlexophorePoint.put(fp, "");
			}
			
			
			List<Linker> liLinker = singleFlexophoreEntityExtractor.getLinkerList();

			for (Linker linker : liLinker) {
				linker.calculateHashCode();
				hmLinker.put(linker, "");
			}
		}
		
	}


}
