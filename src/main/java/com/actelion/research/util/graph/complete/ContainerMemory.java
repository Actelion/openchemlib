package com.actelion.research.util.graph.complete;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.util.datamodel.IntArray;

/**
 * 
 * 
 * ContainerMemory
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Sep 27, 2012 MvK: Start implementation
 */
public class ContainerMemory<S extends AMemorizedObject> {
	
	
	private static int CAPACITY_ADD = 1024;
	
	
	private List<S> li;
	
	private IntArray arrAvailable;
	
	private IFactory<S> factory;
	
	/**
	 * 
	 * @param capacity
	 * @param factory
	 */
	public ContainerMemory(int capacity, IFactory<S> factory) {
		
		this.factory = factory;
		
		arrAvailable = new IntArray(capacity);
		
		li = new ArrayList<S>(capacity);
		
		initResources(capacity);
	}
	
	public void reset(){
		
		arrAvailable.reset();
		for (int i = 0; i < li.size(); i++) {
			arrAvailable.add(i);
		}
	}
	
	
	private void initResources(int capacity) {
		
		int indexStart = li.size();
		
		for (int i = 0; i < capacity; i++) {
			
			int index = indexStart+i;
			
			S s = factory.createObject();
			
			s.setPositionInContainer(index);
			
			li.add(s);
			
			arrAvailable.add(index);
		}
	}
	/**
	 * 
	 * @return a fresh (resetted) instance.
	 */
	public S get(){
		
		if(arrAvailable.length()==0){
			
			initResources(CAPACITY_ADD);
		}
		
		int index = arrAvailable.removeLast();
		
		S s = li.get(index);
		
		s.reset();
		
		return s;
	}
	
	public void back(S s){
		arrAvailable.add(s.getPositionInContainer());
	}
	
	public S getWithCopy(S orign){
		
		S t = get();
		
		t.copyIntoThis(orign);
		
		return t;
	}

	
}
