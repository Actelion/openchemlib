package com.actelion.research.chem.mcs;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.util.datamodel.IntArray;

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
public class ContainerListWithIntVec {
	
	private static int CAPACITY_ADD = 1024;
	
	
	private List<ListWithIntVec> li;
		
	private int sizeVector;
	
	private IntArray arrAvailable;
	
	/**
	 * 
	 * @param sizeVector size of the integer array.
	 * @param capacity
	 */
	public ContainerListWithIntVec(int sizeVector, int capacity) {
		
		this.sizeVector = sizeVector;
		
		arrAvailable = new IntArray(capacity);
		
		li = new ArrayList<ListWithIntVec>(capacity);
		
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
			li.add(new ListWithIntVec(sizeVector, index));
			
			arrAvailable.add(index);
		}
	}
	/**
	 * 
	 * @return a fresh (reseted) instance.
	 */
	public ListWithIntVec get(){
		
		if(arrAvailable.length()==0){
			
			initResources(CAPACITY_ADD);
		}
		
		int index = arrAvailable.removeLast();
		
		ListWithIntVec liv = li.get(index);
		
		liv.reset();
		
		return liv;
	}
	
	public void back(ListWithIntVec liv){
		arrAvailable.add(liv.getPositionInContainer());
	}
	
	public ListWithIntVec getWithCopy(ListWithIntVec orign){
		
		ListWithIntVec liv = get();
		
		liv.copyIntoThis(orign);
		
		return liv;
	}

}
