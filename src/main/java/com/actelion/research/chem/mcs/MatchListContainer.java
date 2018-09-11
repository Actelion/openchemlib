package com.actelion.research.chem.mcs;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

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
class MatchListContainer {
	
	private static final double fractionCapacity = 0.25; 
	
	private MatchList [] arrMatchList;
	
	private HashSet<MatchList> hsMatchList;
	
	private LinkedList<Integer> liIndexAvailable;
	
	private LinkedList<Integer> liIndexUsed;
	
	private List<int []> liMatches; 
	
	public MatchListContainer(int maxNumAtomsMolecule) {
		
		int maximumCapacity = (int)(maxNumAtomsMolecule * fractionCapacity);
		
		arrMatchList = new MatchList[maximumCapacity];
		
		liIndexAvailable = new LinkedList<Integer>();
		
		liIndexUsed = new LinkedList<Integer>();
		
		int nInteger = (maxNumAtomsMolecule+Integer.SIZE-1) / Integer.SIZE;
		
		for (int i = 0; i < arrMatchList.length; i++) {
			arrMatchList[i]=new MatchList(nInteger);
			liIndexAvailable.add(i);
		}
		
		hsMatchList = new HashSet<MatchList>();
		
		liMatches = new ArrayList<int[]>();
	}
	
	public void reset() {
		
		hsMatchList.clear();
		
		for(int index : liIndexUsed) {
			arrMatchList[index].reset();
			liIndexAvailable.add(index);
		}
		
		liIndexUsed.clear();
		
		liMatches = new ArrayList<int[]>();
	}
	
	/**
	 * 
	 * @param arrMatchListFragment match list with length=number of atoms in fragment. 
	 * The fields contain the matching atom indices of the molecule.
	 */
	public void add(int [] arrMatchListFragment) {
		
//		if(liIndexAvailable.isEmpty()){
//			throw new RuntimeException("Maximum capacity ("+arrMatchList.length+")of match list exceeded.");
//		}
//		
//		int index = liIndexAvailable.poll();
//		
//		arrMatchList[index].set(arrMatchListFragment);
//		
//		if(hsMatchList.add(arrMatchList[index])){
//			
//			liIndexUsed.add(index);
//			
//			liMatches.add(arrMatchListFragment);
//			
//		} else {
//			
//			arrMatchList[index].reset();
//			
//			liIndexAvailable.add(index);
//		}
		
		
		liMatches.add(arrMatchListFragment);
	}
	
	public List<int[]> getMatchList(){
		return liMatches;
	}
	
	public boolean isMaxCapacityReached(){
		return liIndexAvailable.isEmpty();
	}
}
