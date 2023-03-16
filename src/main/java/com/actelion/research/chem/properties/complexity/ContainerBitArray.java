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

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.util.Formatter;
import com.actelion.research.util.datamodel.IntArray;

/**
 * Contains a list of reusable bit arrays.
 */
public class ContainerBitArray {
	
	private static boolean ELUSIVE = false;

	private static final int LIMIT2FULL = 1500 * 1000 * 1000;

	/**
	 * If we reach Integer.MAX_VALUE an error will be thrown.
	 */
	private static final int CAPACITY_FULL = Integer.MAX_VALUE-1;
	
	
	private List<IBitArray> li;
			
	private IntArray arrAvailable;

	private int bits;
	
	private IBitArrayFactory<? extends IBitArray> bitArrayCreator;
	
	/**
	 * 
	 * @param bits is the number of available bits in the binary vector. 
	 * @param capacity
	 */
	public ContainerBitArray(int bits, int capacity) {
		
		if(bits == BitArray128.MAX_NUM_BITS){
			bitArrayCreator = new BitArray128Factory();
		} else {
			throw new RuntimeException("Do not know a factory to construct " + bits + " bits array.");
		}
		
		this.bits = bits;
		
		if(capacity > LIMIT2FULL) {
			capacity = CAPACITY_FULL;
		}

		if(isELUSIVE()){
			System.out.println("ContainerBitArray(...) capacity " + Formatter.group(capacity) + ".");
		}

		arrAvailable = new IntArray(capacity);
		
		li = new ArrayList<>(capacity);

		addResources(capacity);
	}
	
	@SuppressWarnings("unchecked")
	public void calculateHash(IBitArray f){
		((IBitArrayFactory<IBitArray>)bitArrayCreator).calculateHash(f);
	}

	public int getSizeBinaryArray(){
		return bits;
	}
	
	public void reset(){
		
		arrAvailable.reset();
		for (int i = 0; i < li.size(); i++) {
			arrAvailable.add(i);
		}
	}
	
	
	private void addResources(int capacity) {
		
		if(li.size() == Integer.MAX_VALUE){
			new RuntimeException("Maximum capacity reached").printStackTrace();
			return;
		}
		
		int indexStart = li.size();
		
		for (int i = 0; i < capacity; i++) {
			
			int index = indexStart+i;
			
			li.add(bitArrayCreator.getNew(index));
			
			arrAvailable.add(index);
			
			if(li.size() == Integer.MAX_VALUE){
				new RuntimeException("Maximum capacity reached").printStackTrace();
				break;
			}
		}
	}
	/**
	 * 
	 * @return a fresh (reset) instance.
	 */
	public IBitArray get(){
		
		if(arrAvailable.length()==0){
			throw new CapacityReachedError("Maximum capacity " + Formatter.group(arrAvailable.getCapacity()) + " reached!");
		}
		
		int index = arrAvailable.removeLast();
		
		IBitArray bitArray = li.get(index);
		
		bitArray.reset();
		
		return bitArray;
	}

	/**
	 * The garbage collector. Put BitArray here if not used anymore.
	 * @param liv
	 */
	public void receycle(IBitArray liv){
		arrAvailable.add(liv.getIndex());
	}
	
	public IBitArray getWithCopy(BitArray128 orign){
		
		IBitArray liv = get();
		
		liv.copyIntoThis(orign);
		
		return liv;
	}
	
	public int getCapacity(){
		return li.size();
	}
	
	public int getAvailable(){
		return arrAvailable.length();
	}
	
	/**
	 * @return the eLUSIVE
	 */
	public static boolean isELUSIVE() {
		return ELUSIVE;
	}

	/**
	 * @param elusive the eLUSIVE to set
	 */
	public static void setELUSIVE(boolean elusive) {
		ELUSIVE = elusive;
	}

	public static class CapacityReachedError extends RuntimeException {
		public CapacityReachedError(String msg) {
			super(msg);
		}
	}

}
