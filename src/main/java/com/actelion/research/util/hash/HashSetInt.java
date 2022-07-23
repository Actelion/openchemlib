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
 */

package com.actelion.research.util.hash;

import com.actelion.research.calc.ArrayUtilsCalc;

import java.util.List;

/**
 * 
 * HashSetInt
 * @author Modest von Korff
 * @version 1.0
 * 4 Mar 2010 MvK: Start implementation
 */
public class HashSetInt {
	
	
	private int [][] data;
	
    /**
     * The default initial capacity - MUST be a power of two.
     */
	private static final int DEFAULT_INITIAL_CAPACITY = 256;
    
    /**
     * The maximum capacity, used if a higher value is implicitly specified
     * by either of the constructors with arguments.
     * MUST be a power of two <= 1<<30.
     */
    private static final int MAXIMUM_CAPACITY = 1 << 30;

    /**
     * The load factor used when none specified in constructor.
     **/
    private static final float DEFAULT_LOAD_FACTOR = 0.75f;
    
    
	private int threshold;
	
	private int size;
	
	private double loadFactor;
	
	public HashSetInt() {
		
		data =  new int [DEFAULT_INITIAL_CAPACITY][];
 		
		loadFactor = DEFAULT_LOAD_FACTOR;
		
		threshold = (int)(DEFAULT_INITIAL_CAPACITY * loadFactor);
		
	}
	
	public HashSetInt(int capacity) {
		double log2 = Math.log10(capacity)/Math.log10(2);
		int capPowOf2 = (int)(Math.pow(2, (int)(log2+1)));
		data =  new int [capPowOf2][];
		loadFactor = DEFAULT_LOAD_FACTOR;
		threshold = (int)(capPowOf2 * loadFactor);
	}
	
	public HashSetInt(List<Integer> li) {
		this(li.size());
		for (Integer v : li) {
			add(v);
		}
	}
	
	public HashSetInt(int [] a) {
		this(a.length);
		for (int v : a) {
			add(v);
		}
	}
	
    /**
     * Returns index for hash code h. 
     */
    static int indexFor(int h, int length) {
        return h & (length-1);
    }
    
    public void clear(){
    	for (int i = 0; i < data.length; i++) {
    		data[i]=null;
		}
    	size=0;
    }
    
    /**
     * If the key is already present nothing will happen.
     * @param v
     * @return false if key is already present.
     */
    public boolean add(int v) {
    	
    	if(isEntry(v)){
    		return false;
    	}
    	
    	int hash = v;
    	
    	int indexMap = indexFor(hash, data.length);
    	
    	int [] arrRowNew = null;
    	if(data[indexMap]==null) {
    		arrRowNew = new int [1];
    	} else {
    		arrRowNew = new int [data[indexMap].length+1];
    		System.arraycopy(data[indexMap], 0, arrRowNew, 0, data[indexMap].length);
    	}
    	
    	arrRowNew[arrRowNew.length-1]=v;
    	
    	data[indexMap] = arrRowNew;
		
		size++;

		if(size>threshold){
			resize(data.length * 2);
		}
		
		return true;
    }
    
    public boolean add(int [] a) {
    	
    	boolean allAdded = true;
    	
    	for (int i = 0; i < a.length; i++) {
			if(!add(a[i])){
				allAdded = false;
			}
		}
    	
    	return allAdded;
    	
    }
    
    public boolean add(List<Integer> li) {
    	
    	boolean allAdded = true;
    	
    	for (int i = 0; i < li.size(); i++) {
			if(!add(li.get(i))){
				allAdded = false;
			}
		}
    	
    	return allAdded;
    	
    }

   
    
    
    void resize(int newCapacity) {
    	
    	// System.err.println("Resize");
    	
    	int [][] dataNew = new int [newCapacity][];
    	transfer(dataNew);
    	
        int oldCapacity = data.length;
        
        if (oldCapacity == MAXIMUM_CAPACITY) {
        	
            threshold = Integer.MAX_VALUE;
            
            return;
        }
        
        threshold = (int)(newCapacity * loadFactor);
    }
    
    public void remove(int v) {
    	int index = indexFor(v, data.length);
    	
    	int [] arr = data[index];
    	
    	if(arr==null){
    		return;
    	}
    	
    	int indRow = -1;
    	for (int i = 0; i < arr.length; i++) {
			if(arr[i]==v){
				indRow=i;
				break;
			}
		}
    	
    	if(indRow==-1)
    		return;
    	
    	if(arr.length==1){
        	data[index]=null;
    	} else {
    	
	    	int [] arrNew = new int [arr.length-1];
	    	
	    	int cc=0;
	    	for (int i = 0; i < arr.length; i++) {
				if(arr[i]!=v){
					arrNew[cc++]=arr[i];
				}
			}
	    	
	    	data[index]=arrNew;
	    	
    	}
    	
    	size--;
    	
    }
    
    
    void transfer(int [][] dataNew) {
    	
    	int size = data.length;
    	
    	for (int i = 0; i < size; i++) {
			int [] arr = data[i];
			
			if(arr==null)
				continue;
			
			for (int j = 0; j < arr.length; j++) {
				
				int indexNew = indexFor(arr[j], dataNew.length);
				
				int [] arrRowNew = null;
				if(dataNew[indexNew]==null) {
					arrRowNew = new int [1];
				} else {
				
					arrRowNew = new int [dataNew[indexNew].length+1];
		    	
		    		System.arraycopy(dataNew[indexNew], 0, arrRowNew, 0, dataNew[indexNew].length);
				}
			
		    	arrRowNew[arrRowNew.length-1]=arr[j];
				
				dataNew[indexNew]= arrRowNew;
			}
			
		}
    	
    	data = dataNew;
    	
    }
    
    public boolean contains(int v) {
    	return isEntry(v);
    }
    
    private boolean isEntry(int v) {
    	int indexMap = indexFor(v, data.length);
    	
    	int [] arr = data[indexMap];
    	
    	if(arr==null)
    		return false;
    	
    	boolean keyfound=false;
    	
    	for (int i = 0; i < arr.length; i++) {
			
			if(arr[i] == v) {
				keyfound=true;
				break;
			}
		}
    	
    	return keyfound;
    }

    public int size(){
    	return size;
    }
    
    /**
     * Deep copy.
     * @return
     */
    public int [] getValues(){
    	int [] a = new int [size];
    	
    	int cc=0;
    	for (int i = 0; i < data.length; i++) {
    		if(data[i]==null)
    			continue;
    		
			for (int j = 0; j < data[i].length; j++) {
				a[cc++]=data[i][j];
			}
		}
    	
    	return a;
    }
    
    public List<Integer> toList(){
    	return ArrayUtilsCalc.toList(getValues());
    }
    

}
