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

package com.actelion.research.util.datamodel;

import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.*;

import com.actelion.research.util.ArrayUtils;
import com.actelion.research.util.BurtleHasher;

public class IntArray implements Serializable {
	
	private static final int START_CAPACITY = 32;
	
	private static final int MAX_DELTA_CAPACITY = (int)Math.pow(2, 20);
	
	private int [] data;
	
	private int size;
	
	private int delta_capacity;
	
	private int hash;
	
	public IntArray() {
		init(START_CAPACITY);
	}
	
	public IntArray(int capacity) {
		init(capacity);
	}

	/**
	 * Shallow constructor.
	 * @param a
	 */
	public IntArray(int [] a) {
		data = a;
		
		size = data.length;
		
		delta_capacity = size/2;
		
		calculateHashCode();
		
	}

	/**
	 * Deep constructor.
	 * @param ia
	 */
	public IntArray(IntArray ia) {
		init(ia.data.length);
		
		System.arraycopy(ia.data, 0, data, 0, ia.size);

		size = ia.size;
		
		delta_capacity = ia.delta_capacity;
		
		hash = ia.hash;
	}
	
	private void init(int capacity){
		data = new int[capacity];
		delta_capacity = capacity/2;
		size = 0;
		hash = -1;
	}
	
    public void calculateHashCode(){
    	int h = BurtleHasher.hashlittle(data, 13, size);
    	
    	hash = h;
    }


    public int getCapacity(){
		return data.length;
	}
    
    /**
     * Computational expensive operation!
     * @param value
     */
    public void removeValue(int value){
    	
    	int ccFound=0;
    	for (int i = 0; i < size; i++) {
			if(data[i]==value){
				ccFound++;
			}
		}
    	
    	if(ccFound==0){
    		return;
    	}
    	
    	int newlen = size-ccFound;
    	
    	int [] arr = new int [newlen];
    	
    	int cc=0;
    	for (int i = 0; i < size; i++) {
			if(data[i]!=value){
				arr[cc++]=data[i];
			}
		}
    	
    	data=arr;
    	
    	size = newlen;
    	
    }
    
    public boolean equals(Object o) {
    	
    	IntArray ia = (IntArray)o;
    	
        boolean eq = true;
        if(length() != ia.length()){
        	return false;
        }
        
        for (int i = 0; i < size; i++) {
          if(data[i] != ia.data[i]) {
              eq = false;
              break;
          }
        }
        
        return eq;
    }

	public int hashCode() {
		return hash;
	}

	public int get(int i){
		return data[i];
	}

	/**
	 *
	 * @return shallow copy.
	 */
	public int [] get(){
		resize(size);
		return data;
	}
	
	public int add(int v){
		
		facultativeResize();
		
		data[size]=v;
		
		int index = size;
		
		size++;

		hash = -1;
		
		return index;
	}

	public int max(){

		int max = Integer.MIN_VALUE;
		for (int i = 0; i < size; i++) {
			if(data[i]>max){
				max=data[i];
			}
		}

		return max;
	}

	private void facultativeResize(){
		
		if(size == data.length){
			
			long newsize = (long)data.length + (long)delta_capacity;
						
			resize(newsize);
			
			if(delta_capacity<MAX_DELTA_CAPACITY){
				delta_capacity *= 2;
			}
		}
	}
	
	public void add(int [] a){

		int newsize = size + a.length;
		
		if(newsize > data.length){
			resize(newsize);
		}
		
		System.arraycopy(a, 0, data, size, a.length);

		size = newsize;

		hash = -1;
		
	}
	
	public void add(List<Integer> li){
		for (int v : li) {
			add(v);
		}
		calculateHashCode();
	}
	
	public void add(byte [] a){
				
		for (int i = 0; i < a.length; i++) {
			add(a[i]);
		}
		
	}

	/**
	 *
	 * @return number of occupied fields.
	 */
	public int length(){
		return size;
	}
	
	/**
	 * 
	 * @return last value in the array and removes it.
	 */
	public int removeLast(){
		int last = get(size-1);
		size--;
		return last;
	}



	private void resize(long newlen){
		
		if(data.length == newlen){
			return;
		}
		
		int intNewlen = 0;
		
		long max = Integer.MAX_VALUE;
		
		if(newlen >= max) {
						
			intNewlen = Integer.MAX_VALUE;
			
			new RuntimeException("Warning! Maximum length of integer array reached.").printStackTrace();
			
		} else {
			intNewlen = (int)newlen;
		}
		
		int [] arr = new int [intNewlen];
		
		System.arraycopy(data, 0, arr, 0, Math.min(data.length, intNewlen));
		
		data = arr;
		
	}
	
	public void set(int index, int value){
		data[index]=value;
	}
	
	public List<Integer> toList() {
		List<Integer> li = new ArrayList<Integer>(length());

		for (int i = 0; i < length(); i++) {
			li.add(get(i));
		}
		
		return li;
	}

	public void clear(){
		reset();
	}
	
	public void reset(){
		size=0;
	}
	
    public String toString() {
    	
    	StringBuilder sb = new StringBuilder();
    	
        DecimalFormat nf = new DecimalFormat("0");

        for (int i = 0; i < size; i++) {
            sb.append(nf.format(data[i]) + " ");
        }

        return sb.toString();
    }
    
    public String toString(String seperator) {
    	
    	int types = length();
		
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < types; i++) {
			
			sb.append(get(i));
			
			if(i < types-1){
				sb.append(seperator);
			}
			
		}
					
		return sb.toString();
    }
    
    public static IntArray read(InputStream s) throws IOException{
    	
    	
    	int size = parseInteger(s);
    	
    	int delta_capacity = parseInteger(s);
    	
    	int hash = parseInteger(s);
    	
    	int [] a = new int [size];
    	
    	for (int i = 0; i < size; i++) {
    		a[i] = parseInteger(s);
		}
    	

    	IntArray ia = new IntArray();
    	
    	ia.data = a;
    	
    	ia.size = size;
    	
    	ia.delta_capacity = delta_capacity;
    	
    	ia.hash = hash;
    	
    	return ia;
    	
    }
    
    public String write2String() throws IOException{
    	
    	StringBuilder sb = new StringBuilder();

    	sb.append(size);
    	sb.append(" ");
    	sb.append(delta_capacity);
    	sb.append(" ");
    	sb.append(hash);
    	sb.append(" ");
    	sb.append(toString());
    	
    	return sb.toString();
    	
    }
    
    public static int parseInteger(InputStream s) throws IOException{
    	
    	int i = -1;
    	StringBuilder sb = new StringBuilder();
    	while(' ' != (i=s.read())){
    		
    		if(i==-1){
    			break;
    		}
    		
    		sb.append((char)i); 
    	}
    	
    	int val = Integer.parseInt(sb.toString());
    	
    	return val;
    }

    public void sort(){
        resize(size);
        Arrays.sort(data);
    }

    public static void shuffle(IntArray arr){
    
    	Random rnd = new Random();
    	
    	int cycles = 7;
    	
    	int size = arr.length();
    	
    	for (int i = 0; i < cycles; i++) {
			
    		for (int j = 0; j < size; j++) {
				int dest = rnd.nextInt(size);
				
				if(dest==j){
					continue;
				}
				
				int v = arr.get(j);
				
				arr.set(j, arr.get(dest));
				
				arr.set(dest, v);
			}
    		
		}
    	
    }
    
    public static IntArray read(String l){
    	IntArray ia = new IntArray();
    	StringTokenizer st = new StringTokenizer(l, ", ;");
    	while(st.hasMoreTokens()){
    		String t = st.nextToken().trim();
    		int i = Integer.parseInt(t);
    		ia.add(i);
    	}
    	return ia;
    }

    public static boolean equals(int [] a, int [] b){

		boolean eq = true;
		if(a.length != b.length){
			return false;
		}

		for (int i = 0; i < a.length; i++) {
			if(a[i] != b[i]) {
				eq = false;
				break;
			}
		}

		return eq;

	}

	public static List<Integer> toList(int [] a) {
		List<Integer> li = new ArrayList<>(a.length);

		for (int i = 0; i < a.length; i++) {
			li.add(a[i]);
		}
		
		return li;
	}

	public static int [] resize(int [] data, int newlen){
		int [] arr = null;

		if(data.length == newlen){
			return arr;
		}

		int intNewlen = 0;

		long max = Integer.MAX_VALUE;

		if(newlen >= max) {

			intNewlen = Integer.MAX_VALUE;

			new RuntimeException("Warning! Maximum length of integer array reached.").printStackTrace();

		} else {
			intNewlen = newlen;
		}

		arr = new int [intNewlen];

		System.arraycopy(data, 0, arr, 0, Math.min(data.length, intNewlen));

		return arr;

	}

}
