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

import com.actelion.research.calc.INumericalDataColumn;

/**
 * 
 * DoubleArray
 * 26 Jun 2010 MvK: Start implementation
 */
public class DoubleArray implements INumericalDataColumn {
	
	private static final int START_CAPACITY = 32;
	
	private static final int MAX_DELTA_CAPACITY = (int)Math.pow(2, 20);
	
	private double [] data;
	
	private int size;
	
	private int delta_capacity;
	
	public DoubleArray() {
		init(START_CAPACITY);
	}
	
	public DoubleArray(int capacity) {
		init(capacity);
	}

	/**
	 * Deep constructor
	 * @param a
	 */
	public DoubleArray(double[] a) {
		init(a.length);
		System.arraycopy(a,0, data, 0, a.length);
		size = a.length;
	}



	private void init(int capacity){
		data = new double[capacity];
		delta_capacity = capacity/2;
		size = 0;
	}
	
	public double get(int i){
		return data[i];
	}
	
	public double [] get(){

		if(size != data.length) {
			resize(size);
		}
		return data;
	}
	
	public int add(double v){
		data[size]=v;
		
		int index = size;
		
		size++;
		
		if(size==data.length){
			resize(data.length + delta_capacity);
			if(delta_capacity<MAX_DELTA_CAPACITY){
				delta_capacity *= 2;
			}
		}
		
		return index;
	}
	
	public double avr(){
		
		double avr = 0;
		
		for (int i = 0; i < size; i++) {
			avr += data[i];
		}
		
		return avr/size;
	}
	
	public double max(){
		
		double max = Double.MAX_VALUE * -1;
		
		for (int i = 0; i < size; i++) {
			if(data[i]>max)
				max = data[i];
		}
		
		return max;
	}
	
	public double min(){
		
		double min = Double.MAX_VALUE;
		
		for (int i = 0; i < size; i++) {
			if(data[i] < min)
				min = data[i];
		}
		
		return min;
	}
	
	private void resize(int newlen){

		double [] arr = new double [newlen];
		
		System.arraycopy(data, 0, arr, 0, Math.min(data.length, newlen));
		
		data = arr;
		
	}
	
	public int size(){
		return size;
	}

	@Override
	public int getValueCount() {
		return size;
	}

	@Override
	public double getValueAt(int i) {
		return data[i];
	}
}
