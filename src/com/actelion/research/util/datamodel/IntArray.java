package com.actelion.research.util.datamodel;

import java.io.IOException;
import java.io.InputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.StringTokenizer;

import com.actelion.research.util.BurtleHasher;

/**
 * 
 * IntArray
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * 7 Apr 2010 MvK: Start implementation
 */
public class IntArray {
	
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
	
	public IntArray(int [] a) {
		data = a;
		
		size = data.length;
		
		delta_capacity = size/2;
		
		calculateHashCode();
		
	}
	
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
	
	public int [] get(){
		resize(size);
		return data;
	}
	
	public int add(int v){
		
		facultativeResize();
		
		data[size]=v;
		
		int index = size;
		
		size++;
		
		facultativeResize();
		
		hash = -1;
		
		return index;
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
		
		int indexStart = size;
		
		size += a.length;
		
		while(size >= data.length){
			
			long newsize = (long)data.length + (long)delta_capacity;
			
			resize(newsize);
			
			if(delta_capacity<MAX_DELTA_CAPACITY){
				delta_capacity *= 2;
			}
		}
		
		System.arraycopy(a, 0, data, indexStart, a.length);
		
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
    		String t = st.nextToken();
    		
    		int i = Integer.parseInt(t);
    		
    		ia.add(i);
    		
    	}
    	
    	return ia;
    }
    
	
}
