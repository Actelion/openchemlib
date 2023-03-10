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

import java.io.IOException;
import java.io.InputStream;

public class BitArray128 extends IndexHash implements IBitArray {

	public static final int HASH_NOT_SET=-1;
	public static final int HASH_SIDE_STEP =-2;

	// Two long
	public static final int MAX_NUM_BITS = 128;
	
	
	// Each bit represents one bond
	long l1;
	long l2;

	// bit index after the last bit set
	protected char sizeAfterLastBitSet;

	/**
	 * 
	 */
	public BitArray128() {
		sizeAfterLastBitSet = 0;
	}
	
	public BitArray128(int index) {
		super(index);
	}
	
	/**
	 * Deep copy
	 * @param f
	 */
	public BitArray128(BitArray128 f) {
		copyIntoThis(f);
	}
	
	/**
	 * Deep copy
	 * The index is not copied
	 * @param orign
	 */
	public void copyIntoThis(IBitArray orign) {
		
		BitArray128 ba = (BitArray128)orign;
		
		setHash(ba.hash);
		
		l1 = ba.l1;
		
		l2 = ba.l2;
		
	}
	
	@Override
	public boolean equals(Object obj) {
		
		if(!(obj instanceof BitArray128)){
			return false;
		}
		
		BitArray128 f = (BitArray128)obj;
		
		if(l1 != f.l1) {
			return false;
		} else if(l2 != f.l2){
			return false;
		}
		
		return true;
		
	}
	
	public void add(IBitArray f){
		int size = f.getSizeAfterLastBitSet();
		
		for (int i = 0; i < size; i++) {
			if(f.isBitSet(i)){
				setBit(i);
			}
		}
	}
	
	public boolean isBitSet(int i) {
    	boolean set=false;
    	
    	int indInInt = Long.SIZE - 1 - (i % Long.SIZE);
    	
        long mask = 1;
        
        mask = mask << indInInt;
        
    	if(i < Long.SIZE) {
    		if((l1 & mask) != 0)
    			set = true;
    	} else {
   		 	if((l2 & mask) != 0)
   		 		set = true;
    	}
        
        return set;
        
    }
	
    public void setBit(int i) {
    	
    	int indInInt = Long.SIZE - 1 - (i % Long.SIZE);
    	
    	long mask = 1;
        
        mask = mask << indInInt;
        
    	if(i < Long.SIZE) {
    		l1 = l1 | mask;
    	} else {
    		l2 = l2 | mask;
    	}
    	
    	setHash(BitArray128.HASH_NOT_SET);
    	
    	if(i >= sizeAfterLastBitSet){
    		sizeAfterLastBitSet = (char)(i + 1);
    	}
    }
    
    public void unsetBit(int i) {
    	
    	int indInInt = Long.SIZE - 1 - (i % Long.SIZE);
    	
    	long mask = 1;
        
        mask = mask << indInInt;
        
    	if(i < Long.SIZE) {
    		l1 = l1 & ~mask;
    	} else {
    		l2 = l2 & ~mask;
    	}
    	
    	setHash(BitArray128.HASH_NOT_SET);
    	
    	if((i+1) == sizeAfterLastBitSet){
    		sizeAfterLastBitSet = 0;
    		for (int j = i; j >= 0; j++) {
    			if(isBitSet(j)){
    				sizeAfterLastBitSet = (char)(j + 1);
    			}
			}
    	}
    	
    }
    
    public boolean isOverlap(IBitArray f) {
		boolean overlap=false;
		
		final int size = f.getSizeAfterLastBitSet();
		
		for (int i = 0; i < size; i++) {
			if(isBitSet(i) && f.isBitSet(i)){
				overlap=true;
				break;
			}
		}
		
		return overlap;
	}
    
    public int getBitsSet(){
    	int bits=0;
    	
    	int size = getSizeAfterLastBitSet();
    	
    	for (int i = 0; i < size; i++) {
    		if(isBitSet(i)){
    			bits++;
    		}
    	}
    	
    	return bits;
    }
    
    /**
	 * @return the sizeAfterLastBitSet
	 */
	public int getSizeAfterLastBitSet() {
		return sizeAfterLastBitSet;
	}

	public void reset() {
    	
		setHash(BitArray128.HASH_NOT_SET);

    	l1 = 0;
    	l2 = 0;
    	
		sizeAfterLastBitSet = 0;
    }
    
    /**
	 * @return the l1
	 */
	public long getL1() {
		return l1;
	}

	/**
	 * @param l1 the l1 to set
	 */
	public void setL1(long l1) {
		this.l1 = l1;
	}

	/**
	 * @return the l2
	 */
	public long getL2() {
		return l2;
	}

	/**
	 * @param l2 the l2 to set
	 */
	public void setL2(long l2) {
		this.l2 = l2;
	}

	@Override
	public String toString() {
		return toStringBinary();
	}
	
	public String toStringBinary() {
        StringBuilder sb = new StringBuilder();
        
        sb.append(toStringBinary(l1));
        sb.append(toStringBinary(l2));
        
        return sb.toString();
    }

	private static String toStringBinary(long v) {
		
		StringBuilder sb = new StringBuilder();
		
		int len = Long.SIZE;

		for (int i = 0; i < len; i++) {
			if ((v & 1) == 1) {
				sb.insert(0, "1");
			} else {
				sb.insert(0, "0");
			}
			
			
			v = v >> 1;
		}

		return sb.toString().trim();
	}
	
    public static BitArray128 read(InputStream s) throws IOException{
    	    	
    	BitArray128 f = new BitArray128();
    	
    	f.setIndex((int)parseLong(s));
    	
    	f.setHash((int)parseLong(s));
    	
    	f.l1 =  parseLong(s);
    	
    	f.l2 =  parseLong(s);
    	
    	f.sizeAfterLastBitSet =  (char)parseLong(s);
    	
    	return f;
    	
    }

    public String write2String() {
    	
    	StringBuilder sb = new StringBuilder();

    	sb.append(getIndex());
    	
    	sb.append(" ");
    	
    	sb.append(hashCode());
    	
    	sb.append(" ");
    	
    	sb.append(l1);
    	
    	sb.append(" ");
    	
    	sb.append(l2);
    	
    	sb.append(" ");
    	
    	sb.append(sizeAfterLastBitSet);
    	
    	return sb.toString();
    	
    }
    
    public static long parseLong(InputStream s) throws IOException{
    	
    	int i = -1;
    	StringBuilder sb = new StringBuilder();
    	while(' ' != (i=s.read())){
    		
    		if(i==-1){
    			break;
    		}
    		
    		sb.append((char)i); 
    	}
    	
    	long val = Long.parseLong(sb.toString());
    	
    	return val;
    }

}
