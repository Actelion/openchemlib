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

public class FragmentDefinedByBondsIdCode implements Comparable<FragmentDefinedByBondsIdCode> {

	
	String idcode;
	
	private int frequency;

	private IBitArray bitArray;
	
	public FragmentDefinedByBondsIdCode() {
		idcode = "";
	}
	
	public FragmentDefinedByBondsIdCode(IBitArray f) {
		this.bitArray = f;
	}
	
	public FragmentDefinedByBondsIdCode(FragmentDefinedByBondsIdCode f) {
		this.bitArray = f.bitArray;
		idcode = f.idcode;
		frequency = f.frequency;
		
	}

	public void addBits(IBitArray ba){
		this.bitArray.add(ba);
	}
	
	/**
	 * @return the idcode
	 */
	public String getIdCode() {
		return idcode;
	}

	public IBitArray getBitArray(){
		return bitArray;
	}
	
	/**
	 * @param idcode the idcode to set
	 */
	public void setIdCode(String idcode) {
		this.idcode = idcode;
	}
	
	public int getFrequency() {
		return frequency;
	}

	public void incrementFrequency() {
		this.frequency++;
	}
	
	public void setFrequency(int frequency) {
		this.frequency = frequency;
	}

	public int getBitsSet() {
		if(bitArray==null){
			return 0;
		}
		
		return bitArray.getBitsSet();
	}
	
	public boolean isOverlappingBits(IBitArray ba){
		return bitArray.isOverlap(ba);
	}
	
    public String write2String() throws IOException{
    	
    	StringBuilder sb = new StringBuilder();
    	
    	sb.append(idcode);
    	
    	sb.append("\t");
    	
    	sb.append(frequency);
    	
    	sb.append(" ");

    	sb.append(bitArray.write2String());
    	
    	return sb.toString();
    	
    }

    public static FragmentDefinedByBondsIdCode read(InputStream s) throws IOException{
    	
    	String idcode = parseString(s);

    	int frequency = (int) BitArray128.parseLong(s);
    	
    	IBitArray f = BitArray128.read(s);
    	
    	FragmentDefinedByBondsIdCode fIdCode = new FragmentDefinedByBondsIdCode(f);
    	
    	fIdCode.idcode = idcode;
    	
    	fIdCode.frequency = frequency;
    	
    	return fIdCode;
    	
    }
    
    public static String parseString(InputStream s) throws IOException{
    	
    	int i = -1;
    	StringBuilder sb = new StringBuilder();
    	while('\t' != (i=s.read())){
    		sb.append((char)i); 
    	}
    	    	
    	return sb.toString();
    }


	@Override
	public int compareTo(FragmentDefinedByBondsIdCode o) {

		if(idcode==null && o.idcode==null)
			return 0;
		else if(o.idcode==null)
			return 1;
		else if(idcode==null)
			return -1;

		return idcode.compareTo(o.idcode);
	}
}
