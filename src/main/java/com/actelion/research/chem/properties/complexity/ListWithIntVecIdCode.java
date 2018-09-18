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

import com.actelion.research.chem.mcs.ListWithIntVec;
import com.actelion.research.util.datamodel.IntArray;

public class ListWithIntVecIdCode extends ListWithIntVec {
	
	private String idcode;

	private int frequency;
		
	public ListWithIntVecIdCode(int lengthIntVector) {
		super(lengthIntVector);
	}
	
	/**
	 * Flat constructor.
	 * @param liv
	 */
	public ListWithIntVecIdCode(ListWithIntVec liv) {
		super(liv);
	}
	
	/**
	 * Deep constructor.
	 * @param liv
	 */
	public ListWithIntVecIdCode(ListWithIntVecIdCode liv) {
		super(liv);
		this.idcode = liv.idcode; 
		this.frequency = liv.frequency; 
	}

	
	public void add(ListWithIntVecIdCode liv){
	
		addAllBits(liv);
		
		frequency += liv.getFrequency();
	}
	
	public String getIdCode() {
		return idcode;
	}

	public void setIdcode(String idcode) {
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

	public String toString(){
		StringBuilder sb = new StringBuilder();
		
		sb.append(frequency);
		sb.append(" ");
		sb.append(idcode);
		sb.append(" ");
		sb.append(super.toString());
		
		return sb.toString();
	}

    public static ListWithIntVecIdCode read(InputStream s) throws IOException{
    	
    	String idcode = parseString(s);
    	
    	int frequency = IntArray.parseInteger(s);
    	
    	ListWithIntVec liv = ListWithIntVec.read(s);
    	    	    	
    	ListWithIntVecIdCode livIdCode = new ListWithIntVecIdCode(liv);
    	
    	livIdCode.idcode = idcode;
    	
    	livIdCode.frequency = frequency;
    	    	
    	return livIdCode;
    	
    }
    
    public String write2String() throws IOException{
    	
    	StringBuilder sb = new StringBuilder();

    	sb.append(idcode);
    	sb.append("\t");
    	sb.append(frequency);
    	sb.append(" ");
    	
    	sb.append(super.write2String());
    	
    	return sb.toString();
    	
    }
    
    public static String parseString(InputStream s) throws IOException{
    	
    	int i = -1;
    	StringBuilder sb = new StringBuilder();
    	while('\t' != (i=s.read())){
    		sb.append((char)i); 
    	}
    	    	
    	return sb.toString();
    }

}