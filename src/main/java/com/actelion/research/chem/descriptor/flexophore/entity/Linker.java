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
 * @author Modest v. Korff
 */

package com.actelion.research.chem.descriptor.flexophore.entity;

import com.actelion.research.util.datamodel.ByteVec;
import com.actelion.research.util.datamodel.IntArray;

/**
 * 
 * Linker
 * @author Modest von Korff
 * @version 1.0
 * Jan 15, 2013 MvK Start implementation
 */
public class Linker {
	
	private static final int CAPACITY_ATOMINDEX = 12;
	
	private int id;
	
	private byte [] arrDistanceHistogram;
	
	private String idcode;
	
	private IntArray iaOriginalAtomIndex;
	
	private int idFlexophorePoint1;
	
	private int idFlexophorePoint2;
	
	private int hash;
	
	public Linker(int id) {
		this.id = id;
		
		hash = -1;
		iaOriginalAtomIndex = new IntArray(CAPACITY_ATOMINDEX);
	}
	
	
	public void calculateHashCode() {
		
		if(arrDistanceHistogram == null || idcode==null){
			hash=-1;
			return;
		}
		
		int hashHist = ByteVec.getHashCode(arrDistanceHistogram);
		
		int hashIdCode = idcode.hashCode();
		
		hash = hashHist ^  hashIdCode;
	}
	
	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}

	public int hashCode() {
		
		return hash;
	}

	public void addOriginalAtomIndex(int atomIndex){
		iaOriginalAtomIndex.add(atomIndex);
	}
	
	public void addOriginalAtomIndex(int [] arrAtomIndex){
		iaOriginalAtomIndex.add(arrAtomIndex);
	}
	
	
	public int [] getOriginalAtomIndex(){
		return iaOriginalAtomIndex.get();
	}

	/**
	 * @return the arrDistanceHistogram
	 */
	public byte[] getDistanceHistogram() {
		return arrDistanceHistogram;
	}


	/**
	 * @param arrDistanceHistogram the arrDistanceHistogram to set
	 */
	public void setDistanceHistogram(byte[] arrDistanceHistogram) {
		this.arrDistanceHistogram = arrDistanceHistogram;
		calculateHashCode();
	}


	/**
	 * @return the idcode
	 */
	public String getIdCode() {
		return idcode;
	}


	/**
	 * @param idcode the idcode to set
	 */
	public void setIdCode(String idcode) {
		this.idcode = idcode;
		calculateHashCode();
	}
	
	
	public boolean equals(Object obj) {
		boolean equal = true;
		
		if(!(obj instanceof Linker)){
			return false;
		}
		
		Linker l = (Linker)obj;
		
		if((arrDistanceHistogram == null) && (l.arrDistanceHistogram != null)){
			return false;
		} else if((arrDistanceHistogram != null) && (l.arrDistanceHistogram == null)){
			return false;
		} else if((idcode == null) && (l.idcode != null)){
			return false;
		} else if((idcode != null) && (l.idcode == null)){
			return false;
		} 
		
		if(arrDistanceHistogram != null){
			
			if(l.arrDistanceHistogram != null){
				
				if(arrDistanceHistogram.length == l.arrDistanceHistogram.length) {
					
					for (int i = 0; i < arrDistanceHistogram.length; i++) {
						if(arrDistanceHistogram[i] != l.arrDistanceHistogram[i]){
							equal = false;
							break;
						}
					}
				}
			}
		}
		
		if(idcode != null){
			if(l.idcode != null){
				if(!idcode.equals(l.idcode)){
					equal = false;
				}
			}
		}
		
		return equal;
	}
	
	/**
	 * @return the idFlexophorePoint1
	 */
	public int getIdFlexophorePoint1() {
		return idFlexophorePoint1;
	}


	/**
	 * @param idFlexophorePoint1 the idFlexophorePoint1 to set
	 */
	public void setIdFlexophorePoint1(int idFlexophorePoint1) {
		this.idFlexophorePoint1 = idFlexophorePoint1;
	}


	/**
	 * @return the idFlexophorePoint2
	 */
	public int getIdFlexophorePoint2() {
		return idFlexophorePoint2;
	}


	/**
	 * @param idFlexophorePoint2 the idFlexophorePoint2 to set
	 */
	public void setIdFlexophorePoint2(int idFlexophorePoint2) {
		this.idFlexophorePoint2 = idFlexophorePoint2;
	}

	

}