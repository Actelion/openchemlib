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

import com.actelion.research.util.ConstantsDWAR;
import com.actelion.research.util.datamodel.IntArray;

import java.util.ArrayList;
import java.util.List;

/**
 * 
 * FlexophorePoint
 * @author Modest von Korff
 * @version 1.0
 * Jan 15, 2013 MvK Start implementation
 */
public class FlexophorePoint {
	
	private static final int CAPACITY_IATYPE = 6;
	
	private static final int CAPACITY_ATOMINDEX = 12;
		
	private IntArray iaInteractionType;
		
	private IntArray iaOriginalAtomIndex;
	
	private List<AtomIndexLinkerId> liAtomIndexLinkerId;
	
	private String idcode;
	
	public int id;
	
	public int hash;
	/**
	 * 
	 */
	public FlexophorePoint(int id) {
		
		this.id = id;
		
		hash = -1;
		
		iaInteractionType = new IntArray(CAPACITY_IATYPE);
		
		iaOriginalAtomIndex = new IntArray(CAPACITY_ATOMINDEX);
		
		liAtomIndexLinkerId = new ArrayList<AtomIndexLinkerId>();
	}
	
	public int getId() {
		return id;
	}

	public boolean equals(Object obj) {
		
		boolean equal = true;
		
		if(!(obj instanceof FlexophorePoint)){
			return false;
		}
		
		FlexophorePoint a = (FlexophorePoint)obj;
		
		if(!idcode.equals(a.idcode)){
			equal = false;
		}else if(!iaInteractionType.equals(a.iaInteractionType)){
			equal = false;
		}
					
		return equal;
	}
	
	public void calculateHashCode() {
		
		if(idcode==null){
			hash=-1;
			return;
		}
		
		int hashIdCode = idcode.hashCode();
		
		int hashArrayInteractionType = iaInteractionType.hashCode();
				
		hash = hashArrayInteractionType ^  hashIdCode;

		
	}
	
	
	public int hashCode() {
		return hash;
	}
		
	public void setInteractionType(byte [] arrInteractionType){
		iaInteractionType.clear();
		iaInteractionType.add(arrInteractionType);
		iaInteractionType.calculateHashCode();
		calculateHashCode();
	}
	
	public IntArray getInteractionTypes(){
		return iaInteractionType;
	}
	
	public List<Integer> getOriginalAtomIndex(){
		return iaOriginalAtomIndex.toList();
	}
	
	public void clearOriginalAtomIndex(){
		iaOriginalAtomIndex.clear();
	}


	public void addOriginalAtomIndex(int originalAtomIndex){
		this.iaOriginalAtomIndex.add(originalAtomIndex);
	}
	
	public List<AtomIndexLinkerId> getAtomIndexLinkerId() {
		return liAtomIndexLinkerId;
	}

	public void addArrOriginalAtomIndexRGroups(AtomIndexLinkerId atomIndexLinkerId) {
		liAtomIndexLinkerId.add(atomIndexLinkerId);
	}


	public String getIdCode() {
		return idcode;
	}

	public void setIdCode(String idcode) {
		this.idcode = idcode;
		calculateHashCode();
	}

	
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("FlexophorePoint [iaInteractionType=");
		sb.append(iaInteractionType);
		sb.append(", iaOriginalAtomIndex=");
		sb.append(iaOriginalAtomIndex);
		sb.append(", liAtomIndexLinkerId=");
		sb.append(liAtomIndexLinkerId);
		sb.append(", idcode=");
		sb.append(idcode);
		sb.append(", id=");
		sb.append(id);
		sb.append(", hash=");
		sb.append(hash);
		sb.append("]");
		return sb.toString();
	}

	public String toStringInteractionTypes() {
		
		int types = iaInteractionType.length();
		
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < types; i++) {
			
			sb.append(iaInteractionType.get(i));
			
			if(i < types-1){
				sb.append(ConstantsDWAR.SEP_VALUE);
			}
			
		}
					
		return sb.toString();
	}
	
	
	
}