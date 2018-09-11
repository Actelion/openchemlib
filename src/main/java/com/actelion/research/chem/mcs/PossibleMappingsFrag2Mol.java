
package com.actelion.research.chem.mcs;

import java.util.ArrayList;
import java.util.List;

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
class PossibleMappingsFrag2Mol {
	
	private PossibleMappingsFrag2Mol parent;
	
	private int indexProcessingArray;
	
	// Index of the fragment atom.
	private int indexAtomFrag;
	
	// Ring bond to an object with a lower index in the array
	private final int [] arrIndexCounterAtomFragRingClosure;
	
	private int ringClosures;
	
	// Array for the mapping molecule atom indices.
	private final int [] arrIndexAtomMolMapping;
	
	private int sizeArrayMappingMolecules;
	
	private final List<PossibleMappingsFrag2Mol> liChild;
	
	public PossibleMappingsFrag2Mol(int indexProcessingArray, int capacityMappings) {
		
		this.indexProcessingArray=indexProcessingArray;
		
		indexAtomFrag=-1;
		
		arrIndexCounterAtomFragRingClosure = new int [SubStructSearchExhaustiveTreeWalker.LEN_RING_CLOSURES];
		
		arrIndexAtomMolMapping = new int [capacityMappings];
		
		sizeArrayMappingMolecules=0;
		
		ringClosures=0;
		
		liChild = new ArrayList<PossibleMappingsFrag2Mol>();
	}
	
	public void reset(){
		
		parent=null;
		
		indexAtomFrag=-1;
		
		ringClosures=0;
		
		sizeArrayMappingMolecules=0;
		
		liChild.clear();
		
	}
	
	public void resetMappingMolecules(){
		sizeArrayMappingMolecules=0;
	}
	
	public PossibleMappingsFrag2Mol(int indexProcessingArray) {
		this(indexProcessingArray, SubStructSearchExhaustiveTreeWalker.LEN_MAPPING_BLOCK);
	}

	public void addChild(PossibleMappingsFrag2Mol child){
		child.setParent(this);
		liChild.add(child);
	}
	
	public void addIndexMappingAtomMolecule(int indexMolAtom){
		arrIndexAtomMolMapping[sizeArrayMappingMolecules]=indexMolAtom;
		sizeArrayMappingMolecules++;
	}
	
	public void addIndexCounterAtomFragmentRingClosure(int indexCounterAtomFragmentRingClosure){
		arrIndexCounterAtomFragRingClosure[ringClosures]=indexCounterAtomFragmentRingClosure;
		ringClosures++;
	}
	
	public int getIndexCounterAtomFragmentRingClosure(int i) {
		return arrIndexCounterAtomFragRingClosure[i];
	}
	
	public int getRingClosures() {
		return ringClosures;
	}
	
	public int getIndexAtomFrag() {
		return indexAtomFrag;
	}

	public int getChilds() {
		return liChild.size();
	}
	
	public PossibleMappingsFrag2Mol getChild(int indexChild) {
		return liChild.get(indexChild);
	}

	public boolean empty(){
		return (sizeArrayMappingMolecules==0) ? true : false;
	}
	
	public int getIndexProcessingArray() {
		return indexProcessingArray;
	}
	

	public int pollIndexMappingAtomMolecule(){
		
		int indexAtomMol = arrIndexAtomMolMapping[sizeArrayMappingMolecules-1];
		
		sizeArrayMappingMolecules--;
		
		return indexAtomMol;
		
	}
	
	public PossibleMappingsFrag2Mol getParent(){
		return parent;
	}
	
	private void setParent(PossibleMappingsFrag2Mol parent){
		this.parent=parent;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		if(parent!=null){
			sb.append("Parent " + parent.indexAtomFrag);
			sb.append(", fragment atom " + indexAtomFrag);
		} else {
			sb.append("No parent, fragment atom " + indexAtomFrag);
		}
		
		sb.append(", mapping mol atms ");
		for (int i = 0; i < sizeArrayMappingMolecules; i++) {
			sb.append(arrIndexAtomMolMapping[i]);
			if(i<sizeArrayMappingMolecules-1){
				sb.append(",");
			}
		}
		
		sb.append(".");
		
		return sb.toString();
	}

	public void setIndexAtomFrag(int indexAtomFrag) {
		this.indexAtomFrag = indexAtomFrag;
	}

}