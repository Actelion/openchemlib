/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * 3. Neither the name of the copyright holder nor the
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

package com.actelion.research.chem.descriptor.flexophore;


import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.util.BurtleHasher;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class PPNode implements Comparable<PPNode>, IPPNode {

	public static final int NUM_BYTES_INTERACTION_TYPE = 3;

	private static final int BYTES_FREQUENCY_INTERACTION = 1;


	
	public static final int DUMMY_INTERACT_ID = 37;
		
	public static final String SEPARATOR_ATOMS = ",";

	public static final String MULT_FREQ = "*";

	public static final int INFO_DEFAULT = -1;

	private static final int MASK1 = 0x000000FF;
	
	private static final int BUFFER_SIZE= 3 * getNumBytesEntry();

	
	private byte [] arrInteractionType;
	
	private byte size;


	transient boolean heteroAtom;

	transient byte modeFlexophore;


	public PPNode(){
		init();
	}

	public PPNode(PPNode node){
		copy(node);
	}

	/**
	 * Caution! After finishing all adding's realize() has to be called!
	 * @param interactionType
	 */
	public void add(int interactionType){
		
		if(interactionType >= ConstantsFlexophoreGenerator.MAX_VAL_INTERACTION_TYPE)
			throw new RuntimeException("Interaction type " + interactionType + " larger than " + ConstantsFlexophoreGenerator.MAX_VAL_INTERACTION_TYPE + ".");

		boolean addedByIncrement = false;
		for (byte i = 0; i < size; i++) {
			int interactionTypeNode = getInteractionType(i);
			if(interactionType == interactionTypeNode){
				incrementAtomTypeCount(i);
				addedByIncrement = true;
				break;
			}

		}

		if(!addedByIncrement){
			initInteractionType(interactionType);
		}
	}


	/**
	 * Only atoms are added that are not yet in the list,
	 * check PPAtom.equals for comparison.
	 * @param node
	 */
	public void addAtoms(PPNode node){
		
		for (int i = 0; i < node.getInteractionTypeCount(); i++) {
			
			if(!containsInteractionID(node.getInteractionType(i))){
				add(node.getInteractionType(i));
			}
		}
		
	}
		
	public int compareTo(PPNode o) {
		int cmp=0;
		
		
		int size1 = getInteractionTypeCount();
		
		int size2 = o.getInteractionTypeCount();
		
		int max1 = getMaximumInteractionType(this);
		
		int max2 = getMaximumInteractionType(o);
		
		if(max1 > max2)
			cmp=1;
		else if(max1 < max2)
			cmp=-1;
		else {
			if(size1 > size2)
				cmp=1;
			else if(size1 < size2)
				cmp=-1;
			else {
				
				for (int i = 0; i < size1; i++) {
					
					int id1 = getInteractionType(i);
					
					int id2 = o.getInteractionType(i);
					
					if(id1 > id2){
						cmp=1;
						
					} else if(id1 < id2){
						
						cmp=-1;
					}
				}
			}
		}
			
		return cmp;
	}
	
	public boolean containsInteractionID(int interactid){
		boolean contains=false;
		
		if(arrInteractionType==null)
			return contains;
		
		int size = getInteractionTypeCount();

		for (int i = 0; i < size; i++) {
			if(getInteractionType(i) == interactid){
				contains=true;
				break;
			}
		}
		
		return contains;
		
	}
	
	/**
	 * Copy of node into this.
	 * @param node
	 */
	public void copy(PPNode node){
		
		arrInteractionType = new byte[node.arrInteractionType.length];
		
		System.arraycopy(node.arrInteractionType, 0, arrInteractionType, 0, node.arrInteractionType.length);
		
		size = node.size;

		modeFlexophore = node.modeFlexophore;
	}
	
	/**
	 * 
	 * node deep copy.
	 */
	public PPNode getCopy(){
		return new PPNode(this);  
	}
	
	public boolean equals(Object o) {
		PPNode n = (PPNode)o;
		return equalAtoms(n);
	}

	@Override
	public int hashCode() {
		return BurtleHasher.hashlittle(arrInteractionType, 13);
	}

	/**
	 * May be called after finishing adding new interaction types. 
	 */
	public void realize(){
		int sizeBytes = size * getNumBytesEntry();
		arrInteractionType = ArrayUtilsCalc.resize(arrInteractionType, sizeBytes);
		sortInteractionTypes();
		calcHasHeteroAtom();
	}
	
	/**
	 * realize() may be called first. 
	 * @param node
	 * @return
	 */
	public boolean equalAtoms(PPNode node) {
		boolean b = true;
		
		int size1 = getInteractionTypeCount();
		
		int size2 = node.getInteractionTypeCount();
		
		if(size1 != size2)
			b = false;
		else {
			
			for (int i = 0; i < size1; i++) {
				
				int id1 = getInteractionType(i);
				
				int id2 = node.getInteractionType(i);
				
				if(id1 != id2){
					b = false;
					break;
				}
			}
		}
			
		return b;
	}
	
	public byte [] get(){
		return arrInteractionType;
	}

	/**
	 * realize() may be called first.
	 * @return
	 */
	public int getInteractionTypeCount(){
		return size;
	}
	
	public int getInteractionType(int i){
		int index = getIndexInteractionTypeInArray(i);
		return getInteractionTypeFromByteArray(arrInteractionType[index], arrInteractionType[index+1], arrInteractionType[index+2]);
	}

	private int getIndexInteractionTypeInArray(int i) {
		return i * getNumBytesEntry();
	}

	private int incrementAtomTypeCount(int i){
		int index = getIndexInteractionTypeInArray(i);
		arrInteractionType[index+NUM_BYTES_INTERACTION_TYPE]++;
		return arrInteractionType[index+NUM_BYTES_INTERACTION_TYPE];
	}
	
	protected static int getInteractionTypeFromByteArray(byte low, byte med, byte high){
		
		// & 0xFF to prevent conservation of -1.
		int v = (high & 0xFF) << 16;

		v = v | (med & 0xFF) << 8;

		v = v | (low & 0xFF);
		
		return v;
	}

	InterActionTypeFreq getInteraction(int i){

		int index = getIndexInteractionTypeInArray(i);

		InterActionTypeFreq iaf =
				new InterActionTypeFreq(getInteractionType(i), arrInteractionType[index+NUM_BYTES_INTERACTION_TYPE]);

		return iaf;
	}
	
	public int getAtomicNo(int i){
		return InteractionAtomTypeCalculator.getAtomicNumber(getInteractionType(i));
	}
	
	public boolean isAromatic(int i){
		return InteractionAtomTypeCalculator.isAromatic(getInteractionType(i));
	}
	
	public static int getAtomicNoFromInteractionType(int interactionType){
		return InteractionAtomTypeCalculator.getAtomicNumber(interactionType);
	}

	public static PPNode getDummy(){
		PPNode node = new PPNode();
		node.add(DUMMY_INTERACT_ID);
		return node;
	}

	public boolean hasHeteroAtom(){
		return heteroAtom;
	}

	private void calcHasHeteroAtom(){
		heteroAtom = false;
		
		int size = getInteractionTypeCount();
		
		for (int i = 0; i < size; i++) {
			if(getAtomicNo(i) !=6) {
				heteroAtom=true;
				break;
			}
		}
	}

	public boolean isCarbonExclusiveNode(){
		boolean carbon = true;
		int size = getInteractionTypeCount();
		for (int i = 0; i < size; i++) {
			if(getAtomicNo(i) !=6) {
				carbon=false;
				break;
			}
		}
		return carbon;
	}
	public boolean containsHetero(){
		return heteroAtom;
	}

	
	private void init(){
		arrInteractionType = new byte [BUFFER_SIZE];
		modeFlexophore = ConstantsFlexophore.MODE_SOFT_PPPOINTS;
	}
	
	/**
	 * Flat copy from node into this.
	 * @param arrInteractionType
	 * @param size
	 */
	public void set(byte [] arrInteractionType, byte size){
		
		this.arrInteractionType = arrInteractionType;
		
		this.size = size;
	}

	public byte getModeFlexophore() {
		return modeFlexophore;
	}

	public void setModeFlexophore(byte modeFlexophore) {
		this.modeFlexophore = modeFlexophore;
	}



	private void initInteractionType(int interactionType){

		if(interactionType > ConstantsFlexophoreGenerator.MAX_VAL_INTERACTION_TYPE){
			throw new RuntimeException("Interaction type to large for PPNode!");
		}

		int index = getIndexInteractionTypeInArray(size);

		setInterActionType(interactionType, index, (byte)1);

		size++;

		int l = size * getNumBytesEntry();

		if(l == arrInteractionType.length){
			arrInteractionType = ArrayUtilsCalc.resize(arrInteractionType, arrInteractionType.length+BUFFER_SIZE);
		}
	}

	private void setInterActionType(int interactionType, int indexInArray, byte frequency){

		byte low = (byte)(interactionType & MASK1);

		byte med = (byte)(interactionType >> 8);

		byte high = (byte)(interactionType >> 16);

		arrInteractionType[indexInArray] = low;

		arrInteractionType[indexInArray+1] = med;

		arrInteractionType[indexInArray+2] = high;

		arrInteractionType[indexInArray+3] = frequency;

	}
	
	public void sortInteractionTypes(){

		List<InterActionTypeFreq> li = new ArrayList<>();

		for (byte i = 0; i < size; i++) {
			li.add(getInteraction(i));
		}

		Comparator<InterActionTypeFreq> c = new Comparator<InterActionTypeFreq>() {
			@Override
			public int compare(InterActionTypeFreq o1, InterActionTypeFreq o2) {

				int cmp = 0;

				if(o1.interactionType > o2.interactionType){
					cmp=1;
				}else if(o1.interactionType < o2.interactionType){
					cmp=-1;
				}

				return cmp;
			}
		};

		Collections.sort(li, c);

		for (int i = 0; i < li.size(); i++) {
			InterActionTypeFreq iaf = li.get(i);

			int index = getIndexInteractionTypeInArray(i);

			setInterActionType(iaf.interactionType, index, iaf.frequency);
		}
	}

	public double getFractionCarbonInteractions(){

		int carbons = 0;

		double sumFreq = 0;
		for (int i = 0; i < getInteractionTypeCount(); i++) {

			InterActionTypeFreq iaf = getInteraction(i);

			int interactionTypeQuery = iaf.interactionType;

			if(InteractionAtomTypeCalculator.isCarbonInteraction(interactionTypeQuery)){
				carbons+=iaf.frequency;
			}

			sumFreq += iaf.frequency;
		}

		double fractionCarbonQuery = carbons / sumFreq;

		return fractionCarbonQuery;

	}
		
	public String toString(){
		StringBuilder sb = new StringBuilder();
		
		sb.append("(");
		
		int size = getInteractionTypeCount();
		
		for (int i = 0; i < size;i++) {
			
			if(i>0){
				sb.append(SEPARATOR_ATOMS);
			}

			InterActionTypeFreq iaf = getInteraction(i);

			sb.append(iaf.interactionType);

			if(iaf.frequency>1){
				sb.append(MULT_FREQ);
				sb.append( iaf.frequency);
			}
		}
		
		sb.append(")");
		return sb.toString();
	}

	public static PPNode getHeteroOnlyNode(PPNode node){

		PPNode nodeHetero = new PPNode();

		int nQuery = node.getInteractionTypeCount();
		for (int i = 0; i < nQuery; i++) {
			int interactionType = node.getInteractionType(i);
			if(InteractionAtomTypeCalculator.getAtomicNumber(interactionType)!=6){
				nodeHetero.add(interactionType);
			}
		}
		nodeHetero.realize();

		return nodeHetero;
	}

	/**
	 * Prints the atom symbols for the types
	 * @return
	 */
	public String toStringText(){

		StringBuilder sb = new StringBuilder();
		sb.append("(");
		for (int i = 0; i < size; i++) {
			if(i>0){
				sb.append(SEPARATOR_ATOMS);
			}
			InterActionTypeFreq iaf = getInteraction(i);
			String s = InteractionAtomTypeCalculator.getString(iaf.interactionType);
			sb.append(s);
			if(iaf.frequency>1){
				sb.append(MULT_FREQ);
				sb.append( iaf.frequency);
			}
		}
		sb.append(")");
		return sb.toString();
	}

	/**
	 * Prints the identifier and the atom symbols for the types. Separated by ':'.
	 * @return
	 */
	public String toStringElusive(){

		StringBuilder sb = new StringBuilder();

		sb.append("(");

		for (int i = 0; i < size; i++) {
			if(i>0){
				sb.append(SEPARATOR_ATOMS);
			}
			InterActionTypeFreq iaf = getInteraction(i);
			String s = InteractionAtomTypeCalculator.getString(iaf.interactionType);
			sb.append(iaf.interactionType);
			sb.append(":");
			sb.append(s);
			if(iaf.frequency>1){
				sb.append(MULT_FREQ);
				sb.append( iaf.frequency);
			}
		}

		sb.append(")");
		return sb.toString();
	}

	public String toStringLongHardPPPoint(){

		StringBuilder sb = new StringBuilder();

		sb.append("(");

		for (int i = 0; i < size; i++) {

			if(i>0){
				sb.append(SEPARATOR_ATOMS);
			}

			InterActionTypeFreq iaf = getInteraction(i);

			String s = ConstantsFlexophoreHardPPPoints.toStringPPPoints(iaf.interactionType);

			sb.append(s);

			if(iaf.frequency>1){
				sb.append(MULT_FREQ);
				sb.append(iaf.frequency);
			}
		}

		sb.append(")");
		return sb.toString();
	}

	/**
	 *
	 * @param strNode i.e. 262,392,4358*2,8582,590088,598407
	 * @return the node with the atom types. If an empty string is given a node without atom types is returned.
	 */
	public static PPNode read(String strNode) {

		if(strNode.length()==0){
			return new PPNode();
		}

		String [] arr = strNode.split(SEPARATOR_ATOMS);

		PPNode n = new PPNode();
		for (int i = 0; i < arr.length; i++) {
			String strAtomType = arr[i];
			if(strAtomType.contains(MULT_FREQ)){
				String [] arrAtType = strAtomType.split("\\"+MULT_FREQ);
				int type = Integer.parseInt(arrAtType[0]);
				int freq = Integer.parseInt(arrAtType[1]);
				for (int j = 0; j < freq; j++) {
					n.add(type);
				}
			} else {
				int type = Integer.parseInt(strAtomType);
				n.add(type);
			}
		}

		n.realize();

		return n;
	}

	private static int getMaximumInteractionType(PPNode n){
		
		int max = 0;
		
		int size = n.getInteractionTypeCount();
		
		for (int i = 0; i < size ; i++) {
			
			int id = n.getInteractionType(i);
			
			if(id < max){
				max=id;
			}
			
		}

		return max;
	}

	public static int getNumBytesEntry(){
		return NUM_BYTES_INTERACTION_TYPE + BYTES_FREQUENCY_INTERACTION;
	}


	private static class InterActionTypeFreq {

		int interactionType;
		byte frequency;

		public InterActionTypeFreq(int interactionType, byte frequency) {
			this.interactionType = interactionType;
			this.frequency = frequency;
		}

		@Override
		public String toString() {

			return "InterActionTypeFreq{" +
					"interactionType=" + interactionType +
					", frequency=" + frequency +
					'}';
		}
	}
	
	
}

