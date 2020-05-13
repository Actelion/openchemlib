/*
 * Copyright (c) 2020.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
 */

package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;

import java.util.*;

public class PPNode implements Comparable<PPNode> {

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
		
		boolean bY = true;
		
		int size = getInteractionTypeCount();
		
		for (int i = 0; i < size; i++) {
			if(getAtomicNo(i) !=6) {
				bY=false;
				break;
			}
		}
		
		return bY;
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

	public String toStringLong(){

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
				sb.append( iaf.frequency);
			}
		}

		sb.append(")");
		return sb.toString();
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

