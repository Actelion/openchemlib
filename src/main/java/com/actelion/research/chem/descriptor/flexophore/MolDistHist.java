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
import com.actelion.research.util.StringFunctions;

import java.awt.*;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class MolDistHist extends DistHist implements Serializable, IMolDistHist  {
	
	private static final long serialVersionUID = 5042020;


	public static final boolean VERBOSE = false;
	
	private static final int SIZE_BUFFER = 16 * getNumBytesEntry();
	
	// The pharmacophore nodes containing 
	// The first entry is the number of atom types in the node, followed by the atom types and the frequency.
	// [nAtomsNode1, byte1_attype1, byte2_attype1, byte3_freq_attype1, nAtomsNode2, byte1_attype2, byte2_attype2, byte3_freq_attype2...]
	
	private byte [] arrNode;
	
	private int posNode;
	
	private boolean finalized;

	private byte modeFlexophore;

	public MolDistHist () {
		initHistogramArray(0);
		init();
	}
	
	public MolDistHist (int nNodes) {
		initHistogramArray(nNodes);
		init();
	}
	
	public MolDistHist (MolDistHist mdh) {
		initHistogramArray(mdh.getNumPPNodes());
		init();
		mdh.copy(this);
	}
	
	public MolDistHist copy(){
		MolDistHist copy = new MolDistHist();
		copy(copy);
		return copy;
	}

	private void init(){
		modeFlexophore = ConstantsFlexophore.MODE_SOFT_PPPOINTS;
	}
	
	public boolean check(){
		boolean bOK = true;
		
		int nodes = getNumPPNodes();
		for (int i = 0; i < nodes; i++) {
			PPNode node = getNode(i);
			int ats = node.getInteractionTypeCount();
			for (int j = 0; j < ats; j++) {
				int inttype = node.getInteractionType(j);

				String s = InteractionAtomTypeCalculator.getString(inttype);

				if(s.length()==0) {
					bOK = false;
					if(VERBOSE)
						System.err.println("Node " + i + " atom " + j + " Interaction type " + inttype + ".");
				}
			}
		}
		return bOK;
	}

	/**
	 * 
	 * @param copy: This object is written into copy.
	 */
	public void copy(MolDistHist copy){
		super.copy(copy);
		
		copy.arrNode = new byte [arrNode.length];
		
		System.arraycopy(arrNode, 0, copy.arrNode, 0, arrNode.length);
		
		copy.posNode = posNode;

		copy.modeFlexophore = modeFlexophore;

		copy.realize();
		
	}
	
	public boolean equals(Object o) {
		boolean bEQ=true;
		
		MolDistHist mdh=null;
		try {
			mdh = (MolDistHist)o;
		} catch (RuntimeException e) {
			return false;
		}
		
		
		if(getNumPPNodes() != mdh.getNumPPNodes())
			return false;
		
		
		for (int i = 0; i < getNumPPNodes(); i++) {
			PPNode n1 = getNode(i);
			PPNode n2 = mdh.getNode(i);
			if(!n1.equals(n2)){
				bEQ = false;
				break;
			}
		}
		
		for (int i = 0; i < getNumPPNodes(); i++) {
			for (int j = i+1; j < getNumPPNodes(); j++) {
				byte [] a1 = getDistHist(i,j);
				byte [] a2 = mdh.getDistHist(i,j);
				for (int k = 0; k < a2.length; k++) {
					if(a1[k]!=a2[k]){
						bEQ = false;
						break;
					}
				}
			}
		}
		
		return bEQ;
	}

	public byte getModeFlexophore() {
		return modeFlexophore;
	}

	protected void initHistogramArray(int nNodes) {
		super.initHistogramArray(nNodes);
		
		arrNode = new byte[nNodes*(PPNode.getNumBytesEntry()+1)];
		
		finalized = false;
	}

	/**
	 *
	 * @param node
	 */
	public void addNode(PPNode node) {

		byte [] arr = node.get();

		byte nInteractionTypeCount = (byte)node.getInteractionTypeCount();

		int newLen = posNode + arr.length + 1;

		if(arrNode.length < newLen){
			int lenBuffer = SIZE_BUFFER;
			while(arrNode.length + lenBuffer < newLen) {
				lenBuffer *= 2;
			}

			resizeNodeArray(arrNode.length + lenBuffer);
		}
		
		arrNode[posNode++] = nInteractionTypeCount;

		System.arraycopy(arr, 0, arrNode, posNode, arr.length);

		posNode += arr.length;

		finalized = false;
	}

	/**
	 * @return the arrNode
	 */
	protected byte[] getArrNode() {
		return arrNode;
	}

	/**
	 * @param arrNode the arrNode to set
	 */
	protected void setArrNode(byte[] arrNode) {
		
		this.arrNode = arrNode;
		
		int pos=0;
		while(arrNode[pos] > 0){
			pos += arrNode[pos] * PPNode.NUM_BYTES_INTERACTION_TYPE + 1;

			if(pos>=arrNode.length)
				break;
		}
		
	}

	/**
	 * Resizes the node array to the needed length.
	 */
	public void realize(){

		int size = getNumPPNodes();
		
		if(size==0){
			throw new RuntimeException("No pharmacophore points in Flexophore.");
		}
		
		int pos = getPositionNode(size-1);
		
		int len = pos + arrNode[pos] * PPNode.getNumBytesEntry() + 1;

		resizeNodeArray(len);
		
		if(getNumPPNodes()==0)
			return;

		finalized = true;
	}
	
	private void resizeNodeArray(int newsize){

		byte [] arr = new byte [newsize];
		
		System.arraycopy(arrNode, 0, arr, 0, Math.min(arrNode.length, newsize));
		
		arrNode = arr;
		
		finalized = false;
	}
	
	
	public int getConnAtom(int at, int index) {
		if(index >= at)
			index++;
		
		return index;
	}

	@Override
	public int hashCode() {
		String s = toString();
		s = s.replace(" ", "");
		return s.hashCode();
	}
	
	public String toString(){

//		Causes errors in debug mode
//		if(!finalized)
//			realize();
		
		StringBuilder sb = new StringBuilder(toStringNodes());

		sb.append(toStringHists());
		
		return sb.toString();
	}

	public String toStringNodes(){

		StringBuilder sb = new StringBuilder();

		sb.append("[");
		for (int i = 0; i < getNumPPNodes(); i++) {
			sb.append(getNode(i).toString());
			if(i<getNumPPNodes()-1){
				sb.append(" ");
			} else {
				sb.append("]");
			}
		}

		return sb.toString();
	}

	public String toStringHists(){

		StringBuilder sb = new StringBuilder();
		
		for (int i = 0; i < getNumPPNodes(); i++) {
			for (int j = i+1; j < getNumPPNodes(); j++) {
				byte [] arrHist = getDistHist(i,j);
				
				if(arrHist!=null)
					sb.append(ArrayUtilsCalc.toString(arrHist));
			}
		}
		
		return sb.toString();
	}

	
	protected boolean isFinalized() {
		return finalized;
	}
	
	
	
	/**
	 * The position of the size of the node. The node index starts with 0 and runs up to size-1 inclusive.
	 * @param index
	 * @return
	 */
	private int getPositionNode(int index){
		int pos = 0;
		for (int i = 0; i < index; i++) {
			pos += arrNode[pos]  * PPNode.getNumBytesEntry() + 1;
		}
		return pos;
	}
	
	/**
	 * !Slow method, it has to iterate through a loop to find the node in the array!
	 * @param index
	 * @return deep copy of the node.
	 *
	 */
	public PPNode getNode(int index){
		PPNode node = new PPNode();
		
		int pos = getPositionNode(index);
		
		int len = arrNode[pos] * PPNode.getNumBytesEntry();

		byte [] arr = new byte[len];

		System.arraycopy(arrNode, pos+1, arr, 0, len);

		node.set(arr, arrNode[pos]);

		node.realize();
		
		return node;
	}
	
	/**
	 * 
	 * @param index
	 * @return number of pharmacophore points at the specified index
	 * slow method because it calls getPositionNode(index).
	 */
	public int getPPPoints(int index){
		
		int pos = getPositionNode(index);
		
		int nPPPoints = arrNode[pos];
		
		return nPPPoints;
	}
	
	public int getSizeBytes(){
		int s = super.getSizeBytes();
		
		s += arrNode.length;	
		
		// bFinalized
		s += (Integer.SIZE / 8)*(1);
		
		return s;
	}

	/**
	 * Only for interface compliance needed.
	 */
	@Override
	public int getNumInevitablePharmacophorePoints() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public boolean isInevitablePharmacophorePoint(int indexNode) {
		// TODO Auto-generated method stub
		return false;
	}

	public static int getNumBytesEntry(){
		return PPNode.getNumBytesEntry()+1;
	}

	/**
	 * reads a MolDistHist from the toString() method.
	 * @param strMolDistHist
	 * @return
	 */
	public static List<PPNode> readNodes(String strMolDistHist){

		int start = strMolDistHist.indexOf('(');
		
		boolean nodesProcessed = false;
		
		List<PPNode> liPPNode = new ArrayList<>();
		while(!nodesProcessed){
			
			int end = StringFunctions.nextClosing(strMolDistHist, start, '(', ')');
			
			String strNode = strMolDistHist.substring(start+1, end);

			PPNode n = PPNode.read(strNode);

			liPPNode.add(n);
			
			start = strMolDistHist.indexOf('(', end);
			
			if(start==-1){
				nodesProcessed = true;
			}
		}

		return liPPNode;
	}

	public static MolDistHist read(String strMolDistHist){

		String pattern = "[0-9]+";

		List<PPNode> liPPNode = readNodes( strMolDistHist);

		int size = liPPNode.size();

		MolDistHist mdh = new MolDistHist(size);

		for (PPNode ppNode : liPPNode) {
			mdh.addNode(ppNode);
		}

		boolean histsProcessed = false;

		List<byte []> liHist = new ArrayList<byte []>();

		int startHist = strMolDistHist.indexOf("][");

		int nHistograms = ((size*size)-size)/2;

		while(!histsProcessed){

			int endHist = StringFunctions.nextClosing(strMolDistHist, startHist, '[', ']');

			String sub = strMolDistHist.substring(startHist, endHist);

			List<Point> li = StringFunctions.match(sub, pattern);

			if(li.size() != ConstantsFlexophoreGenerator.BINS_HISTOGRAM){
				throw new RuntimeException("Error in histogram.");
			}

			byte [] arr = new byte [ConstantsFlexophoreGenerator.BINS_HISTOGRAM];

			int cc=0;
			for (Point p : li) {

				String strCount = sub.substring(p.x, p.y);

				arr[cc++] = (byte)(Integer.parseInt(strCount) & 0xFF);

			}

			liHist.add(arr);

			startHist = strMolDistHist.indexOf('[', endHist);

			if(liHist.size()==nHistograms){
				histsProcessed=true;
			}
		}

		int cc=0;
		for (int i = 0; i < size; i++) {
			for (int j = i+1; j < size; j++) {
				mdh.setDistHist(i,j, liHist.get(cc++));
			}
		}


		return mdh;
	}

}
