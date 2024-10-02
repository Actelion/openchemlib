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

	private int[][] nodeAtoms;  // original atom index list for every node

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
		if(!equalNodes(o)){
			return false;
		}
		MolDistHist mdh=null;
		try {
			mdh = (MolDistHist)o;
		} catch (RuntimeException e) {
			return false;
		}
		boolean bEQ=true;
		
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
	public boolean equalNodes(Object o) {
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

	/**
	 * @return the original atom indexes in node order, provided they have beed added when creating this MolDistHist
	 */
	public int[][] getNodeAtoms() {
		return nodeAtoms;
	}

	/**
	 * Adds the original atom indexes in node order to this MolDistHist
	 * @param nodeAtoms
	 */
	public void setNodeAtoms(int[][] nodeAtoms) {
		this.nodeAtoms = nodeAtoms;
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
	public String toStringNodesElusive(){

		StringBuilder sb = new StringBuilder();

		sb.append("[");
		for (int i = 0; i < getNumPPNodes(); i++) {
			sb.append(getNode(i).toStringElusive());
			if(i<getNumPPNodes()-1){
				sb.append(" ");
			} else {
				sb.append("]");
			}
		}

		return sb.toString();
	}
	public String toStringNodesLineWise(){

		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < getNumPPNodes(); i++) {
			sb.append(getNode(i).toString());
			if(i<getNumPPNodes()-1){
				sb.append("\n");
			}
		}

		return sb.toString();
	}
	public String toStringNodesElusiveLineWise(){

		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < getNumPPNodes(); i++) {
			sb.append(getNode(i).toStringElusive());
			if(i<getNumPPNodes()-1){
				sb.append("\n");
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

	public boolean containsFilledDistHist(){

		if(arrDistHists == null || arrDistHists.length==0){
			return false;
		}

		boolean filled=false;
		for (byte v : arrDistHists) {
			if(v>0){
				filled=true;
				break;
			}
		}

		return filled;
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
	public int getNumMandatoryPharmacophorePoints() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public boolean isMandatoryPharmacophorePoint(int indexNode) {
		return false;
	}
	@Override
	public double getWeightPharmacophorePoint(int indexNode) {
		return 1.0;
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
			if(end==-1){
				throw new RuntimeException("Error for MolDistHist " + strMolDistHist);
			}
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
	public static MolDistHist readNodes2MDH(String strMolDistHist){

		List<PPNode> liPPNode = readNodes(strMolDistHist);
		int size = liPPNode.size();
		MolDistHist mdh = new MolDistHist(size);
		for (PPNode ppNode : liPPNode) {
			mdh.addNode(ppNode);
		}
		return  mdh;
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

		if(nHistograms==0){
			histsProcessed=true;
		}

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
