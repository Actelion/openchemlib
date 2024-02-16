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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class DistHist implements Serializable {
	
	private static final long serialVersionUID = 17042013;
	
	// Array containing the histogram vectors.
	// Only one array, because of memory reasons.
	protected byte [] arrDistHists;
	
	private int numPPNodes;
	
	private int identifier;
	
	public DistHist() {
		arrDistHists = new byte[0];
		numPPNodes = 0;
		identifier=-1;
	}

	/**
	 * Initializes the distance histograms. Filled with 0.
	 * @param nPPNodes
	 */
	public DistHist(int nPPNodes) {
		initHistogramArray(nPPNodes);
	}

	public DistHist(DistHist distHist) {
		distHist.copy(this);
	}
	/**
	 * 
	 * @param copy: This object is written into copy.
	 */
	public void copy(DistHist copy){
		
		copy.identifier = this.identifier;
		
		if(copy.numPPNodes != this.numPPNodes) {
		
			copy.numPPNodes = this.numPPNodes;
			
			copy.arrDistHists = new byte[this.arrDistHists.length];
		}
		
		System.arraycopy(this.arrDistHists, 0, copy.arrDistHists, 0, this.arrDistHists.length);
	}

	public int getBonds() {
		return ((numPPNodes * numPPNodes)-numPPNodes) / 2;
	}
	
	/**
	 * 
	 * @return number of pharmacophore points.
	 */
	public int getNumPPNodes(){
		return numPPNodes;
	}
	
	public int getSizeBytes(){
		int s = 0;
		
		if(getNumPPNodes()>1) {
			s += arrDistHists.length;
		}
		// size, identifier
		s += (Integer.SIZE / 8)*2;
		
		return s;
	}

	protected void initHistogramArray(int nPPNodes) {

		this.numPPNodes = nPPNodes;

		int nBonds = ((numPPNodes * numPPNodes)-numPPNodes) / 2;
		
		arrDistHists = new byte[nBonds* ConstantsFlexophoreGenerator.BINS_HISTOGRAM];
		
	}
	/**
	 * 
	 * @param indexNode1
	 * @param indexNode2
	 * @param size number of nodes
	 * @return
	 */
	public static final int getIndex(int indexNode1, int indexNode2, int size){
		
		int i1 = Math.min(indexNode1,indexNode2);
		
		int i2 = Math.max(indexNode1,indexNode2);
		
		int h = i1+1;
				
		int index = (i1*size+i2) - ((h*h-h)/2) - h;
		
		return index;
	}
	
	public int getMinDist(int indexAt1, int indexAt2) {
		byte [] arr = getDistHist(indexAt1, indexAt2);
		
		int dist=0;
		
		for (int i = 0; i < arr.length; i++) {
			if(arr[i]>0){
				dist = i; 
				break;
			}
		}
		
		return dist;
	}
	
	/**
	 * Calculates for nodes cluster the node with the minimum rmsd to the cluster center.
	 * @param maxDistance maximum distance that a node will become a cluster member
	 * @return List with cluster.
	 */
	public List<ClusterNode> getClusterCenter(int maxDistance){
		
		List<ClusterNode> liCluster = new ArrayList<ClusterNode>();
		
		// Get all cluster
		for (int i = 0; i < getNumPPNodes(); i++) {
			ClusterNode cluster = new ClusterNode(i);
			
			for (int j = 0; j < getNumPPNodes(); j++) {
				if(i!=j){
					byte [] hist = getDistHist(i, j);
					boolean bInCluster=false;
					for (int k = 0; k < maxDistance+1; k++) {
						if(hist[k]>0){
							bInCluster=true;
							break;
						}
					}
					if(bInCluster){
						cluster.add(j);
					}
					
				}
			}
			if(cluster.isCluster())
				liCluster.add(cluster);
		}
		
		// Calculate RMSD
		for (int i = 0; i < liCluster.size(); i++) {
			ClusterNode cl = liCluster.get(i);
			
			int rmsd = ClusterNode.getRMSD(cl, this);
			
			cl.setRMSD(rmsd);
			
		}
		
		// Equal cluster into one list
		// The cluster with the lowest rmsd is the chosen cluster center.
		List<ClusterNode> liClusterNodeMinRMSD = new ArrayList<ClusterNode>();
		for (int i = 0; i < liCluster.size()-1; i++) {
			ClusterNode cl = liCluster.get(i);
			List<ClusterNode> liEqCluster = new ArrayList<ClusterNode>();
			liEqCluster.add(cl);
			
			for (int j = liCluster.size()-1; j > i; j--) {
				if(cl.equals(liCluster.get(j))){
					liEqCluster.add(liCluster.remove(j));
				}
			}
			
			Collections.sort(liEqCluster);
			liClusterNodeMinRMSD.add(liEqCluster.get(0));
			
		}
		
		
		return liClusterNodeMinRMSD;
	}

	/**
	 * The distance histograms are stored in a single array.
	 * @param indexAt1
	 * @param indexAt2
	 * @param arrHist a deep copy is taken.
	 */
	public void setDistHist(int indexAt1, int indexAt2, byte [] arrHist) {
		
		if(indexAt1 >= numPPNodes) {
			throw new ArrayIndexOutOfBoundsException(indexAt1);
		} else if(indexAt2 >= numPPNodes)
			throw new ArrayIndexOutOfBoundsException(indexAt2);
		
		int index = getIndex(indexAt1, indexAt2, numPPNodes);
		
		int posStart = index * ConstantsFlexophoreGenerator.BINS_HISTOGRAM;
		
		System.arraycopy(arrHist, 0, arrDistHists, posStart, ConstantsFlexophoreGenerator.BINS_HISTOGRAM);
		
	}
	public byte [] getDistHists() {
		return arrDistHists;
	}
	
	/**
	 * Flat copy 
	 * @param arr
	 * @param size
	 * @param identifier
	 */
	public void setDistHists(byte [] arr, int size, int identifier) {
		arrDistHists=arr;
		this.numPPNodes = size;
		this.identifier = identifier;
	}
	
	/**
	 * 
	 * @param indexAt1
	 * @param indexAt2
	 * @return deep copy.
	 */
	public byte [] getDistHist(int indexAt1, int indexAt2) {
		
		int index = getIndex(indexAt1, indexAt2, numPPNodes);
		
		int posStart = index * ConstantsFlexophoreGenerator.BINS_HISTOGRAM;
		
		byte [] arr = new byte[ConstantsFlexophoreGenerator.BINS_HISTOGRAM];
		
		System.arraycopy(arrDistHists, posStart, arr, 0, ConstantsFlexophoreGenerator.BINS_HISTOGRAM);
		
		return arr;
	}

	public int getIndexPosStartForDistHist(int indexAt1, int indexAt2) {

		int index = getIndex(indexAt1, indexAt2, numPPNodes);

		int posStart = index * ConstantsFlexophoreGenerator.BINS_HISTOGRAM;

		return posStart;
	}

	public byte getValueAtAbsolutePosition(int indexAbsolutePosition){
		return arrDistHists[indexAbsolutePosition];
	}


	public byte [] getDistHist(int indexAt1, int indexAt2, byte [] arr) {
		
		int index = getIndex(indexAt1, indexAt2, numPPNodes);
		
		int posStart = index * ConstantsFlexophoreGenerator.BINS_HISTOGRAM;
		
		System.arraycopy(arrDistHists, posStart, arr, 0, ConstantsFlexophoreGenerator.BINS_HISTOGRAM);
		
		return arr;
	}

	/**
	 * To get the distance in Angstrom the relative distance has to be multiplied with the histogram range <code>CGMult.RANGE_HISTOGRAM</code>.
	 * @param indexAt1
	 * @param indexAt2
	 * @return relative distance 0: only the first bin is occupied, 
	 * 1: the last bin is occupied.
	 */
	public double getRelMaxDistInHist(int indexAt1, int indexAt2) {
		byte [] arr = getDistHist(indexAt1, indexAt2);

		double max = Integer.MIN_VALUE;
		for (int i = arr.length-1; i >= 0; i--) {
			if(arr[i]>0){
				max = i;
				break;
			}
		}
		return max / arr.length;
	}
	
	public double getMaxDistInHist(int indexAt1, int indexAt2) {
		byte [] arr = getDistHist(indexAt1, indexAt2);

		double max = Integer.MIN_VALUE;
		for (int i = arr.length-1; i >= 0; i--) {
			if(arr[i]>0){
				max = i;
				break;
			}
		}
		return max;
	}
	
	/**
	 * To get the distance in Angstrom the relative distance has to be multiplied with the histogram range <code>CGMult.RANGE_HISTOGRAM</code>.
	 * Maximum relative distance between two nodes in the object.
	 * @return
	 */
	public double getRelMaxDistInHist() {
		int size = getNumPPNodes();
		
		double max = 0;
		for (int i = 0; i < size; i++) {
			for (int j = i+1; j < size; j++) {
				double dist = getRelMaxDistInHist(i, j);
				if(dist>max)
					max=dist;
			}
		}
		
		return max;
	}
	
	public String toStringHistsIndexed(){
		
		StringBuffer b = new StringBuffer();
		
		for (int i = 0; i < numPPNodes; i++) {
			for (int j = i+1; j < numPPNodes; j++) {
				byte [] arrHist = getDistHist(i,j);
				
				if(arrHist!=null)
					b.append(i + "\t" + j + "\t[" + ArrayUtilsCalc.toString(arrHist) + "]\n");
			}
		}
		
		return b.toString();
	}
	
	public int getIdentifier() {
		return identifier;
	}
	
	public void setIdentifier(int identifier) {
		this.identifier = identifier;
	}
	

}
