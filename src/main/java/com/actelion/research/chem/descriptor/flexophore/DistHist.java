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
import com.actelion.research.calc.filter.SlidingWindow;
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
