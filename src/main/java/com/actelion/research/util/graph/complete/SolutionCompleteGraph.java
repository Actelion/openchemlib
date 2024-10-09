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
 */

package com.actelion.research.util.graph.complete;

import com.actelion.research.util.BurtleHasher;

import java.util.Arrays;

/**
 * SolutionCompleteGraph
 * @author Modest von Korff
 * @version 1.0
 * Oct 1, 2012 MvK: Start implementation
 */
public class SolutionCompleteGraph extends AMemorizedObject implements Comparable<SolutionCompleteGraph>{

	//
	// Stores the solution as pairwise node indices
	//

	// Contains the node index in MolDistHist for base
	private byte [] heapIndexBase;

	/**
	 * number of the matching nodes
	 */
	private int sizeHeap;

	// Contains the node index in MolDistHist for query
	private byte [] heapIndexQuery;

	private byte maxIndexNodeQuery;
	
	/**
	 * The index is the index of the node in the query molecule. Not matched query nodes contain a -1.
	 * The value at 'index' is the index of the node in the base molecule.
	 * Contains the same information as heapIndexBase and heapIndexQuery. Used for fast lookup.
	 */
	private byte [] arrSolution;

	private int hash;
	
	private double similarity;
	
	private int nodes;
	
	public SolutionCompleteGraph() {
		
		nodes = -1;
		
		arrSolution = new byte [CompleteGraphMatcher.MAX_NUM_NODES];
		
		heapIndexBase = new byte [CompleteGraphMatcher.MAX_NUM_NODES];
		
		heapIndexQuery = new byte [CompleteGraphMatcher.MAX_NUM_NODES];
		
		for (int i = 0; i < arrSolution.length; i++) {
			arrSolution[i]=CompleteGraphMatcher.DEFAULT_VAL;
			heapIndexBase[i]=CompleteGraphMatcher.DEFAULT_VAL;
			heapIndexQuery[i]=CompleteGraphMatcher.DEFAULT_VAL;
		}
	}
	
	/**
	 * Adds a pair of nodes to this solution.
	 * @param indexNodeBase index of the base node in the complete base graph.
	 * @param indexNodeQuery index of the query node in the complete query graph.
	 */
	public void add(byte indexNodeQuery, byte indexNodeBase) {
		arrSolution[indexNodeQuery]=indexNodeBase;
		if(indexNodeQuery>maxIndexNodeQuery){
			maxIndexNodeQuery=indexNodeQuery;
		}
		heapIndexBase[sizeHeap]=indexNodeBase;
		heapIndexQuery[sizeHeap]=indexNodeQuery;
		sizeHeap++;
		calcHashCode();
	}
	
	public int compareTo(SolutionCompleteGraph s) {
		if(similarity > s.similarity){
			return 1;
		} else if(similarity < s.similarity){
			return -1;
		}
		return 0;
	}

	public int getSizeHeap(){
		return sizeHeap;
	}

	public byte getIndexBaseFromHeap(int indexOnHeapBase){
		return heapIndexBase[indexOnHeapBase];
	}
	
	public byte getIndexQueryFromHeap(int indexOnHeapQuery){
		return heapIndexQuery[indexOnHeapQuery];
	}
	
	public byte getIndexCorrespondingBaseNode(int indexQueryNode) {
		return arrSolution[indexQueryNode];
	}

	/**
	 * The index is the index of the node in the query molecule.
	 * The value at 'index' is the index of the node in the base molecule.
	 * Can contain -1 if a node is not mapping.
	 */
	public byte [] getSolution (){
		return arrSolution;
	}

	public boolean equals(Object obj) {
		
		boolean eq = true;
		SolutionCompleteGraph s = (SolutionCompleteGraph)obj;
		if(sizeHeap != s.sizeHeap){
			return false;
		}
		if(maxIndexNodeQuery != s.maxIndexNodeQuery){
			return false;
		}
		for (int i = 0; i < maxIndexNodeQuery; i++) {
			if(arrSolution[i]!=s.arrSolution[i]){
				eq = false;
				break;
			}
		}
		
		return eq;
	} 
		
	public void copyIntoThis(AMemorizedObject a) {
		SolutionCompleteGraph solution = (SolutionCompleteGraph)a;
		System.arraycopy(solution.heapIndexBase, 0, heapIndexBase, 0, solution.sizeHeap);
		System.arraycopy(solution.heapIndexQuery, 0, heapIndexQuery, 0, solution.sizeHeap);
		System.arraycopy(solution.arrSolution, 0, arrSolution, 0, arrSolution.length);
		sizeHeap = solution.sizeHeap;
		hash = solution.hash;
		similarity = solution.similarity;
		nodes = solution.nodes;
		maxIndexNodeQuery = solution.maxIndexNodeQuery;
	}
	
    private void calcHashCode(){
    	hash = BurtleHasher.hashlittle(arrSolution, 13, maxIndexNodeQuery+1);
    }
    
	public int hashCode() {
		return hash;
	}
	
	public void reset() {
		for (int i = 0; i < arrSolution.length; i++) {
			arrSolution[i] = CompleteGraphMatcher.DEFAULT_VAL;
			heapIndexBase[i] = CompleteGraphMatcher.DEFAULT_VAL;
			heapIndexQuery[i] = CompleteGraphMatcher.DEFAULT_VAL;
		}
		sizeHeap = 0;
		hash = 0;
		similarity = 0;
		nodes = 0;
		maxIndexNodeQuery = 0;
	}

	public double getSimilarity() {
		return similarity;
	}

	public void setSimilarity(double similarity) {
		this.similarity = similarity;
	}

	/**
	 * Only needed for the toString() method.  
	 * @param nodes
	 */
	public void setNodesQuery(int nodes) {
		this.nodes = nodes;
	}
	
	public int getNodesQuery() {
		return nodes;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		if(sizeHeap>0) {
			
			for (int i = 0; i < maxIndexNodeQuery+1; i++) {
				if(sb.length()>0)
					sb.append(" ");
				sb.append(arrSolution[i]);
			}
			
		} else {
			sb.append(Arrays.toString(arrSolution));	
		}
		
		return sb.toString();
	}
}