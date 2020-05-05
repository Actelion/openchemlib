package com.actelion.research.util.graph.complete;

import java.util.Arrays;

import com.actelion.research.util.BurtleHasher;

/**
 * 
 * 
 * SolutionCompleteGraph
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Oct 1, 2012 MvK: Start implementation
 */
public class SolutionCompleteGraph extends AMemorizedObject implements Comparable<SolutionCompleteGraph>{
	
	private byte [] heapIndexBase;
	
	private int sizeHeap;
	
	private byte [] heapIndexQuery;

	private byte maxIndexNodeQuery;
	
	/**
	 * The index is the index of the node in the query molecule.
	 * 
	 * The value at 'index' is the index of the node in the base molecule.
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
		
		if(nodes>-1) {
			
			for (int i = 0; i < nodes; i++) {
				sb.append(arrSolution[i]);
				if(i < nodes -1)
					sb.append(" ");
			}
			
		} else {
			sb.append(Arrays.toString(arrSolution));	
		}
		
		return sb.toString();
	}
}