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

import com.actelion.research.chem.descriptor.flexophore.completegraphmatcher.ObjectiveBlurFlexophoreHardMatchUncovered;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

/**
 * 
 * 
 * CompleteGraphMatcher
 * @author Modest von Korff
 * @version 1.0
 * Sep 26, 2012 MvK: Start implementation
 */
public class CompleteGraphMatcher<T extends ICompleteGraph> {

	public static final boolean DEBUG = false;

	// We need at least two pharmacophore points.
	public static final double TINY_SIM = 0.000001;
	public static final int MIN_NUM_NODES_SIM = 2;

	public static final int MAX_NUM_NODES = 127;
	
	private static final int INIT_CAPACITY_MEMORY = 1000;
	
	private static final int SIZE_LIST_SOLUTION = 1000;
	
	public static final int MAX_NUM_SOLUTIONS = 10000;
	
	public static final byte DEFAULT_VAL = -1;
	
	private IObjectiveCompleteGraph<T> objectiveCompleteGraph;

	private int maxNumSolutions;
	
	private int nodesBase;
	
	private int nodesQuery;
	
	private ContainerMemory<SolutionCompleteGraph> cm;
	
	private List<List<SolutionCompleteGraph>> liliSolution;
	
	private HashSet<SolutionCompleteGraph> hsSolution;
	
	private byte [] arrIndexBaseTmp;
	private byte [] arrIndexQueryTmp;
	
	private SolutionCompleteGraph solutionBest;
	private List<SolutionCompleteGraph> solutionsTop;

	private long validSolutions;
	
	private long createdSolutions;

	private boolean nodeSimilarityWithoutSizeDifference;
	
	public CompleteGraphMatcher(IObjectiveCompleteGraph<T> objectiveCompleteGraph) {
		this.objectiveCompleteGraph = objectiveCompleteGraph;
		init();
	}

	/**
	 * 
	 * @param objectiveCompleteGraph
	 */
	public void setObjective(IObjectiveCompleteGraph<T> objectiveCompleteGraph){
		this.objectiveCompleteGraph = objectiveCompleteGraph;
	}

	public void setVerbose(boolean verbose){
		objectiveCompleteGraph.setVerbose(verbose);
	}

	public IObjectiveCompleteGraph<T> getObjectiveCompleteGraph() {
		return objectiveCompleteGraph;
	}

	private void init(){
		
		maxNumSolutions = MAX_NUM_SOLUTIONS;
		
		cm = new ContainerMemory<>(INIT_CAPACITY_MEMORY, new FactorySolution());
		
		liliSolution = new ArrayList<>(MAX_NUM_NODES);
		
		liliSolution.add(new ArrayList<>());
		for (int i = 1; i < MAX_NUM_NODES; i++) {
			liliSolution.add(new ArrayList<>(SIZE_LIST_SOLUTION));
		}
		
		hsSolution = new HashSet<>(SIZE_LIST_SOLUTION);
		
		arrIndexBaseTmp = new byte [MAX_NUM_NODES];
		arrIndexQueryTmp = new byte [MAX_NUM_NODES];
		
		solutionBest = new SolutionCompleteGraph();
		solutionsTop = new ArrayList<>();
		nodeSimilarityWithoutSizeDifference = false;
	}


	public void setNodeSimilarityWithoutSizeDifference(boolean nodeSimilarityWithoutSizeDifference) {
		this.nodeSimilarityWithoutSizeDifference = nodeSimilarityWithoutSizeDifference;
	}

	public void set(T cgBase, T cgQuery){
		objectiveCompleteGraph.setBase(cgBase);
		objectiveCompleteGraph.setQuery(cgQuery);
	}
	
	
	private void initSearch(){
		cm.reset();
		nodesBase = objectiveCompleteGraph.getBase().getNumPPNodes();
		nodesQuery = objectiveCompleteGraph.getQuery().getNumPPNodes();
		for (List<SolutionCompleteGraph> liSolution : liliSolution) {
			liSolution.clear();
		}
		hsSolution.clear();
		int nodesInSolution = 1;

		// Creates indices for all one node mappings from base and query.
		// This is the start of the mapping process for the complete graph of the Flexophore.
		// It is tested if the nodes are similar for the given threshold.

		List<SolutionCompleteGraph> liSolution = liliSolution.get(nodesInSolution);
		for (byte indexNodeQuery = 0; indexNodeQuery < nodesQuery; indexNodeQuery++) {
			for (byte indexNodeBase = 0; indexNodeBase < nodesBase; indexNodeBase++) {
				if(objectiveCompleteGraph.areNodesMapping(indexNodeQuery, indexNodeBase)){
					// System.out.println(indexNodeQuery + "\t" + indexNodeBase);
					SolutionCompleteGraph solution = cm.get();
					solution.setNodesQuery(nodesQuery);
					solution.add(indexNodeQuery, indexNodeBase);
					liSolution.add(solution);
				}
			}
		}

		solutionBest = new SolutionCompleteGraph();
		solutionsTop.clear();

		// System.out.println(liSolution.size());
	}
	
	public double calculateSimilarity () {
		initSearch();
		if(nodesBase==1 && nodesQuery==1) {
			double sim = objectiveCompleteGraph.getSimilarityNodes(0,0);
			return sim;
		}

		if(objectiveCompleteGraph.isModeFragment() && nodesQuery==1) {

//			double simMax=0;
//			for (byte indexNodeBase = 0; indexNodeBase < nodesBase; indexNodeBase++) {
//				double sim = objectiveCompleteGraph.getSimilarityNodes(0,indexNodeBase);
//				if(sim>simMax){
//					simMax=sim;
//				}
//			}

			List<SolutionCompleteGraph> liSolution = liliSolution.get(1);
			double simMax=0;
			for (SolutionCompleteGraph solutionCompleteGraph : liSolution) {
				double sim = objectiveCompleteGraph.getSimilarity(solutionCompleteGraph);
				solutionCompleteGraph.setSimilarity(sim);
				if(sim>simMax){
					simMax=sim;
					solutionBest=solutionCompleteGraph;
				}
			}

			solutionsTop.clear();
			double simMaxMargin=simMax-TINY_SIM;
			for (SolutionCompleteGraph scg : liSolution) {
				if(scg.getSimilarity()>simMaxMargin){
					solutionsTop.add(scg);
				}
			}

			return simMax;
		}


		int maxNumNodesWithSolution = 0;
		for (int nodesInSolution = 1; nodesInSolution < nodesBase+1; nodesInSolution++) {
			List<SolutionCompleteGraph> liSolution = liliSolution.get(nodesInSolution);
			boolean validSolutionFound = false;
			hsSolution.clear();
			for (SolutionCompleteGraph solution : liSolution) {
				if(getNeighbourSolutions(solution)){
					validSolutionFound = true;
				}
				if(hsSolution.size()> maxNumSolutions){
					break;
				}
			}
			
			if(validSolutionFound){
				maxNumNodesWithSolution = nodesInSolution+1;
				liliSolution.get(maxNumNodesWithSolution).addAll(hsSolution);
				
				// System.out.println("Found " + hsSolution.size() + " valid solutions for " + maxNumNodesWithSolution + " nodes.");
								
				//
				// Remove all smaller solutions
				//
				List<SolutionCompleteGraph> li = liliSolution.get(nodesInSolution);
				for (SolutionCompleteGraph solutionCompleteGraph : li) {
					cm.back(solutionCompleteGraph);
				}
				li.clear();
			} else {
				break;
			}
		}
		
		if(maxNumNodesWithSolution == 0){
			return 0;
		}

		if(!objectiveCompleteGraph.isModeFragment()) {
			if (maxNumNodesWithSolution < MIN_NUM_NODES_SIM) {
				return 0;
			}
		}

		List<SolutionCompleteGraph> hsSolution = liliSolution.get(maxNumNodesWithSolution);
		List<SolutionCompleteGraph> li = new ArrayList<SolutionCompleteGraph>();
		for (SolutionCompleteGraph solution : hsSolution) {
			float similarity = objectiveCompleteGraph.getSimilarity(solution);
			solution.setSimilarity(similarity);
			li.add(solution);
		}
		
		Collections.sort(li);
		solutionBest.copyIntoThis(li.get(li.size()-1));

		if(DEBUG) {
			objectiveCompleteGraph.setVerbose(true);
			objectiveCompleteGraph.getSimilarity(solutionBest);
			objectiveCompleteGraph.setVerbose(false);
		}

		double similarity = solutionBest.getSimilarity();

		solutionsTop.clear();
		double simMaxMargin=similarity-TINY_SIM;
		for (SolutionCompleteGraph scg : li) {
			if(scg.getSimilarity()>simMaxMargin){
				solutionsTop.add(scg);
			}
		}

		return similarity;
	}

	public List<SolutionCompleteGraph> getSolutionsTop() {
		return solutionsTop;
	}

	/**
	 * Be careful! The histogram similarity is still considered if not explicitly set to false in the objective.
	 * @return
	 */
	public double calculateNodeSimilarity () {

		if(!((ObjectiveBlurFlexophoreHardMatchUncovered)objectiveCompleteGraph).isExcludeHistogramSimilarity()){
			throw new RuntimeException("Histogram similarity was not excluded from objective!");
		}

		initSearch();


		if(nodesBase==1 && nodesQuery==1) {
			double sim = objectiveCompleteGraph.getSimilarityNodes(0,0);
			return sim;
		}

		if(objectiveCompleteGraph.isModeFragment() && nodesQuery==1) {
			List<SolutionCompleteGraph> liSolution = liliSolution.get(1);
			double simMax=0;
			for (SolutionCompleteGraph solutionCompleteGraph : liSolution) {
				double similarity = objectiveCompleteGraph.getSimilarityNodes(solutionCompleteGraph);
				if(similarity>simMax){
					simMax=similarity;
					solutionBest=solutionCompleteGraph;
				}
			}
			return simMax;
		}

		int maxNumNodesWithSolution = 0;
		for (int nodesInSolution = 1; nodesInSolution < nodesBase+1; nodesInSolution++) {

			List<SolutionCompleteGraph> liSolution = liliSolution.get(nodesInSolution);
			boolean validSolutionFound = false;
			hsSolution.clear();
			for (SolutionCompleteGraph solution : liSolution) {
				if(getNeighbourSolutions(solution)){
					validSolutionFound = true;
				}
				if(hsSolution.size()> maxNumSolutions){
					break;
				}
			}

			if(validSolutionFound){
				maxNumNodesWithSolution = nodesInSolution+1;
				liliSolution.get(maxNumNodesWithSolution).addAll(hsSolution);
				// System.out.println("Found " + hsSolution.size() + " valid solutions for " + maxNumNodesWithSolution + " nodes.");

				//
				// Remove all smaller solutions
				//
				List<SolutionCompleteGraph> li = liliSolution.get(nodesInSolution);
				for (SolutionCompleteGraph solutionCompleteGraph : li) {
					cm.back(solutionCompleteGraph);
				}
				li.clear();
			} else {
				break;
			}
		}

		if(maxNumNodesWithSolution==0){
			return 0;
		}

		if(!objectiveCompleteGraph.isModeFragment()) {
			if (maxNumNodesWithSolution < MIN_NUM_NODES_SIM) {
				return 0;
			}
		}

		List<SolutionCompleteGraph> hsSolution = liliSolution.get(maxNumNodesWithSolution);

		List<SolutionCompleteGraph> li = new ArrayList<>();

		for (SolutionCompleteGraph solution : hsSolution) {
			double similarity = objectiveCompleteGraph.getSimilarityNodes(solution);
			solution.setSimilarity(similarity);
			li.add(solution);
		}

		Collections.sort(li);

		solutionBest.copyIntoThis(li.get(li.size()-1));

		if(DEBUG) {
			objectiveCompleteGraph.setVerbose(true);
			objectiveCompleteGraph.getSimilarity(solutionBest);
			objectiveCompleteGraph.setVerbose(false);
		}

		double similarity = solutionBest.getSimilarity();

		return similarity;
	}

	public SolutionCompleteGraph getBestMatchingSolution(){
		SolutionCompleteGraph solution = new SolutionCompleteGraph();
		solution.copyIntoThis(solutionBest);
		return solution;
	}
	
	private boolean getNeighbourSolutions(SolutionCompleteGraph solution){
		
		boolean validSolutionFound = false;
		
		int heap = solution.getSizeHeap();

		for (int i = 0; i < arrIndexBaseTmp.length; i++) {
			arrIndexBaseTmp[i]=0;
			arrIndexQueryTmp[i]=0;
		}
		
		for (int i = 0; i < heap; i++) {
			arrIndexBaseTmp[solution.getIndexBaseFromHeap(i)] = 1;
			arrIndexQueryTmp[solution.getIndexQueryFromHeap(i)] = 1;	 
		}
		
		for (byte indexNodeQuery = 0; indexNodeQuery < nodesQuery; indexNodeQuery++) {
			if(arrIndexQueryTmp[indexNodeQuery]==0){
				
				for (byte indexNodeBase = 0; indexNodeBase < nodesBase; indexNodeBase++) {
					if(arrIndexBaseTmp[indexNodeBase] == 0){
						
						if(objectiveCompleteGraph.areNodesMapping(indexNodeQuery, indexNodeBase)){
						
							SolutionCompleteGraph solutionNew = cm.getWithCopy(solution);

							// to do: is this necessary?, no 21.08.2024
							// solutionNew.copyIntoThis(solution);
							
							solutionNew.setNodesQuery(nodesQuery);
							
							solutionNew.add(indexNodeQuery, indexNodeBase);
							
							createdSolutions++;
							
							if(hsSolution.contains(solutionNew)) {
								cm.back(solutionNew);
							} else {
								if(objectiveCompleteGraph.isValidSolution(solutionNew)){
									hsSolution.add(solutionNew);
									validSolutionFound=true;
									validSolutions++;
								} else {
									cm.back(solutionNew);
								}
							}
						} 
					}
				}
			}
		}
		
		return validSolutionFound;
		
	}
	
	
	public long getValidSolutions() {
		return validSolutions;
	}
	
	public long getCreatedSolutions() {
		return createdSolutions;
	}

	public void setMaxNumSolutions(int maxNumSolutions) {
		this.maxNumSolutions = maxNumSolutions;
	}
	
	

}
