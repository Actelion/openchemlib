package com.actelion.research.util.graph.complete;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

/**
 * 
 * 
 * CompleteGraphMatcher
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Sep 26, 2012 MvK: Start implementation
 */
public class CompleteGraphMatcher<T extends ICompleteGraph> {

	public static final boolean DEBUG = false;

	// We need at least two pharmacophore points.
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
	
	private long validSolutions;
	
	private long createdSolutions;
	
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
//	
//	/**
//	 * Only for debugging!!!
//	 */
//	public IObjectiveCompleteGraph<T> getObjective(){
//		return objectiveCompleteGraph;
//	}
	
	private void init(){
		
		maxNumSolutions = MAX_NUM_SOLUTIONS;
		
		cm = new ContainerMemory<SolutionCompleteGraph>(INIT_CAPACITY_MEMORY, new FactorySolution());
		
		liliSolution = new ArrayList<List<SolutionCompleteGraph>>(MAX_NUM_NODES);
		
		liliSolution.add(new ArrayList<SolutionCompleteGraph>());
		for (int i = 1; i < MAX_NUM_NODES; i++) {
			liliSolution.add(new ArrayList<SolutionCompleteGraph>(SIZE_LIST_SOLUTION));
		}
		
		hsSolution = new HashSet<SolutionCompleteGraph>(SIZE_LIST_SOLUTION);
		
		arrIndexBaseTmp = new byte [MAX_NUM_NODES];
		arrIndexQueryTmp = new byte [MAX_NUM_NODES];
		
		solutionBest = new SolutionCompleteGraph();
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
		
		List<SolutionCompleteGraph> liSolution = liliSolution.get(nodesInSolution);
		
		for (byte indexNodeQuery = 0; indexNodeQuery < nodesQuery; indexNodeQuery++) {
			
			for (byte indexNodeBase = 0; indexNodeBase < nodesBase; indexNodeBase++) {
			
				if(objectiveCompleteGraph.areNodesMapping(indexNodeQuery, indexNodeBase)){
					
					SolutionCompleteGraph solution = cm.get();
					
					solution.setNodesQuery(nodesQuery);
					
					solution.add(indexNodeQuery, indexNodeBase);
					
					liSolution.add(solution);
					
				}
			}
		}
	}
	
	public double calculateSimilarity () {
		
		initSearch();
				
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
		
		if(maxNumNodesWithSolution < MIN_NUM_NODES_SIM){
			return 0;
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
							
							solutionNew.copyIntoThis(solution);
							
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
