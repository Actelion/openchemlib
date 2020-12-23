package com.actelion.research.util.graph.complete;

/**
 * 
 * 
 * IObjectiveCompleteGraph
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Oct 1, 2012 MvK: Start implementation
 */
public interface IObjectiveCompleteGraph<T extends ICompleteGraph> {
	
	T getBase();
	
	T getQuery();
	
	void setBase(T cgBase);
	
	void setQuery(T cgQuery);
	
	boolean areNodesMapping(int indexNodeBase, int indexNodeQuery);
	
	boolean isValidSolution(SolutionCompleteGraph solution);
	
	float getSimilarity(SolutionCompleteGraph solution);

	void setVerbose(boolean v);

	void setMatchingInfoInQueryAndBase(SolutionCompleteGraph solution);

}
