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
	
	public T getBase();
	
	public T getQuery();
	
	public void setBase(T cgBase);
	
	public void setQuery(T cgQuery);
	
	public abstract boolean areNodesMapping(int indexNodeBase, int indexNodeQuery);
	
	public abstract boolean isValidSolution(SolutionCompleteGraph solution);
	
	public abstract float getSimilarity(SolutionCompleteGraph solution);

	void setVerbose(boolean v);

}
