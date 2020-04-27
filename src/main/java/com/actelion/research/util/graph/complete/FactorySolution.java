package com.actelion.research.util.graph.complete;


/**
 * 
 * 
 * FactorySolution
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Oct 1, 2012 MvK: Start implementation
 */
public class FactorySolution implements IFactory<SolutionCompleteGraph>  {
	
	public SolutionCompleteGraph createObject() {
		
		SolutionCompleteGraph solution = new SolutionCompleteGraph();
		
		return solution;
	}

}
