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

import com.actelion.research.chem.descriptor.flexophore.IPPNode;

/**
 * IObjectiveCompleteGraph
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
	float getSimilarityNodes(SolutionCompleteGraph solution);

	float getSimilarityHistogram(int indexNode1Query, int indexNode2Query, int indexNode1Base, int indexNode2Base);

	double getSimilarityNodes(IPPNode query, IPPNode base);

	double getSimilarityNodes(int indexNodeQuery, int indexNodeBase);

	void setVerbose(boolean v);

	void setMatchingInfoInQueryAndBase(SolutionCompleteGraph solution);

	boolean isModeFragment();

}
