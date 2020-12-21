package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.util.graph.complete.SolutionCompleteGraph;

public class ModelSolutionSimilarity extends SolutionCompleteGraph {

    private float [] arrSimilarityNodes;

    public ModelSolutionSimilarity(SolutionCompleteGraph scg, float [] arrSimilarityNodes) {
        super.copyIntoThis(scg);
        this.arrSimilarityNodes = arrSimilarityNodes;
    }

    public float getSimilarityNode(int index){
        return arrSimilarityNodes[index];
    }


}
