package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.util.graph.complete.IObjectiveCompleteGraph;

public interface IDescriptorHandlerFlexophore extends DescriptorHandler {

    IObjectiveCompleteGraph getObjectiveCompleteGraph();

    MolDistHistViz createVisualDescriptor(StereoMolecule mol);

    ModelSolutionSimilarity getBestMatch(MolDistHistViz mdhvBase, MolDistHistViz mdhvQuery);

}
