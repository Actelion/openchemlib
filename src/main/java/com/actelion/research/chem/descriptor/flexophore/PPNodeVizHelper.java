package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.flexophore.generator.MultCoordFragIndex;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphIndices;

import java.util.ArrayList;
import java.util.List;

/*
 

 
 Created by Modest von Korff 
 16/02/2024
 
 */
public class PPNodeVizHelper {

    public static PPNodeVizMultCoord createWithoutCoordinates(MultCoordFragIndex multCoordFragIndex, int indexPPPoint, Molecule3D mol){
        PPNodeViz ppNodeViz = new PPNodeViz();
        ppNodeViz.setIndex(indexPPPoint);
        for (int index : multCoordFragIndex.getArrIndexFrag()) {
            int interactionType = mol.getInteractionAtomType(index);
            ppNodeViz.add(interactionType);
            ppNodeViz.addIndexOriginalAtom(index);
        }
        ppNodeViz.realize();

        PPNodeVizMultCoord nodeVizMultCoord = new PPNodeVizMultCoord(ppNodeViz, multCoordFragIndex);

        return nodeVizMultCoord;
    }
    public static PPNodeViz createWithoutCoordinates(int [] arrIndexAtomFrag, int indexPPPoint, Molecule3D mol){

        PPNodeViz ppNodeViz = new PPNodeViz();
        ppNodeViz.setIndex(indexPPPoint);
        for (int index : arrIndexAtomFrag) {
            int interactionType = mol.getInteractionAtomType(index);
            ppNodeViz.add(interactionType);
            ppNodeViz.addIndexOriginalAtom(index);
        }
        ppNodeViz.realize();
        return ppNodeViz;
    }
    public static List<PPNodeViz> createWithoutCoordinates(List<SubGraphIndices> liSubGraphIndices, Molecule3D mol){
        List<PPNodeViz> liPPNodeViz = new ArrayList<>();
        for (int i = 0; i < liSubGraphIndices.size(); i++) {
            SubGraphIndices sgi = liSubGraphIndices.get(i);
            PPNodeViz ppNodeViz = createWithoutCoordinates(sgi.getAtomIndices(), i, mol);
            liPPNodeViz.add(ppNodeViz);
        }
        return liPPNodeViz;
    }

}
