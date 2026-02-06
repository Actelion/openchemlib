package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.flexophore.generator.MultCoordFragIndex;
import com.actelion.research.chem.descriptor.flexophore.redgraph.SubGraphIndices;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Copyright (c) 2026
 * Alipheron AG
 * Hochbergerstrasse 60C
 * CH-4057 Basel
 * Switzerland
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
 * 3. Neither the name of the copyright holder nor the
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
 *
 * @author Modest v. Korff
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

    /**
     * Must already be sorted.
     * @param l0
     * @param l1
     * @return
     */
    public static boolean equals(List<PPNodeViz> l0, List<PPNodeViz> l1){
        if(l0.size()!=l1.size()) return false;
        boolean eq=true;
        for (int i = 0; i < l0.size(); i++) {
            if(!l0.get(i).equals(l1.get(i))){
                eq=false;
                break;
            }
        }
        return eq;
    }
}
