package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.chem.descriptor.flexophore.MolDistHist;
import com.actelion.research.chem.descriptor.flexophore.PPNode;

/*
 
 Copyright (c) 2024 Alipheron AG. All rights reserved.
 
 This file is part of the Alipheron AG software suite.
 
 Licensed under the Alipheron AG Software License Agreement (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at the company's official website or upon request.
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License. 
 
 Created by Modest von Korff 
 30/04/2024
 
 */
public class SubFlexophoreExtractorHelper {

    public static MolDistHist getSubFragment(MolDistHist mdh, int [] arrIndices){

        MolDistHist frag = new MolDistHist(arrIndices.length);

        for (int i = 0; i < arrIndices.length; i++) {
            PPNode node = new PPNode(mdh.getNode(arrIndices[i]));
            frag.addNode(node);
        }

        for (int i = 0; i < arrIndices.length; i++) {
            for (int j = i+1; j < arrIndices.length; j++) {
                byte [] arrHist = mdh.getDistHist(arrIndices[i],arrIndices[j]);
                frag.setDistHist(i,j,arrHist);
            }
        }

        frag.realize();
        return frag;
    }

}
