package com.actelion.research.chem.descriptor.flexophore.completegraphmatcher;

import com.actelion.research.chem.descriptor.flexophore.PPNode;

/**
 * IPPNodeSimilarity
 * <p>Copyright: Actelion Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 02.03.16.
 */
public interface IPPNodeSimilarity {

    double getSimilarity(PPNode query, PPNode base);

    void setVerbose(boolean v);

}
