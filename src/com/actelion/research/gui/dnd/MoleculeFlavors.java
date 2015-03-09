package com.actelion.research.gui.dnd;

/**
 * <p>Title: DD_GUI</p>
 * Copyright 1997-2011 Actelion Ltd., Inc. All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 * @author Christian Rufener
 * @version 1.0
 */

import java.awt.datatransfer.DataFlavor;

@Deprecated
public class MoleculeFlavors
{
    public static final DataFlavor DF_SERIALIZEDOBJECT  = com.actelion.research.chem.dnd.MoleculeFlavors.DF_SERIALIZEDOBJECT;
    public static final DataFlavor DF_MDLMOLFILE        = com.actelion.research.chem.dnd.MoleculeFlavors.DF_MDLMOLFILE;
    public static final DataFlavor DF_MDLMOLFILEV3      = com.actelion.research.chem.dnd.MoleculeFlavors.DF_MDLMOLFILEV3;
    public static final DataFlavor DF_SMILES            = com.actelion.research.chem.dnd.MoleculeFlavors.DF_SMILES;
        public static final DataFlavor[] FLAVORS = { 
//                DF_SERIALIZEDSTRUCTURETRANSFERDATA,
                DF_SERIALIZEDOBJECT,
                 DF_MDLMOLFILE,
                DF_MDLMOLFILEV3,
                 DF_SMILES,
            DataFlavor.stringFlavor
    };

}
