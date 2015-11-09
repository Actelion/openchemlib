/*
 * Project: DD_jfx
 * @(#)MoleculeFlavors.java
 *
 * Copyright (c) 1997- 2015
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * Author: Christian Rufener
 */

package com.actelion.research.gui.dnd;

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
