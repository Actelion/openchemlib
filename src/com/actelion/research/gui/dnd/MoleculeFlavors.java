package com.actelion.research.gui.dnd;

/*
* Copyright (c) 1997 - 2015
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
