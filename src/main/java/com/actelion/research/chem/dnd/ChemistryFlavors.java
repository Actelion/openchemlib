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
*
*/

package com.actelion.research.chem.dnd;

import java.awt.datatransfer.DataFlavor;

public class ChemistryFlavors {
    static class SerializedClassFlavor extends DataFlavor {
        public SerializedClassFlavor(Class representationClass, String humanPresentableName) {
            super(representationClass, humanPresentableName);
            try {
                ClassLoader cl = this.getClass().getClassLoader();
                //System.out.println("MyFlavour Ignoring classloader!");
                Thread.currentThread().setContextClassLoader(cl);
            } catch (Throwable ex) {
                //System.err.println("KLUDGE: This exception is thrown due to the applet context");
            }
        }
    }

    public static final DataFlavor DF_SERIALIZED_MOLECULE = new SerializedClassFlavor(com.actelion.research.chem.StereoMolecule.class, "Native OpenChemLib Molecule");
    public static final DataFlavor DF_MDLMOLFILE = new DataFlavor("chemical/x-mdl-molfile;class=java.lang.String", "MDL Molfile");
    public static final DataFlavor DF_MDLMOLFILEV3 = new DataFlavor("chemical/x-mdl-molfilev3;class=java.lang.String", "MDL Molfile V3");
    public static final DataFlavor DF_SMILES = new DataFlavor("chemical/x-daylight-smiles;class=java.lang.String", "Daylight Smiles");
    public static final DataFlavor DF_IDCODE = new DataFlavor("chemical/x-openmolecules-idcode;class=java.lang.String", "OpenChemLib ID-Code");
    public static final DataFlavor[] MOLECULE_FLAVORS = {
        DF_SERIALIZED_MOLECULE,
        DF_MDLMOLFILE,
        DF_MDLMOLFILEV3,
        DF_SMILES,
        DF_IDCODE,
        DataFlavor.stringFlavor
    };

    public static final DataFlavor DF_SERIALIZED_REACTION = new SerializedClassFlavor(com.actelion.research.chem.reaction.Reaction.class, "Native OpenChemLib Reaction");
    public static final DataFlavor DF_REACTION_SMILES = new DataFlavor("chemical/x-daylight-reactionsmiles;class=java.lang.String", "Daylight Reaction Smiles");
    public static final DataFlavor[] REACTION_FLAVORS = {
        DF_SERIALIZED_REACTION,
        DF_REACTION_SMILES,
    };
}
