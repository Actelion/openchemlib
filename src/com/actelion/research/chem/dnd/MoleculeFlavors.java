/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Christian Rufener
 */

package com.actelion.research.chem.dnd;

import java.awt.datatransfer.DataFlavor;

public class MoleculeFlavors
{
    static class MyFlavor extends DataFlavor
    {
        public MyFlavor(Class representationClass, String humanPresentableName)
        {
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

    public static final DataFlavor DF_SERIALIZEDOBJECT = new MyFlavor(com.actelion.research.chem.Molecule.class, "Actelion Molecule Class");
    public static final DataFlavor DF_MDLMOLFILE = new DataFlavor("chemical/x-mdl-molfile;class=java.lang.String", "MDL Molfile");
    public static final DataFlavor DF_MDLMOLFILEV3 = new DataFlavor("chemical/x-mdl-molfilev3;class=java.lang.String", "MDL Molfile V3");
    public static final DataFlavor DF_SMILES = new DataFlavor("chemical/x-daylight-smiles;class=java.lang.String", "Daylight Smiles");
    public static final DataFlavor[] FLAVORS = {
//                DF_SERIALIZEDSTRUCTURETRANSFERDATA,
        DF_SERIALIZEDOBJECT,
        DF_MDLMOLFILE,
        DF_MDLMOLFILEV3,
        DF_SMILES,
        DataFlavor.stringFlavor
    };
}
