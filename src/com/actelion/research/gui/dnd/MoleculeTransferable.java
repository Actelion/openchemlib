/*
 * Project: DD_jfx
 * @(#)MoleculeTransferable.java
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

import com.actelion.research.chem.*;

import java.awt.datatransfer.*;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class MoleculeTransferable implements Transferable,ClipboardOwner 
{
        private static final List cFlavorList = Arrays.asList( MoleculeFlavors.FLAVORS );
        protected ExtendedMolecule mMol;

        public MoleculeTransferable(ExtendedMolecule mol) {
                mMol = mol;
                }
        public synchronized DataFlavor[] getTransferDataFlavors() {
//            System.out.println("Moleculetransferable getTransferFlavors");
                return MoleculeFlavors.FLAVORS;
                }

        public boolean isDataFlavorSupported( DataFlavor flavor ) {
//            System.out.println("Moleculetransferable  isdataflavor supported : " + flavor);
//                return (cFlavorList.contains(flavor));
            for (int i = 0; i<MoleculeFlavors.FLAVORS.length;i++) {
                if (MoleculeFlavors.FLAVORS.equals(flavor))
                return true;
            }
            return false;
                }

        public synchronized Object getTransferData(DataFlavor flavor)
                                        throws UnsupportedFlavorException,IOException 
    {
//        System.out.println("MoleculeTransferable flavor " + flavor);
        if (flavor.equals(MoleculeFlavors.DF_SERIALIZEDOBJECT)) {
            return new StereoMolecule(mMol);
        } else if (flavor.equals(MoleculeFlavors.DF_MDLMOLFILEV3)) {
            return new MolfileV3Creator(mMol).getMolfile();
        } else if (flavor.equals(MoleculeFlavors.DF_MDLMOLFILE)) {
            return new MolfileCreator(mMol).getMolfile();
        } else if (flavor.equals(MoleculeFlavors.DF_SMILES)) {
            return new SmilesCreator().generateSmiles(mMol);
        } else if (flavor.equals(DataFlavor.stringFlavor)) {
            return new Canonizer(mMol).getIDCode();
        } else
            throw new UnsupportedFlavorException(flavor);
    }

    public String toString()
    {
        return "MoleculeTransferable";
    }

    public void lostOwnership(Clipboard clipboard, Transferable contents)
    {
    }
}

