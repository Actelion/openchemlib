/*
 * Project: Mercury
 * @(#)SDFileMoleculeTransferable.java
 *
 * Copyright (c) 1997-2006 
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
import java.io.IOException;
import com.actelion.research.chem.ExtendedMolecule;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;

import java.awt.datatransfer.DataFlavor;

public class SDFileMoleculeTransferable extends MoleculeTransferable
{
    public static final DataFlavor DF_SERIALIZEDSTRUCTURETRANSFERDATA =new DataFlavor(com.actelion.research.chem.StructureTransferData.class,"Actelion Structure Transfer Data Class");

    private DataFlavor FLAVORS[] = {DF_SERIALIZEDSTRUCTURETRANSFERDATA};
    StructureInfo  si_ = null;
    public SDFileMoleculeTransferable(ExtendedMolecule mol, StructureInfo  si)
    {
        super(mol);
        si_ = si;
    }

    public synchronized DataFlavor[] getTransferDataFlavors()
    {
        DataFlavor[] res = FLAVORS;
        DataFlavor[] d = super.getTransferDataFlavors();
        if (d != null) {
            res = new DataFlavor[d.length+1];
            System.arraycopy(d,0,res,1,d.length);
            res[0] = DF_SERIALIZEDSTRUCTURETRANSFERDATA;
        } 
        return res;
    }

    public boolean isDataFlavorSupported(DataFlavor flavor)
    {
        if (flavor.equals(DF_SERIALIZEDSTRUCTURETRANSFERDATA))
            return true;
        return super.isDataFlavorSupported(flavor);
    }

    public synchronized Object getTransferData(DataFlavor flavor) throws
            UnsupportedFlavorException,IOException
    {
        if(flavor.equals(DF_SERIALIZEDSTRUCTURETRANSFERDATA)) {
            return new StructureTransferData(mMol,si_);
        }

        return super.getTransferData(flavor);

    }

    public String toString()
    {
        return "SDFileMoleculeTransferable";
    }

    public void lostOwnership(Clipboard clipboard,Transferable contents)
    {
    }
}
