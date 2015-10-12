/*
 * Project: DD_jfx
 * @(#)SDFileMoleculeDropAdapter.java
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

import java.awt.dnd.DropTargetDropEvent;
import com.actelion.research.chem.*;

import java.awt.datatransfer.DataFlavor;
import java.awt.dnd.DropTargetDragEvent;

public class SDFileMoleculeDropAdapter extends MoleculeDropAdapter
{
    public SDFileMoleculeDropAdapter()
    {
    }

    protected StereoMolecule createFromDataFlavor(DataFlavor chosen, Object o) throws Exception
    {
        StereoMolecule mol = null;
        if (chosen.equals(SDFileMoleculeTransferable.DF_SERIALIZEDSTRUCTURETRANSFERDATA) && o instanceof StructureTransferData) {
//            System.out.println(this + " createFromDataFlavor " );
            StructureTransferData d =  (StructureTransferData)o;

            SDFileMolecule sdf = new SDFileMolecule(d.getMolecule());
            StructureInfo si = d.getStructureInfo();
			if (si != null) {
				sdf.setFieldData(SDFileMolecule.ID_FIELD,si.getId());
				sdf.setFieldData(SDFileMolecule.NAME_FIELD,si.getName());
				sdf.setFieldData(SDFileMolecule.CMNT_FIELD,si.getComment());
				sdf.setFieldData(SDFileMolecule.DENSITY_FIELD,si.getDensity());
				sdf.setFieldData(SDFileMolecule.PURITY_FIELD,si.getPurity());
			}
            mol = sdf;
//            System.out.println(this + " createFromDataFlavor " + d.getMolecule() );
        } else 
            return super.createFromDataFlavor(chosen,o);
        return mol;
    }
    
    public DataFlavor[] getFlavors()
    {
//        System.out.println("getTransferDataFlavors " + this);
        DataFlavor[] res = super.getFlavors();
        DataFlavor[] d = res;
        if (d != null) {
            res = new DataFlavor[d.length+1];
            System.arraycopy(d,0,res,1,d.length);
            res[0] = SDFileMoleculeTransferable.DF_SERIALIZEDSTRUCTURETRANSFERDATA;
        } 
        return res;
    }
    
    protected boolean isDragFlavorSupported(DropTargetDragEvent e)
    {
        if (e.isDataFlavorSupported(SDFileMoleculeTransferable.DF_SERIALIZEDSTRUCTURETRANSFERDATA))
            return true;
        else
            return super.isDragFlavorSupported(e);
    }

    protected DataFlavor chooseDropFlavor(DropTargetDropEvent e)
    {
        if (e.isDataFlavorSupported(SDFileMoleculeTransferable.DF_SERIALIZEDSTRUCTURETRANSFERDATA))
                return SDFileMoleculeTransferable.DF_SERIALIZEDSTRUCTURETRANSFERDATA;
        else
            return super.chooseDropFlavor(e);
    }
}
