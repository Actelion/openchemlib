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
