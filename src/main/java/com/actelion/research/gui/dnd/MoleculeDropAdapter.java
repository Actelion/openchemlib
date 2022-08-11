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

import com.actelion.research.chem.*;
import com.actelion.research.chem.dnd.ChemistryFlavors;
import com.actelion.research.chem.name.StructureNameResolver;

import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;

public class MoleculeDropAdapter implements DropTargetListener
{
    public static final boolean debug = false;
    private boolean active_ = true;

    private final int mAllowedDropAction = DnDConstants.ACTION_COPY_OR_MOVE;	// changed TLS 9-Jul-04
    public MoleculeDropAdapter()
    {
    }

    public void onDropMolecule(StereoMolecule m, Point pt)
    {
        DEBUG("MoleculeDropAdapter.onDropMolecule(). Override this! " + m);
    }

    @Override
    public void dragEnter(DropTargetDragEvent e)
    {
        DEBUG("DragEnter");
    }

    @Override
    public void dragOver(DropTargetDragEvent e)
    {
       DEBUG("DragOver");
    }

    public void setActive(boolean active)
    {
        active_ = active;
    }

    public boolean isActive()
    {
        return active_;
    }

    @Override
    public void dropActionChanged(DropTargetDragEvent e)
    {
        DEBUG("dropActionChanged");
    }

    @Override
    public void dragExit(DropTargetEvent e)
    {
        DEBUG("dragExit");
    }

    @Override
    public void drop(DropTargetDropEvent e)
    {
        if (active_) {
            // This is neccesary to make sure the correct classloader tries to load the Transferable
            ClassLoader cl = this.getClass().getClassLoader();
            DEBUG("MoleculeDropAdapter   ClassLoader " + cl);
            DEBUG("MoleculeDropAdapter   Ignoring setContextclassloader!!!");
//            Thread.currentThread().setContextClassLoader(cl);
            DEBUG("MoleculeDropAdapter " + e);
            try {
                Transferable tr = e.getTransferable();
                DEBUG("Transferable is " + tr);
                DataFlavor chosen = chooseDropFlavor(e);
                StereoMolecule mol = null;
                if (chosen != null) {
                    e.acceptDrop(DnDConstants.ACTION_COPY_OR_MOVE);
                    DEBUG("Chose is " + chosen);
                    Object o = tr.getTransferData(chosen);
                    DEBUG("Object is " + o);
                    mol = createFromDataFlavor(chosen,o);
                    if (mol != null) {
                        onDropMolecule(mol,e.getLocation());
                        e.dropComplete(true);
                    } else {
                        System.err.println("Drop failed: " + e);
                        e.dropComplete(false);
                    }
                    return;
                } else {
                    System.err.println("Drop failed: " + e);
                    e.rejectDrop();
                }
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }   // active_
    }

    public DataFlavor[] getFlavors()
    {
        return ChemistryFlavors.MOLECULE_FLAVORS;
    }

    protected StereoMolecule createFromDataFlavor(DataFlavor chosen, Object o) throws Exception
    {
        StereoMolecule mol = null;
        if (chosen.equals(ChemistryFlavors.DF_SERIALIZED_MOLECULE) && o instanceof Molecule) {
            mol = new StereoMolecule((Molecule)o);
        } else if (chosen.equals(ChemistryFlavors.DF_MDLMOLFILE)
                   && o instanceof String) {
            mol = new StereoMolecule();
            new MolfileParser().parse(mol, (String)o);
        } else if (chosen.equals(ChemistryFlavors.DF_SMILES) && o instanceof String) {
            mol = new StereoMolecule();
            new SmilesParser().parse(mol, ((String)o).getBytes());
        } else if (chosen.equals(ChemistryFlavors.DF_IDCODE) && o instanceof String) {
            mol = new StereoMolecule();
            new IDCodeParser(true).parse(mol, ((String) o).getBytes());
        } else if (chosen.equals(DataFlavor.stringFlavor) && o instanceof String) {
            try {
                mol = new StereoMolecule();
               new IDCodeParser(true).parse(mol, ((String) o).getBytes());
            } catch(Throwable t) {
                mol = StructureNameResolver.resolve((String) o);
            }
            if (mol == null)
                System.err.println("Unable to instantiate from text flavor: " + o);
        } else {
            System.err.println("Unable to instantiate flavor: " + chosen);
//            throw new InstantiationException("Unable to instantiate flavor " + chosen);
        }
        return mol;
    }
    
    protected boolean isDragFlavorSupported(DropTargetDragEvent e)
    {
        for (int i=0; i<ChemistryFlavors.MOLECULE_FLAVORS.length; i++) {
            if (e.isDataFlavorSupported(ChemistryFlavors.MOLECULE_FLAVORS[i])) {
                return true;
            }
        }
        return false;
    }

    protected DataFlavor chooseDropFlavor(DropTargetDropEvent e)
    {
        for (int i=0; i<ChemistryFlavors.MOLECULE_FLAVORS.length; i++) {
            if (e.isDataFlavorSupported(ChemistryFlavors.MOLECULE_FLAVORS[i])) {
                return ChemistryFlavors.MOLECULE_FLAVORS[i];
            }
        }
        return null;
    }

    public boolean isDropOK(DropTargetDragEvent e)
    {
        if (!isDragFlavorSupported(e)) {
            return false;
        }
        return ((e.getDropAction() & mAllowedDropAction) != 0);
    }

    private void DEBUG(String s)
    {
        if (debug) {
            System.err.println(s);
            System.err.flush();
        }
    }
}
/////
