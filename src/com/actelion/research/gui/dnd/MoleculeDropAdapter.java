package com.actelion.research.gui.dnd;

import com.actelion.research.chem.*;

import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
/**
 * <p>Title: Mercury</p>
 * <p>Description: Actelion Electronic Lab Notebook</p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author Christian Rufener
 * @version 1.0
 */

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


    public void dragEnter(DropTargetDragEvent e)
    {
        DEBUG("DragEnter");
    }

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
    public void dropActionChanged(DropTargetDragEvent e)
    {
        DEBUG("dropActionChanged");
    }

    public void dragExit(DropTargetEvent e)
    {
        DEBUG("dragExit");
    }


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
        return MoleculeFlavors.FLAVORS;
    }


    protected StereoMolecule createFromDataFlavor(DataFlavor chosen, Object o) throws Exception
    {
        StereoMolecule mol = null;
        if (chosen.equals(MoleculeFlavors.DF_SERIALIZEDOBJECT) && o instanceof Molecule) {
            mol = new StereoMolecule((Molecule)o);
        } else if (chosen.equals(MoleculeFlavors.DF_MDLMOLFILE)
                   && o instanceof String) {
            mol = new StereoMolecule();
            new MolfileParser().parse(mol, (String)o);
        } else if (chosen.equals(MoleculeFlavors.DF_SMILES) && o instanceof String) {
            mol = new StereoMolecule();
            new SmilesParser().parse(mol, ((String)o).getBytes());
        } else if (chosen.equals(DataFlavor.stringFlavor) && o instanceof String) {
            try {
                mol = new StereoMolecule();
               new IDCodeParser(true).parse(mol, ((String) o).getBytes());
            } catch(Throwable t) {
                System.err.println("Unable to instantiate from text flavor " + o);
                mol = null;
            }
        } else {
            System.err.println("Unable to instantiate flavor " + chosen);
//            throw new InstantiationException("Unable to instantiate flavor " + chosen);
        }
        return mol;
    }
    
            
    protected boolean isDragFlavorSupported(DropTargetDragEvent e)
    {
        for (int i=0; i<MoleculeFlavors.FLAVORS.length; i++) {
            if (e.isDataFlavorSupported(MoleculeFlavors.FLAVORS[i])) {
                return true;
            }
        }
        return false;
    }

    protected DataFlavor chooseDropFlavor(DropTargetDropEvent e)
    {
        for (int i=0; i<MoleculeFlavors.FLAVORS.length; i++) {
            if (e.isDataFlavorSupported(MoleculeFlavors.FLAVORS[i])) {
                return MoleculeFlavors.FLAVORS[i];
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
