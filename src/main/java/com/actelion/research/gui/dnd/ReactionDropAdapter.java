/*
 * Project: DD_jfx
 * @(#)MoleculeDropAdapter.java
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
import com.actelion.research.chem.dnd.ChemistryFlavors;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;

import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
public class ReactionDropAdapter implements DropTargetListener
{
	public static final boolean debug = false;

	private boolean active_ = true;

	private final int mAllowedDropAction = DnDConstants.ACTION_COPY_OR_MOVE;

	public ReactionDropAdapter()
	{
	}

	public void onDropReaction(Reaction rxn, Point pt)
	{
		DEBUG("ReactionDropAdapter.onDropReaction(). Override this! " + rxn);
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
			DEBUG("ReactionDropAdapter   ClassLoader " + cl);
			DEBUG("ReactionDropAdapter   Ignoring setContextclassloader!!!");
//            Thread.currentThread().setContextClassLoader(cl);
			DEBUG("ReactionDropAdapter " + e);
			try {
				Transferable tr = e.getTransferable();
				DEBUG("Transferable is " + tr);
				DataFlavor chosen = chooseDropFlavor(e);
				Reaction rxn = null;
				if (chosen != null) {
					e.acceptDrop(DnDConstants.ACTION_COPY_OR_MOVE);
					DEBUG("Chose is " + chosen);
					Object o = tr.getTransferData(chosen);
					DEBUG("Object is " + o);
					rxn = createFromDataFlavor(chosen,o);
					if (rxn != null) {
						onDropReaction(rxn, e.getLocation());
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
		}
	}

	public DataFlavor[] getFlavors()
	{
		return ChemistryFlavors.REACTION_FLAVORS;
	}


	protected Reaction createFromDataFlavor(DataFlavor chosen, Object o) throws Exception {
		Reaction rxn = null;
		if (chosen.equals(ChemistryFlavors.DF_SERIALIZED_REACTION) && o instanceof Reaction) {
			rxn = new Reaction((Reaction)o);
		} else if (chosen.equals(ChemistryFlavors.DF_REACTION_SMILES) && o instanceof String) {
			rxn = new SmilesParser().parseReaction(((String)o).getBytes());
		} else if (chosen.equals(DataFlavor.stringFlavor) && o instanceof String) {
			try {
				rxn = ReactionEncoder.decode((String)o, true);
			} catch(Throwable t) {
				System.err.println("Unable to instantiate reaction from text flavor " + o);
				rxn = null;
			}
		} else {
			System.err.println("Unable to instantiate flavor " + chosen);
//            throw new InstantiationException("Unable to instantiate flavor " + chosen);
		}
		return rxn;
	}


	protected boolean isDragFlavorSupported(DropTargetDragEvent e) {
		for (int i=0; i<ChemistryFlavors.REACTION_FLAVORS.length; i++) {
			if (e.isDataFlavorSupported(ChemistryFlavors.REACTION_FLAVORS[i])) {
				return true;
			}
		}
		return false;
	}

	protected DataFlavor chooseDropFlavor(DropTargetDropEvent e) {
		for (int i=0; i<ChemistryFlavors.REACTION_FLAVORS.length; i++) {
			if (e.isDataFlavorSupported(ChemistryFlavors.REACTION_FLAVORS[i])) {
				return ChemistryFlavors.REACTION_FLAVORS[i];
			}
		}
		return null;
	}

	public boolean isDropOK(DropTargetDragEvent e) {
		if (!isDragFlavorSupported(e)) {
			return false;
		}
		return ((e.getDropAction() & mAllowedDropAction) != 0);
	}

	private void DEBUG(String s) {
		if (debug) {
			System.err.println(s);
			System.err.flush();
		}
	}
}
