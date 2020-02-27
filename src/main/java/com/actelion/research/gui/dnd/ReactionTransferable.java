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

import com.actelion.research.chem.IsomericSmilesCreator;
import com.actelion.research.chem.dnd.ChemistryFlavors;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;

import java.awt.datatransfer.*;
import java.io.IOException;

public class ReactionTransferable implements Transferable,ClipboardOwner {
	protected Reaction mReaction;

	public ReactionTransferable(Reaction rxn) {
		mReaction = rxn;
	}

	public synchronized DataFlavor[] getTransferDataFlavors() {
//            System.out.println("ReactionTransferable getTransferFlavors");
		return ChemistryFlavors.REACTION_FLAVORS;
	}

	public boolean isDataFlavorSupported(DataFlavor flavor) {
		if (flavor.equals(ChemistryFlavors.DF_REACTION_SMILES))		// currently no support for reaction smiles
			return false;

		for (DataFlavor f:ChemistryFlavors.REACTION_FLAVORS) {
			if (f.equals(flavor))
				return true;
		}
		return false;
	}

	public synchronized Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException,IOException {
//        System.out.println("ReactionTransferable flavor " + flavor);
		if (flavor.equals(ChemistryFlavors.DF_SERIALIZED_REACTION)) {
			return new Reaction(mReaction);
		} else if (flavor.equals(ChemistryFlavors.DF_REACTION_SMILES)) {
			return IsomericSmilesCreator.createReactionSmiles(mReaction);
		} else if (flavor.equals(DataFlavor.stringFlavor)) {
			return ReactionEncoder.encode(mReaction, true);
		} else
			throw new UnsupportedFlavorException(flavor);
	}

	public String toString()
	{
		return "ReactionTransferable";
	}

	public void lostOwnership(Clipboard clipboard, Transferable contents) {
	}
}

