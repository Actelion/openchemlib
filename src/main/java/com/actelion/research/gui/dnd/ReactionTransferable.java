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
 * @author Christian Rufener
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
			Reaction rxn = new Reaction(mReaction);
			rxn.removeDrawingObjects(); // to include drawing objects make sure they are all serializable
			return rxn;
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

