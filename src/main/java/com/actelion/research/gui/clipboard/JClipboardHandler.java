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
package com.actelion.research.gui.clipboard;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.dnd.ChemistryFlavors;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.dnd.MoleculeTransferable;
import com.actelion.research.gui.dnd.ReactionTransferable;

import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;

/**
 * Clipboard handler for Swing/AWT applications. If available, it uses a native clipboard handler,
 * otherwise it uses the AWT system clipboard. JavaFX applications should use FXClipboardHandler.
 */
public class JClipboardHandler extends AbstractClipboardHandler {
    @Override
    protected String getClipboardText() {
        try {
            Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
            return (String) t.getTransferData(DataFlavor.stringFlavor);
        } catch (Exception e) {
            // ignored!
        }
        return null;
    }

    @Override
    protected StereoMolecule getSerializedMolecule() {
        try {
            Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
            return (StereoMolecule) t.getTransferData(ChemistryFlavors.DF_SERIALIZED_MOLECULE);
        } catch (Exception e) {
            // ignored!
        }
        return null;
    }

    @Override
    protected Reaction getSerializedReaction() {
        try {
            Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
            return (Reaction) t.getTransferData(ChemistryFlavors.DF_SERIALIZED_REACTION);
        } catch (Exception e) {
            // ignored!
        }
        return null;
    }

    @Override
    protected void putMoleculeContent(StereoMolecule mol) {
        MoleculeTransferable transferable = new MoleculeTransferable(mol);
        Toolkit.getDefaultToolkit().getSystemClipboard().setContents(transferable, transferable);
    }

    @Override
    protected void putReactionContent(Reaction rxn) {
        ReactionTransferable transferable = new ReactionTransferable(rxn);
        Toolkit.getDefaultToolkit().getSystemClipboard().setContents(transferable, transferable);
    }

    @Override
    protected void clearClipboard() {
        emptyClipboard();
    }

    public static void emptyClipboard() {
        Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new Transferable() {
            public DataFlavor[] getTransferDataFlavors() {
                return new DataFlavor[0];
            }

            public boolean isDataFlavorSupported(DataFlavor flavor) {
                return false;
            }

            public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException {
                throw new UnsupportedFlavorException(flavor);
            }
        }, null);
    }

    /**
     * Copies an Image to the clipboard
     *
     * @param img Image to be copied
     * @return true on success
     */
    @Override
    public boolean copyImage(java.awt.Image img) {
        return ImageClipboardHandler.copyImage(img);
    }

    @Override
    public java.awt.Image pasteImage() {
        return ImageClipboardHandler.pasteImage();
    }

    /**
     * @deprecated Use ImageClipboardHandler.pasteImage for consistency reasons
     */
    @Deprecated
    public static Image getImage() {
        return ImageClipboardHandler.pasteImage();
    }

    /**
     * @deprecated You may use ImageClipboardHandler.copyImage for consistency reasons
     */
    @Deprecated
    public static void putImage(Image image) {
        ImageClipboardHandler.copyImage(image);
    }
}
