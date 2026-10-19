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

import com.actelion.research.chem.IsomericSmilesCreator;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;
import com.actelion.research.gui.fx.FXMoleculeTransfer;
import javafx.application.Platform;
import javafx.embed.swing.SwingFXUtils;
import javafx.scene.input.Clipboard;
import javafx.scene.input.ClipboardContent;
import javafx.scene.input.DataFormat;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.FutureTask;

/**
 * Clipboard handler for JavaFX applications. If available, it uses a native clipboard handler,
 * otherwise it uses the JavaFX system clipboard. Unlike JClipboardHandler it doesn't need AWT
 * for clipboard access. Molecules are exchanged with Swing applications via idcode text and molfiles.
 * Methods may be called from any thread; JavaFX clipboard access is delegated to the FX application thread.
 */
public class FXClipboardHandler extends AbstractClipboardHandler {
    private static final String MIME_REACTION_SMILES = "chemical/x-daylight-reactionsmiles";

    @Override
    protected String getClipboardText() {
        return runOnFXThread(() -> {
            Clipboard clipboard = Clipboard.getSystemClipboard();
            return clipboard.hasString() ? clipboard.getString() : null;
        });
    }

    @Override
    protected StereoMolecule getSerializedMolecule() {
        // text is interpreted by AbstractClipboardHandler; here we only look at chemistry specific formats
        return runOnFXThread(() -> FXMoleculeTransfer.readMolecule(Clipboard.getSystemClipboard(), false));
    }

    @Override
    protected Reaction getSerializedReaction() {
        return null;    // reaction idcodes, rxn-files and reaction SMILES are interpreted from text by AbstractClipboardHandler
    }

    @Override
    protected void putMoleculeContent(StereoMolecule mol) {
        ClipboardContent content = FXMoleculeTransfer.createContent(mol);
        runOnFXThread(() -> Clipboard.getSystemClipboard().setContent(content));
    }

    @Override
    protected void putReactionContent(Reaction rxn) {
        HashMap<DataFormat, Object> content = new HashMap<>();
        content.put(DataFormat.PLAIN_TEXT, ReactionEncoder.encode(rxn, true));
        try {
            DataFormat format = DataFormat.lookupMimeType(MIME_REACTION_SMILES);
            content.put(format != null ? format : new DataFormat(MIME_REACTION_SMILES), IsomericSmilesCreator.createReactionSmiles(rxn));
        } catch (Exception e) {
            // ignored!
        }
        runOnFXThread(() -> Clipboard.getSystemClipboard().setContent(content));
    }

    @Override
    protected void clearClipboard() {
        runOnFXThread(() -> { Clipboard.getSystemClipboard().clear(); return null; });
    }

    @Override
    public boolean copyImage(Image img) {
        if (img == null)
            return false;

        BufferedImage bi;
        if (img instanceof BufferedImage) {
            bi = (BufferedImage)img;
        } else {
            bi = new BufferedImage(img.getWidth(null), img.getHeight(null), BufferedImage.TYPE_INT_ARGB);
            Graphics2D g = bi.createGraphics();
            g.drawImage(img, 0, 0, null);
            g.dispose();
        }

        final javafx.scene.image.Image fxImage = SwingFXUtils.toFXImage(bi, null);
        Boolean ok = runOnFXThread(() -> {
            ClipboardContent content = new ClipboardContent();
            content.putImage(fxImage);
            return Clipboard.getSystemClipboard().setContent(content);
        });
        return ok != null && ok;
    }

    @Override
    public Image pasteImage() {
        javafx.scene.image.Image fxImage = runOnFXThread(() -> Clipboard.getSystemClipboard().getImage());
        return (fxImage == null) ? null : SwingFXUtils.fromFXImage(fxImage, null);
    }

    private static <T> T runOnFXThread(Callable<T> callable) {
        try {
            if (Platform.isFxApplicationThread())
                return callable.call();

            FutureTask<T> task = new FutureTask<>(callable);
            Platform.runLater(task);
            return task.get();
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
}
