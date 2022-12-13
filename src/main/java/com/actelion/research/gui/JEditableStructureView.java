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
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.editor.GenericEditorArea;
import com.actelion.research.gui.editor.SwingEditorDialog;
import com.actelion.research.gui.hidpi.HiDPIHelper;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;

public class JEditableStructureView extends JStructureView {
    static final long serialVersionUID = 0x20090727;

    private static final String EDIT_MESSAGE = "<double click or drag & drop>";
    private boolean mAllowQueryFeatures;
    private int mAllowedPseudoAtoms;

    public JEditableStructureView() {
        this(null);
		}

	public JEditableStructureView(StereoMolecule mol) {
        super(mol);
		setEditable(true);
		mAllowedPseudoAtoms = GenericEditorArea.DEFAULT_ALLOWED_PSEUDO_ATOMS;
		mAllowQueryFeatures = true;
	    }

	public JEditableStructureView(int dragAction, int dropAction) {
        this(null, dragAction, dropAction);
	    }

	public JEditableStructureView(StereoMolecule mol, int dragAction, int dropAction) {
        super(mol, dragAction, dropAction);
		setEditable(true);
		mAllowedPseudoAtoms = GenericEditorArea.DEFAULT_ALLOWED_PSEUDO_ATOMS;
		mAllowQueryFeatures = true;
	    }

	@Override
	public void paintComponent(Graphics g) {
        Dimension theSize = getSize();
        Insets insets = getInsets();
        theSize.width -= insets.left + insets.right;
        theSize.height -= insets.top + insets.bottom;

		super.paintComponent(g);

		if (isEnabled() && isEditable() && getMolecule().getAllAtoms() == 0) {
	        g.setFont(g.getFont().deriveFont(Font.PLAIN, HiDPIHelper.scale(10)));
	        FontMetrics metrics = g.getFontMetrics();
	        Rectangle2D bounds = metrics.getStringBounds(EDIT_MESSAGE, g);
	        g.drawString(EDIT_MESSAGE, (int)(insets.left+theSize.width-bounds.getWidth())/2,
	                                   (insets.top+theSize.height-metrics.getHeight())/2+metrics.getAscent());
	        }
	    }

    public void mouseClicked(MouseEvent e) {
        if (e.getClickCount() == 2 && isEnabled() && isEditable()) {
            SwingEditorDialog theDialog = createDrawDialog();
            theDialog.getDrawArea().setAllowedPseudoAtoms(mAllowedPseudoAtoms);
	        theDialog.getDrawArea().setAllowQueryFeatures(mAllowQueryFeatures);
            theDialog.getDrawArea().setDisplayMode(getDisplayMode());
            theDialog.addStructureListener(this);
            theDialog.setVisible(true);
            }
        }

    protected SwingEditorDialog createDrawDialog() {
		Component c = this;
		while (!(c instanceof Frame || c instanceof Dialog))
			c = c.getParent();
		return (c instanceof Frame) ? new SwingEditorDialog((Frame) c, getMolecule(), Dialog.ModalityType.DOCUMENT_MODAL) : new SwingEditorDialog((Dialog) c, getMolecule(), Dialog.ModalityType.DOCUMENT_MODAL);
		}

	public void setAllowedPseudoAtoms(int apa) {
		mAllowedPseudoAtoms = apa;
		}

	public void setAllowQueryFeatures(boolean allow) {
		if (mAllowQueryFeatures != allow) {
			mAllowQueryFeatures = allow;
			if (!allow && getMolecule().removeQueryFeatures())
				structureChanged();
			}
		}
    }
