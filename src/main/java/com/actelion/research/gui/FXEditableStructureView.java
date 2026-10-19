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

package com.actelion.research.gui;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.editor.FXEditorDialog;
import com.actelion.research.gui.editor.GenericEditorArea;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.input.MouseButton;
import javafx.scene.input.MouseEvent;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.stage.Modality;
import javafx.stage.Window;

/**
 * JavaFX equivalent of JEditableStructureView: an editable FXStructureView, which opens
 * a structure editor dialog on double click.
 */
public class FXEditableStructureView extends FXStructureView {
	private static final String EDIT_MESSAGE = "<double click or drag & drop>";
	private static final double EDIT_MESSAGE_FONT_SIZE = 10;

	private boolean mAllowQueryFeatures;
	private int mAllowedPseudoAtoms;

	public FXEditableStructureView() {
		this(null);
		}

	public FXEditableStructureView(StereoMolecule mol) {
		super(mol);
		initialize();
		}

	public FXEditableStructureView(int dragAction, int dropAction) {
		this(null, dragAction, dropAction);
		}

	public FXEditableStructureView(StereoMolecule mol, int dragAction, int dropAction) {
		super(mol, dragAction, dropAction);
		initialize();
		}

	private void initialize() {
		setEditable(true);
		mAllowedPseudoAtoms = GenericEditorArea.DEFAULT_ALLOWED_PSEUDO_ATOMS;
		mAllowQueryFeatures = true;

		addEventHandler(MouseEvent.MOUSE_CLICKED, e -> {
			if (e.getButton() == MouseButton.PRIMARY && e.getClickCount() == 2 && !isDisabled() && isEditable()) {
				FXEditorDialog theDialog = createDrawDialog();
				theDialog.getDrawArea().setAllowedPseudoAtoms(mAllowedPseudoAtoms);
				theDialog.getDrawArea().setAllowQueryFeatures(mAllowQueryFeatures);
				theDialog.getDrawArea().setDisplayMode(getDisplayMode());
				theDialog.addStructureListener(this);
				theDialog.show();
				}
			} );
		}

	@Override
	protected void drawOverlay(GraphicsContext gc, double width, double height) {
		if (!isDisabled() && isEditable() && getMolecule().getAllAtoms() == 0) {
			Font font = Font.font(Font.getDefault().getFamily(), FontWeight.NORMAL, EDIT_MESSAGE_FONT_SIZE);
			gc.setFont(font);
			gc.setFill(toFXColor(getForegroundRGB()));
			double textWidth = getTextWidth(EDIT_MESSAGE, font);
			gc.fillText(EDIT_MESSAGE, (width-textWidth)/2, (height+EDIT_MESSAGE_FONT_SIZE)/2);
			}
		}

	protected FXEditorDialog createDrawDialog() {
		Window owner = (getScene() == null) ? null : getScene().getWindow();
		return new FXEditorDialog(owner, getMolecule(), owner == null ? Modality.APPLICATION_MODAL : Modality.WINDOW_MODAL);
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
