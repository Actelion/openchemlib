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

package com.actelion.research.gui.editor;

import com.actelion.research.chem.TextDrawingObject;
import com.actelion.research.gui.generic.*;

import java.awt.*;

public class TextDrawingObjectDialogBuilder implements GenericEventListener<GenericActionEvent> {
    static final long serialVersionUID = 0x20110325;

    private static final int[] TEXT_STYLE = { Font.PLAIN, Font.ITALIC, Font.BOLD, Font.ITALIC | Font.BOLD};
	private static final String[] TEXT_STYLE_LIST = { "plain", "italic", "bold", "italics & bold" };
	private static final String[] TEXT_SIZE_LIST = { "8", "9", "10", "12", "14", "18", "24", "32" };

	private GenericDialog       mDialog;
	private GenericTextField    mTextArea;
	private TextDrawingObject	mTextObject;
	private GenericComboBox     mComboBoxTextSize,mComboBoxStyle;
	private boolean             mOKSelected;

	public TextDrawingObjectDialogBuilder(GenericUIHelper dialogHelper, TextDrawingObject textObject) {
		mDialog = dialogHelper.createDialog("Edit Text", this);
		mTextObject = textObject;
		build();
		}

	/**
	 * @return true if OK was pressed and potential change was applied to molecule
	 */
	public boolean showDialog() {
		mOKSelected = false;
		mDialog.showDialog();
		return mOKSelected;
		}

	private void build() {
        mComboBoxTextSize = mDialog.createComboBox();
        for (String item:TEXT_SIZE_LIST)
        	mComboBoxTextSize.addItem(item);
		mComboBoxTextSize.setEditable(true);
		mComboBoxTextSize.setSelectedItem(""+(int)mTextObject.getSize());

        mComboBoxStyle = mDialog.createComboBox();
        for (String item:TEXT_STYLE_LIST)
        	mComboBoxStyle.addItem(item);
		int styleIndex = 0;
		for (int i=0; i<TEXT_STYLE.length; i++) {
			if (mTextObject.getStyle() == TEXT_STYLE[i]) {
				styleIndex = i;
				break;
				}
			}
		mComboBoxStyle.setSelectedIndex(styleIndex);

		int[] hLayout = { 8, GenericDialog.PREFERRED, 4, GenericDialog.PREFERRED, 8 };
		int[] vLayout = { 8, GenericDialog.PREFERRED, 4, GenericDialog.PREFERRED, 8, GenericDialog.PREFERRED, 8 };
		mDialog.setLayout(hLayout, vLayout);

		mDialog.add(mDialog.createLabel("Text size:"), 1,1);
		mDialog.add(mComboBoxTextSize, 3,1);
		mDialog.add(mDialog.createLabel("Text style:"), 1,3);
		mDialog.add(mComboBoxStyle, 3,3);

		mTextArea = mDialog.createTextField(20, 3);
		mTextArea.setText(mTextObject.getText());
		mDialog.add(mTextArea, 1,5,3,5);
		}

	@Override
	public void eventHappened(GenericActionEvent e) {
		if (e.getWhat() == GenericActionEvent.WHAT_OK) {
			int textSize;
			try {
				textSize = Integer.parseInt(mComboBoxTextSize.getSelectedItem());
				}
			catch (NumberFormatException nfe) {
				mDialog.showMessage("Illegal text size.");
				return;
				}
			int textStyle = TEXT_STYLE[mComboBoxStyle.getSelectedIndex()];
			mTextObject.setValues(mTextArea.getText(), textSize, textStyle);
			mOKSelected = true;
			}

		mDialog.disposeDialog();
		}
	}
