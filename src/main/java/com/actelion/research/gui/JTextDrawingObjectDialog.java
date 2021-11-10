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

import com.actelion.research.chem.TextDrawingObject;
import info.clearthought.layout.TableLayout;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

@Deprecated
public class JTextDrawingObjectDialog extends JDialog implements ActionListener {
    static final long serialVersionUID = 0x20110325;

    private static final int[] TEXT_STYLE = { Font.PLAIN, Font.ITALIC, Font.BOLD, Font.ITALIC | Font.BOLD};
	private static final String[] TEXT_STYLE_LIST = { "plain", "italic", "bold", "italics & bold" };
	private static final String[] TEXT_SIZE_LIST = { "8", "9", "10", "12", "14", "18", "24", "32" };

	private JTextArea			mTextArea;
	private TextDrawingObject	mTextObject;
	private JComboBox			mComboBoxTextSize,mComboBoxStyle;

	protected JTextDrawingObjectDialog(Dialog owner, TextDrawingObject textObject) {
		super(owner, true);
		mTextObject = textObject;
		init(owner);
		}


	protected JTextDrawingObjectDialog(Frame owner, TextDrawingObject textObject) {
		super(owner, true);
		mTextObject = textObject;
		init(owner);
		}


	private void init(Component owner) {
        mComboBoxTextSize = new JComboBox(TEXT_SIZE_LIST);
		mComboBoxTextSize.setEditable(true);
		mComboBoxTextSize.setSelectedItem(""+(int)mTextObject.getSize());

        mComboBoxStyle = new JComboBox(TEXT_STYLE_LIST);
		int styleIndex = 0;
		for (int i=0; i<TEXT_STYLE.length; i++) {
			if (mTextObject.getStyle() == TEXT_STYLE[i]) {
				styleIndex = i;
				break;
				}
			}
		mComboBoxStyle.setSelectedIndex(styleIndex);

		JComponent menuPanel = new JPanel();
		double[][] size = { { 8, TableLayout.PREFERRED, 4, TableLayout.PREFERRED, 8 },
							{ 8, TableLayout.PREFERRED, 4, TableLayout.PREFERRED, 8 } };
		menuPanel.setLayout(new TableLayout(size));

		menuPanel.add(new JLabel("Text size:"), "1,1");
		menuPanel.add(mComboBoxTextSize, "3,1");
		menuPanel.add(new JLabel("Text style:"), "1,3");
		menuPanel.add(mComboBoxStyle, "3,3");

		mTextArea = new JTextArea(mTextObject.getText(), 3, 20);
        mTextArea.setBorder(BorderFactory.createEtchedBorder());
		JPanel textPanel = new JPanel();
        textPanel.setBorder(BorderFactory.createEmptyBorder(0,4,0,4));
		textPanel.add(mTextArea);

		JButton buttonOK = new JButton("OK");
		buttonOK.addActionListener(this);
		JButton buttonCancel = new JButton("Cancel");
		buttonCancel.addActionListener(this);
		JPanel innerButtonPanel = new JPanel();
        innerButtonPanel.setLayout(new GridLayout(1,2,8,0));
        innerButtonPanel.add(buttonCancel);
        innerButtonPanel.add(buttonOK);
        innerButtonPanel.setBorder(BorderFactory.createEmptyBorder(8,8,8,8));
		JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new BorderLayout());
        buttonPanel.add(innerButtonPanel, BorderLayout.EAST);

		getContentPane().setLayout(new BorderLayout());
		getContentPane().add(menuPanel,  BorderLayout.NORTH);
		getContentPane().add(textPanel,  BorderLayout.CENTER);
		getContentPane().add(buttonPanel,  BorderLayout.SOUTH);

		getRootPane().setDefaultButton(buttonOK);

		pack();
		setLocationRelativeTo(owner);

		mTextArea.requestFocus();
		setVisible(true);
		}


	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand() == "OK") {
			int textSize = 0;
			try {
				textSize = Integer.parseInt((String)mComboBoxTextSize.getSelectedItem());
				}
			catch (NumberFormatException nfe) {
				JOptionPane.showMessageDialog(this, "Illegal text size.");
				return;
				}
			int textStyle = TEXT_STYLE[mComboBoxStyle.getSelectedIndex()];
			mTextObject.setValues(mTextArea.getText(), textSize, textStyle);
			}

		dispose();
		}
	}
