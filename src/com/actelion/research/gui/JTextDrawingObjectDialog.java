/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import info.clearthought.layout.TableLayout;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dialog;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;

import com.actelion.research.chem.TextDrawingObject;

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
