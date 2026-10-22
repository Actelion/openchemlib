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

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.clipboard.ClipboardHandler;
import com.actelion.research.gui.hidpi.HiDPIHelper;

import javax.swing.*;
import java.awt.*;

/**
 * Opens a Swing window with JStructureViews to test drag & drop and copy/paste.
 * The window contains a read-only view (drag/copy source only), an editable molecule view and an
 * editable fragment view with query features. Below every view the idcode of its current molecule is shown.
 * To test the exchange with JavaFX, run FXStructureViewTestApp in parallel.
 */
public class JStructureViewTestApp {
	private static final String MOLECULE_SMILES = "CN1CCC[C@H]1c1cccnc1";
	private static final String FRAGMENT_SMARTS = "c1ccc([N;!H0])cc1[F,Cl,Br]";
	private static final int VIEW_WIDTH = 320;
	private static final int VIEW_HEIGHT = 180;

	public static void main(String[] args) {
		SwingUtilities.invokeLater(() -> new JStructureViewTestApp().createWindow());
		}

	private void createWindow() {
		ClipboardHandler clipboardHandler = new ClipboardHandler();

		JStructureView sourceView = new JStructureView(createMolecule());
		sourceView.setClipboardHandler(clipboardHandler);
		sourceView.setSelectable(true);

		JEditableStructureView editableView = new JEditableStructureView();
		editableView.setClipboardHandler(clipboardHandler);

		JEditableStructureView fragmentView = new JEditableStructureView(createFragment());
		fragmentView.setClipboardHandler(clipboardHandler);

		JPanel content = new JPanel();
		content.setLayout(new BoxLayout(content, BoxLayout.Y_AXIS));
		int gap = HiDPIHelper.scale(8);
		content.setBorder(BorderFactory.createEmptyBorder(gap, gap, gap, gap));
		addView(content, "Read-only molecule (drag/copy source, lasso selectable)", sourceView);
		addView(content, "Editable molecule (drop/paste target, double click to edit)", editableView);
		addView(content, "Editable fragment with query features", fragmentView);

		JFrame frame = new JFrame("Swing Structure Views");
		frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		frame.getContentPane().add(content);
		frame.pack();
		frame.setLocation(HiDPIHelper.scale(VIEW_WIDTH + 120), 40);	// next to the window of FXStructureViewTestApp
		frame.setVisible(true);
		}

	private void addView(JPanel content, String title, JStructureView view) {
		Dimension size = new Dimension(HiDPIHelper.scale(VIEW_WIDTH), HiDPIHelper.scale(VIEW_HEIGHT));
		view.setPreferredSize(size);
		view.setMaximumSize(size);
		view.setBorder(BorderFactory.createLineBorder(Color.GRAY));
		JLabel idcodeLabel = new JLabel(getIDCodeText(view.getMolecule()));
		idcodeLabel.setPreferredSize(new Dimension(size.width, idcodeLabel.getPreferredSize().height));
		view.addStructureListener(mol -> idcodeLabel.setText(getIDCodeText(mol)));
		JLabel titleLabel = new JLabel(title);
		titleLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
		view.setAlignmentX(Component.LEFT_ALIGNMENT);
		idcodeLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
		content.add(titleLabel);
		content.add(view);
		content.add(idcodeLabel);
		content.add(Box.createVerticalStrut(HiDPIHelper.scale(4)));
		}

	private static StereoMolecule createMolecule() {
		return new SmilesParser().parseMolecule(MOLECULE_SMILES);
		}

	private static StereoMolecule createFragment() {
		StereoMolecule fragment = new StereoMolecule();
		try {
			new SmilesParser(SmilesParser.SMARTS_MODE_IS_SMARTS).parse(fragment, FRAGMENT_SMARTS);
			}
		catch (Exception e) {
			e.printStackTrace();
			}
		fragment.setFragment(true);
		return fragment;
		}

	private static String getIDCodeText(StereoMolecule mol) {
		if (mol == null || mol.getAllAtoms() == 0)
			return "<empty>";
		return (mol.isFragment() ? "fragment: " : "molecule: ") + mol.getIDCode();
		}
	}
