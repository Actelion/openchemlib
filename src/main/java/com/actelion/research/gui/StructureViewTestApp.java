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
import com.actelion.research.gui.clipboard.FXClipboardHandler;
import com.actelion.research.gui.clipboard.JClipboardHandler;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.control.Label;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

/**
 * Opens a JavaFX window with FXStructureViews and a Swing window with JStructureViews
 * to test drag & drop and copy/paste within and between both toolkits.
 * Each window contains a read-only view (drag/copy source only), an editable molecule view and an
 * editable fragment view with query features. Below every view the idcode of its current molecule is shown.
 * This class intentionally doesn't extend Application: the Java launcher refuses to start Application
 * subclasses with "JavaFX runtime components are missing", if JavaFX is on the classpath rather than the module path.
 */
public class StructureViewTestApp {
	private static final String MOLECULE_SMILES = "CN1CCC[C@H]1c1cccnc1";
	private static final String FRAGMENT_SMARTS = "c1ccc([N;!H0])cc1[F,Cl,Br]";
	private static final int VIEW_WIDTH = 320;
	private static final int VIEW_HEIGHT = 180;

	public static void main(String[] args) {
		Application.launch(FXApp.class, args);
		}

	public static class FXApp extends Application {
		@Override
		public void start(Stage stage) {
			Platform.setImplicitExit(false);
			StructureViewTestApp app = new StructureViewTestApp();
			app.createFXWindow(stage);
			SwingUtilities.invokeLater(() -> app.createSwingWindow());
			}
		}

	private void createFXWindow(Stage stage) {
		FXClipboardHandler clipboardHandler = new FXClipboardHandler();

		FXStructureView sourceView = new FXStructureView(createMolecule());
		sourceView.setClipboardHandler(clipboardHandler);
		sourceView.setSelectable(true);

		FXEditableStructureView editableView = new FXEditableStructureView();
		editableView.setClipboardHandler(clipboardHandler);

		FXEditableStructureView fragmentView = new FXEditableStructureView(createFragment());
		fragmentView.setClipboardHandler(clipboardHandler);

		VBox content = new VBox(4);
		content.setPadding(new Insets(8));
		addFXView(content, "Read-only molecule (drag/copy source, lasso selectable)", sourceView);
		addFXView(content, "Editable molecule (drop/paste target, double click to edit)", editableView);
		addFXView(content, "Editable fragment with query features", fragmentView);

		stage.setTitle("JavaFX Structure Views");
		stage.setScene(new Scene(content));
		stage.setX(40);
		stage.setY(40);
		stage.setOnCloseRequest(e -> exit());
		stage.show();
		}

	private void addFXView(VBox content, String title, FXStructureView view) {
		view.setPrefSize(VIEW_WIDTH, VIEW_HEIGHT);
		view.setStyle("-fx-border-color: #A0A0A0;");
		Label idcodeLabel = new Label();
		idcodeLabel.setMaxWidth(VIEW_WIDTH);
		idcodeLabel.setText(getIDCodeText(view.getMolecule()));
		view.addStructureListener(mol -> idcodeLabel.setText(getIDCodeText(mol)));
		content.getChildren().addAll(new Label(title), view, idcodeLabel);
		}

	private void createSwingWindow() {
		JClipboardHandler clipboardHandler = new JClipboardHandler();

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
		addSwingView(content, "Read-only molecule (drag/copy source, lasso selectable)", sourceView);
		addSwingView(content, "Editable molecule (drop/paste target, double click to edit)", editableView);
		addSwingView(content, "Editable fragment with query features", fragmentView);

		JFrame frame = new JFrame("Swing Structure Views");
		frame.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
		frame.addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				exit();
				}
			} );
		frame.getContentPane().add(content);
		frame.pack();
		frame.setLocation(HiDPIHelper.scale(VIEW_WIDTH + 120), 40);
		frame.setVisible(true);
		}

	private void addSwingView(JPanel content, String title, JStructureView view) {
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

	private static void exit() {
		Platform.exit();
		System.exit(0);
		}
	}
