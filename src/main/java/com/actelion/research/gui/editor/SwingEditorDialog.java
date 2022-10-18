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

package com.actelion.research.gui.editor;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.StructureListener;
import com.actelion.research.gui.hidpi.HiDPIHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

public class SwingEditorDialog extends JDialog implements ActionListener {
	static final long serialVersionUID = 0x20211103;

	private static final String DEFAULT_MOLECULE_TITLE = "Structure Editor";
	private static final String DEFAULT_REACTION_TITLE = "Reaction Editor";

	private StereoMolecule mMolecule;
	private SwingEditorToolbar mToolBar;
	private SwingEditorArea mArea;
	private JPanel         mButtonPanel;
	private boolean        mIsCancelled;
	private ArrayList<StructureListener> mListener;

	/**
	 * Creates a modal chemical editor dialog to edit a single molecule,
	 * which may, of course, consist of multiple disconnected fragments.
	 * Query features can be edited, if the passed mol's fragment bit is true.
	 * @param owner
	 * @param mol
	 */
	public SwingEditorDialog(Dialog owner, StereoMolecule mol) {
		this(owner, mol, Dialog.DEFAULT_MODALITY_TYPE);
	}

	/**
	 * Creates a modal chemical editor dialog to edit a chemical reaction.
	 * Query features can be edited, if the passed rxn's fragment bit is true.
	 * @param owner
	 * @param rxn
	 */
	public SwingEditorDialog(Dialog owner, Reaction rxn) {
		this(owner, rxn, Dialog.DEFAULT_MODALITY_TYPE);
	}

	/**
	 * Creates a chemical editor dialog to edit a single molecule,
	 * which may, of course, consist of multiple disconnected fragments.
	 * Query features can be edited, if the passed mol's fragment bit is true.
	 * @param owner
	 * @param mol
	 * @param modalityType
	 */
	public SwingEditorDialog(Dialog owner, StereoMolecule mol, ModalityType modalityType) {
		super(owner, DEFAULT_MOLECULE_TITLE, modalityType);
		mMolecule = (mol == null) ? new StereoMolecule() : new StereoMolecule(mol);
		initialize(owner, 0);
	}

	/**
	 * Creates a modal chemical editor dialog to edit multiple molecules in one editor pane.
	 * Each of these molecule may consist of multiple disconnected fragments.
	 * Atoms connected by bonds or being in close vicinity are recognized to belong to the
	 * same molecule, while more distant fragments are perceived as separated molecules.
	 * Query features can be edited, if the passed mols' fragment bits are true.
	 * @param owner
	 * @param mol
	 * @param modalityType
	 */
	public SwingEditorDialog(Dialog owner, StereoMolecule[] mol, ModalityType modalityType) {
		super(owner, DEFAULT_MOLECULE_TITLE, modalityType);
		mMolecule = new StereoMolecule();
		initialize(owner, GenericEditorArea.MODE_MULTIPLE_FRAGMENTS);
		if (mol != null)
			mArea.getGenericDrawArea().setFragments(mol);
	}

	/**
	 * Creates a modal chemical editor dialog to edit a chemical reaction.
	 * Query features can be edited, if the passed rxn's fragment bit is true.
	 * @param owner
	 * @param rxn
	 * @param modalityType
	 */
	public SwingEditorDialog(Dialog owner, Reaction rxn, ModalityType modalityType) {
		super(owner, DEFAULT_REACTION_TITLE, modalityType);
		mMolecule = new StereoMolecule();
		initialize(owner, GenericEditorArea.MODE_REACTION);
		if (rxn != null)
			mArea.getGenericDrawArea().setReaction(rxn);
	}

	/**
	 * Creates a modal chemical editor dialog to edit a single molecule,
	 * which may, of course, consist of multiple disconnected fragments.
	 * Query features can be edited, if the passed mol's fragment bit is true.
	 * @param owner
	 * @param mol
	 */
	public SwingEditorDialog(Frame owner, StereoMolecule mol) {
		this(owner, mol, Dialog.DEFAULT_MODALITY_TYPE);
	}

	/**
	 * Creates a modal chemical editor dialog to edit a chemical reaction.
	 * Query features can be edited, if the passed rxn's fragment bit is true.
	 * @param owner
	 * @param rxn
	 */
	public SwingEditorDialog(Frame owner, Reaction rxn) {
		this(owner, rxn, Dialog.DEFAULT_MODALITY_TYPE);
	}

	/**
	 * Creates a chemical editor dialog to edit a single molecule,
	 * which may, of course, consist of multiple disconnected fragments.
	 * Query features can be edited, if the passed mol's fragment bit is true.
	 * @param owner
	 * @param mol
	 * @param modalityType
	 */
	public SwingEditorDialog(Frame owner, StereoMolecule mol, ModalityType modalityType) {
		super(owner, DEFAULT_MOLECULE_TITLE, modalityType);
		mMolecule = (mol == null) ? new StereoMolecule() : new StereoMolecule(mol);
		initialize(owner, 0);
	}

	/**
	 * Creates a modal chemical editor dialog to edit multiple molecules in one editor pane.
	 * Each of these molecule may consist of multiple disconnected fragments.
	 * Atoms connected by bonds or being in close vicinity are recognized to belong to the
	 * same molecule, while more distant fragments are perceived as separated molecules.
	 * Query features can be edited, if the passed mols' fragment bits are true.
	 * @param owner
	 * @param mol
	 * @param modalityType
	 */
	public SwingEditorDialog(Frame owner, StereoMolecule[] mol, ModalityType modalityType) {
		super(owner, DEFAULT_MOLECULE_TITLE, modalityType);
		mMolecule = new StereoMolecule();
		initialize(owner, GenericEditorArea.MODE_MULTIPLE_FRAGMENTS);
		if (mol != null)
			mArea.getGenericDrawArea().setFragments(mol);
	}

	/**
	 * Creates a modal chemical editor dialog to edit a chemical reaction.
	 * Query features can be edited, if the passed rxn's fragment bit is true.
	 * @param owner
	 * @param rxn
	 * @param modalityType
	 */
	public SwingEditorDialog(Frame owner, Reaction rxn, ModalityType modalityType) {
		super(owner, DEFAULT_REACTION_TITLE, modalityType);
		mMolecule = new StereoMolecule();
		initialize(owner, GenericEditorArea.MODE_REACTION);
		if (rxn != null)
			mArea.getGenericDrawArea().setReaction(rxn);
	}

	private void initialize(Component owner, int mode) {
		mArea = new SwingEditorArea(mMolecule, mode);
		mArea.setPreferredSize(new Dimension(HiDPIHelper.scale(mode == GenericEditorArea.MODE_REACTION ? 800 : 480), HiDPIHelper.scale(300)));
		getContentPane().add(mArea, BorderLayout.CENTER);

		mToolBar = new SwingEditorToolbar(mArea);
		getContentPane().add(mToolBar, BorderLayout.WEST);

		mButtonPanel = new JPanel();
		mButtonPanel.setLayout(new BorderLayout());
		JPanel ibp = new JPanel();
		ibp.setLayout(new GridLayout(1, 2, 8, 0));
		JButton bcancel = new JButton("Cancel");
		bcancel.addActionListener(this);
		ibp.add(bcancel);
		JButton bok = new JButton("OK");
		bok.addActionListener(this);
		ibp.add(bok);
		mButtonPanel.add(ibp, BorderLayout.EAST);
		JButton bhelp = new JButton("Help");
		bhelp.addActionListener(this);
		mButtonPanel.add(bhelp, BorderLayout.WEST);
		mButtonPanel.setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));
		getContentPane().add(mButtonPanel, BorderLayout.SOUTH);

		mListener = new ArrayList<>();

		addKeyListener(mArea.getKeyHandler());

		pack();
		setLocationRelativeTo(owner);

		mIsCancelled = true;
	}

	public void addStructureListener(StructureListener listener) {
		mListener.add(listener);
	}

	public GenericEditorArea getDrawArea() {
		return mArea.getGenericDrawArea();
	}

	public void setAccessory(Component accessory) {
		mButtonPanel.add(accessory, BorderLayout.NORTH);
		pack();
	}

	public StereoMolecule getStructure() {
		return mMolecule;
	}

	/**
	 * @return mapped reaction with absolute coordinates, but without drawing objects
	 */
	public Reaction getReaction() {
		return mArea.getGenericDrawArea().getReaction();
	}

	/**
	 * @return mapped reaction with absolute coordinates and drawing objects
	 */
	public Reaction getReactionAndDrawings() {
		return mArea.getGenericDrawArea().getReactionAndDrawings();
	}

	public boolean isCancelled() {
		return mIsCancelled;
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("Help")) {
			mArea.getGenericDrawArea().showHelpDialog();
			return;
		}

		if (e.getActionCommand().equals("OK")) {
			for (StructureListener sl : mListener)
				sl.structureChanged(mMolecule);
			mIsCancelled = false;
		}

		dispose();
	}
}
