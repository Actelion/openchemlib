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
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.clipboard.ClipboardHandler;
import com.actelion.research.gui.hidpi.HiDPIHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.util.ArrayList;

@Deprecated
public class JDrawDialog extends JDialog implements ActionListener,KeyListener {
    static final long serialVersionUID = 0x20061019;

    private static final String DEFAULT_TITLE = "OSIRIS Structure Editor";

	private StereoMolecule mMolecule;
	protected JDrawToolbar   mToolBar;
	protected JDrawArea      mArea;
	private JPanel         mButtonPanel;
	private boolean        mIsCancelled;
	private ArrayList<StructureListener> mListener;

	public JDrawDialog(Dialog owner, StereoMolecule mol) {
		this(owner, mol, DEFAULT_TITLE);
		}

	public JDrawDialog(Dialog owner, StereoMolecule mol, ModalityType modalityType) {
		this(owner, mol, DEFAULT_TITLE, modalityType);
		}

	public JDrawDialog(Dialog owner, StereoMolecule mol, String title) {
		this(owner, mol, title, Dialog.DEFAULT_MODALITY_TYPE);
		}

	public JDrawDialog(Dialog owner, StereoMolecule mol, String title, ModalityType modalityType) {
		super(owner, title, modalityType);
		mMolecule = (mol == null) ? new StereoMolecule() : new StereoMolecule(mol);
		initialize(owner, 0);
		}

	public JDrawDialog(Dialog owner, StereoMolecule[] mol, String title, ModalityType modalityType) {
		super(owner, title, modalityType);
		mMolecule = new StereoMolecule();
		initialize(owner, JDrawArea.MODE_MULTIPLE_FRAGMENTS);
		if (mol != null)
			mArea.setFragments(mol);
		}

	public JDrawDialog(Dialog owner, Reaction rxn, String title, ModalityType modalityType) {
		super(owner, title, modalityType);
		mMolecule = new StereoMolecule();
		initialize(owner, JDrawArea.MODE_REACTION);
		if (rxn != null)
			mArea.setReaction(rxn);
		}

	public JDrawDialog(Frame owner) {
		this(owner, false, DEFAULT_TITLE);
		}

	public JDrawDialog(Frame owner, boolean isFragment) {
		this(owner, isFragment, DEFAULT_TITLE);
		}

	public JDrawDialog(Frame owner, boolean isFragment, String title) {
		super(owner, title, true);
		mMolecule = new StereoMolecule();
		mMolecule.setFragment(isFragment);
		initialize(owner, 0);
		}

	public JDrawDialog(Frame owner, StereoMolecule mol) {
		this(owner, mol, DEFAULT_TITLE);
		}

	public JDrawDialog(Frame owner, StereoMolecule mol, ModalityType modalityType) {
		this(owner, mol, DEFAULT_TITLE, modalityType);
		}

	public JDrawDialog(Frame owner, StereoMolecule mol, String title) {
		this(owner, mol, title, Dialog.DEFAULT_MODALITY_TYPE);
		}

	public JDrawDialog(Frame owner, StereoMolecule mol, String title, ModalityType modalityType) {
		super(owner, title, modalityType);
		mMolecule = (mol == null) ? new StereoMolecule() : new StereoMolecule(mol);
		initialize(owner, 0);
		}

	public JDrawDialog(Frame owner, StereoMolecule[] mol) {
		this(owner, mol, DEFAULT_TITLE);
		}

	public JDrawDialog(Frame owner, StereoMolecule[] mol, ModalityType modalityType) {
		this(owner, mol, DEFAULT_TITLE, modalityType);
		}

	public JDrawDialog(Frame owner, StereoMolecule[] mol, String title) {
		this(owner, mol, title, Dialog.DEFAULT_MODALITY_TYPE);
		}

	public JDrawDialog(Frame owner, StereoMolecule[] mol, String title, ModalityType modalityType) {
		super(owner, title, modalityType);
		mMolecule = new StereoMolecule();
		initialize(owner, JDrawArea.MODE_MULTIPLE_FRAGMENTS);
		if (mol != null)
			mArea.setFragments(mol);
		}

	public JDrawDialog(Frame owner, Reaction rxn) {
		this(owner, rxn, DEFAULT_TITLE);
		}

	public JDrawDialog(Frame owner, Reaction rxn, ModalityType modalityType) {
		this(owner, rxn, DEFAULT_TITLE, modalityType);
		}

	public JDrawDialog(Frame owner, Reaction rxn, String title) {
		this(owner, rxn, title, Dialog.DEFAULT_MODALITY_TYPE);
		}

	public JDrawDialog(Frame owner, Reaction rxn, String title, ModalityType modalityType) {
		super(owner, title, modalityType);
		mMolecule = new StereoMolecule();
		initialize(owner, JDrawArea.MODE_REACTION);
		if (rxn != null)
			mArea.setReaction(rxn);
		}

	private void initialize(Component owner, int mode) {
        addKeyListener(this);

		mArea = new JDrawArea(mMolecule, mode);
		mArea.setClipboardHandler(new ClipboardHandler());
		mArea.setPreferredSize(new Dimension(HiDPIHelper.scale(mode == JDrawArea.MODE_REACTION ? 800 : 480), HiDPIHelper.scale(300)));
		getContentPane().add(mArea, BorderLayout.CENTER);

		mToolBar = new JDrawToolbar(mArea, mode);
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

		pack();
		setLocationRelativeTo(owner);

		mIsCancelled = true;
		}

	public void addStructureListener(StructureListener listener) {
		mListener.add(listener);
		}

	public JDrawArea getDrawArea() {
		return mArea;
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
		return mArea.getReaction();
		}

    /**
     * @return mapped reaction with absolute coordinates and drawing objects
     */
	public Reaction getReactionAndDrawings() {
		return mArea.getReactionAndDrawings();
		}

	public boolean isCancelled() {
	    return mIsCancelled;
	    }

	public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand().equals("Help")) {
            mArea.showHelpDialog();
            return;
            }

        if (e.getActionCommand().equals("OK")) {
	        for (StructureListener sl : mListener)
		        sl.structureChanged(mMolecule);
	        mIsCancelled = false;
            }

	    dispose();
	    }

	public void keyPressed(KeyEvent e) {
        mArea.keyPressed(e);
        }

	public void keyReleased(KeyEvent e) {
		mArea.keyReleased(e);
		}

	public void keyTyped(KeyEvent e) {}
	}
