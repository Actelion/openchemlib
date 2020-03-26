/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
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

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.IReactionMapper;
import com.actelion.research.chem.reaction.MCSReactionMapper;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.clipboard.ClipboardHandler;
import com.actelion.research.gui.hidpi.HiDPIHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;

public class JDrawDialog extends JDialog implements ActionListener,KeyListener {
    static final long serialVersionUID = 0x20061019;

    private static final String DEFAULT_TITLE = "OSIRIS Structure Editor";

	private StereoMolecule mMolecule;
	private JDrawToolbar   mToolBar;
	private JDrawArea      mArea;
	private JPanel         mButtonPanel;
	private boolean        mIsCancelled;
	private ArrayList<StructureListener> mListener;

	public JDrawDialog(Dialog owner, StereoMolecule mol) {
		this(owner, mol, DEFAULT_TITLE);
		}

	public JDrawDialog(Dialog owner, StereoMolecule mol, String title) {
		super(owner, title, true);
		mMolecule = (mol == null) ? new StereoMolecule() : new StereoMolecule(mol);
		initialize(owner, 0);
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

	public JDrawDialog(Frame owner, StereoMolecule mol, String title) {
		super(owner, title, true);
		mMolecule = (mol == null) ? new StereoMolecule() : new StereoMolecule(mol);
		initialize(owner, 0);
		}

	public JDrawDialog(Frame owner, StereoMolecule[] mol) {
		this(owner, mol, DEFAULT_TITLE);
	}

	public JDrawDialog(Frame owner, StereoMolecule[] mol, String title) {
		super(owner, title, true);
		mMolecule = new StereoMolecule();
		initialize(owner, JDrawArea.MODE_REACTION);
		if (mol != null)
			mArea.setFragments(mol);
		}

	public JDrawDialog(Frame owner, Reaction rxn) {
		this(owner, rxn, DEFAULT_TITLE);
		}

	public JDrawDialog(Frame owner, Reaction rxn, String title) {
		super(owner, title, true);
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

		mListener = new ArrayList<StructureListener>();

        mArea.setReactionMapper(new MCSReactionMapper());
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

	public void setReactionMapper(IReactionMapper mapper) {
		mArea.setReactionMapper(mapper);
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
