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

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import com.actelion.research.chem.*;

public class JBondQueryFeatureDialog extends JDialog implements ActionListener {
    static final long serialVersionUID = 0x20070822;

    private ExtendedMolecule	mMol;
	private int					mBond,mFirstSpanItem;
	private JCheckBox			mCBSingle,mCBDouble,mCBTriple,mCBDelocalized,mCBIsBridge,mCBMatchStereo;
	private JComboBox			mComboBoxRing,mComboBoxRingSize,mComboBoxMinAtoms,mComboBoxMaxAtoms;

	protected JBondQueryFeatureDialog(Frame parent, ExtendedMolecule mol, int bond) {
		super(parent, (mol.isSelectedAtom(mol.getBondAtom(0, bond))
                    && mol.isSelectedAtom(mol.getBondAtom(1, bond))) ?
                        "Multiple Bond Properties" : "Bond Properties", true);
		mMol = mol;
		mBond = bond;

		JPanel p1 = new JPanel();
        double[][] size = { {8, TableLayout.PREFERRED, 16, TableLayout.FILL, TableLayout.PREFERRED, TableLayout.PREFERRED, 8},
                            {8, TableLayout.PREFERRED, TableLayout.PREFERRED, TableLayout.PREFERRED, TableLayout.PREFERRED, 8, TableLayout.PREFERRED, 8, TableLayout.PREFERRED, 4} };
        p1.setLayout(new TableLayout(size));

		mCBSingle = new JCheckBox("Single");
		p1.add(mCBSingle,"1,1");

		mCBDouble = new JCheckBox("Double");
		p1.add(mCBDouble,"1,2");

		mCBTriple = new JCheckBox("Triple");
		p1.add(mCBTriple,"1,3");

		mCBDelocalized = new JCheckBox("Delocalized");
		p1.add(mCBDelocalized,"1,4");

		mComboBoxRing = new JComboBox();
		mComboBoxRing.addItem("any ring state");
		mComboBoxRing.addItem("is not in a ring");
		mComboBoxRing.addItem("is any ring bond");
		mComboBoxRing.addItem("is non-aromatic ring bond");
		mComboBoxRing.addItem("is aromatic bond");
		p1.add(mComboBoxRing,"1,6");

        mCBIsBridge = new JCheckBox("Is atom bridge between");
        mCBIsBridge.addActionListener(this);
        p1.add(mCBIsBridge,"3,1,5,1");

        mComboBoxMinAtoms = new JComboBox();
        int itemCount = (1 << Molecule.cBondQFBridgeMinBits);
        for (int i=0; i<itemCount; i++)
            mComboBoxMinAtoms.addItem(""+i);
        p1.add(mComboBoxMinAtoms,"4,2");
        p1.add(new JLabel(" and"),"5,2");
        mComboBoxMinAtoms.addActionListener(this);

        mComboBoxMaxAtoms = new JComboBox();
        populateComboBoxMaxAtoms(0);
        p1.add(mComboBoxMaxAtoms,"4,3");
        p1.add(new JLabel(" atoms"),"5,3");

        mComboBoxRingSize = new JComboBox();
        mComboBoxRingSize.addItem("any ring size");
        mComboBoxRingSize.addItem("is in 3-membered ring");
        mComboBoxRingSize.addItem("is in 4-membered ring");
        mComboBoxRingSize.addItem("is in 5-membered ring");
        mComboBoxRingSize.addItem("is in 6-membered ring");
        mComboBoxRingSize.addItem("is in 7-membered ring");
		p1.add(mComboBoxRingSize, "3,6,5,6");

		mCBMatchStereo = new JCheckBox("Match Stereo Configuration", (mol.getBondQueryFeatures(bond) & Molecule.cBondQFMatchStereo) != 0);
		mCBMatchStereo.addActionListener(this);
		p1.add(mCBMatchStereo, "1,8,5,8");

		JPanel bp = new JPanel();
        bp.setBorder(BorderFactory.createEmptyBorder(12, 8, 8, 8));
        bp.setLayout(new BorderLayout());
        JPanel ibp = new JPanel();
        ibp.setLayout(new GridLayout(1, 6, 8, 0));
        JButton bcancel = new JButton("Cancel");
        bcancel.addActionListener(this);
        ibp.add(bcancel);
        JButton bok = new JButton("OK");
        bok.addActionListener(this);
        ibp.add(bok);
        bp.add(ibp, BorderLayout.EAST);

		getContentPane().add(p1,BorderLayout.CENTER);
		getContentPane().add(bp,BorderLayout.SOUTH);
		getRootPane().setDefaultButton(bok);

        mMol.ensureHelperArrays(Molecule.cHelperRings);
        setInitialStates();

        pack();
        setLocationRelativeTo(parent);
        setVisible(true);
		}


	public void actionPerformed(ActionEvent e) {
        if (e.getSource() == mCBIsBridge) {
            enableItems();
            return;
            }
        if (e.getSource() == mComboBoxMinAtoms) {
            int minAtoms = mComboBoxMinAtoms.getSelectedIndex();
            if (mFirstSpanItem != minAtoms)  {
                int maxAtoms = mFirstSpanItem + mComboBoxMaxAtoms.getSelectedIndex();
                int itemCount = populateComboBoxMaxAtoms(minAtoms);
                if (maxAtoms < minAtoms)
                    mComboBoxMaxAtoms.setSelectedIndex(0);
                else if (maxAtoms < minAtoms+itemCount)
                    mComboBoxMaxAtoms.setSelectedIndex(maxAtoms-minAtoms);
                else
                    mComboBoxMaxAtoms.setSelectedIndex(itemCount-1);

                mFirstSpanItem = minAtoms;
                }
            }

        if (e.getActionCommand() == "Cancel")
			dispose();
		else if (e.getActionCommand() == "OK") {
			setQueryFeatures();
			dispose();
			}
		}


	private void setInitialStates() {
		int queryFeatures = mMol.getBondQueryFeatures(mBond);
		int bondType = (mMol.getBondType(mBond) == Molecule.cBondTypeDelocalized
		             || mMol.isDelocalizedBond(mBond)) ?
						0 : mMol.getBondOrder(mBond);

		if ((queryFeatures & Molecule.cBondQFSingle) != 0 || bondType == 1)
			mCBSingle.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFDouble) != 0 || bondType == 2)
			mCBDouble.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFTriple) != 0 || bondType == 3)
			mCBTriple.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFDelocalized) != 0 || bondType == 0)
			mCBDelocalized.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFMatchStereo) != 0)
			mCBMatchStereo.setSelected(true);

		if ((queryFeatures & Molecule.cBondQFNotRing) != 0)
			mComboBoxRing.setSelectedIndex(1);
		else if ((queryFeatures & Molecule.cBondQFRing) != 0) {
			if ((queryFeatures & Molecule.cBondQFNotAromatic) != 0)
				mComboBoxRing.setSelectedIndex(3);
			else if ((queryFeatures & Molecule.cBondQFAromatic) != 0)
				mComboBoxRing.setSelectedIndex(4);
			else
				mComboBoxRing.setSelectedIndex(2);
			}

		int ringSize = (queryFeatures & Molecule.cBondQFRingSize) >> Molecule.cBondQFRingSizeShift;
		mComboBoxRingSize.setSelectedIndex((ringSize == 0) ? 0 : ringSize-2);

        if ((queryFeatures & Molecule.cBondQFBridge) != 0) {
            mCBIsBridge.setSelected(true);
            int minAtoms = (queryFeatures & Molecule.cBondQFBridgeMin) >> Molecule.cBondQFBridgeMinShift;
            int atomSpan = (queryFeatures & Molecule.cBondQFBridgeSpan) >> Molecule.cBondQFBridgeSpanShift;
            mComboBoxMinAtoms.setSelectedIndex(minAtoms);
            populateComboBoxMaxAtoms(minAtoms);
            mComboBoxMaxAtoms.setSelectedIndex(atomSpan);
            }
        
        enableItems();
        }


    private int populateComboBoxMaxAtoms(int minAtoms) {
        mComboBoxMaxAtoms.removeAllItems();
        int itemCount = (1 << Molecule.cBondQFBridgeSpanBits);
        for (int i=0; i<itemCount; i++)
            mComboBoxMaxAtoms.addItem(""+(minAtoms+i));
        return itemCount;
        }

    private void enableItems() {
        boolean bridgeIsSelected = mCBIsBridge.isSelected();
        mCBSingle.setEnabled(!bridgeIsSelected);
        mCBDouble.setEnabled(!bridgeIsSelected);
        mCBTriple.setEnabled(!bridgeIsSelected);
        mCBDelocalized.setEnabled(!bridgeIsSelected);
        mCBMatchStereo.setEnabled(!bridgeIsSelected
        						&& mMol.getBondOrder(mBond) == 2	// exclude BINAP-type stereo bonds for now
        						&& mMol.getBondParity(mBond) != Molecule.cBondParityNone
        						&& mMol.getBondParity(mBond) != Molecule.cBondParityUnknown);
        mComboBoxRing.setEnabled(!bridgeIsSelected);
        mComboBoxRingSize.setEnabled(!bridgeIsSelected);
        mComboBoxMinAtoms.setEnabled(bridgeIsSelected);
        mComboBoxMaxAtoms.setEnabled(bridgeIsSelected);
        }


    private void setQueryFeatures() {
        if (isSelectedBond(mBond)) {
            for (int bond=0; bond<mMol.getAllBonds(); bond++)
                if (isSelectedBond(bond))
                    setQueryFeatures(bond);
            }
        else {
            setQueryFeatures(mBond);
            }
        }


    private void setQueryFeatures(int bond) {
		int queryFeatures = 0;

        if (mCBIsBridge.isSelected()) {
            int minAtoms = mComboBoxMinAtoms.getSelectedIndex();
            int atomSpan = mComboBoxMaxAtoms.getSelectedIndex();
            queryFeatures |= (minAtoms << Molecule.cBondQFBridgeMinShift);
            queryFeatures |= (atomSpan << Molecule.cBondQFBridgeSpanShift);
            queryFeatures &= ~Molecule.cBondQFBondTypes;
            }
        else {
            int bondType = -1;
            if (mCBSingle.isSelected()) {
                mMol.setBondType(bond, Molecule.cBondTypeSingle);
                bondType = 1;
                }
            else if (mCBDouble.isSelected()) {
                mMol.setBondType(bond, Molecule.cBondTypeDouble);
                bondType = 2;
                }
            else if (mCBTriple.isSelected()) {
                mMol.setBondType(bond, Molecule.cBondTypeTriple);
                bondType = 3;
                }
            else if (mCBDelocalized.isSelected()) {
                if (!mMol.isDelocalizedBond(bond))
                    mMol.setBondType(bond, Molecule.cBondTypeDelocalized);
                bondType = 0;
                }

            if (mCBSingle.isSelected() && bondType != 1)
    			queryFeatures |= Molecule.cBondQFSingle;
    		if (mCBDouble.isSelected() && bondType != 2)
    			queryFeatures |= Molecule.cBondQFDouble;
    		if (mCBTriple.isSelected() && bondType != 3)
    			queryFeatures |= Molecule.cBondQFTriple;
    		if (mCBDelocalized.isSelected() && bondType != 0)
    			queryFeatures |= Molecule.cBondQFDelocalized;
    		if (mCBMatchStereo.isSelected())
    			queryFeatures |= Molecule.cBondQFMatchStereo;

    		if (!mMol.isAromaticBond(bond)) {
    			if (mComboBoxRing.getSelectedIndex() == 4)
    				queryFeatures |= Molecule.cBondQFAromatic;
    			else if (mComboBoxRing.getSelectedIndex() == 3)
    				queryFeatures |= Molecule.cBondQFNotAromatic | Molecule.cBondQFRing;

	    		if (!mMol.isRingBond(bond)) {
	    			if (mComboBoxRing.getSelectedIndex() == 2)
	    				queryFeatures |= Molecule.cBondQFRing;
	    			else if (mComboBoxRing.getSelectedIndex() == 1)
	    				queryFeatures |= Molecule.cBondQFNotRing;
	    			}
    			}
            }

        if (mComboBoxRingSize.getSelectedIndex() != 0)
            queryFeatures |= ((mComboBoxRingSize.getSelectedIndex()+2) << Molecule.cBondQFRingSizeShift);

		mMol.setBondQueryFeature(bond, Molecule.cBondQFAllFeatures, false);
		mMol.setBondQueryFeature(bond, queryFeatures, true);
		}

    private boolean isSelectedBond(int bond) {
        return mMol.isSelectedAtom(mMol.getBondAtom(0, bond))
            && mMol.isSelectedAtom(mMol.getBondAtom(1, bond));
        }
	}
