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

import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.Molecule;
import info.clearthought.layout.TableLayout;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

@Deprecated
public class JBondQueryFeatureDialog extends JDialog implements ActionListener {
    static final long serialVersionUID = 0x20070822;

    private ExtendedMolecule	mMol;
	private int					mBond,mFirstSpanItem;
	private JCheckBox			mCBSingle,mCBDouble,mCBTriple,mCBDelocalized,mCBMetalLigand,mCBIsBridge,mCBMatchStereo;
	private JComboBox			mComboBoxRing,mComboBoxRingSize,mComboBoxMinAtoms,mComboBoxMaxAtoms;

	protected JBondQueryFeatureDialog(Dialog parent, ExtendedMolecule mol, int bond) {
		super(parent, (mol.isSelectedAtom(mol.getBondAtom(0, bond))
				&& mol.isSelectedAtom(mol.getBondAtom(1, bond))) ?
				"Multiple Bond Properties" : "Bond Properties", true);
		init(parent, mol, bond);
		}

	protected JBondQueryFeatureDialog(Frame parent, ExtendedMolecule mol, int bond) {
		super(parent, (mol.isSelectedAtom(mol.getBondAtom(0, bond))
				&& mol.isSelectedAtom(mol.getBondAtom(1, bond))) ?
				"Multiple Bond Properties" : "Bond Properties", true);
		init(parent, mol, bond);
		}

	private void init(Component parent, ExtendedMolecule mol, int bond) {
		mMol = mol;
		mBond = bond;

		JPanel p1 = new JPanel();
        double[][] size = { {8, TableLayout.FILL, TableLayout.PREFERRED, TableLayout.PREFERRED, 8},
                            {8, TableLayout.PREFERRED, TableLayout.PREFERRED, TableLayout.PREFERRED,
								TableLayout.PREFERRED, TableLayout.PREFERRED, 8,
								TableLayout.PREFERRED, 8, TableLayout.PREFERRED, 8, TableLayout.PREFERRED, 16,
								TableLayout.PREFERRED, 4, TableLayout.PREFERRED, 4, TableLayout.PREFERRED, 16} };
        p1.setLayout(new TableLayout(size));

		mCBSingle = new JCheckBox("Single");
		p1.add(mCBSingle,"1,1,3,1");

		mCBDouble = new JCheckBox("Double");
		p1.add(mCBDouble,"1,2,3,2");

		mCBTriple = new JCheckBox("Triple");
		p1.add(mCBTriple,"1,3,3,3");

		mCBDelocalized = new JCheckBox("Delocalized");
		p1.add(mCBDelocalized,"1,4,3,4");

		mCBMetalLigand = new JCheckBox("Coordinate (0-order)");
		p1.add(mCBMetalLigand,"1,5,3,5");

		mComboBoxRing = new JComboBox();
		mComboBoxRing.addItem("any ring state");
		mComboBoxRing.addItem("is not in a ring");
		mComboBoxRing.addItem("is any ring bond");
		mComboBoxRing.addItem("is non-aromatic ring bond");
		mComboBoxRing.addItem("is aromatic bond");
		mComboBoxRing.addActionListener(this);
		p1.add(mComboBoxRing,"1,7,3,7");

		mComboBoxRingSize = new JComboBox();
		mComboBoxRingSize.addItem("any ring size");
		mComboBoxRingSize.addItem("is in 3-membered ring");
		mComboBoxRingSize.addItem("is in 4-membered ring");
		mComboBoxRingSize.addItem("is in 5-membered ring");
		mComboBoxRingSize.addItem("is in 6-membered ring");
		mComboBoxRingSize.addItem("is in 7-membered ring");
		p1.add(mComboBoxRingSize, "1,9,3,9");

		mCBMatchStereo = new JCheckBox("Match Stereo Configuration", (mol.getBondQueryFeatures(bond) & Molecule.cBondQFMatchStereo) != 0);
		mCBMatchStereo.addActionListener(this);
		p1.add(mCBMatchStereo, "1,11,3,11");

		mCBIsBridge = new JCheckBox("Is atom bridge between");
        mCBIsBridge.addActionListener(this);
        p1.add(mCBIsBridge,"1,13,3,13");

        mComboBoxMinAtoms = new JComboBox();
        int itemCount = (1 << Molecule.cBondQFBridgeMinBits);
        for (int i=0; i<itemCount; i++)
            mComboBoxMinAtoms.addItem(""+i);
        p1.add(mComboBoxMinAtoms,"2,15");
        p1.add(new JLabel(" and"),"3,15");
        mComboBoxMinAtoms.addActionListener(this);

        mComboBoxMaxAtoms = new JComboBox();
        populateComboBoxMaxAtoms(0);
        p1.add(mComboBoxMaxAtoms,"2,17");
        p1.add(new JLabel(" atoms"),"3,17");

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
        if (e.getSource() == mCBIsBridge
		 || e.getSource() == mComboBoxRing) {
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
		int bondOrder = (mMol.getBondType(mBond) == Molecule.cBondTypeDelocalized
					  || mMol.isDelocalizedBond(mBond)) ?
						4 : mMol.getBondOrder(mBond);

		if ((queryFeatures & Molecule.cBondQFSingle) != 0 || bondOrder == 1)
			mCBSingle.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFDouble) != 0 || bondOrder == 2)
			mCBDouble.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFTriple) != 0 || bondOrder == 3)
			mCBTriple.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFDelocalized) != 0 || bondOrder == 4)
			mCBDelocalized.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFMetalLigand) != 0 || bondOrder == 0)
			mCBMetalLigand.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFMatchStereo) != 0)
			mCBMatchStereo.setSelected(true);

		int ringState = queryFeatures & Molecule.cBondQFRingState;
		int aromState = queryFeatures & Molecule.cBondQFAromState;
		if (ringState == Molecule.cBondQFNotRing)
			mComboBoxRing.setSelectedIndex(1);
		else if (aromState == Molecule.cBondQFAromatic)
			mComboBoxRing.setSelectedIndex(4);
		else if (ringState == Molecule.cBondQFRing) {
			if (aromState == 0)
				mComboBoxRing.setSelectedIndex(2);
			else if (aromState == Molecule.cBondQFNotAromatic)
				mComboBoxRing.setSelectedIndex(3);
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
		mCBMetalLigand.setEnabled(!bridgeIsSelected);
        mCBMatchStereo.setEnabled(!bridgeIsSelected
        						&& mMol.getBondOrder(mBond) == 2	// exclude BINAP-type stereo bonds for now
        						&& mMol.getBondParity(mBond) != Molecule.cBondParityNone
        						&& mMol.getBondParity(mBond) != Molecule.cBondParityUnknown);
        mComboBoxRing.setEnabled(!bridgeIsSelected);
        mComboBoxRingSize.setEnabled(!bridgeIsSelected && mComboBoxRing.getSelectedIndex() != 1);
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
        	// priority in order of bond orders
            int bondOrder = -1;
            if (mCBSingle.isSelected()) {
                mMol.setBondType(bond, Molecule.cBondTypeSingle);
				bondOrder = 1;
                }
            else if (mCBDelocalized.isSelected() && !mMol.isDelocalizedBond(bond)) {
	            mMol.setBondType(bond, Molecule.cBondTypeDelocalized);
	            bondOrder = 4;
	            }
            else if (mCBDouble.isSelected()) {
                mMol.setBondType(bond, Molecule.cBondTypeDouble);
				bondOrder = 2;
                }
            else if (mCBTriple.isSelected()) {
                mMol.setBondType(bond, Molecule.cBondTypeTriple);
				bondOrder = 3;
                }
 			else if (mCBMetalLigand.isSelected()) {
				mMol.setBondType(bond, Molecule.cBondTypeMetalLigand);
				bondOrder = 0;
				}

            if (mCBSingle.isSelected() && bondOrder != 1)
    			queryFeatures |= Molecule.cBondQFSingle;
    		if (mCBDouble.isSelected() && bondOrder != 2)
    			queryFeatures |= Molecule.cBondQFDouble;
    		if (mCBTriple.isSelected() && bondOrder != 3)
    			queryFeatures |= Molecule.cBondQFTriple;
    		if (mCBDelocalized.isSelected() && bondOrder != 4)
    			queryFeatures |= Molecule.cBondQFDelocalized;
			if (mCBMetalLigand.isSelected() && bondOrder != 0)
				queryFeatures |= Molecule.cBondQFMetalLigand;
    		if (mCBMatchStereo.isSelected())
    			queryFeatures |= Molecule.cBondQFMatchStereo;

			if (mComboBoxRing.getSelectedIndex() != 0) {
				if (mComboBoxRing.getSelectedIndex() == 1) {
					if (!mMol.isRingBond(bond))
						queryFeatures |= Molecule.cBondQFNotRing;
					}
				else if (mComboBoxRing.getSelectedIndex() == 2) {
					if (!mMol.isRingBond(bond))
						queryFeatures |= Molecule.cBondQFRing;
					}
				else if (mComboBoxRing.getSelectedIndex() == 3) {
					if (!mMol.isAromaticBond(bond))
						queryFeatures |= Molecule.cBondQFNotAromatic | Molecule.cBondQFRing;
					}
				else if (mComboBoxRing.getSelectedIndex() == 4) {
					if (!mMol.isAromaticBond(bond))
						queryFeatures |= Molecule.cBondQFAromatic;
					}
    			}

			if (mComboBoxRingSize.getSelectedIndex() != 0) {
				int ringSize = mComboBoxRingSize.getSelectedIndex() + 2;
				int implicitSize = mMol.getBondRingSize(bond);
				if (ringSize != implicitSize)
					queryFeatures |= (ringSize << Molecule.cBondQFRingSizeShift);
				}
            }

		mMol.setBondQueryFeature(bond, Molecule.cBondQFAllFeatures, false);
		mMol.setBondQueryFeature(bond, queryFeatures, true);
		}

    private boolean isSelectedBond(int bond) {
        return mMol.isSelectedAtom(mMol.getBondAtom(0, bond))
            && mMol.isSelectedAtom(mMol.getBondAtom(1, bond));
        }
	}
