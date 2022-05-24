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
import com.actelion.research.gui.hidpi.HiDPIHelper;
import info.clearthought.layout.TableLayout;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Arrays;

@Deprecated
public class JAtomQueryFeatureDialog extends JDialog
					 implements ActionListener,ItemListener {
    static final long serialVersionUID = 0x20060720;

	private Component           mParent;
    private ExtendedMolecule	mMol;
	private int					mAtom;
	private JCheckBox			mCBAny,mCBBlocked,mCBSubstituted,mCBMatchStereo,mCBExcludeGroup;
	private JComboBox			mChoiceArom,mChoiceRingState,mChoiceRingSize,mChoiceCharge,
	                            mChoiceNeighbours,mChoiceHydrogen,mChoicePi,mChoiceReactionParityHint;
	private JTextField			mTFAtomList;
	private JLabel				mLabelAtomList;

	protected JAtomQueryFeatureDialog(Dialog parent, ExtendedMolecule mol, int atom, boolean includeReactionHints) {
		super(parent, (mol.isSelectedAtom(atom)) ? "Multiple Atom Properties" : "Atom Properties", true);
		init(parent, mol, atom, includeReactionHints);
		}

	protected JAtomQueryFeatureDialog(Frame parent, ExtendedMolecule mol, int atom, boolean includeReactionHints) {
		super(parent, (mol.isSelectedAtom(atom)) ? "Multiple Atom Properties" : "Atom Properties", true);
		init(parent, mol, atom, includeReactionHints);
		}

	private void init(Component parent, ExtendedMolecule mol, int atom, boolean includeReactionHints) {
		mParent = parent;
		mMol = mol;
		mAtom = atom;
		boolean opaque = false;
		
		JPanel p1 = new JPanel();
		p1.setOpaque(opaque);

		int gap = HiDPIHelper.scale(8);
		int[] gap1 = { gap, gap/2, gap*3/2, gap/2, gap/2, gap/2, gap/2, gap/2, gap/2, gap*3/2, 0, 0, 0 };
		int[] gap2 = { gap*3/2, gap/2 };
		double[] size2 = new double[1 + 2*gap1.length + (includeReactionHints ? 2*gap2.length : 0)];
		int index = 0;
		for (int g:gap1) {
			size2[index++] = g;
			size2[index++] = TableLayout.PREFERRED;
			}
		if (includeReactionHints)
			for (int g:gap2) {
				size2[index++] = g;
				size2[index++] = TableLayout.PREFERRED;
				}
		size2[index++] = gap;
        double[][] size = { {gap, TableLayout.PREFERRED, gap, TableLayout.PREFERRED, gap}, size2 };
        p1.setLayout(new TableLayout(size));
		
		mCBAny = new JCheckBox("any atomic number");
		mCBAny.setOpaque(opaque);
		mCBAny.addItemListener(this);
		p1.add(mCBAny, "1,1,3,1");

		mLabelAtomList = new JLabel("excluded atoms:");
		mLabelAtomList.setOpaque(opaque);
		mTFAtomList = new JTextField(8);

		p1.add(mLabelAtomList, "1,3");
		p1.add(mTFAtomList, "3,3");

		mChoiceArom = new JComboBox();
		mChoiceArom.setOpaque(opaque);
		mChoiceArom.addItem("any aromatic state");
		mChoiceArom.addItem("is aromatic");
		mChoiceArom.addItem("is not aromatic");
		p1.add(mChoiceArom, "1,5,3,5");

		mChoiceRingState = new JComboBox();
		mChoiceRingState.setOpaque(opaque);
		mChoiceRingState.addItem("any ring state");
		mChoiceRingState.addItem("is chain atom");
		mChoiceRingState.addItem("is any ring atom");
		mChoiceRingState.addItem("has 2 ring bonds");
		mChoiceRingState.addItem("has 3 ring bonds");
		mChoiceRingState.addItem("has >3 ring bonds");
		p1.add(mChoiceRingState, "1,7,3,7");

		mChoiceRingSize = new JComboBox();
		mChoiceRingSize.setOpaque(opaque);
		mChoiceRingSize.addItem("any ring size");
		mChoiceRingSize.addItem("is in 3-membered ring");
        mChoiceRingSize.addItem("is in 4-membered ring");
        mChoiceRingSize.addItem("is in 5-membered ring");
        mChoiceRingSize.addItem("is in 6-membered ring");
        mChoiceRingSize.addItem("is in 7-membered ring");
		p1.add(mChoiceRingSize, "1,9,3,9");

		mChoiceCharge = new JComboBox();
		mChoiceCharge.setOpaque(opaque);
		mChoiceCharge.addItem("any atom charge");
		mChoiceCharge.addItem("has no charge");
		mChoiceCharge.addItem("has negative charge");
		mChoiceCharge.addItem("has positive charge");
		p1.add(mChoiceCharge, "1,11,3,11");

		mChoiceNeighbours = new JComboBox();
		mChoiceNeighbours.setOpaque(opaque);
		mChoiceNeighbours.addItem("any non-H neighbour count");
		mChoiceNeighbours.addItem("has exactly 1 neighbour");
        mChoiceNeighbours.addItem("has exactly 2 neighbours");
        mChoiceNeighbours.addItem("has exactly 3 neighbours");
        mChoiceNeighbours.addItem("has less than 3 neighbours");
        mChoiceNeighbours.addItem("has less than 4 neighbours");
        mChoiceNeighbours.addItem("has more than 1 neighbour");
        mChoiceNeighbours.addItem("has more than 2 neighbours");
        mChoiceNeighbours.addItem("has more than 3 neighbours");
		p1.add(mChoiceNeighbours, "1,13,3,13");

		mChoiceHydrogen = new JComboBox();
		mChoiceHydrogen.setOpaque(opaque);
		mChoiceHydrogen.addItem("any hydrogen count");
		mChoiceHydrogen.addItem("no hydrogen");
		mChoiceHydrogen.addItem("exactly 1 hydrogen");
        mChoiceHydrogen.addItem("exactly 2 hydrogens");
		mChoiceHydrogen.addItem("at least 1 hydrogen");
		mChoiceHydrogen.addItem("at least 2 hydrogens");
		mChoiceHydrogen.addItem("at least 3 hydrogens");
        mChoiceHydrogen.addItem("less than 2 hydrogens");
        mChoiceHydrogen.addItem("less than 3 hydrogens");
		p1.add(mChoiceHydrogen, "1,15,3,15");

        mChoicePi = new JComboBox();
        mChoicePi.setOpaque(opaque);
        mChoicePi.addItem("any pi electron count");
        mChoicePi.addItem("no pi electrons");
        mChoicePi.addItem("exactly 1 pi electron");
        mChoicePi.addItem("exactly 2 pi electrons");
        mChoicePi.addItem("at least 1 pi electron");
        p1.add(mChoicePi, "1,17,3,17");

		mCBBlocked = new JCheckBox("prohibit further substitution");
		mCBBlocked.setOpaque(opaque);
		mCBBlocked.addItemListener(this);
		p1.add(mCBBlocked, "1,19,3,19");

		mCBSubstituted = new JCheckBox("require further substitution");
		mCBSubstituted.setOpaque(opaque);
		mCBSubstituted.addItemListener(this);
		p1.add(mCBSubstituted, "1,21,3,21");

		mCBMatchStereo = new JCheckBox("match stereo center");
		mCBMatchStereo.setOpaque(opaque);
		p1.add(mCBMatchStereo, "1,23,3,23");

		mCBExcludeGroup = new JCheckBox("is part of exclude group");
		mCBExcludeGroup.setOpaque(opaque);
		p1.add(mCBExcludeGroup, "1,25,3,25");

		if (includeReactionHints) {
			p1.add(new JLabel("Stereo center hint for product:"), "1,27,3,27");
			mChoiceReactionParityHint = new JComboBox();
			mChoiceReactionParityHint.setOpaque(opaque);
			mChoiceReactionParityHint.addItem("Make unknown in product");
			mChoiceReactionParityHint.addItem("Keep reactant configuration");
			mChoiceReactionParityHint.addItem("Invert reactant configuration");
			mChoiceReactionParityHint.addItem("Racemise configuration");
			p1.add(mChoiceReactionParityHint, "1,29,3,29");
			}

		JPanel buttonpanel = new JPanel();
        buttonpanel.setBorder(BorderFactory.createEmptyBorder(12, 8, 8, 8));
        buttonpanel.setLayout(new BorderLayout());
        JPanel ibp = new JPanel();
        ibp.setLayout(new GridLayout(1, 2, 8, 0));
        JButton bcancel = new JButton("Cancel");
        bcancel.addActionListener(this);
        ibp.add(bcancel);
        JButton bok = new JButton("OK");
        bok.addActionListener(this);
        ibp.add(bok);
        buttonpanel.add(ibp, BorderLayout.EAST);

		getContentPane().add(p1,BorderLayout.CENTER);
		getContentPane().add(buttonpanel,BorderLayout.SOUTH);
		getRootPane().setDefaultButton(bok);

		pack();

        setLocationRelativeTo(parent);

		mMol.ensureHelperArrays(Molecule.cHelperCIP);
		setInitialStates();
        setVisible(true);
		}


	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand() == "Cancel")
			dispose();
		else if (e.getActionCommand() == "OK") {
			setQueryFeatures();
			dispose();
			}
		}


	public void itemStateChanged(ItemEvent e) {
		if (e.getItemSelectable() == mCBAny) {
			if (mCBAny.isSelected()) {
				mTFAtomList.setText("");
				mLabelAtomList.setText("excluded atoms:");
				}
			else {
				mTFAtomList.setText(mMol.getAtomLabel(mAtom));
				mLabelAtomList.setText("allowed atoms:");
				}
			}
		else if (e.getItemSelectable() == mCBBlocked && e.getStateChange() == ItemEvent.SELECTED) {
			mCBSubstituted.setSelected(false);
			mChoiceNeighbours.setSelectedIndex(0);
		    }
		else if (e.getItemSelectable() == mCBSubstituted && e.getStateChange() == ItemEvent.SELECTED)
			mCBBlocked.setSelected(false);
		}


	private void setInitialStates() {
		long queryFeatures = mMol.getAtomQueryFeatures(mAtom);

		if ((queryFeatures & Molecule.cAtomQFAny) != 0) {
			mCBAny.setSelected(true);
			mLabelAtomList.setText("excluded atoms:");
			}
		else
			mLabelAtomList.setText("allowed atoms:");

		mTFAtomList.setText(mMol.getAtomList(mAtom) == null ? "" : mMol.getAtomListString(mAtom));

		long aromState = queryFeatures & Molecule.cAtomQFAromState;
		if (aromState == Molecule.cAtomQFAromatic)
			mChoiceArom.setSelectedIndex(1);
		else if (aromState == Molecule.cAtomQFNotAromatic)
			mChoiceArom.setSelectedIndex(2);

		long ringState = queryFeatures & Molecule.cAtomQFRingState;
		if (ringState == (Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds))
			mChoiceRingState.setSelectedIndex(1);
		else if (ringState == Molecule.cAtomQFNotChain)
			mChoiceRingState.setSelectedIndex(2);
		else if (ringState == (Molecule.cAtomQFNotChain | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds))
			mChoiceRingState.setSelectedIndex(3);
		else if (ringState == (Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot4RingBonds))
			mChoiceRingState.setSelectedIndex(4);
		else if (ringState == (Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds))
			mChoiceRingState.setSelectedIndex(5);

		long ringSize = (queryFeatures & Molecule.cAtomQFSmallRingSize) >> Molecule.cAtomQFSmallRingSizeShift;
		mChoiceRingSize.setSelectedIndex((ringSize == 0) ? 0 : (int)ringSize-2);

		long neighbourFeatures = queryFeatures & Molecule.cAtomQFNeighbours;
		if (neighbourFeatures == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot1Neighbour))
		    mChoiceNeighbours.setSelectedIndex(1);
        else if (neighbourFeatures == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot2Neighbours))
            mChoiceNeighbours.setSelectedIndex(2);
        else if (neighbourFeatures == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot3Neighbours))
            mChoiceNeighbours.setSelectedIndex(3);
        else if (neighbourFeatures == (Molecule.cAtomQFNot3Neighbours | Molecule.cAtomQFNot4Neighbours))
            mChoiceNeighbours.setSelectedIndex(4);
        else if (neighbourFeatures == Molecule.cAtomQFNot4Neighbours)
            mChoiceNeighbours.setSelectedIndex(5);
        else if (neighbourFeatures == (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour))
            mChoiceNeighbours.setSelectedIndex(6);
        else if (neighbourFeatures == (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour | Molecule.cAtomQFNot2Neighbours))
            mChoiceNeighbours.setSelectedIndex(7);
        else if (neighbourFeatures == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot4Neighbours))
            mChoiceNeighbours.setSelectedIndex(8);

		long chargeFeatures = queryFeatures & Molecule.cAtomQFCharge;
		if (chargeFeatures == (Molecule.cAtomQFNotChargeNeg | Molecule.cAtomQFNotChargePos))
			mChoiceCharge.setSelectedIndex(1);
		else if (chargeFeatures == (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargePos))
			mChoiceCharge.setSelectedIndex(2);
		else if (chargeFeatures == (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargeNeg))
			mChoiceCharge.setSelectedIndex(3);

		long hydrogenFeatures = queryFeatures & Molecule.cAtomQFHydrogen;
		if (hydrogenFeatures == (Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen))
			mChoiceHydrogen.setSelectedIndex(1);
		else if (hydrogenFeatures == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen))
			mChoiceHydrogen.setSelectedIndex(2);
        else if (hydrogenFeatures == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot3Hydrogen))
            mChoiceHydrogen.setSelectedIndex(3);
		else if (hydrogenFeatures == Molecule.cAtomQFNot0Hydrogen)
			mChoiceHydrogen.setSelectedIndex(4);
		else if (hydrogenFeatures == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen))
			mChoiceHydrogen.setSelectedIndex(5);
		else if (hydrogenFeatures == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot2Hydrogen))
			mChoiceHydrogen.setSelectedIndex(6);
        else if (hydrogenFeatures == (Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen))
            mChoiceHydrogen.setSelectedIndex(7);
        else if (hydrogenFeatures == Molecule.cAtomQFNot3Hydrogen)
            mChoiceHydrogen.setSelectedIndex(8);

        long piFeatures = queryFeatures & Molecule.cAtomQFPiElectrons;
        if (piFeatures == (Molecule.cAtomQFNot1PiElectron | Molecule.cAtomQFNot2PiElectrons))
            mChoicePi.setSelectedIndex(1);
        else if (piFeatures == (Molecule.cAtomQFNot0PiElectrons | Molecule.cAtomQFNot2PiElectrons))
            mChoicePi.setSelectedIndex(2);
        else if (piFeatures == (Molecule.cAtomQFNot0PiElectrons | Molecule.cAtomQFNot1PiElectron))
            mChoicePi.setSelectedIndex(3);
        else if (piFeatures == Molecule.cAtomQFNot0PiElectrons)
            mChoicePi.setSelectedIndex(4);

		if ((queryFeatures & Molecule.cAtomQFNoMoreNeighbours) != 0)
			mCBBlocked.setSelected(true);

		if ((queryFeatures & Molecule.cAtomQFMoreNeighbours) != 0)
			mCBSubstituted.setSelected(true);

		if ((queryFeatures & Molecule.cAtomQFMatchStereo) != 0)
			mCBMatchStereo.setSelected(true);

		if ((queryFeatures & Molecule.cAtomQFExcludeGroup) != 0)
			mCBExcludeGroup.setSelected(true);

		if (mChoiceReactionParityHint != null) {
			long rxnStereo = queryFeatures & Molecule.cAtomQFRxnParityHint;
			if (rxnStereo == Molecule.cAtomQFRxnParityRetain)
				mChoiceReactionParityHint.setSelectedIndex(1);
			else if (rxnStereo == Molecule.cAtomQFRxnParityInvert)
				mChoiceReactionParityHint.setSelectedIndex(2);
			else if (rxnStereo == Molecule.cAtomQFRxnParityRacemize)
				mChoiceReactionParityHint.setSelectedIndex(3);
			}
		}


	private void setQueryFeatures() {
        int[] atomList = createAtomList();
        if (mMol.isSelectedAtom(mAtom)) {
            for (int atom=0; atom<mMol.getAllAtoms(); atom++)
                if (mMol.isSelectedAtom(atom))
                    setQueryFeatures(atom, atomList);
            }
        else {
            setQueryFeatures(mAtom, atomList);
            }
        }


    private void setQueryFeatures(int atom, int[] atomList) {
		int queryFeatures = 0;

		if (mCBAny.isSelected()) {
			queryFeatures |= Molecule.cAtomQFAny;
			mMol.setAtomList(atom, atomList, true);
			}
		else
			mMol.setAtomList(atom, atomList, false);

		if (!mMol.isAromaticAtom(atom)) {
			if (mChoiceArom.getSelectedIndex() == 1)
				queryFeatures |= Molecule.cAtomQFAromatic;
			else if (mChoiceArom.getSelectedIndex() == 2)
				queryFeatures |= Molecule.cAtomQFNotAromatic;
			}

		int ringBonds = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (mMol.isRingBond(mMol.getConnBond(atom, i)))
				ringBonds++;
		switch (mChoiceRingState.getSelectedIndex()) {
		case 1:
			if (ringBonds == 0)
				queryFeatures |= (Molecule.cAtomQFNot2RingBonds
								| Molecule.cAtomQFNot3RingBonds
								| Molecule.cAtomQFNot4RingBonds);
			break;
		case 2:
			queryFeatures |= Molecule.cAtomQFNotChain;
			break;
		case 3:
			if (ringBonds < 3)
				queryFeatures |= (Molecule.cAtomQFNotChain
								| Molecule.cAtomQFNot3RingBonds
								| Molecule.cAtomQFNot4RingBonds);
			break;
		case 4:
			if (ringBonds < 4)
				queryFeatures |= (Molecule.cAtomQFNotChain
								| Molecule.cAtomQFNot2RingBonds
								| Molecule.cAtomQFNot4RingBonds);
			break;
		case 5:
			queryFeatures |= (Molecule.cAtomQFNotChain
							| Molecule.cAtomQFNot2RingBonds
							| Molecule.cAtomQFNot3RingBonds);
			break;
			}

        if (mChoiceRingSize.getSelectedIndex() != 0)
            queryFeatures |= ((mChoiceRingSize.getSelectedIndex()+2) << Molecule.cAtomQFSmallRingSizeShift);

        switch (mChoiceCharge.getSelectedIndex()) {
        case 1:
            queryFeatures |= (Molecule.cAtomQFCharge & ~Molecule.cAtomQFNotCharge0);
        	break;
        case 2:
            queryFeatures |= (Molecule.cAtomQFCharge & ~Molecule.cAtomQFNotChargeNeg);
        	break;
        case 3:
            queryFeatures |= (Molecule.cAtomQFCharge & ~Molecule.cAtomQFNotChargePos);
        	break;
        	}

        switch (mChoiceNeighbours.getSelectedIndex()) {
        case 1:
            if (mMol.getConnAtoms(atom) == 1)
                queryFeatures |= Molecule.cAtomQFNoMoreNeighbours;
            else if (mMol.getConnAtoms(atom) < 1)
                queryFeatures |= (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot1Neighbour);
            break;
        case 2:
            if (mMol.getConnAtoms(atom) == 2)
                queryFeatures |= Molecule.cAtomQFNoMoreNeighbours;
            else if (mMol.getConnAtoms(atom) < 2)
                queryFeatures |= (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot2Neighbours);
            break;
        case 3:
            if (mMol.getConnAtoms(atom) == 3)
                queryFeatures |= Molecule.cAtomQFNoMoreNeighbours;
            else if (mMol.getConnAtoms(atom) < 3)
                queryFeatures |= (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot3Neighbours);
            break;
        case 4: // less than 3 non-H neighbours
            if (mMol.getConnAtoms(atom) == 2)
                queryFeatures |= Molecule.cAtomQFNoMoreNeighbours;
            else if (mMol.getConnAtoms(atom) < 2)
                queryFeatures |= (Molecule.cAtomQFNot3Neighbours | Molecule.cAtomQFNot4Neighbours);
            break;
        case 5: // less than 4 non-H neighbours
            if (mMol.getConnAtoms(atom) == 3)
                queryFeatures |= Molecule.cAtomQFNoMoreNeighbours;
            else if (mMol.getConnAtoms(atom) < 3)
                queryFeatures |= Molecule.cAtomQFNot4Neighbours;
            break;
        case 6: // more than 1 non-H neighbour
            if (mMol.getConnAtoms(atom) == 1)
                queryFeatures |= Molecule.cAtomQFMoreNeighbours;
            else if (mMol.getConnAtoms(atom) < 1)
                queryFeatures |= (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour);
            break;
        case 7: // more than 2 non-H neighbours
            if (mMol.getConnAtoms(atom) == 2)
                queryFeatures |= Molecule.cAtomQFMoreNeighbours;
            else if (mMol.getConnAtoms(atom) < 2)
                queryFeatures |= (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour | Molecule.cAtomQFNot2Neighbours);
            break;
        case 8: // more than 3 non-H neighbours
            if (mMol.getConnAtoms(atom) == 3)
                queryFeatures |= Molecule.cAtomQFMoreNeighbours;
            else if (mMol.getConnAtoms(atom) < 3)
                queryFeatures |= (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot4Neighbours);
            break;
            }

		switch (mChoiceHydrogen.getSelectedIndex()) {
		case 1:	// no hydrogens
			queryFeatures |= (Molecule.cAtomQFNot1Hydrogen
							| Molecule.cAtomQFNot2Hydrogen
			                | Molecule.cAtomQFNot3Hydrogen);
			break;
		case 2:	// exactly 1 hydrogen
			queryFeatures |= (Molecule.cAtomQFNot0Hydrogen
			                | Molecule.cAtomQFNot2Hydrogen
			                | Molecule.cAtomQFNot3Hydrogen);
			break;
        case 3: // exactly 2 hydrogen
            queryFeatures |= (Molecule.cAtomQFNot0Hydrogen
                            | Molecule.cAtomQFNot1Hydrogen
                            | Molecule.cAtomQFNot3Hydrogen);
            break;
		case 4:	// at least 1 hydrogen
			queryFeatures |= Molecule.cAtomQFNot0Hydrogen;
			break;
		case 5:	// at least 2 hydrogens
			queryFeatures |= (Molecule.cAtomQFNot0Hydrogen
							| Molecule.cAtomQFNot1Hydrogen);
			break;
		case 6:	// at least 3 hydrogens
			queryFeatures |= (Molecule.cAtomQFNot0Hydrogen
							| Molecule.cAtomQFNot1Hydrogen
							| Molecule.cAtomQFNot2Hydrogen);
			break;
        case 7: // less than 2 hydrogens
            queryFeatures |= (Molecule.cAtomQFNot2Hydrogen
                            | Molecule.cAtomQFNot3Hydrogen);
            break;
        case 8: // less than 3 hydrogens
            queryFeatures |= (Molecule.cAtomQFNot3Hydrogen);
            break;
			}

        switch (mChoicePi.getSelectedIndex()) {
        case 1: // no pi electrons
            queryFeatures |= (Molecule.cAtomQFNot1PiElectron
                            | Molecule.cAtomQFNot2PiElectrons);
            break;
        case 2: // exactly 1 pi electron
            queryFeatures |= (Molecule.cAtomQFNot0PiElectrons
                            | Molecule.cAtomQFNot2PiElectrons);
            break;
        case 3: // exactly 2 pi electrons
            queryFeatures |= (Molecule.cAtomQFNot0PiElectrons
                            | Molecule.cAtomQFNot1PiElectron);
            break;
        case 4: // at least 1 pi electron
            queryFeatures |= Molecule.cAtomQFNot0PiElectrons;
            break;
            }

		if (mCBBlocked.isSelected()
		 && (mMol.getFreeValence(atom) > 0
		  || (mMol.getAtomCharge(atom)==0 && (mMol.getAtomicNo(atom)==5 || mMol.isNitrogenFamily(atom) || mMol.isChalcogene(atom)))))
			queryFeatures |= Molecule.cAtomQFNoMoreNeighbours;

		if (mCBSubstituted.isSelected()
		 && (mMol.getFreeValence(atom) > 0
		  || (mMol.getAtomCharge(atom)==0 && (mMol.getAtomicNo(atom)==5 || mMol.isNitrogenFamily(atom) || mMol.isChalcogene(atom)))))
			queryFeatures |= Molecule.cAtomQFMoreNeighbours;

		if (mCBMatchStereo.isSelected())
			queryFeatures |= Molecule.cAtomQFMatchStereo;

		if (mCBExcludeGroup.isSelected())
			queryFeatures |= Molecule.cAtomQFExcludeGroup;

		if (mChoiceReactionParityHint != null) {
		    switch (mChoiceReactionParityHint.getSelectedIndex()) {
			    case 1:
				    queryFeatures |= Molecule.cAtomQFRxnParityRetain;
				    break;
			    case 2:
				    queryFeatures |= Molecule.cAtomQFRxnParityInvert;
				    break;
			    case 3:
				    queryFeatures |= Molecule.cAtomQFRxnParityRacemize;
				    break;
			    }
			}

	    mMol.setAtomQueryFeature(atom, 0xFFFFFFFF, false);
		mMol.setAtomQueryFeature(atom, queryFeatures, true);
		}


	private int[] createAtomList() {
		String listString = mTFAtomList.getText();
		if (listString.length() == 0)
			return null;

        int[] list = null;
		int delimiterIndex;
		do {
			String label;
			delimiterIndex = listString.indexOf(',');
			if (delimiterIndex == -1)
				label = listString;
			else {
				label = listString.substring(0, delimiterIndex);
				if (delimiterIndex == listString.length() - 1)
					listString = "";
				else
					listString = listString.substring(delimiterIndex+1);
				}

			int atomicNo = Molecule.getAtomicNoFromLabel(label);
			if (atomicNo != 0) {
				if (atomicNo == 1) {
					JOptionPane.showMessageDialog(mParent, "'H' cannot be part of an atom list and is removed.");
					}
                else if (list == null) {
                    list = new int[1];
                    list[0] = atomicNo;
                    }
                else {
    				boolean found = false;
    				for (int i=0; i<list.length; i++) {
    					if (atomicNo == list[i]) {
    						found = true;
    						break;
    						}
    					}
    				if (!found) {
    					int[] newList = new int[list.length+1];
                        for (int i=0; i<list.length; i++)
                            newList[i] = list[i];
                        newList[list.length] = atomicNo;
                        list = newList;
                        }
                    }
				}
			} while (delimiterIndex != -1);

		if (list != null)
		    Arrays.sort(list);

		return list;
		}
	}
