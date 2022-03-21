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

package com.actelion.research.gui.editor;

import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.Molecule;
import com.actelion.research.gui.generic.*;
import com.actelion.research.gui.hidpi.HiDPIHelper;

import java.util.Arrays;

public class AtomQueryFeatureDialogBuilder implements DialogEventConsumer {
    static final long serialVersionUID = 0x20060720;

	private GenericDialog       mDialog;
	private GenericLabel        mLabelAtomList;
	private GenericTextField    mTFAtomList;
	private GenericCheckBox     mCBAny,mCBBlocked,mCBSubstituted,mCBMatchStereo,mCBExcludeGroup;
	private GenericComboBox     mChoiceArom,mChoiceRingState,mChoiceRingSize,mChoiceCharge,
								mChoiceNeighbours,mChoiceHydrogen,mChoicePi,mChoiceReactionParityHint;
    private ExtendedMolecule	mMol;
	private int					mAtom;

	public AtomQueryFeatureDialogBuilder(GenericDialogHelper dialogHelper, ExtendedMolecule mol, int atom, boolean includeReactionHints) {
		mDialog = dialogHelper.createDialog(mol.isSelectedAtom(atom) ? "Atom Query Features (Multiple)" : "Atom Query Features", this);
		build(mol, atom, includeReactionHints);
		}

	private void build(ExtendedMolecule mol, int atom, boolean includeReactionHints) {
		mMol = mol;
		mAtom = atom;

		int gap = HiDPIHelper.scale(8);
		int[] gap1 = { gap, gap/2, gap*3/2, gap/2, gap/2, gap/2, gap/2, gap/2, gap/2, gap*3/2, 0, 0, 0 };
		int[] gap2 = { gap*3/2, gap/2 };
		int[] hLayout = {gap, GenericDialog.PREFERRED, gap, GenericDialog.PREFERRED, gap};
		int[] vLayout = new int[1 + 2*gap1.length + (includeReactionHints ? 2*gap2.length : 0)];
		int index = 0;
		for (int g:gap1) {
			vLayout[index++] = g;
			vLayout[index++] = GenericDialog.PREFERRED;
			}
		if (includeReactionHints)
			for (int g:gap2) {
				vLayout[index++] = g;
				vLayout[index++] = GenericDialog.PREFERRED;
				}
		vLayout[index++] = gap;
        mDialog.setLayout(hLayout, vLayout);
		
		mCBAny = mDialog.createCheckBox("any atomic number");
		mCBAny.setEventConsumer(this);
		mDialog.add(mCBAny, 1,1,3,1);

		mLabelAtomList = mDialog.createLabel("excluded atoms:");
		mTFAtomList = mDialog.createTextField(16, 1);

		mDialog.add(mLabelAtomList, 1,3);
		mDialog.add(mTFAtomList, 3,3);

		mChoiceArom = mDialog.createComboBox();
		mChoiceArom.addItem("any aromatic state");
		mChoiceArom.addItem("is aromatic");
		mChoiceArom.addItem("is not aromatic");
		mDialog.add(mChoiceArom, 1,5,3,5);

		mChoiceRingState = mDialog.createComboBox();
		mChoiceRingState.addItem("any ring state");
		mChoiceRingState.addItem("is chain atom");
		mChoiceRingState.addItem("is any ring atom");
		mChoiceRingState.addItem("has 2 ring bonds");
		mChoiceRingState.addItem("has 3 ring bonds");
		mChoiceRingState.addItem("has >3 ring bonds");
		mDialog.add(mChoiceRingState, 1,7,3,7);

		mChoiceRingSize = mDialog.createComboBox();
		mChoiceRingSize.addItem("any ring size");
		mChoiceRingSize.addItem("is in 3-membered ring");
        mChoiceRingSize.addItem("is in 4-membered ring");
        mChoiceRingSize.addItem("is in 5-membered ring");
        mChoiceRingSize.addItem("is in 6-membered ring");
        mChoiceRingSize.addItem("is in 7-membered ring");
		mDialog.add(mChoiceRingSize, 1,9,3,9);

		mChoiceCharge = mDialog.createComboBox();
		mChoiceCharge.addItem("any atom charge");
		mChoiceCharge.addItem("has no charge");
		mChoiceCharge.addItem("has negative charge");
		mChoiceCharge.addItem("has positive charge");
		mDialog.add(mChoiceCharge, 1,11,3,11);

		mChoiceNeighbours = mDialog.createComboBox();
		mChoiceNeighbours.addItem("any non-H neighbour count");
		mChoiceNeighbours.addItem("has exactly 1 neighbour");
        mChoiceNeighbours.addItem("has exactly 2 neighbours");
        mChoiceNeighbours.addItem("has exactly 3 neighbours");
        mChoiceNeighbours.addItem("has less than 3 neighbours");
        mChoiceNeighbours.addItem("has less than 4 neighbours");
        mChoiceNeighbours.addItem("has more than 1 neighbour");
        mChoiceNeighbours.addItem("has more than 2 neighbours");
        mChoiceNeighbours.addItem("has more than 3 neighbours");
		mDialog.add(mChoiceNeighbours, 1,13,3,13);

		mChoiceHydrogen = mDialog.createComboBox();
		mChoiceHydrogen.addItem("any hydrogen count");
		mChoiceHydrogen.addItem("no hydrogen");
		mChoiceHydrogen.addItem("exactly 1 hydrogen");
        mChoiceHydrogen.addItem("exactly 2 hydrogens");
		mChoiceHydrogen.addItem("at least 1 hydrogen");
		mChoiceHydrogen.addItem("at least 2 hydrogens");
		mChoiceHydrogen.addItem("at least 3 hydrogens");
        mChoiceHydrogen.addItem("less than 2 hydrogens");
        mChoiceHydrogen.addItem("less than 3 hydrogens");
		mDialog.add(mChoiceHydrogen, 1,15,3,15);

        mChoicePi = mDialog.createComboBox();
        mChoicePi.addItem("any pi electron count");
        mChoicePi.addItem("no pi electrons");
        mChoicePi.addItem("exactly 1 pi electron");
        mChoicePi.addItem("exactly 2 pi electrons");
        mChoicePi.addItem("at least 1 pi electron");
		mDialog.add(mChoicePi, 1,17,3,17);

		mCBBlocked = mDialog.createCheckBox("prohibit further substitution");
		mCBBlocked.setEventConsumer(this);
		mDialog.add(mCBBlocked, 1,19,3,19);

		mCBSubstituted = mDialog.createCheckBox("require further substitution");
		mCBSubstituted.setEventConsumer(this);
		mDialog.add(mCBSubstituted, 1,21,3,21);

		mCBMatchStereo = mDialog.createCheckBox("match stereo center");
		mDialog.add(mCBMatchStereo, 1,23,3,23);

		mCBExcludeGroup = mDialog.createCheckBox("is part of exclude group");
		mDialog.add(mCBExcludeGroup, 1,25,3,25);

		if (includeReactionHints) {
			mDialog.add(mDialog.createLabel("Stereo center hint for product:"), 1,27,3,27);
			mChoiceReactionParityHint = mDialog.createComboBox();
			mChoiceReactionParityHint.addItem("Copy from generic product");
			mChoiceReactionParityHint.addItem("Keep reactant configuration");
			mChoiceReactionParityHint.addItem("Invert reactant configuration");
			mChoiceReactionParityHint.addItem("Racemise configuration");
			mDialog.add(mChoiceReactionParityHint, 1,29,3,29);
			}

		mMol.ensureHelperArrays(Molecule.cHelperCIP);
		setInitialStates();

		mDialog.showDialog();
		}

	@Override
	public void dialogEventHappened(DialogEvent e) {
		if (e.getWhat() == DialogEvent.WHAT_OK) {
			setQueryFeatures();
			mDialog.disposeDialog();
			}
		else if (e.getWhat() == DialogEvent.WHAT_CANCEL) {
			mDialog.disposeDialog();
			}
		else if (e.getSource() == mCBAny) {
			if (e.getValue() == 1) {
				mTFAtomList.setText("");
				mLabelAtomList.setText("excluded atoms:");
				}
			else {
				mTFAtomList.setText(mMol.getAtomLabel(mAtom));
				mLabelAtomList.setText("allowed atoms:");
				}
			}
		else if (e.getSource() == mCBBlocked) {
			mCBSubstituted.setSelected(false);
			mChoiceNeighbours.setSelectedIndex(0);
		    }
		else if (e.getSource() == mCBSubstituted)
			mCBBlocked.setSelected(false);
		}

	private void setInitialStates() {
		int queryFeatures = mMol.getAtomQueryFeatures(mAtom);

		if ((queryFeatures & Molecule.cAtomQFAny) != 0) {
			mCBAny.setSelected(true);
			mLabelAtomList.setText("excluded atoms:");
			}
		else
			mLabelAtomList.setText("allowed atoms:");

		mTFAtomList.setText(mMol.getAtomList(mAtom) == null ? "" : mMol.getAtomListString(mAtom));

		int aromState = queryFeatures & Molecule.cAtomQFAromState;
		if (aromState == Molecule.cAtomQFAromatic)
			mChoiceArom.setSelectedIndex(1);
		else if (aromState == Molecule.cAtomQFNotAromatic)
			mChoiceArom.setSelectedIndex(2);

		int ringState = queryFeatures & Molecule.cAtomQFRingState;
		switch (ringState) {
			case Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds:
				mChoiceRingState.setSelectedIndex(1);
				break;
			case Molecule.cAtomQFNotChain:
				mChoiceRingState.setSelectedIndex(2);
				break;
			case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds:
				mChoiceRingState.setSelectedIndex(3);
				break;
			case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot4RingBonds:
				mChoiceRingState.setSelectedIndex(4);
				break;
			case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds:
				mChoiceRingState.setSelectedIndex(5);
				break;
			}

		int ringSize = (queryFeatures & Molecule.cAtomQFRingSize) >> Molecule.cAtomQFRingSizeShift;
		mChoiceRingSize.setSelectedIndex((ringSize == 0) ? 0 : ringSize-2);

		int neighbourFeatures = queryFeatures & Molecule.cAtomQFNeighbours;
		switch (neighbourFeatures) {
		case Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot1Neighbour:
		    mChoiceNeighbours.setSelectedIndex(1);
		    break;
        case Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot2Neighbours:
            mChoiceNeighbours.setSelectedIndex(2);
            break;
        case Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot3Neighbours:
            mChoiceNeighbours.setSelectedIndex(3);
            break;
        case Molecule.cAtomQFNot3Neighbours | Molecule.cAtomQFNot4Neighbours:
            mChoiceNeighbours.setSelectedIndex(4);
            break;
        case Molecule.cAtomQFNot4Neighbours:
            mChoiceNeighbours.setSelectedIndex(5);
            break;
        case Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour:
            mChoiceNeighbours.setSelectedIndex(6);
            break;
        case Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour | Molecule.cAtomQFNot2Neighbours:
            mChoiceNeighbours.setSelectedIndex(7);
            break;
        case Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot4Neighbours:
            mChoiceNeighbours.setSelectedIndex(8);
            break;
		    }

		int chargeFeatures = queryFeatures & Molecule.cAtomQFCharge;
		switch (chargeFeatures) {
		case Molecule.cAtomQFNotChargeNeg | Molecule.cAtomQFNotChargePos:
			mChoiceCharge.setSelectedIndex(1);
			break;
		case Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargePos:
			mChoiceCharge.setSelectedIndex(2);
			break;
		case Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargeNeg:
			mChoiceCharge.setSelectedIndex(3);
			break;
			}

		int hydrogenFeatures = queryFeatures & Molecule.cAtomQFHydrogen;
		switch (hydrogenFeatures) {
			case Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen:
				mChoiceHydrogen.setSelectedIndex(1);
				break;
			case Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen:
				mChoiceHydrogen.setSelectedIndex(2);
				break;
            case Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot3Hydrogen:
                mChoiceHydrogen.setSelectedIndex(3);
                break;
			case Molecule.cAtomQFNot0Hydrogen:
				mChoiceHydrogen.setSelectedIndex(4);
				break;
			case Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen:
				mChoiceHydrogen.setSelectedIndex(5);
				break;
			case Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot2Hydrogen:
				mChoiceHydrogen.setSelectedIndex(6);
				break;
            case Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen:
                mChoiceHydrogen.setSelectedIndex(7);
                break;
            case Molecule.cAtomQFNot3Hydrogen:
                mChoiceHydrogen.setSelectedIndex(8);
                break;
			}

        int piFeatures = queryFeatures & Molecule.cAtomQFPiElectrons;
        switch (piFeatures) {
            case Molecule.cAtomQFNot1PiElectron | Molecule.cAtomQFNot2PiElectrons:
                mChoicePi.setSelectedIndex(1);
                break;
            case Molecule.cAtomQFNot0PiElectrons | Molecule.cAtomQFNot2PiElectrons:
                mChoicePi.setSelectedIndex(2);
                break;
            case Molecule.cAtomQFNot0PiElectrons | Molecule.cAtomQFNot1PiElectron:
                mChoicePi.setSelectedIndex(3);
                break;
            case Molecule.cAtomQFNot0PiElectrons:
                mChoicePi.setSelectedIndex(4);
                break;
            }

		if ((queryFeatures & Molecule.cAtomQFNoMoreNeighbours) != 0)
			mCBBlocked.setSelected(true);

		if ((queryFeatures & Molecule.cAtomQFMoreNeighbours) != 0)
			mCBSubstituted.setSelected(true);

		if ((queryFeatures & Molecule.cAtomQFMatchStereo) != 0)
			mCBMatchStereo.setSelected(true);

		if ((queryFeatures & Molecule.cAtomQFExcludeGroup) != 0)
			mCBExcludeGroup.setSelected(true);

		if (mChoiceReactionParityHint != null) {
			int rxnStereo = queryFeatures & Molecule.cAtomQFRxnParityHint;
			switch (rxnStereo) {
				case Molecule.cAtomQFRxnParityRetain:
					mChoiceReactionParityHint.setSelectedIndex(1);
					break;
				case Molecule.cAtomQFRxnParityInvert:
					mChoiceReactionParityHint.setSelectedIndex(2);
					break;
				case Molecule.cAtomQFRxnParityRacemize:
					mChoiceReactionParityHint.setSelectedIndex(3);
					break;
				}
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
            queryFeatures |= ((mChoiceRingSize.getSelectedIndex()+2) << Molecule.cAtomQFRingSizeShift);

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
					mDialog.showMessage("'H' cannot be part of an atom list and is removed.");
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
