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

public class AtomQueryFeatureDialogBuilder implements GenericEventListener<GenericActionEvent> {
	public static final String[] RING_SIZE_SHORT_TEXT = {
			"",
			"r0",
			"r",
			"r3",
			"r3-4",
			"r3-5",
			"r3-6",
			"r3-7",
			"r4",
			"r4-5",
			"r4-6",
			"r4-7",
			"r>3",
			"r5",
			"r5-6",
			"r5-7",
			"r>4",
			"r6",
			"r6-7",
			"r>5",
			"r7",
			"r>6",
			"r>7"
	};
	private static final String[] RING_SIZE_OPTIONS = {
			"any",
			"is not in a ring",
			"is in any ring",
			"3 members",
			"3-4 members",
			"3-5 members",
			"3-6 members",
			"3-7 members",
			"4 members",
			"4-5 members",
			"4-6 members",
			"4-7 members",
			"> 3 members",
			"5 members",
			"5-6 members",
			"5-7 members",
			"> 4 members",
			"6 members",
			"6-7 members",
			"> 5 members",
			"7 members",
			"> 6 members",
			"> 7 members"
			};
	public static final long[] RING_SIZE_VALUES = {
			0,    // special handling: 0 means all allowed
			Molecule.cAtomQFRingSize0,
			Molecule.cAtomQFNewRingSize ^ Molecule.cAtomQFRingSize0,
			Molecule.cAtomQFRingSize3,
			Molecule.cAtomQFRingSize3 | Molecule.cAtomQFRingSize4,
			Molecule.cAtomQFRingSize3 | Molecule.cAtomQFRingSize4 | Molecule.cAtomQFRingSize5,
			Molecule.cAtomQFRingSize3 | Molecule.cAtomQFRingSize4 | Molecule.cAtomQFRingSize5 | Molecule.cAtomQFRingSize6,
			Molecule.cAtomQFNewRingSize ^ (Molecule.cAtomQFRingSize0 | Molecule.cAtomQFRingSizeLarge),
			Molecule.cAtomQFRingSize4,
			Molecule.cAtomQFRingSize4 | Molecule.cAtomQFRingSize5,
			Molecule.cAtomQFRingSize4 | Molecule.cAtomQFRingSize5 | Molecule.cAtomQFRingSize6,
			Molecule.cAtomQFRingSize4 | Molecule.cAtomQFRingSize5 | Molecule.cAtomQFRingSize6 | Molecule.cAtomQFRingSize7,
			Molecule.cAtomQFNewRingSize ^ (Molecule.cAtomQFRingSize0 | Molecule.cAtomQFRingSize3),
			Molecule.cAtomQFRingSize5,
			Molecule.cAtomQFRingSize5 | Molecule.cAtomQFRingSize6,
			Molecule.cAtomQFRingSize5 | Molecule.cAtomQFRingSize6 | Molecule.cAtomQFRingSize7,
			Molecule.cAtomQFRingSize5 | Molecule.cAtomQFRingSize6 | Molecule.cAtomQFRingSize7 | Molecule.cAtomQFRingSizeLarge,
			Molecule.cAtomQFRingSize6,
			Molecule.cAtomQFRingSize6 | Molecule.cAtomQFRingSize7,
			Molecule.cAtomQFRingSize6 | Molecule.cAtomQFRingSize7 | Molecule.cAtomQFRingSizeLarge,
			Molecule.cAtomQFRingSize7,
			Molecule.cAtomQFRingSize7 | Molecule.cAtomQFRingSizeLarge,
			Molecule.cAtomQFRingSizeLarge,
			};

	private GenericDialog       mDialog;
	private GenericLabel        mLabelAtomList;
	private GenericTextField    mTFAtomList;
	private GenericCheckBox     mCBAny,mCBBlocked,mCBSubstituted,mCBMatchStereo,mCBExcludeGroup;
	private GenericComboBox     mChoiceArom,mChoiceRingState,mChoiceSmallRingSize,mChoiceRingSize,mChoiceCharge,
								mChoiceNeighbours,mChoiceHydrogen,mChoicePi,mChoiceReactionParityHint;
    private ExtendedMolecule	mMol;
	private int					mAtom;
	private long                mRingSizeCustomValue;
	private boolean             mOKSelected;

	public AtomQueryFeatureDialogBuilder(GenericUIHelper dialogHelper, ExtendedMolecule mol, int atom, boolean includeReactionHints) {
		mDialog = dialogHelper.createDialog(mol.isSelectedAtom(atom) ? "Atom Query Features (Multiple)" : "Atom Query Features", this);
		build(mol, atom, includeReactionHints);
		}

	/**
	 * @return true if OK was pressed and potential change was applied to molecule
	 */
	public boolean showDialog() {
		mOKSelected = false;
		mDialog.showDialog();

		return mOKSelected;
		}

	private void build(ExtendedMolecule mol, int atom, boolean includeReactionHints) {
		mMol = mol;
		mAtom = atom;

		int gap = HiDPIHelper.scale(8);
		int[] gap1 = { gap, gap/2, gap*3/2, gap/2, gap/2, gap/2, gap/2, gap/2, gap/2, gap/2, gap*3/2, gap/4, gap/4, gap/4 };
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
		
		mCBAny = mDialog.createCheckBox("Any atomic number");
		mCBAny.addEventConsumer(this);
		mDialog.add(mCBAny, 1,1,3,1);

		mLabelAtomList = mDialog.createLabel("Excluded atoms:");
		mTFAtomList = mDialog.createTextField(16, 1);

		mDialog.add(mLabelAtomList, 1,3);
		mDialog.add(mTFAtomList, 3,3);

		mChoiceArom = mDialog.createComboBox();
		mChoiceArom.addItem("any");
		mChoiceArom.addItem("is aromatic");
		mChoiceArom.addItem("is not aromatic");
		mDialog.add(mDialog.createLabel("Aromaticity:"), 1,5);
		mDialog.add(mChoiceArom, 3,5);

		mChoiceRingState = mDialog.createComboBox();
		mChoiceRingState.addItem("any");
		mChoiceRingState.addItem("none (not in a ring)");
		mChoiceRingState.addItem("at least 2");
		mChoiceRingState.addItem("exactly 2");
		mChoiceRingState.addItem("exactly 3");
		mChoiceRingState.addItem("more than 3");
		mDialog.add(mDialog.createLabel("Ring bonds:"), 1,7);
		mDialog.add(mChoiceRingState, 3,7);

		mChoiceSmallRingSize = mDialog.createComboBox();
		mChoiceSmallRingSize.addItem("any");
		mChoiceSmallRingSize.addItem("3 members");
		mChoiceSmallRingSize.addItem("4 members");
		mChoiceSmallRingSize.addItem("5 members");
		mChoiceSmallRingSize.addItem("6 members");
		mChoiceSmallRingSize.addItem("7 members");
		mDialog.add(mDialog.createLabel("Smallest ring size:"), 1,9);
		mDialog.add(mChoiceSmallRingSize, 3,9);

		mChoiceRingSize = mDialog.createComboBox();
		for (String option:RING_SIZE_OPTIONS)
			mChoiceRingSize.addItem(option);
		mDialog.add(mDialog.createLabel("Any ring size:"), 1,11);
		mDialog.add(mChoiceRingSize, 3,11);

		mChoiceCharge = mDialog.createComboBox();
		mChoiceCharge.addItem("any");
		mChoiceCharge.addItem("not charged");
		mChoiceCharge.addItem("has negative charge");
		mChoiceCharge.addItem("has positive charge");
		mDialog.add(mDialog.createLabel("Charge:"), 1,13);
		mDialog.add(mChoiceCharge, 3,13);

		mChoiceNeighbours = mDialog.createComboBox();
		mChoiceNeighbours.addItem("any");
		mChoiceNeighbours.addItem("exactly 1");
        mChoiceNeighbours.addItem("exactly 2");
        mChoiceNeighbours.addItem("exactly 3");
        mChoiceNeighbours.addItem("less than 3");
        mChoiceNeighbours.addItem("less than 4");
		mChoiceNeighbours.addItem("at least 1");
        mChoiceNeighbours.addItem("at least 2");
        mChoiceNeighbours.addItem("at least 3");
        mChoiceNeighbours.addItem("at least 4");
		mDialog.add(mDialog.createLabel("Non-H neighbours:"), 1,15);
		mDialog.add(mChoiceNeighbours, 3,15);

		mChoiceHydrogen = mDialog.createComboBox();
		mChoiceHydrogen.addItem("any");
		mChoiceHydrogen.addItem("none");
		mChoiceHydrogen.addItem("exactly 1");
        mChoiceHydrogen.addItem("exactly 2");
		mChoiceHydrogen.addItem("at least 1");
		mChoiceHydrogen.addItem("at least 2");
		mChoiceHydrogen.addItem("at least 3");
        mChoiceHydrogen.addItem("less than 2");
        mChoiceHydrogen.addItem("less than 3");
		mDialog.add(mDialog.createLabel("Hydrogen count:"), 1,17);
		mDialog.add(mChoiceHydrogen, 3,17);

        mChoicePi = mDialog.createComboBox();
        mChoicePi.addItem("any");
        mChoicePi.addItem("none");
        mChoicePi.addItem("exactly 1");
        mChoicePi.addItem("exactly 2");
        mChoicePi.addItem("at least 1");
		mDialog.add(mDialog.createLabel("Pi-electron count:"), 1,19);
		mDialog.add(mChoicePi,3,19);

		mCBBlocked = mDialog.createCheckBox("prohibit further substitution");
		mCBBlocked.addEventConsumer(this);
		mDialog.add(mCBBlocked, 1,21,3,21);

		mCBSubstituted = mDialog.createCheckBox("require further substitution");
		mCBSubstituted.addEventConsumer(this);
		mDialog.add(mCBSubstituted, 1,23,3,23);

		mCBMatchStereo = mDialog.createCheckBox("match stereo center");
		mDialog.add(mCBMatchStereo, 1,25,3,25);

		mCBExcludeGroup = mDialog.createCheckBox("is part of exclude group");
		mDialog.add(mCBExcludeGroup, 1,27,3,27);

		if (includeReactionHints) {
			mDialog.add(mDialog.createLabel("Stereo center hint for product:"), 1,29,3,29);
			mChoiceReactionParityHint = mDialog.createComboBox();
			mChoiceReactionParityHint.addItem("Copy from generic product");
			mChoiceReactionParityHint.addItem("Keep reactant configuration");
			mChoiceReactionParityHint.addItem("Invert reactant configuration");
			mChoiceReactionParityHint.addItem("Racemise configuration");
			mDialog.add(mChoiceReactionParityHint, 1,31,3,31);
			}

		mMol.ensureHelperArrays(Molecule.cHelperCIP);
		setInitialStates();
		}

	@Override
	public void eventHappened(GenericActionEvent e) {
		if (e.getWhat() == GenericActionEvent.WHAT_OK) {
			setQueryFeatures();
			mOKSelected = true;
			mDialog.disposeDialog();
			}
		else if (e.getWhat() == GenericActionEvent.WHAT_CANCEL) {
			mDialog.disposeDialog();
			}
		else if (e.getSource() == mCBAny) {
			if (e.getValue() == 1) {
				mTFAtomList.setText("");
				mLabelAtomList.setText("Excluded atoms:");
				}
			else {
				mTFAtomList.setText(mMol.getAtomLabel(mAtom));
				mLabelAtomList.setText("Allowed atoms:");
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
		long queryFeatures = mMol.getAtomQueryFeatures(mAtom);

		if ((queryFeatures & Molecule.cAtomQFAny) != 0) {
			mCBAny.setSelected(true);
			mLabelAtomList.setText("Excluded atoms:");
			}
		else
			mLabelAtomList.setText("Allowed atoms:");

		mTFAtomList.setText(mMol.getAtomList(mAtom) == null ? "" : mMol.getAtomListString(mAtom));

		long aromState = queryFeatures & Molecule.cAtomQFAromState;
		if (aromState == Molecule.cAtomQFAromatic)
			mChoiceArom.setSelectedIndex(1);
		else if (aromState == Molecule.cAtomQFNotAromatic)
			mChoiceArom.setSelectedIndex(2);
		else
			mChoiceArom.setSelectedIndex(0);

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
		else
			mChoiceRingState.setSelectedIndex(0);

		int smallRingSize = (int)((queryFeatures & Molecule.cAtomQFSmallRingSize) >> Molecule.cAtomQFSmallRingSizeShift);
		mChoiceSmallRingSize.setSelectedIndex((smallRingSize == 0) ? 0 : smallRingSize-2);

		long ringSize = queryFeatures & Molecule.cAtomQFNewRingSize;
		int index = -1;
		for (int i=0; i<RING_SIZE_VALUES.length; i++) {
			if (ringSize == RING_SIZE_VALUES[i]) {
				index = i;
				break;
				}
			}
		if (index != -1) {
			mChoiceRingSize.setSelectedIndex(index);
			}
		else {
			StringBuilder customOption = new StringBuilder("Custom:");
			if ((ringSize & Molecule.cAtomQFRingSize0) != 0)
				customOption.append(" 0");
			if ((ringSize & Molecule.cAtomQFRingSize3) != 0)
				customOption.append(" 3");
			if ((ringSize & Molecule.cAtomQFRingSize4) != 0)
				customOption.append(" 4");
			if ((ringSize & Molecule.cAtomQFRingSize5) != 0)
				customOption.append(" 5");
			if ((ringSize & Molecule.cAtomQFRingSize6) != 0)
				customOption.append(" 6");
			if ((ringSize & Molecule.cAtomQFRingSize7) != 0)
				customOption.append(" 7");
			if ((ringSize & Molecule.cAtomQFRingSizeLarge) != 0)
				customOption.append(" >=8");
			mRingSizeCustomValue = ringSize;
			mChoiceRingSize.addItem(customOption.toString());
			mChoiceRingSize.setSelectedIndex(RING_SIZE_VALUES.length);
			}

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
		else if (neighbourFeatures == Molecule.cAtomQFNot0Neighbours)
			mChoiceNeighbours.setSelectedIndex(6);
        else if (neighbourFeatures == (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour))
            mChoiceNeighbours.setSelectedIndex(7);
		else if (neighbourFeatures == (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour | Molecule.cAtomQFNot2Neighbours))
            mChoiceNeighbours.setSelectedIndex(8);
        else if (neighbourFeatures == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot4Neighbours))
            mChoiceNeighbours.setSelectedIndex(9);
		else
			mChoiceNeighbours.setSelectedIndex(0);

		long chargeFeatures = queryFeatures & Molecule.cAtomQFCharge;
		if (chargeFeatures == (Molecule.cAtomQFNotChargeNeg | Molecule.cAtomQFNotChargePos))
			mChoiceCharge.setSelectedIndex(1);
		else if (chargeFeatures == (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargePos))
			mChoiceCharge.setSelectedIndex(2);
		else if (chargeFeatures == (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargeNeg))
			mChoiceCharge.setSelectedIndex(3);
		else
			mChoiceCharge.setSelectedIndex(0);

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
        else if (hydrogenFeatures == (Molecule.cAtomQFNot3Hydrogen))
            mChoiceHydrogen.setSelectedIndex(8);
		else
			mChoiceHydrogen.setSelectedIndex(0);

		long piFeatures = queryFeatures & Molecule.cAtomQFPiElectrons;
        if (piFeatures == (Molecule.cAtomQFNot1PiElectron | Molecule.cAtomQFNot2PiElectrons))
            mChoicePi.setSelectedIndex(1);
        else if (piFeatures == (Molecule.cAtomQFNot0PiElectrons | Molecule.cAtomQFNot2PiElectrons))
            mChoicePi.setSelectedIndex(2);
        else if (piFeatures == (Molecule.cAtomQFNot0PiElectrons | Molecule.cAtomQFNot1PiElectron))
            mChoicePi.setSelectedIndex(3);
        else if (piFeatures == Molecule.cAtomQFNot0PiElectrons)
            mChoicePi.setSelectedIndex(4);
	    else
		    mChoicePi.setSelectedIndex(0);

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
			else
				mChoiceReactionParityHint.setSelectedIndex(0);
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
		long queryFeatures = 0;

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

	    if (mChoiceSmallRingSize.getSelectedIndex() != 0)
		    queryFeatures |= (long)(mChoiceSmallRingSize.getSelectedIndex()+2) << Molecule.cAtomQFSmallRingSizeShift;

	    int ringSizeIndex = mChoiceRingSize.getSelectedIndex();
        if (ringSizeIndex == RING_SIZE_VALUES.length)
	        queryFeatures |= mRingSizeCustomValue;
        else if (ringSizeIndex != 0)
	        queryFeatures |= RING_SIZE_VALUES[ringSizeIndex];

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
        case 6: // more than 0 non-H neighbour
	        if (mMol.getConnAtoms(atom) == 0)
		        queryFeatures |= Molecule.cAtomQFMoreNeighbours;
	        break;
        case 7: // more than 1 non-H neighbour
            if (mMol.getConnAtoms(atom) == 1)
                queryFeatures |= Molecule.cAtomQFMoreNeighbours;
            else if (mMol.getConnAtoms(atom) < 1)
                queryFeatures |= (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour);
            break;
        case 8: // more than 2 non-H neighbours
            if (mMol.getConnAtoms(atom) == 2)
                queryFeatures |= Molecule.cAtomQFMoreNeighbours;
            else if (mMol.getConnAtoms(atom) < 2)
                queryFeatures |= (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour | Molecule.cAtomQFNot2Neighbours);
            break;
        case 9: // more than 3 non-H neighbours
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
