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

public class BondQueryFeatureDialogBuilder implements GenericEventListener<GenericActionEvent> {
	private GenericDialog       mDialog;
    private ExtendedMolecule	mMol;
	private int					mBond,mFirstSpanItem;
	private GenericCheckBox     mCBSingle,mCBDouble,mCBTriple,mCBQuadruple,mCBQuintuple,mCBDelocalized,
								mCBMetalLigand,mCBIsBridge,mCBMatchFormalOrder,mCBMatchStereo;
	private GenericComboBox     mComboBoxRing,mComboBoxRingSize,mComboBoxMinAtoms,mComboBoxMaxAtoms;
	private boolean             mOKSelected;

	public BondQueryFeatureDialogBuilder(GenericUIHelper dialogHelper, ExtendedMolecule mol, int bond) {
		mDialog = dialogHelper.createDialog((mol.isSelectedAtom(mol.getBondAtom(0, bond))
				&& mol.isSelectedAtom(mol.getBondAtom(1, bond))) ?
				"Bond Query Features (Multiple)" : "Bond Query Features", this);
		build(mol, bond);
		}

	/**
	 * @return true if OK was pressed and potential change was applied to molecule
	 */
	public boolean showDialog() {
		mOKSelected = false;
		mDialog.showDialog();
		return mOKSelected;
		}

	private void build(ExtendedMolecule mol, int bond) {
		mMol = mol;
		mBond = bond;

		int gap = HiDPIHelper.scale(8);
		int[] hLayout = {gap, GenericDialog.FILL, GenericDialog.PREFERRED, GenericDialog.PREFERRED, gap};
        int[] vLayout = {gap, GenericDialog.PREFERRED, gap, GenericDialog.PREFERRED, GenericDialog.PREFERRED,
					        GenericDialog.PREFERRED, GenericDialog.PREFERRED,
		                    GenericDialog.PREFERRED, GenericDialog.PREFERRED, GenericDialog.PREFERRED, gap,
					        GenericDialog.PREFERRED, gap, GenericDialog.PREFERRED, gap, GenericDialog.PREFERRED, gap,
		                    GenericDialog.PREFERRED, 2*gap,
					        GenericDialog.PREFERRED, gap/2, GenericDialog.PREFERRED, gap/2, GenericDialog.PREFERRED, 2*gap};
        mDialog.setLayout(hLayout, vLayout);

		mDialog.add(mDialog.createLabel("Desired Bond type(s):"), 1, 1, 3, 1);

		mCBSingle = mDialog.createCheckBox("Single");
		mDialog.add(mCBSingle,1,3,3,3);

		mCBDouble = mDialog.createCheckBox("Double");
		mDialog.add(mCBDouble,1,4,3,4);

		mCBTriple = mDialog.createCheckBox("Triple");
		mDialog.add(mCBTriple,1,5,3,5);

		mCBQuadruple = mDialog.createCheckBox("Quadruple");
		mDialog.add(mCBQuadruple,1,6,3,6);

		mCBQuintuple = mDialog.createCheckBox("Quintuple");
		mDialog.add(mCBQuintuple,1,7,3,7);

		mCBDelocalized = mDialog.createCheckBox("Delocalized");
		mDialog.add(mCBDelocalized,1,8,3,8);

		mCBMetalLigand = mDialog.createCheckBox("Coordinate (0-order)");
		mDialog.add(mCBMetalLigand,1,9,3,9);

		mComboBoxRing = mDialog.createComboBox();
		mComboBoxRing.addItem("any ring state");
		mComboBoxRing.addItem("is not in a ring");
		mComboBoxRing.addItem("is any ring bond");
		mComboBoxRing.addItem("is non-aromatic ring bond");
		mComboBoxRing.addItem("is aromatic bond");
		mComboBoxRing.addEventConsumer(this);
		mDialog.add(mComboBoxRing,1,11,3,11);

		mComboBoxRingSize = mDialog.createComboBox();
		mComboBoxRingSize.addItem("any ring size");
		mComboBoxRingSize.addItem("is in 3-membered ring");
		mComboBoxRingSize.addItem("is in 4-membered ring");
		mComboBoxRingSize.addItem("is in 5-membered ring");
		mComboBoxRingSize.addItem("is in 6-membered ring");
		mComboBoxRingSize.addItem("is in 7-membered ring");
		mDialog.add(mComboBoxRingSize, 1,13,3,13);

		mCBMatchFormalOrder = mDialog.createCheckBox("Match formal bond order");
		mCBMatchFormalOrder.setSelected((mol.getBondQueryFeatures(bond) & Molecule.cBondQFMatchFormalOrder) != 0);
		mCBMatchFormalOrder.addEventConsumer(this);
		mDialog.add(mCBMatchFormalOrder, 1,15,3,15);

		mCBMatchStereo = mDialog.createCheckBox("Match Stereo Configuration");
		mCBMatchStereo.setSelected((mol.getBondQueryFeatures(bond) & Molecule.cBondQFMatchStereo) != 0);
		mCBMatchStereo.addEventConsumer(this);
		mDialog.add(mCBMatchStereo, 1,17,3,17);

		mCBIsBridge = mDialog.createCheckBox("Is atom bridge between");
        mCBIsBridge.addEventConsumer(this);
		mDialog.add(mCBIsBridge,1,19,3,19);

        mComboBoxMinAtoms = mDialog.createComboBox();
        int itemCount = (1 << Molecule.cBondQFBridgeMinBits);
        for (int i=0; i<itemCount; i++)
            mComboBoxMinAtoms.addItem(""+i);
		mDialog.add(mComboBoxMinAtoms,2,21);
		mDialog.add(mDialog.createLabel(" and"),3,21);
        mComboBoxMinAtoms.addEventConsumer(this);

        mComboBoxMaxAtoms = mDialog.createComboBox();
        populateComboBoxMaxAtoms(0);
		mDialog.add(mComboBoxMaxAtoms,2,23);
		mDialog.add(mDialog.createLabel(" atoms"),3,23);

        mMol.ensureHelperArrays(Molecule.cHelperRings);
        setInitialStates();
		}

	@Override
	public void eventHappened(GenericActionEvent e) {
		if (e.getWhat() == GenericActionEvent.WHAT_CANCEL) {
			mDialog.disposeDialog();
			}
		else if (e.getWhat() == GenericActionEvent.WHAT_OK) {
			setQueryFeatures();
			mOKSelected = true;
			mDialog.disposeDialog();
			}
        else if (e.getSource() == mCBIsBridge || e.getSource() == mComboBoxRing) {
            enableItems();
            }
        else if (e.getSource() == mComboBoxMinAtoms) {
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
		}


	private void setInitialStates() {
		int queryFeatures = mMol.getBondQueryFeatures(mBond);
		int bondOrder = (mMol.getBondType(mBond) == Molecule.cBondTypeDelocalized
					  || mMol.isDelocalizedBond(mBond)) ?
						6 : mMol.getBondOrder(mBond);

		if ((queryFeatures & Molecule.cBondQFSingle) != 0 || bondOrder == 1)
			mCBSingle.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFDouble) != 0 || bondOrder == 2)
			mCBDouble.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFTriple) != 0 || bondOrder == 3)
			mCBTriple.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFQuadruple) != 0 || bondOrder == 4)
			mCBQuadruple.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFQuintuple) != 0 || bondOrder == 5)
			mCBQuintuple.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFDelocalized) != 0 || bondOrder == 6)
			mCBDelocalized.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFMetalLigand) != 0 || bondOrder == 0)
			mCBMetalLigand.setSelected(true);
		if ((queryFeatures & Molecule.cBondQFMatchFormalOrder) != 0)
			mCBMatchFormalOrder.setSelected(true);
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
			else
				mComboBoxRing.setSelectedIndex(0);
			}
		else
			mComboBoxRing.setSelectedIndex(0);

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
        else
	        mComboBoxMaxAtoms.setSelectedIndex(0);

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
	    mCBQuadruple.setEnabled(!bridgeIsSelected);
	    mCBQuintuple.setEnabled(!bridgeIsSelected);
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

	    mMol.validateBondQueryFeatures();
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
	        if (mCBQuadruple.isSelected() && bondOrder != 4)
		        queryFeatures |= Molecule.cBondQFQuadruple;
	        if (mCBQuintuple.isSelected() && bondOrder != 5)
		        queryFeatures |= Molecule.cBondQFQuintuple;
    		if (mCBDelocalized.isSelected() && !mMol.isDelocalizedBond(bond) && bondOrder != 4)
    			queryFeatures |= Molecule.cBondQFDelocalized;
			if (mCBMetalLigand.isSelected() && bondOrder != 0)
				queryFeatures |= Molecule.cBondQFMetalLigand;
    		if (mCBMatchFormalOrder.isSelected())
    			queryFeatures |= Molecule.cBondQFMatchFormalOrder;
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
