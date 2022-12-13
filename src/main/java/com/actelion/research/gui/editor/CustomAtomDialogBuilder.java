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

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.generic.*;

import javax.swing.*;

public class CustomAtomDialogBuilder implements GenericEventListener<GenericActionEvent> {
    private static final String[] RADICAL_STATES = { "None", "One electron (duplet)", "Two electrons (triplet)", "Two electrons (singulet)" };

	private GenericEditorArea mEditorArea;
	private GenericDialog mDialog;
    private StereoMolecule mMol;
    private int mAtom,mOldAtomicNo,mOldAtomMass,mOldAtomValence,mOldAtomRadical;
    private String mOldCustomLabel;
    private GenericTextField mTextFieldLabel,mTextFieldMass,mTextFieldValence;
    private GenericComboBox mComboBoxRadical;
    private boolean mOKSelected;

	public CustomAtomDialogBuilder(GenericUIHelper dialogHelper, GenericEditorArea editorArea,
	                               int atomicNo, int mass, int valence, int radical, String label) {
		mDialog = dialogHelper.createDialog("Atom Properties", this);
		mEditorArea = editorArea;
		mAtom = -1;
		mOldAtomicNo = atomicNo;
		mOldAtomMass = mass;
		mOldAtomValence = valence;
		mOldAtomRadical = radical;
		mOldCustomLabel = label;
		build();

	}

	public CustomAtomDialogBuilder(GenericUIHelper dialogHelper, GenericEditorArea editorArea, StereoMolecule mol, int atom) {
		mDialog = dialogHelper.createDialog("Atom Properties", this);
		mEditorArea = editorArea;
		mMol = mol;
		mAtom = atom;
		mOldAtomicNo = mMol.getAtomicNo(atom);
		mOldAtomMass = mMol.getAtomMass(atom);
		mOldAtomValence = mMol.getAtomAbnormalValence(atom);
		mOldAtomRadical = mMol.getAtomRadical(atom);
		mOldCustomLabel = mMol.getAtomCustomLabel(atom);
		build();
		}

	/**
	 * @return true if OK was pressed and potential change was applied to molecule
	 */
	public boolean showDialog() {
		mOKSelected = false;
		mDialog.showDialog();
		return mOKSelected;
		}

	private void build() {
        int[] hLayout = {8, GenericDialog.PREFERRED, 8, GenericDialog.PREFERRED, 8 };
		int[] vLayout = {8, GenericDialog.PREFERRED, 4, GenericDialog.PREFERRED,
                        12, GenericDialog.PREFERRED, 4, GenericDialog.PREFERRED,
                        12, GenericDialog.PREFERRED, 4, GenericDialog.PREFERRED,
                        12, GenericDialog.PREFERRED, 8};
        mDialog.setLayout(hLayout, vLayout);

        mTextFieldLabel = mDialog.createTextField(1,1);
        mTextFieldLabel.addEventConsumer(this);
		mDialog.add(mDialog.createLabel("Atom Label:"), 1,1);
		mDialog.add(mTextFieldLabel, 3,1);
		mDialog.add(mDialog.createLabel("(examples: 'D', 'Li', 'Cys', 'R12', 'R3@C')"), 1,3,3,3);

        mTextFieldMass = mDialog.createTextField(1,1);
        mTextFieldMass.addEventConsumer(this);
		mDialog.add(mDialog.createLabel("Atom Mass:"), 1,5);
		mDialog.add(mTextFieldMass, 3,5);
		mDialog.add(mDialog.createLabel("(empty for natural abundance)"), 1,7,3,7);

        mTextFieldValence = mDialog.createTextField(1,1);
        mTextFieldValence.addEventConsumer(this);
		mDialog.add(mDialog.createLabel("Abnormal Valence:"), 1,9);
		mDialog.add(mTextFieldValence, 3,9);
		mDialog.add(mDialog.createLabel("(empty for default valence)"), 1,11,3,11);

        if (mAtom == -1) {
	        String label = Molecule.cAtomLabel[mOldAtomicNo];
	        mTextFieldLabel.setText(mOldCustomLabel == null ? label : mOldCustomLabel+"@"+label);
	        if (mOldAtomMass != 0)
		        mTextFieldMass.setText(""+mOldAtomMass);
	        if (mOldAtomValence != -1)
		        mTextFieldValence.setText(""+mOldAtomValence);
	        }
        else {
        	String label = mMol.getAtomLabel(mAtom);
        	String customLabel = mMol.getAtomCustomLabel(mAtom);
            mTextFieldLabel.setText(customLabel == null ? label : customLabel+"@"+label);
            if (mMol.getAtomMass(mAtom) != 0)
                mTextFieldMass.setText(""+mMol.getAtomMass(mAtom));
            if (mMol.getAtomAbnormalValence(mAtom) != -1)
                mTextFieldValence.setText(""+mMol.getAtomAbnormalValence(mAtom));
            }

		mComboBoxRadical = mDialog.createComboBox();
        for (String s:RADICAL_STATES)
	        mComboBoxRadical.addItem(s);
        int state = (mAtom == -1) ? mOldAtomRadical : mMol.getAtomRadical(mAtom);
        mComboBoxRadical.setSelectedIndex(state == Molecule.cAtomRadicalStateD ? 1 :
                                          state == Molecule.cAtomRadicalStateT ? 2 :
                                          state == Molecule.cAtomRadicalStateS ? 3 : 0);
		mDialog.add(mDialog.createLabel("Radical State:"), 1,13);
		mDialog.add(mComboBoxRadical, 3,13);
		}

	@Override
	public void eventHappened(GenericActionEvent e) {
		if (e.getSource() instanceof JTextField) {
			processAtomLabel(false);
			}
		else if (e.getWhat() == GenericActionEvent.WHAT_CANCEL) {
			if (mAtom != -1) {
				mMol.setAtomicNo(mAtom, mOldAtomicNo);
				mMol.setAtomMass(mAtom, mOldAtomMass);
				mMol.setAtomAbnormalValence(mAtom, mOldAtomValence);
				mMol.setAtomRadical(mAtom, mOldAtomRadical);
				mMol.setAtomCustomLabel(mAtom, mOldCustomLabel);
				}
			mDialog.disposeDialog();
			}
		else if (e.getWhat() == GenericActionEvent.WHAT_OK) {
			processAtomLabel(true);
			mDialog.disposeDialog();
			}
		}

	private void processAtomLabel(boolean updateDefault) {
		String text = mTextFieldLabel.getText();
		String customLabel = null;

		if (text.length() != 0) {
			int index = text.indexOf('@');
			if (index != -1) {
				customLabel = text.substring(0, index);
				text = text.substring(index+1);
				}
			}

		if (text.length() != 0) {
			int atomicNo = Molecule.getAtomicNoFromLabel(text, mEditorArea.getAllowedPseudoAtoms());
			if (atomicNo != 0 || text.equals("?")) {
			    int mass = 0;
			    if (mTextFieldMass.getText().length() != 0) {
    			    try {
    			        mass = Integer.parseInt(mTextFieldMass.getText());
    	                if (mass < Molecule.cRoundedMass[atomicNo] - 18
                         || mass > Molecule.cRoundedMass[atomicNo] + 12) {
    	                	mDialog.showMessage("Your mass is out of range!");
    	                    return;
    	                    }
    			        }
    			    catch (NumberFormatException nfe) {
				        mDialog.showMessage("Your mass is not a number!");
    			        return;
    			        }
			        }

                int valence = -1;
			    if (mTextFieldValence.getText().length() != 0) {
                    try {
                        valence = Integer.parseInt(mTextFieldValence.getText());
                        if (valence < 0 || valence > 15) {
	                        mDialog.showMessage("Your valence is out of range!");
                            return;
                            }
                        }
                    catch (NumberFormatException nfe) {
	                    mDialog.showMessage("Your valence is not a number!");
                        return;
                        }
                    }

                int	radical = mComboBoxRadical.getSelectedIndex() == 1 ? Molecule.cAtomRadicalStateD :
                			  mComboBoxRadical.getSelectedIndex() == 2 ? Molecule.cAtomRadicalStateT :
                   			  mComboBoxRadical.getSelectedIndex() == 3 ? Molecule.cAtomRadicalStateS : 0;

			    // set the current property set for the custom atom
				if (updateDefault)
				    mEditorArea.setCustomAtom(atomicNo, mass, valence, radical,customLabel);

			    if (mAtom != -1) {
				    mMol.changeAtom(mAtom, atomicNo, mass, valence, radical);
				    mMol.setAtomCustomLabel(mAtom, customLabel);
			        }

                mOKSelected = true;
                mDialog.disposeDialog();
			    }
			}
		}
	}
