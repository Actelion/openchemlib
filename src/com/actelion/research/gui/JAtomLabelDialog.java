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

import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.Molecule;

public class JAtomLabelDialog extends JDialog implements ActionListener {
    public static final long serialVersionUID = 0x20080213;

    private static final String[] RADICAL_STATES = { "None", "One electron (duplet)", "Two electrons (triplet)", "Two electrons (singulet)" };

    private Frame mOwner;
    private ExtendedMolecule mMol;
    private int mAtom;
    private JTextField mTextFieldLabel,mTextFieldMass,mTextFieldValence;
    private JComboBox mComboBoxRadical;

	protected JAtomLabelDialog(Frame owner, ExtendedMolecule mol, int atom) {
		super(owner, true);
		mOwner = owner;
		mMol = mol;
		mAtom = atom;
		init();
		}


	private void init() {
	    JPanel p1 = new JPanel();
        double[][] size = { {8, TableLayout.PREFERRED, 8, TableLayout.PREFERRED, 8 },
                            {8, TableLayout.PREFERRED, 4, TableLayout.PREFERRED,
                             12, TableLayout.PREFERRED, 4, TableLayout.PREFERRED,
                             12, TableLayout.PREFERRED, 4, TableLayout.PREFERRED,
                             12, TableLayout.PREFERRED, 8} };
        p1.setLayout(new TableLayout(size));

        mTextFieldLabel = new JTextField(1);
        mTextFieldLabel.addActionListener(this);
		p1.add(new JLabel("Atom Label:", JLabel.RIGHT), "1,1");
		p1.add(mTextFieldLabel, "3,1");
        p1.add(new JLabel("(examples: 'D', 'Li', 'Cys', 'R12', 'R3@C')", JLabel.RIGHT), "1,3,3,3");

        mTextFieldMass = new JTextField(1);
        mTextFieldMass.addActionListener(this);
        p1.add(new JLabel("Atom Mass:", JLabel.RIGHT), "1,5");
        p1.add(mTextFieldMass, "3,5");
        p1.add(new JLabel("(empty for natural abundance)", JLabel.RIGHT), "1,7,3,7");

        mTextFieldValence = new JTextField(1);
        mTextFieldValence.addActionListener(this);
        p1.add(new JLabel("Abnormal Valence:", JLabel.RIGHT), "1,9");
        p1.add(mTextFieldValence, "3,9");
        p1.add(new JLabel("(empty for default valence)", JLabel.RIGHT), "1,11,3,11");

        if (mAtom != -1) {
        	String label = mMol.getAtomLabel(mAtom);
        	String customLabel = mMol.getAtomCustomLabel(mAtom);
            mTextFieldLabel.setText(customLabel == null ? label : customLabel+"@"+label);
            if (mMol.getAtomMass(mAtom) != 0)
                mTextFieldMass.setText(""+mMol.getAtomMass(mAtom));
            if (mMol.getAtomAbnormalValence(mAtom) != -1)
                mTextFieldValence.setText(""+mMol.getAtomAbnormalValence(mAtom));
            }

        mComboBoxRadical = new JComboBox(RADICAL_STATES);
        if (mAtom != -1) {
        	int state = mMol.getAtomRadical(mAtom);
        	mComboBoxRadical.setSelectedIndex(state == Molecule.cAtomRadicalStateD ? 1 :
        									  state == Molecule.cAtomRadicalStateT ? 2 :
        									  state == Molecule.cAtomRadicalStateS ? 3 : 0);
        	}
        p1.add(new JLabel("Radical State:", JLabel.RIGHT), "1,13");
        p1.add(mComboBoxRadical, "3,13");

        JPanel p2 = new JPanel();
        p2.setBorder(BorderFactory.createEmptyBorder(12, 8, 8, 8));
        p2.setLayout(new BorderLayout());
        JPanel ibp = new JPanel();
        ibp.setLayout(new GridLayout(1, 2, 8, 0));
        JButton bcancel = new JButton("Cancel");
        bcancel.addActionListener(this);
        ibp.add(bcancel);
        JButton bok = new JButton("OK");
        bok.addActionListener(this);
        ibp.add(bok);
        p2.add(ibp, BorderLayout.EAST);

		getContentPane().add(p1, BorderLayout.CENTER);
		getContentPane().add(p2, BorderLayout.SOUTH);
		getRootPane().setDefaultButton(bok);

		pack();
		setLocationRelativeTo(mOwner);
        setVisible(true);
		}


	public void actionPerformed(ActionEvent e) {
		if (e.getSource() instanceof JTextField)
			processAtomLabel();
		else if (e.getActionCommand() == "Cancel")
			dispose();
		else if (e.getActionCommand() == "OK")
			processAtomLabel();
		}

	private void processAtomLabel() {
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
			int atomicNo = Molecule.getAtomicNoFromLabel(text);
			if (atomicNo != 0 || text.equals("?")) {
			    int mass = 0;
			    if (mTextFieldMass.getText().length() != 0) {
    			    try {
    			        mass = Integer.parseInt(mTextFieldMass.getText());
    	                if (mass < Molecule.cRoundedMass[atomicNo] - 18
                         || mass > Molecule.cRoundedMass[atomicNo] + 12) {
    	                    JOptionPane.showMessageDialog(mOwner, "Your mass is out of range!");
    	                    return;
    	                    }
    			        }
    			    catch (NumberFormatException nfe) {
                        JOptionPane.showMessageDialog(mOwner, "Your mass is not a number!");
    			        return;
    			        }
			        }

                int valence = -1;
			    if (mTextFieldValence.getText().length() != 0) {
                    try {
                        valence = Integer.parseInt(mTextFieldValence.getText());
                        if (valence < 0 || valence > 15) {
                            JOptionPane.showMessageDialog(mOwner, "Your valence is out of range!");
                            return;
                            }
                        }
                    catch (NumberFormatException nfe) {
                        JOptionPane.showMessageDialog(mOwner, "Your valence is not a number!");
                        return;
                        }
                    }

                int	radical = mComboBoxRadical.getSelectedIndex() == 1 ? Molecule.cAtomRadicalStateD :
                			  mComboBoxRadical.getSelectedIndex() == 2 ? Molecule.cAtomRadicalStateT :
                   			  mComboBoxRadical.getSelectedIndex() == 3 ? Molecule.cAtomRadicalStateS : 0;

                mMol.changeAtom(mAtom, atomicNo, mass, valence, radical);
                if (customLabel != null)
                	mMol.setAtomCustomLabel(mAtom, customLabel);

                dispose();
			    }
			}
		}
	}
