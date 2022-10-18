package com.actelion.research.gui.editor;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;

import javax.swing.*;
import java.awt.*;

public class SwingEditorPanel extends JPanel {
	static final long serialVersionUID = 0x20211103;

	private SwingEditorToolbar mToolBar;
	private SwingEditorArea mArea;

	/**
	 * Creates a chemical editor panel to edit the given molecule.
	 * Query features can be edited, if the passed mol's fragment bit is true.
	 * @param mol
	 */
	public SwingEditorPanel(StereoMolecule mol) {
		this(mol, 0);
	}

	/**
	 * Creates a chemical editor panel to draw a molecule, reaction, or set of molecules,
	 * depending on the mode.
	 * Query features can be edited, if the passed mol's fragment bit is true.
	 * @param mol
	 * @param mode GenericDrawArea.MODE_...
	 */
	public SwingEditorPanel(StereoMolecule mol, int mode) {
		initialize(mol, mode);
	}

	/**
	 * Creates a chemical editor panel to edit multiple molecules.
	 * Each of these molecule may consist of multiple disconnected fragments.
	 * Atoms connected by bonds or being in close vicinity are recognized to belong to the
	 * same molecule, while more distant fragments are perceived as separated molecules.
	 * Query features can be edited, if the passed mols' fragment bits are true.
	 * @param mol
	 */
	public SwingEditorPanel(StereoMolecule[] mol) {
		initialize(null, GenericEditorArea.MODE_MULTIPLE_FRAGMENTS);
		mArea.getGenericDrawArea().setFragments(mol);
	}

	/**
	 * Creates a chemical editor panel to edit the given chemical reaction.
	 * Query features can be edited, if the passed rxn's fragment bit is true.
	 * @param rxn
	 */
	public SwingEditorPanel(Reaction rxn) {
		initialize(null, GenericEditorArea.MODE_REACTION);
		mArea.getGenericDrawArea().setReaction(rxn);
	}

	private void initialize(StereoMolecule mol, int mode) {
		setLayout(new BorderLayout());

		mArea = new SwingEditorArea(mol != null ? mol : new StereoMolecule(), mode);
		add(mArea, BorderLayout.CENTER);

		mToolBar = new SwingEditorToolbar(mArea);
		add(mToolBar, BorderLayout.WEST);
		}

	public GenericEditorArea getDrawArea() {
		return mArea.getGenericDrawArea();
	}

	public SwingEditorToolbar getSwingDrawToolbar() {
		return mToolBar;
	}

	public SwingEditorArea getSwingDrawArea() {
		return mArea;
	}

	public void cleanStructure(){
		mArea.getGenericDrawArea().toolChanged(GenericEditorToolbar.cButtonCleanStructure);
	}
}
