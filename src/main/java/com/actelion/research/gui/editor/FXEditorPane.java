package com.actelion.research.gui.editor;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.Pane;

public class FXEditorPane extends BorderPane {
	private FXEditorToolbar mToolBar;
	private FXEditorArea mArea;

	/**
	 * Creates a chemical editor panel to edit the given molecule.
	 * Query features can be edited, if the passed mol's fragment bit is true.
	 * @param mol
	 */
	public FXEditorPane(StereoMolecule mol) {
		this(mol, 0);
	}

	/**
	 * Creates a chemical editor panel to draw a molecule, reaction, or set of molecules,
	 * depending on the mode.
	 * Query features can be edited, if the passed mol's fragment bit is true.
	 * @param mol
	 * @param mode GenericDrawArea.MODE_...
	 */
	public FXEditorPane(StereoMolecule mol, int mode) {
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
	public FXEditorPane(StereoMolecule[] mol) {
		initialize(null, GenericEditorArea.MODE_MULTIPLE_FRAGMENTS);
		mArea.getGenericDrawArea().setFragments(mol);
	}

	/**
	 * Creates a chemical editor panel to edit the given chemical reaction.
	 * Query features can be edited, if the passed rxn's fragment bit is true.
	 * @param rxn
	 */
	public FXEditorPane(Reaction rxn) {
		initialize(null, GenericEditorArea.MODE_REACTION);
		mArea.getGenericDrawArea().setReaction(rxn);
	}

	private void initialize(StereoMolecule mol, int mode) {
		mArea = new FXEditorArea(mol != null ? mol : new StereoMolecule(), mode);
		mToolBar = new FXEditorToolbar(mArea, mode);
		widthProperty().addListener((observable, oldValue, newValue) -> mArea.setWidth((double) newValue-mToolBar.getWidth()));
		heightProperty().addListener((observable, oldValue, newValue) -> mArea.setHeight((double) newValue));
		setLeft(mToolBar);
		setCenter(new Pane(mArea));
	}

	public GenericEditorArea getDrawArea() {
		return mArea.getGenericDrawArea();
	}

	public FXEditorArea getFXDrawArea() {
		return mArea;
	}

	public void cleanStructure(){
		mArea.getGenericDrawArea().toolChanged(GenericEditorToolbar.cButtonCleanStructure);
	}
}
