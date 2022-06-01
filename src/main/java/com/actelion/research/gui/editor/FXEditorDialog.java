package com.actelion.research.gui.editor;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.StructureListener;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import javafx.event.ActionEvent;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.SceneAntialiasing;
import javafx.scene.control.Button;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.GridPane;
import javafx.stage.Modality;
import javafx.stage.Stage;
import javafx.stage.StageStyle;
import javafx.stage.Window;

import java.util.ArrayList;

public class FXEditorDialog extends Stage {
	private static final String DEFAULT_MOLECULE_TITLE = "Structure Editor";
	private static final String DEFAULT_REACTION_TITLE = "Reaction Editor";

	private StereoMolecule mMolecule;
	private FXEditorArea mArea;
	private boolean        mIsCancelled;
	private ArrayList<StructureListener> mListener;

	/**
	 * Creates a modal chemical editor dialog to edit a single molecule,
	 * which may, of course, consist of multiple disconnected fragments.
	 * Query features can be edited, if the passed mol's fragment bit is true.
	 * @param owner
	 * @param mol
	 */
	public FXEditorDialog(Window owner, StereoMolecule mol) {
		this(owner, mol, Modality.APPLICATION_MODAL);
	}

	/**
	 * Creates a modal chemical editor dialog to edit a chemical reaction.
	 * Query features can be edited, if the passed rxn's fragment bit is true.
	 * @param owner
	 * @param rxn
	 */
	public FXEditorDialog(Window owner, Reaction rxn) {
		this(owner, rxn, Modality.APPLICATION_MODAL);
	}

	/**
	 * Creates a chemical editor dialog to edit a single molecule,
	 * which may, of course, consist of multiple disconnected fragments.
	 * Query features can be edited, if the passed mol's fragment bit is true.
	 * @param owner
	 * @param mol
	 * @param modality
	 */
	public FXEditorDialog(Window owner, StereoMolecule mol, Modality modality) {
		super(StageStyle.DECORATED);
		mMolecule = (mol == null) ? new StereoMolecule() : new StereoMolecule(mol);
		initialize(owner, modality, DEFAULT_MOLECULE_TITLE, 0);
	}

	/**
	 * Creates a modal chemical editor dialog to edit multiple molecules in one editor pane.
	 * Each of these molecule may consist of multiple disconnected fragments.
	 * Atoms connected by bonds or being in close vicinity are recognized to belong to the
	 * same molecule, while more distant fragments are perceived as separated molecules.
	 * Query features can be edited, if the passed mols' fragment bits are true.
	 * @param owner
	 * @param mol
	 * @param modality
	 */
	public FXEditorDialog(Window owner, StereoMolecule[] mol, Modality modality) {
		super(StageStyle.DECORATED);
		mMolecule = new StereoMolecule();
		initialize(owner, modality, DEFAULT_MOLECULE_TITLE, GenericEditorArea.MODE_MULTIPLE_FRAGMENTS);
		if (mol != null)
			mArea.getGenericDrawArea().setFragments(mol);
	}

	/**
	 * Creates a modal chemical editor dialog to edit a chemical reaction.
	 * Query features can be edited, if the passed rxn's fragment bit is true.
	 * @param owner
	 * @param rxn
	 * @param modality
	 */
	public FXEditorDialog(Window owner, Reaction rxn, Modality modality) {
		super(StageStyle.DECORATED);
		mMolecule = new StereoMolecule();
		initialize(owner, modality, DEFAULT_REACTION_TITLE, GenericEditorArea.MODE_REACTION);
		if (rxn != null)
			mArea.getGenericDrawArea().setReaction(rxn);
	}

	public void initialize(Window owner, Modality modality, String title, int mode) {
		initOwner(owner);
		initModality(modality);
		setTitle(title);

		FXEditorPane editorPane =  new FXEditorPane(mMolecule, mode);
		mArea = editorPane.getFXDrawArea();

		Button buttonHelp = new Button("Help");
		Button buttonCancel = new Button("Cancel");
		Button buttonOK = new Button("OK");
		buttonHelp.addEventHandler(ActionEvent.ACTION, e -> mArea.getGenericDrawArea().showHelpDialog() );
		buttonCancel.addEventHandler(ActionEvent.ACTION, e -> close() );
		buttonOK.addEventHandler(ActionEvent.ACTION, e -> {
				for (StructureListener sl : mListener)
					sl.structureChanged(mMolecule);
				mIsCancelled = false;
				close();
				} );

		int gap = HiDPIHelper.scale(8);
		GridPane okPane = new GridPane();
		okPane.setPadding(new Insets(gap, gap, gap, gap));
		okPane.setHgap(gap);
		okPane.add(buttonCancel, 0, 0);
		okPane.add(buttonOK, 1, 0);

		GridPane helpPane = new GridPane();
		helpPane.setPadding(new Insets(gap, gap, gap, gap));
		helpPane.add(buttonHelp, 0, 0);

		BorderPane buttonPane = new BorderPane();
		buttonPane.setLeft(helpPane);
		buttonPane.setRight(okPane);

		BorderPane overallPane = new BorderPane();
		overallPane.setCenter(editorPane);
		overallPane.setBottom(buttonPane);

//		String css = getClass().getResource("/resources/fxeditor.css").toExternalForm();
		Scene scene = new Scene(overallPane, HiDPIHelper.scale(mode == GenericEditorArea.MODE_REACTION ? 800 : 480), HiDPIHelper.scale(400), true, SceneAntialiasing.BALANCED);
//		scene.getStylesheets().add(css);
//		editorPane.getScene3D().widthProperty().bind(scene.widthProperty());
//		editorPane.getScene3D().heightProperty().bind(scene.heightProperty());
		setScene(scene);

		mListener = new ArrayList<>();
		mIsCancelled = true;
	}

	public boolean isCancelled() {
		return mIsCancelled;
	}

	public StereoMolecule getStructure() {
		return mMolecule;
	}

	/**
	 * @return mapped reaction with absolute coordinates, but without drawing objects
	 */
	public Reaction getReaction() {
		return mArea.getGenericDrawArea().getReaction();
	}

	/**
	 * @return mapped reaction with absolute coordinates and drawing objects
	 */
	public Reaction getReactionAndDrawings() {
		return mArea.getGenericDrawArea().getReactionAndDrawings();
	}

	public GenericEditorArea getDrawArea() {
		return mArea.getGenericDrawArea();
	}

	public void addStructureListener(StructureListener listener) {
		mListener.add(listener);
	}
}
