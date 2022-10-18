package com.actelion.research.gui.editor;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.MolfileParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.generic.GenericPoint;
import javafx.scene.input.DataFormat;
import javafx.scene.input.Dragboard;
import javafx.scene.input.TransferMode;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.Pane;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

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
		mToolBar = new FXEditorToolbar(mArea);
		widthProperty().addListener((observable, oldValue, newValue) -> mArea.setWidth((double) newValue-mToolBar.getWidth()));
		heightProperty().addListener((observable, oldValue, newValue) -> mArea.setHeight((double) newValue));
		setLeft(mToolBar);
		setCenter(new Pane(mArea));

		initializeDragAndDrop();
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

	private void initializeDragAndDrop() {
		setOnDragEntered(event -> {
				/* the drag-and-drop gesture entered the target */
				System.out.println("onDragEntered");
				/* show to the user that it is an actual gesture target */
				if (event.getGestureSource() != this &&
						event.getDragboard().hasString()) {
//                            target.setFill(Color.GREEN);
				}
				event.consume();
		});

		setOnDragOver(event -> {
			/* data is dragged over the target */
			// System.out.println("onDragOver");

			/* accept it only if it is  not dragged from the same node
			 * and if it has a string data */
			Dragboard db = event.getDragboard();
			if (event.getGestureSource() != this && /*db.hasString() */ getAcceptedFormats(db).size() > 0) {
				/* allow for both copying and moving, whatever user chooses */
				event.acceptTransferModes(TransferMode.COPY_OR_MOVE);
			}

			event.consume();
			});

		setOnDragExited(event -> event.consume() );

		setOnDragDropped(event -> {
			/* data dropped */
			System.out.println("onDragDropped");
			/* if there is a string data on dragboard, read it and use it */
			Dragboard db = event.getDragboard();
//                Set<DataFormat> formats = db.getContentTypes();
			boolean success = false;
//                for (DataFormat f : formats) {
//                    System.out.printf("D&D Format %s\n", f);
//                }

			List<DataFormat> formats = getAcceptedFormats(db);

			StereoMolecule m = new StereoMolecule();
			for (DataFormat format : formats) {
				if (isJavaFormat(format, MOLFILE_FORMAT)) {
					MolfileParser p = new MolfileParser();
					p.parse(m, db.getContent(format).toString());
					success = true;
					break;
				}
			}
			if (!success && db.hasString()) {
				try {
					IDCodeParser p = new IDCodeParser(true);
					p.parse(m, db.getString());
					success = true;
				} catch (Exception e) {
					System.err.println("Cannot parse data for molecules ");
				}
			}
			if (success) {
				mArea.getGenericDrawArea().addPastedOrDropped(m, new GenericPoint(event.getX(), event.getY()));
			}

			/* let the source know whether the string was successfully
			 * transferred and used */
			event.setDropCompleted(success);
			event.consume();
		});

//        setOnDragDetected(event -> {
//                        /* drag was detected, start drag-and-drop gesture*/
//                        System.out.println("onDragDetected");
//
//                        /* allow any transfer mode */
//                        Dragboard db = source.startDragAndDrop(TransferMode.ANY);
//
//                        /* put a string on dragboard */
//                        ClipboardContent content = new ClipboardContent();
//                        content.putString(source.getText());
//                        db.setContent(content);
//
//                        event.consume();
//                });
	}

	private static String MOLFILE_FORMAT = "JAVA_DATAFLAVOR:chemical/x-mdl-molfilev3; class=java.lang.String";

	private List<DataFormat> getAcceptedFormats(Dragboard db) {
		Set<DataFormat> formats = db.getContentTypes();
		List<DataFormat> res = new ArrayList<>();
		for (DataFormat f : formats) {
			if (f.equals(DataFormat.PLAIN_TEXT))
				res.add(f);
			else if (isJavaFormat(f, MOLFILE_FORMAT)) {
				res.add(f);
			}
		}
		return res;
	}

	private boolean isJavaFormat(DataFormat f, String identifier) {
		Set<String> s = f.getIdentifiers();
		for (String a : s) {
			if (a.equals(identifier)) {
				return true;
			}
		}
		return false;
	}
}
