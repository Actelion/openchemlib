package com.actelion.research.gui.editor;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;
import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.SceneAntialiasing;
import javafx.stage.Stage;

public class FXEditorTestApp extends Application {
	private static final boolean IS_REACTION_MODE = true;
	private static final boolean IS_TEST_DIALOG = true;

	public static void main(String[] args) {
		launch(args);
	}

	@Override
	public void start(Stage primaryStage) {
//		Molecule.setDefaultAverageBondLength(12);

		FXEditorPane editorPane;
		if (IS_REACTION_MODE) {
			Reaction rxn = ReactionEncoder.decode("gOYDGaDDHRTve`H!gKXHL@aJWFe`H#qB`ip qiV`#!B_vq?Dw}lL{y?[G|S !BTqa`FbpX?`@##" , true);

			if (IS_TEST_DIALOG) {
				new FXEditorDialog(null, rxn).show();
				return;
			}
			else {
				editorPane = new FXEditorPane(rxn);
			}
		}
		else {
			StereoMolecule mol = new SmilesParser().parseMolecule("Nc1cc(OCCO)cc(N)c1");
			mol.setFragment(true);

			if (IS_TEST_DIALOG) {
				new FXEditorDialog(null, mol).show();
				return;
			}
			else {
				editorPane = new FXEditorPane(mol);
			}
		}

		Scene scene = new Scene(editorPane, 800, 600, true, SceneAntialiasing.BALANCED);
		primaryStage.setTitle("JavaFX Molecule Editor");
		primaryStage.setScene(scene);
		primaryStage.show();
	}
}
