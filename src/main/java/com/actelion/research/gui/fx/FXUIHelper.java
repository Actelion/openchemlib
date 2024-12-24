package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.*;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import javafx.application.Platform;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.web.WebView;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javafx.stage.Window;

import java.io.File;
import java.io.IOException;
import java.net.URL;

public class FXUIHelper implements GenericUIHelper {
	private final Node mParentNode;
	private Stage mHelpDialog;

	public FXUIHelper(Node parent) {
		mParentNode = parent;
	}

	@Override
	public void showMessage(String message) {
		new Alert(Alert.AlertType.INFORMATION, message).showAndWait();
	}

	@Override
	public GenericDialog createDialog(String title, GenericEventListener<GenericActionEvent> consumer) {
		Window window = mParentNode.getScene().getWindow();
		GenericDialog dialog = new FXDialog(window, title);
		dialog.setEventConsumer(consumer);
		return dialog;
	}

	@Override
	public GenericPopupMenu createPopupMenu(GenericEventListener<GenericActionEvent> consumer) {
		return new FXPopupMenu(mParentNode, consumer);
	}

	@Override
	public void grabFocus() {
		mParentNode.requestFocus();
	}

	@Override
	public void setCursor(int cursor) {
		mParentNode.setCursor(FXCursorHelper.getCursor(cursor));
	}

	@Override
	public GenericImage createImage(String name) {
		return new FXImage(name);
	}

	@Override
	public GenericImage createImage(int width, int height) {
		return new FXImage(width, height);
	}

	@Override
	public void runLater(Runnable r) {
		Platform.runLater(r);
	}

	@Override
	public File openChemistryFile(boolean isReaction) {
		FileChooser fileChooser = new FileChooser();
		if (isReaction) {
			fileChooser.setTitle("Please select a reaction file");
			fileChooser.getExtensionFilters().addAll(
					new FileChooser.ExtensionFilter("Reaction Files", "*.rxn"));
		}
		else {
			fileChooser.setTitle("Please select a molecule file");
			fileChooser.getExtensionFilters().addAll(
					new FileChooser.ExtensionFilter("MDL Molfiles", "*.mol"),
					new FileChooser.ExtensionFilter("Mol2-Files", "*.mol2"));
		}
		return fileChooser.showOpenDialog(mParentNode.getScene().getWindow());
	}

	@Override
	public void showHelpDialog(String url, String title) {
		if (mHelpDialog == null) {
			WebView view = new WebView();
			view.setZoom(HiDPIHelper.getUIScaleFactor());
			view.getEngine().load(createURL(url).toExternalForm());
			Scene scene = new Scene(view);

			mHelpDialog = new Stage();
			mHelpDialog.setScene(scene);
			mHelpDialog.setMinWidth(scene.getRoot().prefWidth(640));
			mHelpDialog.setMinHeight(scene.getRoot().prefHeight(480));
			mHelpDialog.show();
		}
		else {
			mHelpDialog.toFront();
		}
	}

	public static URL createURL(String urlText) {
		String ref = null;
		int index = urlText.indexOf('#');
		if (index != -1) {
			ref = urlText.substring(index);
			urlText = urlText.substring(0, index);
		}
		URL theURL = FXUIHelper.class.getResource(urlText);
		if (ref != null) {
			try {
				theURL = new URL(theURL, ref);
			}
			catch (IOException e) {
				return null;
			}
		}
		return theURL;
	}
}
