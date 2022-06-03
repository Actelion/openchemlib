package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.*;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.gui.swing.SwingCursorHelper;
import javafx.application.Platform;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.web.WebView;
import javafx.stage.Stage;
import javafx.stage.Window;

import java.io.File;
import java.io.IOException;
import java.net.URL;

public class FXUIHelper implements GenericUIHelper {
	private Node mParentNode;
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
		return null;
/*		return isReaction ?
				FileHelper.getFile(mDrawArea, "Please select a reaction file",
						FileHelper.cFileTypeRXN | CompoundFileHelper.cFileTypeRD)
				: FileHelper.getFile(mDrawArea, "Please select a molecule file",
				FileHelper.cFileTypeMOL | CompoundFileHelper.cFileTypeMOL2);
*/	}

	@Override
	public void showHelpDialog(String url, String title) {
		if (mHelpDialog == null) {
			mHelpDialog = new Stage();
			WebView view = new WebView();
			view.setZoom(HiDPIHelper.getUIScaleFactor());
			view.getEngine().load(createURL(url).toExternalForm());
			Scene scene = new Scene(view, HiDPIHelper.scale(640), HiDPIHelper.scale(480));
			mHelpDialog.setScene(scene);
			mHelpDialog.show();
		}


/*		if (mHelpDialog == null || !mHelpDialog.isVisible()) {
			JEditorPane helpPane = new JEditorPane();
			helpPane.setEditorKit(HiDPIHelper.getUIScaleFactor() == 1f ? new HTMLEditorKit() : new ScaledEditorKit());
			helpPane.setEditable(false);
			try {
				helpPane.setPage(getClass().getResource(url));
			}
			catch (Exception ex) {
				helpPane.setText(ex.toString());
			}

			Component c = getParent();
			if (c instanceof Frame)
				mHelpDialog = new JDialog((Frame)c, title, false);
			else
				mHelpDialog = new JDialog((Dialog)c, title, false);

			mHelpDialog.setSize(HiDPIHelper.scale(520), HiDPIHelper.scale(440));
			mHelpDialog.getContentPane().add(new JScrollPane(helpPane,
					JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
					JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));
			int x = (c.getX()>=8 + mHelpDialog.getWidth()) ? c.getX() - 8 - mHelpDialog.getWidth() : c.getX() + 8 + c.getWidth();
			mHelpDialog.setLocation(x, c.getY());
			mHelpDialog.setVisible(true);
		}
		else {
			Component c = getParent();
			int x = (mHelpDialog.getX() + mHelpDialog.getWidth() / 2>=c.getX() + c.getWidth() / 2) ?
					c.getX() - 8 - mHelpDialog.getWidth() : c.getX() + 8 + c.getWidth();
			mHelpDialog.setLocation(x, c.getY());
		}*/
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
