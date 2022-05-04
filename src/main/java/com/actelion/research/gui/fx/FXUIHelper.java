package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.*;
import javafx.scene.Node;
import javafx.scene.control.Alert;
import javafx.stage.Window;

import java.io.File;

public class FXUIHelper implements GenericUIHelper {
	private Node mParentNode;

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
//		mDrawArea.setCursor(CursorHelper.getCursor(cursor));
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
}
