package com.actelion.research.gui.generic;

import java.io.File;

public interface GenericUIHelper {
	void showMessage(String message);
	void showHelpDialog(String url, String title);
	File openChemistryFile(boolean isReaction);
	GenericDialog createDialog(String title, GenericEventListener<GenericActionEvent> consumer);
	GenericPopupMenu createPopupMenu(GenericEventListener<GenericActionEvent> consumer);
	GenericImage createImage(String name);
	GenericImage createImage(int width, int height);
	void grabFocus();
	void setCursor(int cursor);
	void runLater(Runnable r);
}
