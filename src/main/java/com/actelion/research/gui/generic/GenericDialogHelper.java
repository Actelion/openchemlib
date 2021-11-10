package com.actelion.research.gui.generic;

import com.actelion.research.gui.editor.DialogEventConsumer;

import java.io.File;

public interface GenericDialogHelper {
	void showMessage(String message);
	void showHelpDialog(String url);
	File openChemistryFile(boolean isReaction);
	GenericDialog createDialog(String title, DialogEventConsumer consumer);
	GenericPopupMenu createPopupMenu(DialogEventConsumer consumer);
	void grabFocus();
	void setCursor(int cursor);
}
