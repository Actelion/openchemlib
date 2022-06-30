package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericEventHandler;
import com.actelion.research.gui.generic.GenericKeyEvent;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyEvent;


public class FXKeyHandler extends GenericEventHandler<GenericKeyEvent> {
	public FXKeyHandler(Object source) {
		super(source);
	}

	public void fireEvent(KeyEvent ke, int what) {
		int key = ke.getCode() == KeyCode.ALT ? GenericKeyEvent.KEY_ALT
				: ke.getCode() == KeyCode.CONTROL ? GenericKeyEvent.KEY_CTRL
				: ke.getCode() == KeyCode.SHIFT ? GenericKeyEvent.KEY_SHIFT
				: ke.getCode() == KeyCode.DELETE ? GenericKeyEvent.KEY_DELETE
				: ke.getCode() == KeyCode.BACK_SPACE ? GenericKeyEvent.KEY_BACK_SPACE
				: ke.getCode() == KeyCode.HELP ? GenericKeyEvent.KEY_HELP
				: ke.getCode() == KeyCode.ESCAPE ? GenericKeyEvent.KEY_ESCAPE
				: ke.getText().length() != 1 ? 0 : ke.getText().charAt(0);
		if (key >= 1 && key <= 26)  // strangely, if Ctrl is pressed, letters are encoded from 1-26
			key = 'a' + key - 1;
		GenericKeyEvent gke = new GenericKeyEvent(what, key, ke.isAltDown(), ke.isControlDown(), ke.isShiftDown(), ke.isShortcutDown(), ke.getSource());
		fireEvent(gke);
	}
}
