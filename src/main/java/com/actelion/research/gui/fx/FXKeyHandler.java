package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericKeyEvent;
import com.actelion.research.gui.generic.GenericKeyListener;
import javafx.scene.input.KeyEvent;

import java.util.ArrayList;


public class FXKeyHandler {
	private ArrayList<GenericKeyListener> mListeners;

	public FXKeyHandler() {
		mListeners = new ArrayList<>();
	}

	public void addListener(GenericKeyListener gkl) {
		mListeners.add(gkl);
	}

	public void removeListener(GenericKeyListener gkl) {
		mListeners.remove(gkl);
	}

	public void fireEvent(KeyEvent ke, int what) {
		// TODO
		GenericKeyEvent gke = new GenericKeyEvent(what, ke.hashCode(), ke.getCharacter().charAt(0), ke.isShortcutDown());

		for (GenericKeyListener gkl:mListeners)
			gkl.keyActionHappened(gke);
	}
}
