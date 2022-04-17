package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericMouseEvent;
import com.actelion.research.gui.generic.GenericMouseListener;
import javafx.scene.input.MouseButton;
import javafx.scene.input.MouseEvent;

import java.util.ArrayList;

public class FXMouseHandler {
	private ArrayList<GenericMouseListener> mListeners;

	public FXMouseHandler() {
		mListeners = new ArrayList<>();
	}

	public void addListener(GenericMouseListener gml) {
		mListeners.add(gml);
	}

	public void removeListener(GenericMouseListener gml) {
		mListeners.remove(gml);
	}

	public void fireEvent(MouseEvent me, int what) {
		int button = (me.getButton() == MouseButton.PRIMARY) ? 1
				   : (me.getButton() == MouseButton.MIDDLE) ? 2 : 3;

		GenericMouseEvent gme = new GenericMouseEvent(what, button, me.getClickCount(), (int)me.getX(), (int)me.getY(),
				me.isShiftDown(), me.isControlDown(), me.isAltDown(), me.isPopupTrigger());

		for (GenericMouseListener gml:mListeners)
			gml.mouseActionHappened(gme);
	}
}
