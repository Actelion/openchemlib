package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericEventHandler;
import com.actelion.research.gui.generic.GenericMouseEvent;
import javafx.scene.input.MouseButton;
import javafx.scene.input.MouseEvent;

public class FXMouseHandler extends GenericEventHandler<GenericMouseEvent> {

	public FXMouseHandler(Object source) {
		super(source);
	}

	public void fireEvent(MouseEvent me, int what) {
		int button = (me.getButton() == MouseButton.PRIMARY) ? 1
				   : (me.getButton() == MouseButton.MIDDLE) ? 2 : 3;

		GenericMouseEvent gme = new GenericMouseEvent(what, button, me.getClickCount(), (int)me.getX(), (int)me.getY(),
				me.isShiftDown(), me.isControlDown(), me.isAltDown(), me.isPopupTrigger(), getSource());

		super.fireEvent(gme);
	}
}
