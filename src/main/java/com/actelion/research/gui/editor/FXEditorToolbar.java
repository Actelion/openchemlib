package com.actelion.research.gui.editor;

import com.actelion.research.gui.fx.FXDrawContext;
import com.actelion.research.gui.fx.FXMouseHandler;
import com.actelion.research.gui.generic.GenericCanvas;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericMouseEvent;
import javafx.scene.canvas.Canvas;


public class FXEditorToolbar extends Canvas implements GenericCanvas {
	private GenericEditorToolbar mGenericToolbar;

	public FXEditorToolbar(FXEditorArea fxEditorArea, int mode) {
		mGenericToolbar = new GenericEditorToolbar(this, fxEditorArea.getGenericDrawArea(), mode);

		FXMouseHandler mouseHandler = new FXMouseHandler();
		mouseHandler.addListener(mGenericToolbar);

		setOnMousePressed(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_PRESSED) );
		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_RELEASED) );
		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_DRAGGED) );

		setWidth(mGenericToolbar.getWidth());
		setHeight(mGenericToolbar.getHeight());
	}

	@Override
	public double getCanvasWidth() {
		return getWidth();
	}

	@Override
	public double getCanvasHeight() {
		return getHeight();
	}

	@Override
	public void repaint() {
		// probably not needed!!!
	}

	@Override
	public GenericDrawContext getDrawContext() {
		return new FXDrawContext(getGraphicsContext2D());
	}
}
