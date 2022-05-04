package com.actelion.research.gui.editor;

import com.actelion.research.gui.fx.FXDrawContext;
import com.actelion.research.gui.fx.FXMouseHandler;
import com.actelion.research.gui.generic.GenericCanvas;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericMouseEvent;
import javafx.application.Platform;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;


public class FXEditorToolbar extends Canvas implements GenericCanvas {
	private GenericEditorToolbar mGenericToolbar;
	private volatile boolean mDrawPending;

	public FXEditorToolbar(FXEditorArea fxEditorArea, int mode) {
		mGenericToolbar = new GenericEditorToolbar(this, fxEditorArea.getGenericDrawArea(), mode);

		FXMouseHandler mouseHandler = new FXMouseHandler(mGenericToolbar);
		mouseHandler.addListener(mGenericToolbar);

		setOnMousePressed(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_PRESSED) );
		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_RELEASED) );
		setOnMouseDragged(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_DRAGGED) );

		setWidth(mGenericToolbar.getWidth());
		setHeight(mGenericToolbar.getHeight());

		GraphicsContext gc = getGraphicsContext2D();
		mGenericToolbar.paintContent(new FXDrawContext(gc));
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
		if (!mDrawPending) {
			mDrawPending = true;
			Platform.runLater(() -> draw() );
		}
	}

	private void draw() {
		mDrawPending = false;
		mGenericToolbar.paintContent(new FXDrawContext(getGraphicsContext2D()));
	}

	@Override
	public GenericDrawContext getDrawContext() {
		return new FXDrawContext(getGraphicsContext2D());
	}
}
