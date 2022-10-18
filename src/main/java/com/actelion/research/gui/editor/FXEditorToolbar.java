package com.actelion.research.gui.editor;

import com.actelion.research.gui.fx.FXDrawContext;
import com.actelion.research.gui.fx.FXMouseHandler;
import com.actelion.research.gui.generic.GenericCanvas;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericMouseEvent;
import javafx.application.Platform;
import javafx.scene.Node;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundFill;
import javafx.scene.layout.Pane;
import javafx.scene.paint.Color;
import javafx.scene.paint.Paint;


public class FXEditorToolbar extends Canvas implements GenericCanvas {
	private GenericEditorToolbar mGenericToolbar;
	private volatile boolean mDrawPending;

	public FXEditorToolbar(FXEditorArea fxEditorArea) {
		mGenericToolbar = new GenericEditorToolbar(this, fxEditorArea.getGenericDrawArea());

		FXMouseHandler mouseHandler = new FXMouseHandler(mGenericToolbar);
		mouseHandler.addListener(mGenericToolbar);

		setOnMousePressed(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_PRESSED) );
		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_RELEASED) );
		setOnMouseMoved(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_MOVED) );
		setOnMouseDragged(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_DRAGGED) );
		setOnMouseExited(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_EXITED) );

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
	public int getBackgroundRGB() {
		Node parent = getParent();
		while (parent != null && !(parent instanceof Pane))
			parent = parent.getParent();
		if (parent != null && parent instanceof Pane) {
			Background bg = ((Pane)parent).getBackground();
			if (bg != null) {
				BackgroundFill bgf = bg.getFills().get(0);
				if (bgf != null) {
					Paint paint = bgf.getFill();
					if (paint instanceof Color) {
						Color c = (Color)paint;
						return (Math.round(255f*(float)c.getRed()) << 16) + (Math.round(255f*(float)c.getGreen()) << 8) + Math.round(255f*(float)c.getBlue());
					}
				}
			}
		}
		return 0x00E0E0E0;  // we just assume light gray
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
		GraphicsContext gc = getGraphicsContext2D();
		gc.clearRect(0, 0, getWidth(), getHeight());
		mGenericToolbar.paintContent(new FXDrawContext(gc));
	}

	@Override
	public GenericDrawContext getDrawContext() {
		return new FXDrawContext(getGraphicsContext2D());
	}
}
