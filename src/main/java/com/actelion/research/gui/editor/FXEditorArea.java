package com.actelion.research.gui.editor;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.clipboard.ClipboardHandler;
import com.actelion.research.gui.dnd.MoleculeDropAdapter;
import com.actelion.research.gui.fx.FXDrawContext;
import com.actelion.research.gui.fx.FXKeyHandler;
import com.actelion.research.gui.fx.FXMouseHandler;
import com.actelion.research.gui.fx.FXUIHelper;
import com.actelion.research.gui.generic.*;
import javafx.application.Platform;
import javafx.scene.Node;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundFill;
import javafx.scene.layout.Pane;
import javafx.scene.paint.Color;
import javafx.scene.paint.Paint;

import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.dnd.DnDConstants;


public class FXEditorArea extends Canvas implements GenericCanvas {
	private static final int ALLOWED_DROP_ACTIONS = DnDConstants.ACTION_COPY_OR_MOVE;

	private GenericEditorArea mDrawArea;
	private FXKeyHandler mKeyHandler;
	private volatile boolean mDrawPending;

	public FXEditorArea(StereoMolecule mol, int mode) {
//		setFocusable(true);

		mDrawArea = new GenericEditorArea(mol, mode, new FXUIHelper(this), this);

		widthProperty().addListener(evt -> repaint());
		heightProperty().addListener(evt -> repaint());

		initializeDragAndDrop(ALLOWED_DROP_ACTIONS);

		FXMouseHandler mouseHandler = new FXMouseHandler(mDrawArea);
		mouseHandler.addListener(mDrawArea);

		setOnMousePressed(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_PRESSED) );
		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_RELEASED) );
//		setOnMouseClicked(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_CLICKED) ); not used
		setOnMouseEntered(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_ENTERED) );
//		setOnMouseExited(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_EXITED) ); not used
		setOnMouseMoved(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_MOVED) );
		setOnMouseDragged(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_DRAGGED) );

		mKeyHandler = new FXKeyHandler(mDrawArea);
		mKeyHandler.addListener(mDrawArea);

		setOnKeyPressed(ke -> mKeyHandler.fireEvent(ke, GenericKeyEvent.KEY_PRESSED) );
		setOnKeyReleased(ke -> mKeyHandler.fireEvent(ke, GenericKeyEvent.KEY_RELEASED) );
//		setOnKeyTyped(ke -> mKeyHandler.fireEvent(ke, GenericKeyEvent.KEY_TYPED) );

		getGenericDrawArea().setClipboardHandler(new ClipboardHandler());
	}

	private void draw() {
		mDrawPending = false;

		double width = getWidth();
		double height = getHeight();

		GraphicsContext gc = getGraphicsContext2D();
		gc.clearRect(0, 0, width, height);

		gc.setStroke(Color.RED);
		gc.strokeLine(0, 0, width, height);
		gc.strokeLine(0, height, width, 0);

		mDrawArea.paintContent(new FXDrawContext(gc));
		}

	@Override
	public boolean isResizable() {
		return true;
		}

	@Override
	public double prefWidth(double height) {
		return getWidth();
		}

	@Override
	public double prefHeight(double width) {
		return getHeight();
		}

	public FXKeyHandler getKeyHandler() {
		return mKeyHandler;
	}

	public GenericEditorArea getGenericDrawArea() {
		return mDrawArea;
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
	public GenericDrawContext getDrawContext() {
		return new FXDrawContext(getGraphicsContext2D());
	}

	@Override
	public void repaint() {
		if (!mDrawPending) {
			mDrawPending = true;
			Platform.runLater(() -> draw() );
		}
	}

/*	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		((Graphics2D)g).setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);

		mDrawArea.paintContent(new FXDrawContext(getGraphicsContext2D()));
	}*/

	private void initializeDragAndDrop(int dropAction) {
		if (dropAction != DnDConstants.ACTION_NONE) {
			MoleculeDropAdapter d = new MoleculeDropAdapter() {
				public void onDropMolecule(StereoMolecule mol, Point p) {
					mDrawArea.addPastedOrDropped(mol, new GenericPoint(p.x, p.y));
				}
			};

			// TODO
//			new DropTarget(this, dropAction, d, true, new FXEditorArea.OurFlavorMap());
		}
	}

// This class is needed for inter-jvm drag&drop. Although not neccessary for standard environments, it prevents
// nasty "no native data was transfered" errors. It still might create ClassNotFoundException in the first place by
// the SystemFlavorMap, but as I found it does not hurt, since the context classloader will be installed after
// the first call. I know, that this depends heavely on a specific behaviour of the systemflavormap, but for now
// there's nothing I can do about it.
	static class OurFlavorMap implements java.awt.datatransfer.FlavorMap {
		@Override
		public java.util.Map<DataFlavor, String> getNativesForFlavors(DataFlavor[] dfs) {
			java.awt.datatransfer.FlavorMap m = java.awt.datatransfer.SystemFlavorMap.getDefaultFlavorMap();
			return m.getNativesForFlavors(dfs);
		}

		@Override
		public java.util.Map<String, DataFlavor> getFlavorsForNatives(String[] natives) {
			java.awt.datatransfer.FlavorMap m = java.awt.datatransfer.SystemFlavorMap.getDefaultFlavorMap();
			return m.getFlavorsForNatives(natives);
		}
	}
}
