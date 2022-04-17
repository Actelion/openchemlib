package com.actelion.research.gui.editor;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.clipboard.ClipboardHandler;
import com.actelion.research.gui.dnd.MoleculeDropAdapter;
import com.actelion.research.gui.fx.FXDialogHelper;
import com.actelion.research.gui.fx.FXDrawContext;
import com.actelion.research.gui.fx.FXKeyHandler;
import com.actelion.research.gui.fx.FXMouseHandler;
import com.actelion.research.gui.generic.GenericCanvas;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericKeyEvent;
import com.actelion.research.gui.generic.GenericMouseEvent;
import javafx.scene.canvas.Canvas;

import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.dnd.DnDConstants;


public class FXEditorArea extends Canvas implements GenericCanvas {
	private static final int ALLOWED_DROP_ACTIONS = DnDConstants.ACTION_COPY_OR_MOVE;

	private GenericDrawArea mDrawArea;
	private FXKeyHandler mKeyHandler;

	public FXEditorArea(StereoMolecule mol, int mode) {
//		setFocusable(true);

		mDrawArea = new GenericDrawArea(mol, mode, new FXDialogHelper(this), this);

		initializeDragAndDrop(ALLOWED_DROP_ACTIONS);

		FXMouseHandler mouseHandler = new FXMouseHandler();
		mouseHandler.addListener(mDrawArea);

		setOnMousePressed(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_PRESSED) );
		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_RELEASED) );
//		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_CLICKED) ); not used
		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_ENTERED) );
//		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_EXITED) ); not used
		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_MOVED) );
		setOnMouseReleased(me -> mouseHandler.fireEvent(me, GenericMouseEvent.MOUSE_DRAGGED) );

		mKeyHandler = new FXKeyHandler();
		mKeyHandler.addListener(mDrawArea);

		setOnKeyPressed(ke -> mKeyHandler.fireEvent(ke, GenericKeyEvent.KEY_PRESSED) );
		setOnKeyReleased(ke -> mKeyHandler.fireEvent(ke, GenericKeyEvent.KEY_RELEASED) );
		setOnKeyTyped(ke -> mKeyHandler.fireEvent(ke, GenericKeyEvent.KEY_TYPED) );

		getGenericDrawArea().setClipboardHandler(new ClipboardHandler());
	}

	public FXKeyHandler getKeyHandler() {
		return mKeyHandler;
	}

	public GenericDrawArea getGenericDrawArea() {
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
	public GenericDrawContext getDrawContext() {
		return new FXDrawContext(getGraphicsContext2D());
	}

	@Override
	public void repaint() {}

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
					mDrawArea.addPastedOrDropped(mol, p);
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
