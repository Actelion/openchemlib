/*
 * Copyright (c) 1997 - 2016
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the the copyright holder nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

package com.actelion.research.gui;

import com.actelion.research.chem.*;
import com.actelion.research.chem.name.StructureNameResolver;
import com.actelion.research.gui.clipboard.IClipboardHandler;
import com.actelion.research.gui.fx.FXDrawContext;
import com.actelion.research.gui.fx.FXMoleculeTransfer;
import com.actelion.research.gui.generic.GenericDepictor;
import com.actelion.research.gui.generic.GenericPolygon;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.gui.generic.GenericShape;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.util.ColorHelper;
import javafx.animation.PauseTransition;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.scene.Cursor;
import javafx.scene.SnapshotParameters;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.ContextMenu;
import javafx.scene.control.MenuItem;
import javafx.scene.control.SeparatorMenuItem;
import javafx.scene.input.*;
import javafx.scene.layout.Pane;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.stage.FileChooser;
import javafx.stage.Window;
import javafx.util.Duration;

import java.awt.dnd.DnDConstants;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

/**
 * JavaFX equivalent of JStructureView. Molecules can be exchanged with JStructureView and other
 * Swing components via drag & drop and the clipboard.
 * Except for the constructors and setClipboardHandler(), all methods should be called on the FX application thread.
 * Like JStructureView, sizes are scaled with HiDPIHelper.getUIScaleFactor(), which is 1.0 on Windows and macOS,
 * where JavaFX scales automatically, but reflects Xft.dpi on Linux, which JavaFX doesn't consider for scaling.
 */
public class FXStructureView extends Pane implements StructureListener {
	private static final int DEFAULT_SELECTION_BACKGROUND = 0x80A4C0;
	private static final int DEFAULT_EDITABLE_BACKGROUND = 0xFFFFFF;
	private static final int DEFAULT_INACTIVE_BACKGROUND = 0xF2F2F2;
	private static final int DEFAULT_FOREGROUND = 0x000000;

	private static final String ITEM_COPY = "Copy Structure";
	private static final String ITEM_COPY_SMILES = "Copy Structure As SMILES-String";
	private static final String ITEM_PASTE= "Paste Structure";
	private static final String ITEM_PASTE_WITH_NAME = ITEM_PASTE+" or Name";
	private static final String ITEM_CLEAR = "Clear Structure";
	private static final String ITEM_OPEN_MOLFILE = "Open Molfile...";
	private static final String ITEM_SAVE_MOLFILE_V2 = "Save As Molfile V2...";
	private static final String ITEM_SAVE_MOLFILE_V3 = "Save As Molfile V3...";

	private static final long WARNING_MILLIS = 1200;
	private static final int WARNING_FONT_SIZE = 12;

	private static final int DRAG_MARGIN = 12;
	private static final int BORDER_ARC = 10;

	private final Canvas mCanvas;
	private ArrayList<StructureListener> mListener;
	private String mIDCode;
	private StereoMolecule mMol,mDisplayMol;
	private GenericDepictor mDepictor;
	private boolean mShowBorder,mAllowFragmentStatusChangeOnPasteOrDrop,mIsDraggingThis,mIsEditable,mDisableBorder,
			mIsSelectable,mIsSelecting,mIsLassoSelect;
	private volatile boolean mDrawPending;
	private int mChiralTextPosition,mDisplayMode;
	private String[] mAtomText;
	private IClipboardHandler mClipboardHandler;
	protected int mAllowedDragAction;
	protected int mAllowedDropAction;
	private String mWarningMessage;
	private PauseTransition mWarningTimer;
	private int[] mAtomHiliteColor;
	private float[] mAtomHiliteRadius;
	private double mTextSizeFactor,mX1,mY1,mX2,mY2;
	private GenericPolygon mLassoRegion;
	private Color mBackground,mForeground;
	private ContextMenu mPopup;

	public FXStructureView() {
		this(null);
		}

	/**
	 * This creates a standard structure view where the displayed molecule is
	 * used for D&D and clipboard transfer after removing atom colors and bond highlights.
	 * The default will support copy/paste and drag&drop from this view only,
	 * but dropping anything onto this view doesn't have an effect.
	 * Call setEditable(true) to allow changes through drag&drop and pasting.
	 * @param mol used for display, clipboard copy and d&d
	 */
	public FXStructureView(StereoMolecule mol) {
		this(mol, DnDConstants.ACTION_COPY_OR_MOVE, DnDConstants.ACTION_COPY_OR_MOVE);
		}

	/**
	 * This creates a structure view that distinguishes between displayed molecule
	 * and the one being used for D&D and clipboard transfer. Use this if the displayed
	 * molecule is structurally different, e.g. uses custom atom labels or additional
	 * illustrative atoms or bonds, which shall not be copied.
	 * Custom atom colors or highlighted bonds don't require a displayMol.
	 * The default will support copy/paste and drag&drop FROM this view only,
	 * but dropping anything onto this view doesn't have an effect.
	 * Call setEditable(true) to allow changes through drag&drop and pasting.
	 * @param mol used for clipboard copy and d&d; used for display if displayMol is null
	 * @param displayMol null if mol shall be displayed
	 */
	public FXStructureView(StereoMolecule mol, StereoMolecule displayMol) {
		this(mol, displayMol, DnDConstants.ACTION_COPY_OR_MOVE, DnDConstants.ACTION_COPY_OR_MOVE);
		}

	/**
	 * @param dragAction DnDConstants.ACTION_xxx
	 * @param dropAction DnDConstants.ACTION_xxx
	 */
	public FXStructureView(int dragAction, int dropAction) {
		this(null, dragAction, dropAction);
		}

	/**
	 * This creates a standard structure view where the displayed molecule is
	 * used for D&D and clipboard transfer after removing atom colors and bond highlights.
	 * The default will support copy/paste and drag&drop from this view only,
	 * but dropping anything onto this view doesn't have an effect.
	 * Call setEditable(true) to allow changes through drag&drop and pasting.
	 * @param mol used for display, clipboard copy and d&d
	 * @param dragAction DnDConstants.ACTION_xxx
	 * @param dropAction DnDConstants.ACTION_xxx
	 */
	public FXStructureView(StereoMolecule mol, int dragAction, int dropAction) {
		this(mol, null, dragAction, dropAction);
		}

	/**
	 * This creates a structure view that distinguishes between displayed molecule
	 * and the one being used for D&D and clipboard transfer. Use this if the displayed
	 * molecule is structurally different, e.g. uses custom atom labels or additional
	 * illustrative atoms or bonds, which shall not be copied.
	 * Custom atom colors or highlighted bonds don't require a displayMol.
	 * The default will support copy/paste and drag&drop from this view only,
	 * but dropping anything onto this view doesn't have an effect.
	 * Call setEditable(true) to allow changes through drag&drop and pasting.
	 * @param mol used for clipboard copy and d&d; used for display if displayMol is null
	 * @param displayMol null if mol shall be displayed
	 * @param dragAction DnDConstants.ACTION_xxx
	 * @param dropAction DnDConstants.ACTION_xxx
	 */
	public FXStructureView(StereoMolecule mol, StereoMolecule displayMol, int dragAction, int dropAction) {
		mMol = (mol == null) ? new StereoMolecule() : new StereoMolecule(mol);
		mDisplayMol = (displayMol == null) ? mMol : displayMol;
		mDisplayMode = AbstractDepictor.cDModeHiliteAllQueryFeatures;
		mTextSizeFactor = 1.0;
		mIsEditable = false;

		// The canvas is not managed to not let its size influence the preferred size of this pane.
		mCanvas = new Canvas();
		mCanvas.setManaged(false);
		mCanvas.widthProperty().addListener(e -> repaint());
		mCanvas.heightProperty().addListener(e -> repaint());
		getChildren().add(mCanvas);

		disabledProperty().addListener(e -> repaint());

		initializeMouseHandling();
		initializeDragAndDrop(dragAction, dropAction);
		}

	@Override
	protected void layoutChildren() {
		super.layoutChildren();
		Insets insets = getInsets();
		mCanvas.relocate(insets.getLeft(), insets.getTop());
		mCanvas.setWidth(Math.max(0, getWidth() - insets.getLeft() - insets.getRight()));
		mCanvas.setHeight(Math.max(0, getHeight() - insets.getTop() - insets.getBottom()));
		}

	public void setSelectable(boolean s) {
		mIsSelectable = s;
		}

	/**
	 * Call this in order to get clipboard support:
	 * setClipboardHandler(new ClipboardHandler());
	 */
	public void setClipboardHandler(IClipboardHandler h) {
		mClipboardHandler = h;
		}

	public IClipboardHandler getClipboardHandler() {
		return mClipboardHandler;
		}

	public int getDisplayMode() {
		return mDisplayMode;
		}

	/**
	 * Sets the display mode for the Depictor. The default is
	 * AbstractDepictor.cDModeHiliteAllQueryFeatures.
	 * @param mode
	 */
	public void setDisplayMode(int mode) {
		if (mDisplayMode != mode) {
			mDisplayMode = mode;
			repaint();
			}
		}

	/**
	 * Sets a multiplication factor to the text size of all labels. The default is 1.0.
	 * @param factor text size factor
	 */
	public void setFactorTextSize(double factor) {
		mTextSizeFactor = factor;
		repaint();
		}

	public void setDisableBorder(boolean b) {
		mDisableBorder = b;
		}

	/**
	 * Defines additional atom text to be displayed in top right
	 * position of some/all atom labels. If the atom is charged, then
	 * the atom text is drawn right of the atom charge.
	 * If using atom text make sure to update it accordingly, if atom
	 * indexes change due to molecule changes.
	 * @param atomText null or String array matching atom indexes (may contain null entries)
	 */
	public void setAtomText(String[] atomText) {
		mAtomText = atomText;
		repaint();
		}

	/**
	 * Defines the background color. If null (default), then the background color depends
	 * on whether the view is editable and enabled.
	 * @param c null or background color
	 */
	public void setBackgroundColor(Color c) {
		mBackground = c;
		repaint();
		}

	/**
	 * @return explicitly set background color or null
	 */
	public Color getBackgroundColor() {
		return mBackground;
		}

	/**
	 * Defines the foreground color, i.e. the color of bonds and carbon atoms.
	 * If null (default), then black or a contrast color to a dark background is used.
	 * @param c null or foreground color
	 */
	public void setForegroundColor(Color c) {
		mForeground = c;
		repaint();
		}

	/**
	 * @return explicitly set foreground color or null
	 */
	public Color getForegroundColor() {
		return mForeground;
		}

	public boolean isEditable() {
		return mIsEditable;
		}

	public void setEditable(boolean b) {
		if (mIsEditable != b) {
			mIsEditable = b;
			repaint();
			}
		}

	/**
	 * When fragment status change on drop is allowed then dropping a fragment (molecule)
	 * on a molecule (fragment) inverts the status of the view's chemical object.
	 * As default status changes are prohibited.
	 * @param allow
	 */
	public void setAllowFragmentStatusChangeOnPasteOrDrop(boolean allow) {
		mAllowFragmentStatusChangeOnPasteOrDrop = allow;
		}

	/**
	 * If you want this structure view to also draw an atom background with specific colors for every atom,
	 * then you need to call this method before or just after one of the structureChanged() calls.
	 * @param argb if alpha < 1 then the background is mixed in accordingly
	 * @param radius <= 1.0; if null, then a default of 0.5 of the average bond length is used
	 */
	public void setAtomHighlightColors(int[] argb, float[] radius) {
		mAtomHiliteColor = argb;
		mAtomHiliteRadius = radius;
		repaint();
		}

	public boolean canDrop() {
		return mIsEditable && !isDisabled() && !mIsDraggingThis;
		}

	/**
	 * Schedules a redraw of the view on the FX application thread. Multiple calls are combined.
	 */
	public void repaint() {
		if (!mDrawPending) {
			mDrawPending = true;
			Platform.runLater(() -> draw() );
			}
		}

	private synchronized void draw() {
		mDrawPending = false;

		double width = mCanvas.getWidth();
		double height = mCanvas.getHeight();

		GraphicsContext gc = mCanvas.getGraphicsContext2D();
		gc.clearRect(0, 0, width, height);

		if (width <= 0 || height <= 0)
			return;

		int bgRGB = getBackgroundRGB();
		gc.setFill(toFXColor(bgRGB));
		gc.fillRect(0, 0, width, height);

		if (mDisplayMol != null && mDisplayMol.getAllAtoms() != 0) {
			mDepictor = new GenericDepictor(mDisplayMol);
			mDepictor.setDisplayMode(mDisplayMode);
			mDepictor.setFactorTextSize(mTextSizeFactor);
			mDepictor.setAtomText(mAtomText);
			mDepictor.setAtomHighlightColors(mAtomHiliteColor, mAtomHiliteRadius);

			if (isDisabled())
				mDepictor.setOverruleColor(ColorHelper.getContrastColor(0x808080, bgRGB), bgRGB);
			else
				mDepictor.setForegroundColor(getForegroundRGB(), bgRGB);

			gc.save();
			FXDrawContext context = new FXDrawContext(gc);
			mDepictor.validateView(context, new GenericRectangle(0, 0, width, height),
								   AbstractDepictor.cModeInflateToMaxAVBL | mChiralTextPosition | HiDPIHelper.scale(AbstractDepictor.cOptAvBondLen));

			if (mShowBorder && !mDisableBorder) {
				GenericRectangle rect = mDepictor.getBoundingRect();
				if (rect != null) {
					gc.setFill(toFXColor(ColorHelper.perceivedBrightness(bgRGB) < 0.5f ?
							ColorHelper.brighter(bgRGB, 0.85f) : ColorHelper.darker(bgRGB, 0.85f)));
					double arc = Math.min(rect.height/4, Math.min(rect.width/4, HiDPIHelper.scale(BORDER_ARC)));
					gc.fillRoundRect(rect.x, rect.y, rect.width, rect.height, arc, arc);
					}
				}

			mDepictor.paint(context);
			gc.restore();
			}
		else {
			mDepictor = null;
			}

		if (mWarningMessage != null) {
			double fontSize = HiDPIHelper.scale(WARNING_FONT_SIZE);
			gc.setFont(Font.font(Font.getDefault().getFamily(), FontWeight.BOLD, fontSize));
			gc.setFill(Color.RED);
			double textWidth = getTextWidth(mWarningMessage, gc.getFont());
			gc.fillText(mWarningMessage, (width-textWidth)/2, 1.4*fontSize);
			}

		if (mIsSelecting) {
			gc.setStroke(toFXColor(lassoRGB(bgRGB)));
			gc.setLineWidth(1.0);
			if (mIsLassoSelect) {
				int count = mLassoRegion.getSize();
				double[] x = new double[count];
				double[] y = new double[count];
				for (int i=0; i<count; i++) {
					x[i] = mLassoRegion.getX(i);
					y[i] = mLassoRegion.getY(i);
					}
				gc.strokePolygon(x, y, count);
				}
			else {
				gc.strokeRect(Math.min(mX1, mX2), Math.min(mY1, mY2), Math.abs(mX2 - mX1), Math.abs(mY2 - mY1));
				}
			}

		drawOverlay(gc, width, height);
		}

	/**
	 * Override this to draw additional content on top of the structure.
	 * @param gc graphics context of the view's canvas
	 * @param width canvas width
	 * @param height canvas height
	 */
	protected void drawOverlay(GraphicsContext gc, double width, double height) {}

	protected static double getTextWidth(String text, Font font) {
		Text t = new Text(text);
		t.setFont(font);
		return t.getLayoutBounds().getWidth();
		}

	protected int getBackgroundRGB() {
		if (mBackground != null)
			return toRGB(mBackground);
		return (isEditable() && !isDisabled()) ? DEFAULT_EDITABLE_BACKGROUND : DEFAULT_INACTIVE_BACKGROUND;
		}

	protected int getForegroundRGB() {
		if (mForeground != null)
			return toRGB(mForeground);
		return ColorHelper.getContrastColor(DEFAULT_FOREGROUND, getBackgroundRGB());
		}

	private int lassoRGB(int bgRGB) {
		return ColorHelper.createColor(DEFAULT_SELECTION_BACKGROUND, ColorHelper.perceivedBrightness(bgRGB) < 0.5f ? 0.65f : 0.35f);
		}

	protected static Color toFXColor(int rgb) {
		return Color.rgb((rgb >> 16) & 0xFF, (rgb >> 8) & 0xFF, rgb & 0xFF);
		}

	private static int toRGB(Color c) {
		return ((int)Math.round(255*c.getRed()) << 16) + ((int)Math.round(255*c.getGreen()) << 8) + (int)Math.round(255*c.getBlue());
		}

	public void setIDCode(String idcode) {
		int index = (idcode == null) ? -1 : idcode.indexOf(' ');
		if (index == -1)
			setIDCode(idcode, null);
		else
			setIDCode(idcode.substring(0, index), idcode.substring(index+1));
		}

	public synchronized void setIDCode(String idcode, String coordinates) {
		if (idcode != null && idcode.isEmpty())
			idcode = null;

		if (mIDCode == null && idcode == null)
			return;

		if (mIDCode != null && mIDCode.equals(idcode))
			return;

		new IDCodeParser(true).parse(mMol, idcode, coordinates);
		mDisplayMol = mMol;

		mIDCode = idcode;
		repaint();
		informListeners();
		}

	/**
	 * Updates the molecule used for display, drag & drop and clipboard transfer.
	 * Also triggers a repaint().
	 * @param mol new molecule used for display, clipboard copy and d&d; may be null
	 */
	@Override
	public synchronized void structureChanged(StereoMolecule mol) {
		if (mol == null) {
			mMol.clear();
			}
		else {
			mol.copyMolecule(mMol);
			}

		mDisplayMol = mMol;
		structureChanged();
		}

	/**
	 * Updates both molecules used for display and for drag & drop/clipboard transfer.
	 * Also triggers a repaint().
	 * @param mol new molecule used for clipboard copy and d&d; may be null
	 * @param displayMol new molecule used for display, may be null
	 */
	public synchronized void structureChanged(StereoMolecule mol, StereoMolecule displayMol) {
		if (mol == null) {
			mMol.clear();
			}
		else {
			mol.copyMolecule(mMol);
			}

		mDisplayMol = displayMol;
		structureChanged();
		}

	/**
	 * Should only be called if FXStructureView's internal Molecule is changed
	 * from outside as: theStructureView.getMolecule().setFragment(false);
	 * The caller is responsible to update displayMol also, if it is different from
	 * the molecule.
	 */
	public synchronized void structureChanged() {
		mIDCode = null;
		repaint();
		informListeners();
		}

	public StereoMolecule getMolecule() {
		return mMol;
		}

	public StereoMolecule getDisplayMolecule() {
		return mDisplayMol;
		}

	public AbstractDepictor getDepictor() {
		return mDepictor;
		}

	public void addStructureListener(StructureListener l) {
		if (mListener == null)
			mListener = new ArrayList<>();

		mListener.add(l);
		}

	public void removeStructureListener(StructureListener l) {
		if (mListener != null)
			mListener.remove(l);
		}

	public void setChiralDrawPosition(int p) {
		mChiralTextPosition = p;
		repaint();
		}

	private void initializeMouseHandling() {
		mCanvas.addEventHandler(MouseEvent.MOUSE_MOVED, e -> {
			boolean isInRect = isInDragRect(e.getX(), e.getY());
			updateBorder(isInRect);
			setCursor(isInRect ? Cursor.OPEN_HAND : Cursor.DEFAULT);
			} );

		mCanvas.addEventHandler(MouseEvent.MOUSE_EXITED, e -> {
			if (!mIsSelecting)
				updateBorder(false);
			} );

		mCanvas.addEventHandler(MouseEvent.MOUSE_PRESSED, e -> {
			if (mIsSelectable
			 && e.getButton() == MouseButton.PRIMARY
			 && !e.isPopupTrigger()
			 && !isInDragRect(e.getX(), e.getY())) {
				mX1 = e.getX();
				mY1 = e.getY();
				mIsSelecting = true;
				mIsLassoSelect = !e.isAltDown();
				if (mIsLassoSelect) {
					mLassoRegion = new GenericPolygon();
					mLassoRegion.addPoint(mX1, mY1);
					mLassoRegion.addPoint(mX1, mY1);
					}
				}
			} );

		mCanvas.addEventHandler(MouseEvent.MOUSE_RELEASED, e -> {
			if (mIsSelecting) {
				mIsSelecting = false;
				repaint();
				}
			} );

		mCanvas.addEventHandler(MouseEvent.MOUSE_DRAGGED, e -> {
			if (mIsSelecting)
				updateSelection(e.getX(), e.getY());
			} );

		mCanvas.setOnContextMenuRequested(e -> {
			showPopupMenu(e.getScreenX(), e.getScreenY());
			e.consume();
			} );
		}

	/**
	 * @return true if a drag is allowed and x,y are within the molecule's area, where a drag would be started
	 */
	private boolean isInDragRect(double x, double y) {
		return mDepictor != null
			&& (mAllowedDragAction & DnDConstants.ACTION_COPY) != 0
			&& shrink(mDepictor.getBoundingRect()).contains(x, y);
		}

	private GenericRectangle shrink(GenericRectangle rect) {
		int margin = HiDPIHelper.scale(DRAG_MARGIN);
		double marginX = Math.min(margin, rect.width / 6);
		double marginY = Math.min(margin, rect.height / 6);
		return new GenericRectangle(rect.x+marginX, rect.y+marginY, rect.width-2*marginX, rect.height-2*marginY);
		}

	private void updateSelection(double x, double y) {
		mX2 = x;
		mY2 = y;

		GenericShape selectedShape = null;

		if (mIsLassoSelect) {
			if ((Math.abs(mX2 - mLassoRegion.getX(mLassoRegion.getSize() - 1)) >= 3)
			 || (Math.abs(mY2 - mLassoRegion.getY(mLassoRegion.getSize() - 1)) >= 3)) {
				mLassoRegion.removeLastPoint();
				mLassoRegion.addPoint(mX2, mY2);
				mLassoRegion.addPoint(mX1, mY1);

				selectedShape = mLassoRegion;
				}
			}
		else {
			selectedShape = new GenericRectangle(Math.min(mX1, mX2), Math.min(mY1, mY2), Math.abs(mX2 - mX1), Math.abs(mY2 - mY1));
			}

		if (selectedShape != null && mDepictor != null) {
			for (int atom=0; atom<mDisplayMol.getAllAtoms(); atom++) {
				boolean isSelected = selectedShape.contains((int)mDepictor.getAtomX(atom), (int)mDepictor.getAtomY(atom));
				if (isSelected != mDisplayMol.isSelectedAtom(atom))
					mDisplayMol.setAtomSelection(atom, isSelected);
				}
			}

		repaint();
		}

	private void showPopupMenu(double screenX, double screenY) {
		if (mMol == null)
			return;

		if (mPopup != null)
			mPopup.hide();

		mPopup = new ContextMenu();
		boolean hasAtoms = (mMol.getAllAtoms() != 0);

		if (mClipboardHandler != null) {
			addMenuItem(ITEM_COPY, hasAtoms, () -> mClipboardHandler.copyMolecule(mMol));
			addMenuItem(ITEM_COPY_SMILES, hasAtoms, () -> copySmiles());

			if (mIsEditable) {
				String itemText = StructureNameResolver.getInstance() == null ? ITEM_PASTE : ITEM_PASTE_WITH_NAME;
				addMenuItem(itemText, true, () -> paste());
				}

			mPopup.getItems().add(new SeparatorMenuItem());
			}

		if (mIsEditable) {
			addMenuItem(ITEM_CLEAR, hasAtoms, () -> {
				mMol.clear();
				mDisplayMol = mMol;
				structureChanged();
				} );

			mPopup.getItems().add(new SeparatorMenuItem());

			addMenuItem(ITEM_OPEN_MOLFILE, true, () -> openMolfile());

			mPopup.getItems().add(new SeparatorMenuItem());
			}

		addMenuItem(ITEM_SAVE_MOLFILE_V2, hasAtoms, () -> saveMolfile(false));
		addMenuItem(ITEM_SAVE_MOLFILE_V3, hasAtoms, () -> saveMolfile(true));

		mPopup.show(this, screenX, screenY);
		}

	private void addMenuItem(String text, boolean enabled, Runnable action) {
		MenuItem item = new MenuItem(text);
		item.setDisable(!enabled);
		item.setOnAction(e -> action.run());
		mPopup.getItems().add(item);
		}

	private void copySmiles() {
		ClipboardContent content = new ClipboardContent();
		content.putString(new IsomericSmilesCreator(mMol).getSmiles());
		Clipboard.getSystemClipboard().setContent(content);
		}

	private void paste() {
		if (!mIsEditable)
			return;

		int smartsMode = mAllowFragmentStatusChangeOnPasteOrDrop ? SmilesParser.SMARTS_MODE_GUESS
					   : mMol.isFragment() ? SmilesParser.SMARTS_MODE_IS_SMARTS : SmilesParser.SMARTS_MODE_IS_SMILES;
		StereoMolecule mol = mClipboardHandler.pasteMolecule(true, smartsMode);
		if (mol != null) {
			if (!mAllowFragmentStatusChangeOnPasteOrDrop)
				mol.setFragment(mMol.isFragment());
			mMol = mol;
			mDisplayMol = mol;
			structureChanged();
			}
		else {
			showWarningMessage("No molecule on clipboard!");
			}
		}

	private void openMolfile() {
		FileChooser chooser = createMolfileChooser("Open Structure");
		File file = chooser.showOpenDialog(getWindow());
		if (file == null)
			return;

		StereoMolecule mol = new MolfileParser().getCompactMolecule(file);
		if (mol != null) {
			if (!mAllowFragmentStatusChangeOnPasteOrDrop)
				mol.setFragment(mMol.isFragment());
			mMol = mol;
			mDisplayMol = mol;
			structureChanged();
			}
		else {
			showWarningMessage("Couldn't read molecule!");
			}
		}

	private void saveMolfile(boolean isV3) {
		String name = mMol.getName();
		if (name == null || name.isEmpty())
			name = "Unknown Molecule";

		FileChooser chooser = createMolfileChooser("Save Structure");
		chooser.setInitialFileName(name+".mol");
		File file = chooser.showSaveDialog(getWindow());
		if (file == null)
			return;

		String molfile = isV3 ? new MolfileV3Creator(mMol, false).getMolfile() : new MolfileCreator(mMol, false).getMolfile();
		try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), StandardCharsets.UTF_8))) {
			writer.write(molfile);
			}
		catch (IOException ioe) {
			showWarningMessage("Couldn't save molfile: "+ioe.getMessage());
			}
		}

	private FileChooser createMolfileChooser(String title) {
		FileChooser chooser = new FileChooser();
		chooser.setTitle(title);
		chooser.getExtensionFilters().add(new FileChooser.ExtensionFilter("MDL Molfiles", "*.mol"));
		return chooser;
		}

	private Window getWindow() {
		return (getScene() == null) ? null : getScene().getWindow();
		}

	protected void showWarningMessage(String msg) {
		if (!Platform.isFxApplicationThread()) {
			Platform.runLater(() -> showWarningMessage(msg));
			return;
			}

		mWarningMessage = msg;
		repaint();

		if (mWarningTimer != null)
			mWarningTimer.stop();
		mWarningTimer = new PauseTransition(Duration.millis(WARNING_MILLIS));
		mWarningTimer.setOnFinished(e -> {
			mWarningMessage = null;
			repaint();
			} );
		mWarningTimer.play();
		}

	private void informListeners() {
		if (mListener != null)
			for (int i = 0; i<mListener.size(); i++)
				mListener.get(i).structureChanged(mMol);
		}

	private void initializeDragAndDrop(int dragAction, int dropAction) {
		mAllowedDragAction = dragAction;
		mAllowedDropAction = dropAction;
		mAllowFragmentStatusChangeOnPasteOrDrop = false;

		if (dragAction != DnDConstants.ACTION_NONE) {
			mCanvas.setOnDragDetected(e -> {
				// If selectable, then dragging outside of the molecule's area selects atoms.
				if (!mIsSelecting
				 && mMol.getAllAtoms() != 0
				 && (!mIsSelectable || isInDragRect(e.getX(), e.getY()))) {
					ClipboardContent content = getMoleculeContent(e.getX(), e.getY());
					if (content != null) {
						Dragboard db = mCanvas.startDragAndDrop(toTransferModes(mAllowedDragAction));
						db.setContent(content);
						SnapshotParameters parameters = new SnapshotParameters();
						parameters.setFill(Color.TRANSPARENT);
						db.setDragView(mCanvas.snapshot(parameters, null), e.getX(), e.getY());
						mIsDraggingThis = true;
						}
					}
				e.consume();
				} );

			mCanvas.setOnDragDone(e -> {
				mIsDraggingThis = false;
				e.consume();
				} );
			}

		if (dropAction != DnDConstants.ACTION_NONE) {
			setOnDragEntered(e -> {
				if (isAcceptableDrag(e))
					onDragEnter();
				e.consume();
				} );

			setOnDragOver(e -> {
				if (isAcceptableDrag(e)) {
					e.acceptTransferModes(toTransferModes(mAllowedDropAction));
					onDragOver();
					}
				e.consume();
				} );

			setOnDragExited(e -> {
				if (isAcceptableDrag(e))
					onDragExit();
				e.consume();
				} );

			setOnDragDropped(e -> {
				StereoMolecule m = canDrop() ? FXMoleculeTransfer.readMolecule(e.getDragboard(), true) : null;
				if (m != null) {
					boolean isFragment = mMol.isFragment();
					mMol = new StereoMolecule(m);
					mMol.removeAtomColors();
					mMol.removeBondHiliting();
					if (!mAllowFragmentStatusChangeOnPasteOrDrop)
						mMol.setFragment(isFragment);
					mDisplayMol = mMol;
					mIDCode = null;
					repaint();
					informListeners();
					onDrop();
					}
				updateBorder(false);
				e.setDropCompleted(m != null);
				e.consume();
				} );
			}
		}

	private boolean isAcceptableDrag(DragEvent e) {
		return e.getGestureSource() != mCanvas && canDrop() && FXMoleculeTransfer.hasMolecule(e.getDragboard());
		}

	private static TransferMode[] toTransferModes(int dndAction) {
		ArrayList<TransferMode> modes = new ArrayList<>();
		if ((dndAction & DnDConstants.ACTION_COPY) != 0)
			modes.add(TransferMode.COPY);
		if ((dndAction & DnDConstants.ACTION_MOVE) != 0)
			modes.add(TransferMode.MOVE);
		if ((dndAction & DnDConstants.ACTION_LINK) != 0)
			modes.add(TransferMode.LINK);
		return modes.toArray(new TransferMode[0]);
		}

	/**
	 * Creates the content to be dragged. Override this, if the content shall depend on the drag start position.
	 * @param x drag start x within the canvas
	 * @param y drag start y within the canvas
	 * @return content or null, if no drag shall be started
	 */
	protected ClipboardContent getMoleculeContent(double x, double y) {
		return FXMoleculeTransfer.createContent(mMol);
		}

	// Drag notifications if needed by subclasses; called when a drag, which could be dropped, enters, exits, or moves over this view
	protected void onDragEnter() {}
	protected void onDragExit() {}
	protected void onDragOver() {}
	protected void onDrop() {}

	private void updateBorder(boolean showBorder) {
		if (mShowBorder != showBorder) {
			mShowBorder = showBorder;
			repaint();
			}
		}
	}
