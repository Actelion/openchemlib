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
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.io.CompoundFileHelper;
import com.actelion.research.chem.io.RDFileParser;
import com.actelion.research.chem.io.RXNFileParser;
import com.actelion.research.chem.name.StructureNameResolver;
import com.actelion.research.chem.reaction.IReactionMapper;
import com.actelion.research.chem.reaction.MCSReactionMapper;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionArrow;
import com.actelion.research.gui.clipboard.IClipboardHandler;
import com.actelion.research.gui.dnd.MoleculeDropAdapter;
import com.actelion.research.gui.editor.EditorEvent;
import com.actelion.research.gui.generic.*;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.gui.hidpi.ScaledEditorKit;
import com.actelion.research.gui.swing.SwingDrawContext;
import com.actelion.research.util.ColorHelper;
import com.actelion.research.gui.swing.SwingCursorHelper;

import javax.swing.*;
import javax.swing.text.html.HTMLEditorKit;
import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DropTarget;
import java.awt.event.*;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.TreeMap;

@Deprecated
public class JDrawArea extends JPanel implements ActionListener, KeyListener, MouseListener, MouseMotionListener {
	static final long serialVersionUID = 0x20061019;

	public static final int MODE_MULTIPLE_FRAGMENTS = 1;
	public static final int MODE_MARKUSH_STRUCTURE = 2;
	public static final int MODE_REACTION = 4;
	public static final int MODE_DRAWING_OBJECTS = 8;

	private static final int MAX_CONNATOMS = 8;
	private static final int MIN_BOND_LENGTH_SQUARE = 100;

	private static final int KEY_IS_ATOM_LABEL = 1;
	private static final int KEY_IS_SUBSTITUENT = 2;
	private static final int KEY_IS_VALID_START = 3;
	private static final int KEY_IS_INVALID = 4;

	private static final String ITEM_COPY_STRUCTURE = "Copy Structure";
	private static final String ITEM_COPY_REACTION = "Copy Reaction";
	private static final String ITEM_PASTE_STRUCTURE = "Paste Structure";
	private static final String ITEM_PASTE_REACTION = "Paste Reaction";
	private static final String ITEM_PASTE_WITH_NAME = ITEM_PASTE_STRUCTURE+" or Name";
	private static final String ITEM_LOAD_REACTION = "Open Reaction File...";
	private static final String ITEM_ADD_AUTO_MAPPING = "Auto-Map Reaction";
	private static final String ITEM_REMOVE_MAPPING = "Remove Manual Atom Mapping";
	private static final String ITEM_FLIP_HORIZONTALLY = "Flip Horizontally";
	private static final String ITEM_FLIP_VERTICALLY = "Flip Vertically";

	private static final long WARNING_MILLIS = 1200;

	private static final float FRAGMENT_MAX_CLICK_DISTANCE = 24.0f;
	private static final float FRAGMENT_GROUPING_DISTANCE = 1.4f;	// in average bond lengths
	private static final float FRAGMENT_CLEANUP_DISTANCE = 1.5f;	// in average bond lengths
	private static final float DEFAULT_ARROW_LENGTH = 0.08f;		// relative to panel width

	protected static final int UPDATE_NONE = 0;
	protected static final int UPDATE_REDRAW = 1;
	// redraw molecules and drawing objects with their current coordinates
	protected static final int UPDATE_CHECK_VIEW = 2;
	// redraw with on-the-fly coordinate transformation only if current coords do not fit within view area
	// (new coords are generated for one draw() only; the original coords are not changed)
	protected static final int UPDATE_CHECK_COORDS = 3;
	// redraw with in-place coordinate transformation only if current coords do not fit within view area
	// (the original atom and object coords are replaced by the new ones)
	protected static final int UPDATE_SCALE_COORDS = 4;
	// redraw with in-place coordinate transformation; molecules and objects are scaled to fill
	// the view unless the maximum average bond length reaches the optimum.
	protected static final int UPDATE_SCALE_COORDS_USE_FRAGMENTS = 5;
	// as UPDATE_SCALE_COORDS but uses fragments from mFragment rather than creating a
	// fresh mFragment list from mMol. Used for setting a reaction or fragment list from outside.
	protected static final int UPDATE_INVENT_COORDS = 6;
	// redraw with in-place coordinate transformation; first all molecules' coordinates
	// are generated from scratch, then molecules and objects are scaled to fill
	// the view unless the maximum average bond length reaches the optimum.

	private static final Color DEFAULT_SELECTION_BACKGROUND = new Color(128,164,192);

	private static final int ALLOWED_DROP_ACTIONS = DnDConstants.ACTION_COPY_OR_MOVE;

	private static final int cRequestNone = 0;
	private static final int cRequestNewBond = 1;
	private static final int cRequestNewChain = 2;
	private static final int cRequestMoveSingle = 3;
	private static final int cRequestMoveSelected = 4;
	private static final int cRequestLassoSelect = 5;
	private static final int cRequestSelectRect = 6;
	private static final int cRequestZoomAndRotate = 7;
	private static final int cRequestMapAtoms = 8;
	private static final int cRequestCopySelected = 9;
	private static final int cRequestMoveObject = 10;
	private static final int cRequestCopyObject = 11;

	private static IReactionMapper sMapper;
	private Dimension mSize;
	private int mMode, mChainAtoms, mCurrentTool, mOtherAtom, mOtherMass, mOtherValence, mOtherRadical,
		mCurrentHiliteAtom, mCurrentHiliteBond, mPendingRequest,
		mCurrentCursor, mReactantCount, mUpdateMode, mDisplayMode, mAtom1, mAtom2;
	private int[] mChainAtom, mFragmentNo, mHiliteBondSet;
	private double mX1, mY1, mX2, mY2;
	private double[] mX, mY, mChainAtomX, mChainAtomY;
	private boolean mShiftIsDown, mAltIsDown, mControlIsDown, mMouseIsDown,
					mIsAddingToSelection, mAtomColorSupported, mAllowQueryFeatures;
	private boolean[] mIsSelectedAtom, mIsSelectedObject;
	private String mOtherLabel,mWarningMessage;
	private String[] mAtomText;
	private ExtendedDepictor mDepictor;
	private StereoMolecule mMol;	    // molecule being modified directly by the drawing editor
	private Molecule mUndoMol;          // molecule in undo buffer
	private StereoMolecule[] mFragment;	// in case of MODE_MULTIPLE_FRAGMENTS contains valid stereo fragments
										// for internal and external read-only-access (reconstructed at any change)
	private DrawingObjectList mDrawingObjectList, mUndoDrawingObjectList;
	private AbstractDrawingObject mCurrentHiliteObject;
	private GenericPolygon mLassoRegion;
	private ArrayList<GenericEventListener<EditorEvent>> mListeners;
	private IClipboardHandler mClipboardHandler;
	private JDialog mHelpDialog;
	private StringBuilder mAtomKeyStrokeBuffer;

	/**
	 * @param mol  an empty or valid stereo molecule
	 * @param mode 0 or a meaningful combination of the mode flags, e.g. MODE_REACTION | MODE_DRAWING_OBJECTS
	 */
	public JDrawArea(StereoMolecule mol, int mode)
	{
		mMol = mol;
		mMode = mode;

		setFocusable(true);
		addKeyListener(this);
		addMouseListener(this);
		addMouseMotionListener(this);

		mListeners = new ArrayList<>();

		mCurrentTool = JDrawToolbar.cToolStdBond;
		mCurrentHiliteAtom = -1;
		mCurrentHiliteBond = -1;
		mCurrentHiliteObject = null;
		mAtom1 = -1;
		mOtherAtom = 6;
		mOtherMass = 0;
		mOtherValence = -1;
		mOtherRadical = 0;
		mOtherLabel = null;
		mAllowQueryFeatures = true;
		mPendingRequest = cRequestNone;
		mCurrentCursor = SwingCursorHelper.cPointerCursor;
		mAtomKeyStrokeBuffer = new StringBuilder();

		if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE)) != 0) {
			mMode |= (MODE_MULTIPLE_FRAGMENTS);
		}

		if ((mMode & (MODE_DRAWING_OBJECTS | MODE_REACTION)) != 0) {
			mDrawingObjectList = new DrawingObjectList();
		}

		mUpdateMode = UPDATE_SCALE_COORDS;

		initializeDragAndDrop(ALLOWED_DROP_ACTIONS);

/*		dumpBytesOfGif("/chainCursor.gif");
		dumpBytesOfGif("/deleteCursor.gif");
		dumpBytesOfGif("/handCursor.gif");
		dumpBytesOfGif("/handPlusCursor.gif");
		dumpBytesOfGif("/fistCursor.gif");
		dumpBytesOfGif("/lassoCursor.gif");
		dumpBytesOfGif("/lassoPlusCursor.gif");
		dumpBytesOfGif("/rectCursor.gif");
		dumpBytesOfGif("/rectPlusCursor.gif");	*/
	}

	/**
	 * Call this after initialization to get clipboard support
	 *
	 * @param h
	 */
	public void setClipboardHandler(IClipboardHandler h)
	{
		mClipboardHandler = h;
	}

	private void update(int mode)
	{
		mUpdateMode = Math.max(mUpdateMode, mode);
		repaint();
	}

	public static void setReactionMapper(IReactionMapper mapper)
	{
		sMapper = mapper;
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		((Graphics2D) g).setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);

		Dimension theSize = getSize();
		if (mSize == null || mSize.width != theSize.width || mSize.height != theSize.height) {
			mSize = theSize;
			if (mUpdateMode < UPDATE_CHECK_COORDS) {
				mUpdateMode = UPDATE_CHECK_COORDS;
			}
		}

		Color background = UIManager.getColor("TextArea.background");
		Color foreground = UIManager.getColor("TextArea.foreground");

		g.setColor(background);
		g.fillRect(0, 0, theSize.width, theSize.height);

		if ((mMode & MODE_REACTION) != 0 && mDrawingObjectList.size() == 0) {
			float mx = 0.5f * (float) theSize.width;
			float my = 0.5f * (float) theSize.height;
			float dx = 0.5f * DEFAULT_ARROW_LENGTH * (float) theSize.width;
			ReactionArrow arrow = new ReactionArrow();
			arrow.setCoordinates(mx - dx, my, mx + dx, my);
			arrow.setDeletable(false);
			mDrawingObjectList.add(arrow);
		}

		SwingDrawContext context = new SwingDrawContext((Graphics2D)g);

		boolean isScaledView = false;
		if (mUpdateMode != UPDATE_NONE) {
			if ((mMode & MODE_MULTIPLE_FRAGMENTS) != 0
			 && mUpdateMode != UPDATE_SCALE_COORDS_USE_FRAGMENTS)
				analyzeFragmentMembership();

			mDepictor = ((mMode & MODE_REACTION) != 0) ?
				new ExtendedDepictor(new Reaction(mFragment, mReactantCount), mDrawingObjectList, false)
				: ((mMode & MODE_MARKUSH_STRUCTURE) != 0) ?
				new ExtendedDepictor(mFragment, mReactantCount, mDrawingObjectList)
				: ((mMode & MODE_MULTIPLE_FRAGMENTS) != 0) ?
				new ExtendedDepictor(mFragment, mDrawingObjectList)
				: new ExtendedDepictor(mMol, mDrawingObjectList);

			mDepictor.setForegroundColor(foreground, background);
			mDepictor.setFragmentNoColor(((mMode & MODE_MULTIPLE_FRAGMENTS) == 0) ? 0
					: LookAndFeelHelper.isDarkLookAndFeel() ? ColorHelper.brighter(background.getRGB(), 0.85f)
															: ColorHelper.darker(background.getRGB(), 0.85f));
			mDepictor.setDisplayMode(mDisplayMode
				| AbstractDepictor.cDModeHiliteAllQueryFeatures
				| ((mCurrentTool == JDrawToolbar.cToolMapper) ?
					AbstractDepictor.cDModeShowMapping
				  | AbstractDepictor.cDModeSuppressCIPParity : 0));

			if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE | MODE_MULTIPLE_FRAGMENTS)) == 0) {
				mDepictor.getMoleculeDepictor(0).setAtomText(mAtomText);
			}

			switch (mUpdateMode) {
				case UPDATE_INVENT_COORDS:
				case UPDATE_SCALE_COORDS:
				case UPDATE_SCALE_COORDS_USE_FRAGMENTS:
					cleanupCoordinates(context, ((Graphics2D) g));
					break;
				case UPDATE_CHECK_COORDS:
					DepictorTransformation t1 = mDepictor.updateCoords(context, new GenericRectangle(0, 0, theSize.width, theSize.height), 0);
					if (t1 != null && (mMode & MODE_MULTIPLE_FRAGMENTS) != 0) {
						// in fragment mode depictor transforms mFragment[] rather than mMol
						t1.applyTo(mMol);
					}
					break;
				case UPDATE_CHECK_VIEW:
					DepictorTransformation t2 = mDepictor.validateView(context, new GenericRectangle(0, 0, theSize.width, theSize.height), 0);
					isScaledView = (t2 != null && !t2.isVoidTransformation());
					break;
			}
			mUpdateMode = UPDATE_NONE;
		}

		if (mDepictor != null) {
			mDepictor.paintFragmentNumbers(context);
		}

		// don't hilite anything when the view is scaled and object coords don't reflect screen coords
		if (!isScaledView) {
			drawHiliting(context, g);
		}

		if (mDepictor != null) {
			mDepictor.paintStructures(context);
			mDepictor.paintDrawingObjects(context);
		}

		if (mCurrentHiliteAtom != -1 && mAtomKeyStrokeBuffer.length() != 0) {
			int x = (int) mMol.getAtomX(mCurrentHiliteAtom);
			int y = (int) mMol.getAtomY(mCurrentHiliteAtom);
			String s = mAtomKeyStrokeBuffer.toString();
			int validity = getAtomKeyStrokeValidity(s);
			g.setColor((validity == KEY_IS_ATOM_LABEL) ? foreground
					 : (validity == KEY_IS_SUBSTITUENT) ? Color.BLUE
					 : (validity == KEY_IS_VALID_START) ? Color.GRAY : Color.RED);
			if (validity == KEY_IS_INVALID)
				s = s + "<unknown>";
			g.setFont(g.getFont().deriveFont(0, 24));
			g.drawString(s, x, y);
		}

		g.setColor(foreground);
		switch (mPendingRequest) {
			case cRequestNewBond:
				int x1, y1, x2, y2, xdiff, ydiff;
				x1 = (int) mX1;
				y1 = (int) mY1;
				if (mCurrentHiliteAtom == -1 || mCurrentHiliteAtom == mAtom1) {
					x2 = (int) mX2;
					y2 = (int) mY2;
				} else {
					x2 = (int) mMol.getAtomX(mCurrentHiliteAtom);
					y2 = (int) mMol.getAtomY(mCurrentHiliteAtom);
				}
				switch (mCurrentTool) {
					case JDrawToolbar.cToolStdBond:
						g.drawLine(x1, y1, x2, y2);
						break;
					case JDrawToolbar.cToolUpBond:
						int[] x = new int[3];
						int[] y = new int[3];
						xdiff = (y1 - y2) / 9;
						ydiff = (x2 - x1) / 9;
						x[0] = x1;
						y[0] = y1;
						x[1] = x2 + xdiff;
						y[1] = y2 + ydiff;
						x[2] = x2 - xdiff;
						y[2] = y2 - ydiff;
						g.fillPolygon(x, y, 3);
						break;
					case JDrawToolbar.cToolDownBond:
						int xx1, xx2, yy1, yy2;
						xdiff = x2 - x1;
						ydiff = y2 - y1;
						for (int i = 2; i < 17; i += 2) {
							xx1 = x1 + i * xdiff / 17 - i * ydiff / 128;
							yy1 = y1 + i * ydiff / 17 + i * xdiff / 128;
							xx2 = x1 + i * xdiff / 17 + i * ydiff / 128;
							yy2 = y1 + i * ydiff / 17 - i * xdiff / 128;
							g.drawLine(xx1, yy1, xx2, yy2);
						}
						break;
				}
				break;
			case cRequestNewChain:
				if (mChainAtoms > 0) {
					g.drawLine((int) mX1, (int) mY1, (int) mChainAtomX[0], (int) mChainAtomY[0]);
				}
				if (mChainAtoms > 1) {
					for (int i = 1; i < mChainAtoms; i++) {
						g.drawLine((int) mChainAtomX[i - 1], (int) mChainAtomY[i - 1],
							(int) mChainAtomX[i], (int) mChainAtomY[i]);
					}
				}
				break;
			case cRequestLassoSelect:
				g.setColor(lassoColor());
				java.awt.Polygon p = new java.awt.Polygon();
				for (int i=0; i<mLassoRegion.getSize(); i++)
					p.addPoint(Math.round((float)mLassoRegion.getX(i)), Math.round((float)mLassoRegion.getY(i)));
				g.drawPolygon(p);
				g.setColor(foreground);
				break;
			case cRequestSelectRect:
				int x = (mX1 < mX2) ? (int) mX1 : (int) mX2;
				int y = (mY1 < mY2) ? (int) mY1 : (int) mY2;
				int w = (int) Math.abs(mX2 - mX1);
				int h = (int) Math.abs(mY2 - mY1);
				g.setColor(lassoColor());
				g.drawRect(x, y, w, h);
				g.setColor(foreground);
				break;
			case cRequestMapAtoms:
				x1 = (int) mX1;
				y1 = (int) mY1;
				if (mCurrentHiliteAtom == -1 || mCurrentHiliteAtom == mAtom1) {
					x2 = (int) mX2;
					y2 = (int) mY2;
				} else {
					x2 = (int) mMol.getAtomX(mCurrentHiliteAtom);
					y2 = (int) mMol.getAtomY(mCurrentHiliteAtom);
				}
				g.setColor(mapToolColor());
				g.drawLine(x1, y1, x2, y2);
				g.setColor(foreground);
				break;
		}

		if (mWarningMessage != null) {
			int fontSize = HiDPIHelper.scale(12);
			g.setFont(getFont().deriveFont(Font.BOLD, (float)fontSize));
			Color original = g.getColor();
			g.setColor(Color.RED);
			FontMetrics metrics = g.getFontMetrics();
			Rectangle2D bounds = metrics.getStringBounds(mWarningMessage, g);
			g.drawString(mWarningMessage, (int)(theSize.width-bounds.getWidth())/2, metrics.getHeight());
			g.setColor(original);
		}
	}

	public static Color lassoColor() {
		Color selectionColor = selectionColor();
		return ColorHelper.createColor(selectionColor, LookAndFeelHelper.isDarkLookAndFeel() ? 0.65f : 0.35f);
	}

	public static Color selectionColor() {
		Color selectionColor = UIManager.getColor("TextArea.selectionBackground");
		return (selectionColor != null) ? selectionColor : DEFAULT_SELECTION_BACKGROUND;
	}

	public static Color mapToolColor() {
		Color background = UIManager.getColor("TextArea.background");
		return ColorHelper.getContrastColor(new Color(128, 0, 0), background);
	}

	public static Color chainHiliteColor() {
		Color background = UIManager.getColor("TextArea.background");
		Color selectionColor = selectionColor();
		return ColorHelper.intermediateColor(selectionColor, background, 0.5f);
	}

	private void drawHiliting(GenericDrawContext context, Graphics g)
	{
		if (mHiliteBondSet != null) {
			g.setColor(chainHiliteColor());
			for (int i = 0; i < mHiliteBondSet.length; i++) {
				hiliteBond(g, mHiliteBondSet[i]);
			}
		}

		if (mCurrentHiliteAtom != -1) {
			g.setColor(selectionColor());
			hiliteAtom(g, mCurrentHiliteAtom);
			if (mCurrentTool == JDrawToolbar.cToolMapper) {
				int mapNo = mMol.getAtomMapNo(mCurrentHiliteAtom);
				if (mapNo != 0) {
					for (int atom = 0; atom < mMol.getAtoms(); atom++) {
						if (atom != mCurrentHiliteAtom
							&& mMol.getAtomMapNo(atom) == mapNo) {
							hiliteAtom(g, atom);
						}
					}
				}
			}
		}

		if (mCurrentHiliteBond != -1) {
			g.setColor(selectionColor());
			hiliteBond(g, mCurrentHiliteBond);
		}

		if (mCurrentHiliteObject != null) {
			mCurrentHiliteObject.hilite(context);
		}
	}

	private void hiliteAtom(Graphics g, int atom)
	{
		int radius = (int) (0.32f * mMol.getAverageBondLength());

		int x = (int) mMol.getAtomX(atom);
		int y = (int) mMol.getAtomY(atom);
		g.fillOval(x - radius, y - radius, 2 * radius, 2 * radius);
	}

	private void hiliteBond(Graphics g, int bond)
	{
		int width = (int) (0.32f * mMol.getAverageBondLength());

		int x1 = (int) mMol.getAtomX(mMol.getBondAtom(0, bond));
		int y1 = (int) mMol.getAtomY(mMol.getBondAtom(0, bond));
		int x2 = (int) mMol.getAtomX(mMol.getBondAtom(1, bond));
		int y2 = (int) mMol.getAtomY(mMol.getBondAtom(1, bond));

		Stroke oldStroke = ((Graphics2D) g).getStroke();
		((Graphics2D) g).setStroke(new BasicStroke((float) width, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
		((Graphics2D) g).drawLine(x1, y1, x2, y2);
		((Graphics2D) g).setStroke(oldStroke);
	}

	public void addDrawAreaListener(GenericEventListener<EditorEvent> l)
	{
		mListeners.add(l);
	}

	protected void buttonPressed(int button)
	{
		switch (button) {
			case JDrawToolbar.cButtonClear:
				clearAll();
				return;
			case JDrawToolbar.cButtonCleanStructure:
				storeState();
				fireMoleculeChanged();
				update(UPDATE_INVENT_COORDS);
				return;
			case JDrawToolbar.cButtonUndo:
				restoreState();
				fireMoleculeChanged();
				update(UPDATE_CHECK_VIEW);
				return;
		}
	}

	public void clearAll() {
		if (mDrawingObjectList != null) {
			mDrawingObjectList.clear();
			update(UPDATE_REDRAW);
		}
		storeState();
		boolean isFragment = mMol.isFragment();
		mMol.clear();
		mMol.setFragment(isFragment);
		if (mUndoMol.getAllAtoms() != 0) {
			fireMoleculeChanged();
		}
		update(UPDATE_REDRAW);
	}

	public void toolChanged(int newTool)
	{
		if (mCurrentTool != newTool) {
			setOtherAtom(-1, 0, -1, 0, null);
			if (mCurrentTool == JDrawToolbar.cToolMapper
				|| newTool == JDrawToolbar.cToolMapper) {
				update(UPDATE_REDRAW);
			}

			mCurrentTool = newTool;
		}
	}

	private void setOtherAtom(int atomicNo, int mass, int valence, int radical, String customLabel)
	{
		mOtherAtom = atomicNo;
		mOtherMass = mass;
		mOtherValence = valence;
		mOtherRadical = radical;
		mOtherLabel = customLabel;
	}

	public void actionPerformed(ActionEvent e)
	{
		String command = e.getActionCommand();
		if (command.equals(ITEM_COPY_STRUCTURE) || command.equals(ITEM_COPY_REACTION)) {
			copy();
		} else if (command.equals(ITEM_PASTE_REACTION)) {
			pasteReaction();
		} else if (command.startsWith(ITEM_PASTE_STRUCTURE)) {
			pasteMolecule();
		} else if (e.getActionCommand().equals(ITEM_LOAD_REACTION)) {
			openReaction();
		} else if (e.getActionCommand().equals(ITEM_ADD_AUTO_MAPPING)) {
			autoMapReaction();
			fireMoleculeChanged();
			mUpdateMode = Math.max(mUpdateMode, UPDATE_REDRAW);
			repaint();
		} else if (e.getActionCommand().equals(ITEM_REMOVE_MAPPING)) {
			removeManualMapping();
		} else if (e.getActionCommand().equals(ITEM_FLIP_HORIZONTALLY)) {
			flip(true);
		} else if (e.getActionCommand().equals(ITEM_FLIP_VERTICALLY)) {
			flip(false);
		} else if (command.startsWith("atomColor")) {
			int index = command.indexOf(':');
			int atom = Integer.parseInt(command.substring(9, index));
			int color = Integer.parseInt(command.substring(index + 1));
			if (mMol.isSelectedAtom(atom)) {
				for (int i = 0; i < mMol.getAtoms(); i++) {
					if (mMol.isSelectedAtom(i)) {
						mMol.setAtomColor(i, color);
					}
				}
			} else {
				mMol.setAtomColor(atom, color);
			}
		}
	}

	private void removeManualMapping() {
		boolean changed = false;
		for (int atom = 0; atom < mMol.getAtoms(); atom++) {
			if (mMol.getAtomMapNo(atom) != 0 && !mMol.isAutoMappedAtom(atom)) {
				if (!changed) {
					storeState();
					changed = true;
				}

				mMol.setAtomMapNo(atom, 0, false);
				}
			}

		if (changed) {
			autoMapReaction();
			fireMoleculeChanged();
			mUpdateMode = Math.max(mUpdateMode, UPDATE_REDRAW);
			repaint();
		}
	}

	/**
	 * Checks, whether a copy operation would copy a molecule or reaction.
	 * @param doCopy if true, then the chemistry object is copied to the clipboard
	 * @return true, if reaction
	 */
	private boolean analyseCopy(boolean doCopy)
	{
		boolean isReaction = ((mMode & MODE_REACTION) != 0);
		boolean selectionFound = false;
		boolean isBothSideSelection = false;
		boolean isOnProductSide = false;
		ReactionArrow arrow = null;

		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
			if (mMol.isSelectedAtom(atom)) {
				if (!selectionFound) {
					selectionFound = true;
					if (!isReaction) {
						break;
					}

					arrow = (ReactionArrow) mDrawingObjectList.get(0);
					isOnProductSide = arrow.isOnProductSide(mMol.getAtomX(atom), mMol.getAtomY(atom));
				} else {
					if (isOnProductSide != arrow.isOnProductSide(mMol.getAtomX(atom), mMol.getAtomY(atom))) {
						isBothSideSelection = true;
						break;
					}
				}
			}
		}

		if (!doCopy)
			return isReaction && (isBothSideSelection || !selectionFound);

		if (isReaction) {
			if (isBothSideSelection) {
				copyReaction(true);
				return true;
			} else if (selectionFound) {
				copyMolecule(true);
				return false;
			} else {
				copyReaction(false);
				return true;
			}
		}

		copyMolecule(selectionFound);
		return false;
	}

	private void copy()
	{
		analyseCopy(true);
	}

	private boolean copyReaction(boolean selectionOnly)
	{
		Reaction rxn = selectionOnly ? getSelectedReaction() : getReaction();
		if (rxn != null && mClipboardHandler != null) {
			return mClipboardHandler.copyReaction(rxn);
		}

		return false;
	}

	private Reaction getSelectedReaction()
	{
		Reaction rxn = new Reaction();
		for (int i = 0; i < mFragment.length; i++) {
			StereoMolecule selectedMol = getSelectedCopy(mFragment[i]);
			if (selectedMol != null) {
				if (i < mReactantCount) {
					rxn.addReactant(selectedMol);
				} else {
					rxn.addProduct(selectedMol);
				}
			}
		}
		return rxn;
	}

	private StereoMolecule getSelectedCopy(StereoMolecule sourceMol)
	{
		int atomCount = 0;
		for (int atom = 0; atom < sourceMol.getAllAtoms(); atom++) {
			if (sourceMol.isSelectedAtom(atom)) {
				atomCount++;
			}
		}

		if (atomCount == 0) {
			return null;
		}

		int bondCount = 0;
		for (int bond = 0; bond < sourceMol.getAllBonds(); bond++) {
			if (sourceMol.isSelectedBond(bond)) {
				bondCount++;
			}
		}

		boolean[] includeAtom = new boolean[sourceMol.getAllAtoms()];
		for (int atom = 0; atom < sourceMol.getAllAtoms(); atom++) {
			includeAtom[atom] = sourceMol.isSelectedAtom(atom);
		}

		StereoMolecule destMol = new StereoMolecule(atomCount, bondCount);
		sourceMol.copyMoleculeByAtoms(destMol, includeAtom, false, null);
		return destMol;
	}

	private boolean copyMolecule(boolean selectionOnly)
	{
		if (mMol.getAllAtoms() != 0 && mClipboardHandler != null) {
			return mClipboardHandler.copyMolecule(selectionOnly ? getSelectedCopy(mMol) : mMol);
		}

		return false;
	}

	private void paste()
	{
		if ((mMode & MODE_REACTION) != 0) {
			if (pasteReaction()) {
				return;
			}
		}

		pasteMolecule();
	}

	private boolean pasteReaction()
	{
		boolean ret = false;
		if (mClipboardHandler != null) {
			Reaction rxn = mClipboardHandler.pasteReaction();
			if (rxn != null) {
				for (int i = 0; i < rxn.getMolecules(); i++) {
					rxn.getMolecule(i).setFragment(mMol.isFragment());
				}
				storeState();
				setReaction(rxn);
				ret = true;
			}
			else {
				showWarningMessage("No reaction on clipboard!");
			}
		}
		return ret;
	}

	private boolean pasteMolecule()
	{
		boolean ret = false;
		if (mClipboardHandler != null) {
			StereoMolecule mol = mClipboardHandler.pasteMolecule();
			if (mol != null && mol.getAllAtoms() != 0) {
				if (mol.getAllBonds() != 0)
					new Depictor2D(mol).updateCoords((Graphics2D)getGraphics(),
									new GenericRectangle(0, 0, this.getWidth(), this.getHeight()),
									AbstractDepictor.cModeInflateToMaxAVBL + (int)mMol.getAverageBondLength());

				storeState();
				if (mMol.getAllAtoms() == 0) {
					boolean isFragment = mMol.isFragment();
					mol.copyMolecule(mMol);
					mMol.setFragment(isFragment);
					moleculeChanged(true);
				} else {
					int originalAtoms = mMol.getAllAtoms();
					mMol.addMolecule(mol);
					for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
						mMol.setAtomSelection(atom, atom >= originalAtoms);
					}
					moleculeChanged(true);
				}
				ret = true;
			}
			else {
				showWarningMessage("No molecule on clipboard!");
			}
		}
		return ret;
	}

	private void openReaction() {
		File rxnFile = FileHelper.getFile(this, "Please select a reaction file",
				FileHelper.cFileTypeRXN | CompoundFileHelper.cFileTypeRD);
		if (rxnFile != null) {
			try {
				Reaction reaction = null;

				if (FileHelper.getFileType(rxnFile.getName()) == FileHelper.cFileTypeRXN) {
					reaction = new RXNFileParser().getReaction(rxnFile);
				}
				else {
					RDFileParser rdfParser = new RDFileParser(rxnFile);
					if (rdfParser.isReactionNext())
						reaction = rdfParser.getNextReaction();
				}

				if (reaction != null) {
					for (int i = 0; i < reaction.getMolecules(); i++) {
						reaction.getMolecule(i).setFragment(mMol.isFragment());
					}
					storeState();
					setReaction(reaction);
				}
			}
			catch (Exception ex) {}
		}
	}

	private void showWarningMessage(String msg) {
		mWarningMessage = msg;
		repaint();
		new Thread(() -> {
			try { Thread.sleep(WARNING_MILLIS); } catch (InterruptedException ie) {}
			mWarningMessage = null;
			repaint();
		} ).start();
	}

	public void mousePressed(MouseEvent e)
	{
		if (mCurrentHiliteAtom != -1 && mAtomKeyStrokeBuffer.length() != 0)
			expandAtomKeyStrokes(mAtomKeyStrokeBuffer.toString());

		mAtomKeyStrokeBuffer.setLength(0);

		if (handlePopupTrigger(e)) {
			return;
		}

		if ((e.getModifiers() & InputEvent.BUTTON1_MASK) != 0) {
			if (e.getClickCount() == 2) {
				return;
			}

			mMouseIsDown = true;
			updateCursor();
			mousePressedButton1(e);
		}
	}

	public void mouseReleased(MouseEvent e)
	{
		if (handlePopupTrigger(e)) {
			return;
		}

		if ((e.getModifiers() & InputEvent.BUTTON1_MASK) != 0) {
			if (e.getClickCount() == 2) {
				handleDoubleClick(e);
				return;
			}

			mMouseIsDown = false;
			updateCursor();
			mouseReleasedButton1(e);
		}
	}

	public void mouseEntered(MouseEvent e)
	{
		requestFocus();
		updateCursor();
	}

	public void mouseExited(MouseEvent e)
	{
	}

	public void mouseClicked(MouseEvent e)
	{
	}

	public void mouseDragged(MouseEvent e)
	{
		mMouseIsDown = true;
		mX2 = e.getX();
		mY2 = e.getY();

		boolean repaintNeeded = trackHiliting(mX2, mY2, true);

		switch (mPendingRequest) {
			case cRequestNewChain:
				double lastX, lastY;
				if (mChainAtoms > 0) {
					lastX = mChainAtomX[mChainAtoms - 1];
					lastY = mChainAtomY[mChainAtoms - 1];
				} else {
					lastX = 0;
					lastY = 0;
				}
				double avbl = mMol.getAverageBondLength();
				double s0 = (int) avbl;
				double s1 = (int) (0.866 * avbl);
				double s2 = (int) (0.5 * avbl);
				double dx = mX2 - mX1;
				double dy = mY2 - mY1;
				if (Math.abs(dy) > Math.abs(dx)) {
					mChainAtoms = (int) (2 * Math.abs(dy) / (s0 + s2));
					if (Math.abs(dy) % (s0 + s2) > s0) {
						mChainAtoms++;
					}
					mChainAtomX = new double[mChainAtoms];
					mChainAtomY = new double[mChainAtoms];
					if (mX2 < mX1) {
						s1 = -s1;
					}
					if (mY2 < mY1) {
						s0 = -s0;
						s2 = -s2;
					}
					for (int i = 0; i < mChainAtoms; i++) {
						mChainAtomX[i] = mX1 + ((i + 1) / 2) * s1;
						mChainAtomY[i] = mY1 + ((i + 1) / 2) * (s0 + s2);
						if ((i & 1) == 0) {
							mChainAtomY[i] += s0;
						}
					}
				} else {
					mChainAtoms = (int) (Math.abs(dx) / s1);
					mChainAtomX = new double[mChainAtoms];
					mChainAtomY = new double[mChainAtoms];
					if (mX2 < mX1) {
						s1 = -s1;
					}
					if (mY2 < mY1) {
						s2 = -s2;
					}
					for (int i = 0; i < mChainAtoms; i++) {
						mChainAtomX[i] = mX1 + (i + 1) * s1;
						mChainAtomY[i] = mY1;
						if ((i & 1) == 0) {
							mChainAtomY[i] += s2;
						}
					}
				}
				if (mChainAtoms > 0) {
					mChainAtom = new int[mChainAtoms];
					for (int i = 0; i < mChainAtoms; i++) {
						mChainAtom[i] = mMol.findAtom(mChainAtomX[i], mChainAtomY[i]);
						if (mChainAtom[i] != -1) {
							mChainAtomX[i] = mMol.getAtomX(mChainAtom[i]);
							mChainAtomY[i] = mMol.getAtomY(mChainAtom[i]);
						}
					}
					if (mChainAtomX[mChainAtoms - 1] != lastX
						|| mChainAtomY[mChainAtoms - 1] != lastY) {
						repaintNeeded = true;
					}
				} else if (lastX != 0 || lastY != 0) {
					repaintNeeded = true;
				}
				break;
			case cRequestNewBond:
				if ((mX2 - mX1) * (mX2 - mX1) + (mY2 - mY1) * (mY2 - mY1) < MIN_BOND_LENGTH_SQUARE) {
					suggestNewX2AndY2(mAtom1);
				}

				repaintNeeded = true;
				break;
			case cRequestMoveSingle:
				mMol.setAtomX(mAtom1, mX[mAtom1] + mX2 - mX1);
				mMol.setAtomY(mAtom1, mY[mAtom1] + mY2 - mY1);
				if (mAtom2 != -1) {
					mMol.setAtomX(mAtom2, mX[mAtom2] + mX2 - mX1);
					mMol.setAtomY(mAtom2, mY[mAtom2] + mY2 - mY1);
				}

				update(UPDATE_CHECK_VIEW);
				break;
			case cRequestCopySelected:
				duplicateSelected();
				mPendingRequest = cRequestMoveSelected;
			case cRequestMoveSelected:
				if (mDrawingObjectList != null) {
					for (AbstractDrawingObject drawingObject : mDrawingObjectList) {
						if (drawingObject.isSelected()) {
							drawingObject.translate(mX2, mY2);
						}
					}
				}
				for (int i = 0; i < mMol.getAllAtoms(); i++) {
					if (mMol.isSelectedAtom(i)) {
						mMol.setAtomX(i, mX[i] + mX2 - mX1);
						mMol.setAtomY(i, mY[i] + mY2 - mY1);
					}
				}
				update(UPDATE_CHECK_VIEW);
				break;
			case cRequestCopyObject:
				mDrawingObjectList.add(mCurrentHiliteObject.clone());
				mPendingRequest = cRequestMoveObject;
			case cRequestMoveObject:
				mCurrentHiliteObject.translate(mX2, mY2);
				update(UPDATE_CHECK_VIEW);
				break;
			case cRequestZoomAndRotate:
				boolean selectedAtomsFound = false;
				for (int atom = 0; atom < mMol.getAllAtoms() && !selectedAtomsFound; atom++) {
					selectedAtomsFound = mMol.isSelectedAtom(atom);
				}
				boolean selectedObjectsFound = false;
				if (mDrawingObjectList != null) {
					for (int i = 0; i < mDrawingObjectList.size() && !selectedObjectsFound; i++) {
						selectedObjectsFound = mDrawingObjectList.get(i).isSelected();
					}
				}
				double magnification = (Math.abs(mY2 - mY1) < 20) ? 1.0 : Math.exp((mY2 - mY1) / 100);
				double angleChange = (Math.abs(mX2 - mX1) < 20) ? 0.0f : (mX2 - mX1) / 50;
				boolean selectedOnly = (selectedAtomsFound || selectedObjectsFound);
				if (mDrawingObjectList != null && (!selectedOnly || selectedObjectsFound)) {
					for (int i = 0; i < mDrawingObjectList.size(); i++) {
						if (!selectedOnly || mDrawingObjectList.get(i).isSelected()) {
							mDrawingObjectList.get(i).zoomAndRotate(magnification, angleChange);
						}
					}
					update(UPDATE_CHECK_VIEW);
				}
				if (!selectedOnly || selectedAtomsFound) {
					mMol.zoomAndRotate(magnification, angleChange, selectedOnly);
					update(UPDATE_CHECK_VIEW);
				}
				break;
			case cRequestLassoSelect:
			case cRequestSelectRect:
				GenericShape selectedShape = null;
				if (mPendingRequest == cRequestLassoSelect) {
					if ((Math.abs(mX2 - mLassoRegion.getX(mLassoRegion.getSize() - 1))<3)
							&& (Math.abs(mY2 - mLassoRegion.getY(mLassoRegion.getSize() - 1))<3))
						break;

					mLassoRegion.removeLastPoint();
					mLassoRegion.addPoint(mX2, mY2);
					mLassoRegion.addPoint(mX1, mY1);

					selectedShape = mLassoRegion;
				}
				else {
					selectedShape = new GenericRectangle(Math.min(mX1, mX2), Math.min(mY1, mY2), Math.abs(mX2 - mX1), Math.abs(mY2 - mY1));
				}

				if (mDrawingObjectList != null) {
					for (int i = 0; i<mDrawingObjectList.size(); i++) {
						AbstractDrawingObject object = mDrawingObjectList.get(i);
						boolean isSelected = object.isSurroundedBy(selectedShape);
						if ((!mIsAddingToSelection || !mIsSelectedObject[i]) && (isSelected != object.isSelected())) {
							object.setSelected(isSelected);
							mUpdateMode = Math.max(mUpdateMode, UPDATE_REDRAW);
						}
					}
				}
				for (int i = 0; i<mMol.getAllAtoms(); i++) {
					boolean isSelected = selectedShape.contains((int)mMol.getAtomX(i), (int)mMol.getAtomY(i));
					if ((!mIsAddingToSelection || !mIsSelectedAtom[i]) && (isSelected != mMol.isSelectedAtom(i))) {
						mMol.setAtomSelection(i, isSelected);
						mUpdateMode = Math.max(mUpdateMode, UPDATE_REDRAW);
					}
				}
				repaintNeeded = true;
				break;
			case cRequestMapAtoms:
				repaintNeeded = true;
				break;
		}

		if (repaintNeeded) {
			repaint();
		}
	}

	public void mouseMoved(MouseEvent e)
	{
		mMouseIsDown = false;
		int x = e.getX();
		int y = e.getY();

		if (trackHiliting(x, y, false)) {
			repaint();
		}

		updateCursor();
	}

	public void keyPressed(KeyEvent e)
	{
		if (e.getKeyCode() == KeyEvent.VK_SHIFT) {
			mShiftIsDown = true;
			updateCursor();
		}
		if (e.getKeyCode() == KeyEvent.VK_ALT) {
			mAltIsDown = true;
			updateCursor();
		}
		if (e.getKeyCode() == KeyEvent.VK_CONTROL) {
			mControlIsDown = true;
			updateCursor();
		}

		if (mControlIsDown && e.getKeyCode() == KeyEvent.VK_Z) {
			restoreState();
			fireMoleculeChanged();
			update(UPDATE_CHECK_VIEW);
		} else if (e.getKeyCode() == KeyEvent.VK_DELETE) {
			storeState();
			if (mCurrentTool == JDrawToolbar.cToolMapper) {
				boolean found = false;
				for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
					if (mMol.getAtomMapNo(atom) != 0) {
						mMol.setAtomMapNo(atom, 0, false);
						found = true;
					}
				}
				if (found) {
					fireMoleculeChanged();
					update(UPDATE_REDRAW);
				}
			} else if (!deleteHilited()) {
				if (mMol.deleteSelectedAtoms()) {
					fireMoleculeChanged();
					update(UPDATE_REDRAW);
				}
			}
		} else if (e.getKeyCode() == KeyEvent.VK_HELP || (mCurrentHiliteAtom == -1 && e.getKeyChar() == '?')) {
			showHelpDialog();
			return;
		} else if (mCurrentHiliteBond != -1) {
			char ch = e.getKeyChar();
			if (ch == 'q' && mMol.isFragment()) {
				showBondQFDialog(mCurrentHiliteBond);
			} else if (ch == 'v') { // ChemDraw uses the same key
				if (mMol.addRingToBond(mCurrentHiliteBond, 3, false, Molecule.getDefaultAverageBondLength())) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
			} else if (ch >= '4' && ch <= '7') {
				if (mMol.addRingToBond(mCurrentHiliteBond, ch - '0', false, Molecule.getDefaultAverageBondLength())) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
			} else if (ch == 'a' || ch == 'b') {    // ChemDraw uses 'a', we use 'b' since a long time
				if (mMol.addRingToBond(mCurrentHiliteBond, 6, true, Molecule.getDefaultAverageBondLength())) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
			} else {
				boolean bondChanged =
					  (ch == '0') ? changeHighlightedBond(Molecule.cBondTypeMetalLigand)
					: (ch == '1') ? changeHighlightedBond(Molecule.cBondTypeSingle)
					: (ch == '2') ? changeHighlightedBond(Molecule.cBondTypeDouble)
					: (ch == '3') ? changeHighlightedBond(Molecule.cBondTypeTriple)
					: (ch == 'u') ? changeHighlightedBond(Molecule.cBondTypeUp)
					: (ch == 'd') ? changeHighlightedBond(Molecule.cBondTypeDown)
					: (ch == 'c') ? changeHighlightedBond(Molecule.cBondTypeCross)
					: (ch == 'm') ? changeHighlightedBond(Molecule.cBondTypeMetalLigand)
					: false;
				if (bondChanged) {
					fireMoleculeChanged();
					update(UPDATE_REDRAW);
				}
			}
		} else if (mCurrentHiliteAtom != -1) {
			char ch = e.getKeyChar();
			boolean isFirst = (mAtomKeyStrokeBuffer.length() == 0);
			if (isFirst && (ch == '+' || ch == '-')) {
				storeState();
				if (mMol.changeAtomCharge(mCurrentHiliteAtom, ch == '+')) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
			} else if (isFirst && ch == '.') {
				storeState();
				int newRadical = (mMol.getAtomRadical(mCurrentHiliteAtom) == Molecule.cAtomRadicalStateD) ?
					0 : Molecule.cAtomRadicalStateD;
				mMol.setAtomRadical(mCurrentHiliteAtom, newRadical);
				fireMoleculeChanged();
				update(UPDATE_CHECK_COORDS);
			} else if (isFirst && ch == ':') {
				storeState();
				int newRadical = (mMol.getAtomRadical(mCurrentHiliteAtom) == Molecule.cAtomRadicalStateT) ? Molecule.cAtomRadicalStateS
					: (mMol.getAtomRadical(mCurrentHiliteAtom) == Molecule.cAtomRadicalStateS) ? 0 : Molecule.cAtomRadicalStateT;
				mMol.setAtomRadical(mCurrentHiliteAtom, newRadical);
				fireMoleculeChanged();
				update(UPDATE_CHECK_COORDS);
			} else if (isFirst && ch == 'l') {
				mAtomKeyStrokeBuffer.append("Cl");
				update(UPDATE_REDRAW);
			} else if (isFirst && ch == 'q' && mMol.isFragment()) {
				showAtomQFDialog(mCurrentHiliteAtom);
			} else if (isFirst && ch == '?') {
				storeState();
				if (mMol.changeAtom(mCurrentHiliteAtom, 0, 0, -1, 0)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
			} else if (isFirst && ch > 48 && ch <= 57) {
				if (mMol.getFreeValence(mCurrentHiliteAtom) > 0) {
					storeState();
					int chainAtoms = ch - 47;
					int atom1 = mCurrentHiliteAtom;
					int hydrogenCount = mMol.getAllAtoms() - mMol.getAtoms();
					for (int i = 1; i < chainAtoms; i++) {
						suggestNewX2AndY2(atom1);
						int atom2 = mMol.addAtom(mX2, mY2);
						if (atom2 == -1) {
							break;
						}

						mMol.addBond(atom1, atom2);
						atom1 = atom2 - hydrogenCount;	// new atom was added behind all hydrogens and travels now to the front
						mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
					}
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
			} else if (!isFirst && e.getKeyCode() == KeyEvent.VK_ESCAPE) {
				mAtomKeyStrokeBuffer.setLength(0);
				update(UPDATE_REDRAW);
			} else if (!isFirst && e.getKeyCode() == KeyEvent.VK_BACK_SPACE) {
				mAtomKeyStrokeBuffer.setLength(mAtomKeyStrokeBuffer.length() - 1);
				update(UPDATE_REDRAW);
			} else if ((ch >= 65 && ch <= 90)
				|| (ch >= 97 && ch <= 122)
				|| (ch >= 48 && ch <= 57)
				|| (ch == '-')) {
				mAtomKeyStrokeBuffer.append(ch);
				update(UPDATE_REDRAW);
			} else if (ch == '\n' || ch == '\r') {
				expandAtomKeyStrokes(mAtomKeyStrokeBuffer.toString());
			}
		} else if (mCurrentHiliteAtom == -1 && mCurrentHiliteBond == -1) {
			if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE | MODE_MULTIPLE_FRAGMENTS)) == 0) {
				char ch = e.getKeyChar();
				if (ch == 'h') {
					flip(true);
				} if (ch == 'v') {
					flip(false);
				}
			}

		}
	}

	private boolean changeHighlightedBond(int type)
	{
		storeState();
		return mMol.changeBond(mCurrentHiliteBond, type);
	}

	public void showHelpDialog()
	{
		if (mHelpDialog == null || !mHelpDialog.isVisible()) {
			JEditorPane helpPane = new JEditorPane();
			helpPane.setEditorKit(HiDPIHelper.getUIScaleFactor() == 1f ? new HTMLEditorKit() : new ScaledEditorKit());
			helpPane.setEditable(false);
			try {
				helpPane.setPage(getClass().getResource("/html/editor/editor.html"));
			} catch (Exception ex) {
				helpPane.setText(ex.toString());
			}
			Component c = this;
			while (!(c instanceof Frame || c instanceof Dialog)) {
				c = c.getParent();
			}

			if (c instanceof Frame) {
				mHelpDialog = new JDialog((Frame) c, "Idorsia Structure Editor Help", false);
			} else {
				mHelpDialog = new JDialog((Dialog) c, "Idorsia Structure Editor Help", false);
			}

			mHelpDialog.setSize(HiDPIHelper.scale(520), HiDPIHelper.scale(440));
			mHelpDialog.getContentPane().add(new JScrollPane(helpPane,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));
			int x = (c.getX() >= 8 + mHelpDialog.getWidth()) ? c.getX() - 8 - mHelpDialog.getWidth() : c.getX() + 8 + c.getWidth();
			mHelpDialog.setLocation(x, c.getY());
			mHelpDialog.setVisible(true);
		} else {
			Component c = this;
			while (!(c instanceof Frame || c instanceof Dialog)) {
				c = c.getParent();
			}

			int x = (mHelpDialog.getX() + mHelpDialog.getWidth() / 2 >= c.getX() + c.getWidth() / 2) ?
				c.getX() - 8 - mHelpDialog.getWidth() : c.getX() + 8 + c.getWidth();
			mHelpDialog.setLocation(x, c.getY());
		}
	}

	public void keyReleased(KeyEvent e)
	{
		if (e.getKeyCode() == KeyEvent.VK_SHIFT) {
			mShiftIsDown = false;
			updateCursor();
		}
		if (e.getKeyCode() == KeyEvent.VK_ALT) {
			mAltIsDown = false;
			updateCursor();
		}
		if (e.getKeyCode() == KeyEvent.VK_CONTROL) {
			mControlIsDown = false;
			updateCursor();
		}
		if ((e.getModifiers() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()) != 0) {
			if (e.getKeyCode() == KeyEvent.VK_C) {
				copy();
			} else if (e.getKeyCode() == KeyEvent.VK_V) {
				paste();
			}
		}
	}

	public void keyTyped(KeyEvent e)
	{
	}

	private boolean handlePopupTrigger(MouseEvent e)
	{
		if (e.isPopupTrigger()) {
			JPopupMenu popup = null;

			if (mClipboardHandler != null) {
				popup = new JPopupMenu();

				JMenuItem menuItem1 = new JMenuItem(analyseCopy(false) ? ITEM_COPY_REACTION : ITEM_COPY_STRUCTURE);
				menuItem1.addActionListener(this);
				if (mMol.getAllAtoms() == 0)
					menuItem1.setEnabled(false);
				popup.add(menuItem1);

				JMenuItem menuItem2 = null;
				if ((mMode & MODE_REACTION) != 0) {
					menuItem2 = new JMenuItem(ITEM_PASTE_REACTION);
					menuItem2.addActionListener(this);
					popup.add(menuItem2);
				}

				String itemText = (StructureNameResolver.getInstance() == null) ? ITEM_PASTE_STRUCTURE : ITEM_PASTE_WITH_NAME;
				JMenuItem menuItem3 = new JMenuItem(itemText);
				menuItem3.addActionListener(this);
				popup.add(menuItem3);

				JMenuItem menuItem4 = null;
				if ((mMode & MODE_REACTION) != 0) {
					menuItem4 = new JMenuItem(ITEM_LOAD_REACTION);
					menuItem4.addActionListener(this);
					popup.addSeparator();
					popup.add(menuItem4);
				}
			}

			if ((mMode & MODE_REACTION) != 0 && mCurrentTool == JDrawToolbar.cToolMapper) {
				if (popup == null)
					popup = new JPopupMenu();
				else
					popup.addSeparator();

				JMenuItem menuItem1 = new JMenuItem(ITEM_ADD_AUTO_MAPPING);
				menuItem1.addActionListener(this);
				popup.add(menuItem1);

				JMenuItem menuItem2 = new JMenuItem(ITEM_REMOVE_MAPPING);
				menuItem2.addActionListener(this);
				popup.add(menuItem2);
			}

			if (mCurrentTool == JDrawToolbar.cToolZoom) {
				if (popup == null)
					popup = new JPopupMenu();
				else
					popup.addSeparator();

				JMenuItem menuItem1 = new JMenuItem(ITEM_FLIP_HORIZONTALLY);
				menuItem1.addActionListener(this);
				popup.add(menuItem1);
				JMenuItem menuItem2 = new JMenuItem(ITEM_FLIP_VERTICALLY);
				menuItem2.addActionListener(this);
				popup.add(menuItem2);
			}

			if (mAtomColorSupported && mCurrentHiliteAtom != -1) {
				int atomColor = mMol.getAtomColor(mCurrentHiliteAtom);
				if (popup == null) {
					popup = new JPopupMenu();
				} else {
					popup.addSeparator();
				}
				JMenu colorMenu = new JMenu("Set Atom Color");
				addColorToMenu(colorMenu, Color.BLACK, Molecule.cAtomColorNone, atomColor == Molecule.cAtomColorNone);
				addColorToMenu(colorMenu, new Color(AbstractDepictor.COLOR_BLUE), Molecule.cAtomColorBlue, atomColor == Molecule.cAtomColorBlue);
				addColorToMenu(colorMenu, new Color(AbstractDepictor.COLOR_DARK_RED), Molecule.cAtomColorDarkRed, atomColor == Molecule.cAtomColorDarkRed);
				addColorToMenu(colorMenu, new Color(AbstractDepictor.COLOR_RED), Molecule.cAtomColorRed, atomColor == Molecule.cAtomColorRed);
				addColorToMenu(colorMenu, new Color(AbstractDepictor.COLOR_DARK_GREEN), Molecule.cAtomColorDarkGreen, atomColor == Molecule.cAtomColorDarkGreen);
				addColorToMenu(colorMenu, new Color(AbstractDepictor.COLOR_GREEN), Molecule.cAtomColorGreen, atomColor == Molecule.cAtomColorGreen);
				addColorToMenu(colorMenu, new Color(AbstractDepictor.COLOR_MAGENTA), Molecule.cAtomColorMagenta, atomColor == Molecule.cAtomColorMagenta);
				addColorToMenu(colorMenu, new Color(AbstractDepictor.COLOR_ORANGE), Molecule.cAtomColorOrange, atomColor == Molecule.cAtomColorOrange);
				popup.add(colorMenu);
			}

			if (System.getProperty("development") != null) {
				if (popup == null) {
					popup = new JPopupMenu();
				} else {
					popup.addSeparator();
				}
				JMenuItem menuItem1 = new JMenuItem("Show Atom & Bond Numbers");
				menuItem1.addActionListener(ev -> setDisplayMode(AbstractDepictor.cDModeAtomNo | AbstractDepictor.cDModeBondNo));
				popup.add(menuItem1);

				JMenuItem menuItem2 = new JMenuItem("Show Symmetry");
				menuItem2.addActionListener(ev -> setDisplayMode(AbstractDepictor.cDModeShowSymmetrySimple));
				popup.add(menuItem2);

				JMenuItem menuItem3 = new JMenuItem("Show Normal");
				menuItem3.addActionListener(ev -> setDisplayMode(mCurrentTool == JDrawToolbar.cToolMapper ? AbstractDepictor.cDModeShowMapping : 0));
				popup.add(menuItem3);
				}

			if (popup != null) {
				popup.show(this, e.getX(), e.getY());
			}
			return true;
		}
		return false;
	}

	private void addColorToMenu(JMenu menu, Color color, int atomColor, boolean isSelected)
	{
		JRadioButtonMenuItem item = new JRadioButtonMenuItem("	  ", isSelected);
		item.setActionCommand("atomColor" + mCurrentHiliteAtom + ":" + atomColor);
		item.addActionListener(this);
		item.setBackground(color);
		menu.add(item);
	}

	private void handleDoubleClick(MouseEvent e)
	{
		int x = e.getX();
		int y = e.getY();

		int atom = mMol.findAtom(x, y);
		int bond = mMol.findBond(x, y);

		if (mCurrentTool == JDrawToolbar.cToolLassoPointer) {
			if (mMol.isFragment()) {
				if (atom != -1) {
					showAtomQFDialog(atom);
				} else if (bond != -1) {
					showBondQFDialog(bond);
				} else if (mCurrentHiliteObject != null) {
					if (!mShiftIsDown) {
						for (int i = 0; i < mMol.getAllAtoms(); i++) {
							mMol.setAtomSelection(i, false);
						}
						for (int i = 0; i < mDrawingObjectList.size(); i++) {
							((AbstractDrawingObject) mDrawingObjectList.get(i)).setSelected(false);
						}
					}

					mCurrentHiliteObject.setSelected(true);
					update(UPDATE_REDRAW);
				}
			} else {
				int rootAtom = -1;
				if (atom != -1) {
					rootAtom = atom;
				} else if (bond != -1) {
					rootAtom = mMol.getBondAtom(0, bond);
				}

				if (rootAtom != -1 || mCurrentHiliteObject != null) {
					if (!mShiftIsDown) {
						for (int i = 0; i < mMol.getAllAtoms(); i++) {
							mMol.setAtomSelection(i, false);
						}
						if (mDrawingObjectList != null) {
							for (AbstractDrawingObject drawingObject : mDrawingObjectList) {
								drawingObject.setSelected(false);
							}
						}
					}

					if (rootAtom != -1) {
						if ((mMode & MODE_MULTIPLE_FRAGMENTS) != 0) {
							int fragment = mFragmentNo[rootAtom];
							for (int i = 0; i < mMol.getAllAtoms(); i++) {
								if (mFragmentNo[i] == fragment) {
									mMol.setAtomSelection(i, true);
								}
							}
						} else {
							int[] fragmentMember = mMol.getFragmentAtoms(rootAtom);
							for (int i = 0; i < fragmentMember.length; i++) {
								mMol.setAtomSelection(fragmentMember[i], true);
							}
						}
					} else {
						mCurrentHiliteObject.setSelected(true);
					}

					update(UPDATE_REDRAW);
				}
			}
		} else if (mCurrentTool == JDrawToolbar.cToolZoom) {
			int fragment = -2;
			if ((mMode & MODE_MULTIPLE_FRAGMENTS) != 0) {
				fragment = findFragment(x, y);
			}

			if (fragment != -1) {
				double minX = Integer.MAX_VALUE;
				double maxX = Integer.MIN_VALUE;
				for (int i = 0; i < mMol.getAllAtoms(); i++) {
					if (fragment == -2 || mFragmentNo[i] == fragment) {
						if (minX > mMol.getAtomX(i)) {
							minX = mMol.getAtomX(i);
						}
						if (maxX < mMol.getAtomX(i)) {
							maxX = mMol.getAtomX(i);
						}
					}
				}

				if (maxX > minX) {
					double centerX = (maxX + minX) / 2;
					for (int i = 0; i < mMol.getAllAtoms(); i++) {
						if (fragment == -2 || mFragmentNo[i] == fragment) {
							mMol.setAtomX(i, 2 * centerX - mMol.getAtomX(i));
						}
					}
					for (int i = 0; i < mMol.getAllBonds(); i++) {
						if (fragment == -2 || mFragmentNo[mMol.getBondAtom(0, i)] == fragment) {
							switch (mMol.getBondType(i)) {
								case Molecule.cBondTypeUp:
									mMol.setBondType(i, Molecule.cBondTypeDown);
									break;
								case Molecule.cBondTypeDown:
									mMol.setBondType(i, Molecule.cBondTypeUp);
									break;
							}
						}
					}
				}

				fireMoleculeChanged();
				update(UPDATE_REDRAW);
			}
		} else if (mCurrentTool == JDrawToolbar.cToolAtomOther) {
			Component c = this;
			while (c.getParent() != null) {
				c = c.getParent();
			}
			JOptionPane.showMessageDialog((Frame)c, "Please hold 'Ctrl' while pressing the left mouse button\nto open the atom property dialog.");
		}
	}

	private void showAtomQFDialog(int atom)
	{
		if (mAllowQueryFeatures) {
			Component c = this;
			while (!(c instanceof Frame || c instanceof Dialog) && c.getParent() != null) {
				c = c.getParent();
			}
			storeState();
			boolean showReactionHints = ((mMode & MODE_REACTION) != 0);
			if (c instanceof Dialog)
				new JAtomQueryFeatureDialog((Dialog) c, mMol, atom, showReactionHints);
			else
				new JAtomQueryFeatureDialog((Frame) c, mMol, atom, showReactionHints);
			fireMoleculeChanged();
			update(UPDATE_REDRAW);
		}
	}

	private void showBondQFDialog(int bond)
	{
		if (mAllowQueryFeatures) {
			Component c = this;
			while (!(c instanceof Frame || c instanceof Dialog) && c.getParent() != null) {
				c = c.getParent();
			}
			storeState();
			if (c instanceof Dialog)
				new JBondQueryFeatureDialog((Dialog) c, mMol, bond);
			else
				new JBondQueryFeatureDialog((Frame) c, mMol, bond);
			fireMoleculeChanged();
			update(UPDATE_REDRAW);
		}
	}

	private void mousePressedButton1(MouseEvent e)
	{
		mX1 = e.getX();
		mY1 = e.getY();

		switch (mCurrentTool) {
			case JDrawToolbar.cToolZoom:
				// mX1,mY1 define anker for rotation and zooming
				double x = mX1;
				double y = mY1;

				mAtom1 = mMol.findAtom(mX1, mY1);
				if (mAtom1 != -1) {
					mX1 = mMol.getAtomX(mAtom1);
					mY1 = mMol.getAtomY(mAtom1);
				}

				mMol.zoomAndRotateInit(x, y);
				if (mDrawingObjectList != null) {
					for (AbstractDrawingObject drawingObject : mDrawingObjectList) {
						drawingObject.zoomAndRotateInit(x, y);
					}
				}
				storeState();
				mPendingRequest = cRequestZoomAndRotate;
				break;
			case JDrawToolbar.cToolLassoPointer:
				mPendingRequest = cRequestNone;

				// if atom was hit -> move atom (and if atom is selected then all selected stuff)
				mAtom1 = mMol.findAtom(mX1, mY1);
				if (mAtom1 != -1) {
					mAtom2 = -1;
					mX1 = mMol.getAtomX(mAtom1);
					mY1 = mMol.getAtomY(mAtom1);
					if (mMol.isSelectedAtom(mAtom1)) {
						mPendingRequest = mShiftIsDown ? cRequestCopySelected : cRequestMoveSelected;
					} else {
						mPendingRequest = cRequestMoveSingle;
					}
				}

				// if bond was hit -> move bond (and if atom is selected then all selected stuff)
				if (mPendingRequest == cRequestNone) {
					int bondClicked = mMol.findBond(mX1, mY1);
					if (bondClicked != -1) {
						mAtom1 = mMol.getBondAtom(0, bondClicked);
						mAtom2 = mMol.getBondAtom(1, bondClicked);
						if (mMol.isSelectedBond(bondClicked)) {
							mPendingRequest = mShiftIsDown ? cRequestCopySelected : cRequestMoveSelected;
						} else {
							mPendingRequest = cRequestMoveSingle;
						}
					}
				}

				// if object was hit -> move object (and if atom is selected then all selected stuff)
				if (mPendingRequest == cRequestNone) {
					if (mCurrentHiliteObject != null) {
						if (mCurrentHiliteObject.isSelected()) {
							mPendingRequest = mShiftIsDown ? cRequestCopySelected : cRequestMoveSelected;
						} else {
							mPendingRequest = (mShiftIsDown && !(mCurrentHiliteObject instanceof ReactionArrow)) ?
								cRequestCopyObject : cRequestMoveObject;
						}
					}
				}

				if (mPendingRequest != cRequestNone) {
					mX = new double[mMol.getAllAtoms()];
					mY = new double[mMol.getAllAtoms()];
					for (int i = 0; i < mMol.getAllAtoms(); i++) {
						mX[i] = mMol.getAtomX(i);
						mY[i] = mMol.getAtomY(i);
					}
					if (mDrawingObjectList != null) {
						for (AbstractDrawingObject drawingObject : mDrawingObjectList) {
							drawingObject.translateInit(mX1, mY1);
						}
					}

					storeState();
					break;
				}

				// is lasso- or rectangle-selecting
				mIsSelectedAtom = new boolean[mMol.getAllAtoms()];
				if (mDrawingObjectList != null) {
					mIsSelectedObject = new boolean[mDrawingObjectList.size()];
				}

				mIsAddingToSelection = e.isShiftDown();
				for (int i = 0; i < mMol.getAllAtoms(); i++) {
					mIsSelectedAtom[i] = mMol.isSelectedAtom(i);
				}
				if (mDrawingObjectList != null) {
					for (int i = 0; i < mDrawingObjectList.size(); i++) {
						mIsSelectedObject[i] = mDrawingObjectList.get(i).isSelected();
					}
				}

				if (e.isAltDown()) {
					mPendingRequest = cRequestSelectRect;
				} else {
					mLassoRegion = new GenericPolygon();
					mLassoRegion.addPoint(mX1, mY1);
					mLassoRegion.addPoint(mX1, mY1);
					mPendingRequest = cRequestLassoSelect;
				}
				break;
			case JDrawToolbar.cToolDelete:
				storeState();
				deleteAt(mX1, mY1);
				break;
			case JDrawToolbar.cToolUnknownParity:
				int theAtom = mMol.findAtom(mX1, mY1);
				if (theAtom != -1) {
					storeState();
					mMol.setAtomConfigurationUnknown(theAtom, !mMol.isAtomConfigurationUnknown(theAtom));
					fireMoleculeChanged();
					update(UPDATE_REDRAW);
				}
				break;
			case JDrawToolbar.cToolESRAbs:
			case JDrawToolbar.cToolESRAnd:
			case JDrawToolbar.cToolESROr:
				if (mCurrentHiliteBond != -1 && qualifiesForESR(mCurrentHiliteBond)) {
					storeState();
					setESRInfo(mCurrentHiliteBond, (mCurrentTool == JDrawToolbar.cToolESRAbs) ?
						Molecule.cESRTypeAbs :
						(mCurrentTool == JDrawToolbar.cToolESRAnd) ?
							Molecule.cESRTypeAnd : Molecule.cESRTypeOr);
					fireMoleculeChanged();
					update(UPDATE_REDRAW);
				}
				break;
			case JDrawToolbar.cToolStdBond:
			case JDrawToolbar.cToolUpBond:
			case JDrawToolbar.cToolDownBond:
				mAtom1 = mMol.findAtom(mX1, mY1);
				if (mAtom1 == -1) {
					int bond = mMol.findBond(mX1, mY1);
					if (bond != -1) {
						storeState();
						int bondType = Molecule.cBondTypeIncreaseOrder;
						if (mCurrentTool == JDrawToolbar.cToolUpBond) {
							bondType = Molecule.cBondTypeUp;
						} else if (mCurrentTool == JDrawToolbar.cToolDownBond) {
							bondType = Molecule.cBondTypeDown;
						}
						if (mMol.changeBond(bond, bondType)) {
							fireMoleculeChanged();
							update(UPDATE_REDRAW);
						}
						break;
					}
				} else {
					if (mMol.getAllConnAtomsPlusMetalBonds(mAtom1) == MAX_CONNATOMS) {
						return;
					}
					mX1 = mMol.getAtomX(mAtom1);
					mY1 = mMol.getAtomY(mAtom1);
				}
				mPendingRequest = cRequestNewBond;
				suggestNewX2AndY2(mAtom1);
				repaint();
				break;
			case JDrawToolbar.cToolChain:
				mAtom1 = mMol.findAtom(mX1, mY1);
				if (mAtom1 != -1) {
					if (mMol.getAllConnAtomsPlusMetalBonds(mAtom1) == MAX_CONNATOMS) {
						return;
					}
					mX1 = mMol.getAtomX(mAtom1);
					mY1 = mMol.getAtomY(mAtom1);
				}
				mPendingRequest = cRequestNewChain;
				mChainAtoms = 0;
				break;
			case JDrawToolbar.cTool3Ring:
				storeState();
				if (mMol.addRing(mX1, mY1, 3, false, Molecule.getDefaultAverageBondLength())) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cTool4Ring:
				storeState();
				if (mMol.addRing(mX1, mY1, 4, false, Molecule.getDefaultAverageBondLength())) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cTool5Ring:
				storeState();
				if (mMol.addRing(mX1, mY1, 5, false, Molecule.getDefaultAverageBondLength())) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cTool6Ring:
				storeState();
				if (mMol.addRing(mX1, mY1, 6, false, Molecule.getDefaultAverageBondLength())) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cTool7Ring:
				storeState();
				if (mMol.addRing(mX1, mY1, 7, false, Molecule.getDefaultAverageBondLength())) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAromRing:
				storeState();
				if (mMol.addRing(mX1, mY1, 6, true, Molecule.getDefaultAverageBondLength())) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolPosCharge:
				storeState();
				if (mMol.changeAtomCharge(mX1, mY1, true)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolNegCharge:
				storeState();
				if (mMol.changeAtomCharge(mX1, mY1, false)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomH:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 1, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomC:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 6, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomN:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 7, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomO:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 8, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomSi:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 14, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomP:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 15, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomS:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 16, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomF:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 9, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomCl:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 17, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomBr:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 35, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomI:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 53, 0, -1, 0, null)) {
					fireMoleculeChanged();
					update(UPDATE_CHECK_COORDS);
				}
				break;
			case JDrawToolbar.cToolAtomOther:
				if (mOtherAtom == -1 || e.isControlDown()) {
					int atom = mMol.findAtom(mX1, mY1);
					if (atom != -1) {
						Component c = this;
						while (c.getParent() != null) {
							c = c.getParent();
						}
						storeState();
						new JAtomLabelDialog((Frame) c, mMol, atom);
						mOtherAtom = mMol.getAtomicNo(atom);
						mOtherMass = mMol.getAtomMass(atom);
						mOtherValence = mMol.getAtomAbnormalValence(atom);
						mOtherRadical = mMol.getAtomRadical(atom);
						mOtherLabel = mMol.getAtomCustomLabel(atom);
						fireMoleculeChanged();
						update(UPDATE_REDRAW);
					}
				}
				else {
					storeState();
					if (mMol.addOrChangeAtom(mX1, mY1, mOtherAtom, mOtherMass, mOtherValence, mOtherRadical, mOtherLabel)) {
						fireMoleculeChanged();
						update(UPDATE_CHECK_COORDS);
					}
				}
				break;
			case JDrawToolbar.cToolMapper:
				mAtom1 = mMol.findAtom(mX1, mY1);
				if (mAtom1 != -1 && mAtom1 < mMol.getAtoms()) {
					mX1 = mMol.getAtomX(mAtom1);
					mY1 = mMol.getAtomY(mAtom1);
					mPendingRequest = cRequestMapAtoms;
				}
				break;
			case JDrawToolbar.cToolText:
				TextDrawingObject object = null;
				if (mCurrentHiliteObject == null) {
					object = new TextDrawingObject();
					object.setCoordinates(mX1, mY1);
					mDrawingObjectList.add(object);
				} else if (mCurrentHiliteObject instanceof TextDrawingObject) {
					object = (TextDrawingObject) mCurrentHiliteObject;
				}
				editTextObject(object);
				storeState();
				update(UPDATE_CHECK_COORDS);
				break;
		}
	}

	private void mouseReleasedButton1(MouseEvent e)
	{
		int pendingRequest = mPendingRequest;
		mPendingRequest = cRequestNone;
		switch (pendingRequest) {
			case cRequestNewBond:
				int stopAtom;

				stopAtom = mMol.findAtom(mX2, mY2);
				if (stopAtom != -1
					&& mMol.getAllConnAtomsPlusMetalBonds(stopAtom) == MAX_CONNATOMS) {
					return;
				}

				storeState();

				if (mAtom1 == -1) {
					mAtom1 = mMol.addAtom(mX1, mY1);
				}
				if (stopAtom == -1) {
					stopAtom = mMol.addAtom(mX2, mY2);
				}

				if ((mAtom1 != -1 && stopAtom != -1)) {
					int bondType = (mMol.isMetalAtom(mAtom1) || mMol.isMetalAtom(stopAtom)) ?
							Molecule.cBondTypeMetalLigand : Molecule.cBondTypeSingle;
					if (mCurrentTool == JDrawToolbar.cToolUpBond) {
						bondType = Molecule.cBondTypeUp;
					}
					if (mCurrentTool == JDrawToolbar.cToolDownBond) {
						bondType = Molecule.cBondTypeDown;
					}
					mMol.addOrChangeBond(mAtom1, stopAtom, bondType);
				}

				fireMoleculeChanged();
				update(UPDATE_CHECK_COORDS);
				break;
			case cRequestNewChain:
				storeState();
				if (mChainAtoms > 0) {
					if (mAtom1 == -1) {
						mAtom1 = mMol.addAtom(mX1, mY1);
					}

					if (mChainAtom[0] == -1) {
						mChainAtom[0] = mMol.addAtom(mChainAtomX[0],
							mChainAtomY[0]);
					}

					if (mChainAtom[0] != -1) {
						mMol.addBond(mAtom1, mChainAtom[0]);
					}
				}

				if (mChainAtoms > 1) {
					for (int i = 1; i < mChainAtoms; i++) {
						if (mChainAtom[i] == -1) {
							mChainAtom[i] = mMol.addAtom(mChainAtomX[i],
								mChainAtomY[i]);
						}
						if (mChainAtom[i] != -1) {
							mMol.addBond(mChainAtom[i - 1], mChainAtom[i]);
						}
					}
				}
				fireMoleculeChanged();
				update(UPDATE_CHECK_COORDS);
				break;
			case cRequestMoveSingle:
			case cRequestMoveSelected:
			case cRequestZoomAndRotate:
				fireMoleculeChanged();
			case cRequestMoveObject:
				update(UPDATE_CHECK_COORDS);
				break;
			case cRequestLassoSelect:
			case cRequestSelectRect:
				boolean selectionChanged = false;
				for (int i = 0; i < mMol.getAllAtoms(); i++) {
					if (mIsSelectedAtom[i] != mMol.isSelectedAtom(i)) {
						selectionChanged = true;
						break;
					}
				}
				if (selectionChanged) {
					fireEvent(new EditorEvent(this, EditorEvent.WHAT_SELECTION_CHANGED, true));
				}
				repaint();
				break;
			case cRequestMapAtoms:
				boolean mapNoChanged = false;
				int atom2 = mCurrentHiliteAtom;
//				System.out.printf("Map Request Atom %d => %d (%d)\n", mAtom1, mAtom2, atom2);
				int mapNoAtom1 = mMol.getAtomMapNo(mAtom1);
				if (atom2 == -1) {
					storeState();
					if (mapNoAtom1 != 0) {
						mapNoChanged = true;
						for (int atom = 0; atom < mMol.getAtoms(); atom++) {
							if (mMol.getAtomMapNo(atom) == mapNoAtom1) {
								mMol.setAtomMapNo(atom, 0, false);
							}
						}

						autoMapReaction();
					}
				} else {
					storeState();
					mapNoChanged = true;

					// If we only clicked on a mapped atom, then reset the map number
					if (mAtom1 == atom2) {
						int mapNo = mMol.getAtomMapNo(mAtom1);
						for (int atom = 0; atom < mMol.getAtoms(); atom++) {
							if (mMol.getAtomMapNo(atom) == mapNo) {
								mMol.setAtomMapNo(atom, 0, false);
							}
						}
					} else {
						// remove old mapping numbers of atom1 and atom2
						int mapNoAtom2 = mMol.getAtomMapNo(atom2);
						for (int atom = 0; atom < mMol.getAtoms(); atom++) {
							if (mMol.getAtomMapNo(atom) == mapNoAtom1
								|| mMol.getAtomMapNo(atom) == mapNoAtom2) {
								mMol.setAtomMapNo(atom, 0, false);
							}
						}
						int freeMapNo = 1;
						for (int atom = 0; atom < mMol.getAtoms(); atom++) {
							if (mMol.getAtomMapNo(atom) == freeMapNo) {
								freeMapNo++;
								atom = -1;
							}
						}
						mMol.setAtomMapNo(mAtom1, freeMapNo, false);
						mMol.setAtomMapNo(atom2, freeMapNo, false);
					}

					autoMapReaction();
				}

				if (mapNoChanged) {
					fireMoleculeChanged();
					mUpdateMode = Math.max(mUpdateMode, UPDATE_REDRAW);
				}

				repaint();
				break;
		}
	}

	/**
	 * Takes the manually mapped atom mapping numbers from the display molecule,
	 * copies them into the current fragments, creates a reaction from these,
	 * uses the MCS-mapper to map the reaction and returns the rxn's mapping
	 * into the display molecule.
	 */
	private void autoMapReaction()
	{
		if (sMapper == null)
			sMapper = new MCSReactionMapper();

//		if (sMapper == null) {
//			new MoleculeAutoMapper(mMol).autoMap();
//			return;
//		}

		// We assume that we use the MCS-mapper, which doesn't care about manually mapped atoms.
		// Thus, we need a hack to ensure that manually mapped atoms are reliably part of
		// the MCS-mapper's result. Therefore we give every manually mapped atom pair a
		// unique fake atom mass. Original atom masses are copied and later restored.
		// We also provide an updated SSSearcher(), which requires atom mass equivalence for
		// atoms to match.

		SSSearcher sss = new SSSearcher() {
			@Override
			public boolean areAtomsSimilar(int moleculeAtom, int fragmentAtom) {
				if (mMolecule.getAtomicNo(moleculeAtom) == mFragment.getAtomicNo(fragmentAtom)) {
					if (mMolecule.getAtomMass(moleculeAtom) != mFragment.getAtomMass(fragmentAtom))
						return false;

					if (mMolecule.isAromaticAtom(moleculeAtom) || mFragment.isAromaticAtom(fragmentAtom))
						return true;
					}

				return super.areAtomsSimilar(moleculeAtom, fragmentAtom);
				}

			@Override
			public boolean areBondsSimilar(int moleculeBond, int fragmentBond) {
				if (mMolecule.isAromaticBond(moleculeBond)
				 || mMolecule.isDelocalizedBond(moleculeBond)
				 || mFragment.isAromaticBond(fragmentBond)
				 || mFragment.isDelocalizedBond(fragmentBond))
					return true;

				return super.areBondsSimilar(moleculeBond, fragmentBond);
			}
		};

		Reaction rxn = getReaction();

		// manual mapNos a put as negative keys!!!
		TreeMap<Integer,Integer> oldToNewMapNo = new TreeMap<>();
		int nextMapNo = 1;

		final int fakeAtomMassBase = 512;

		// Mark the manually mapped atoms such that the mapper uses them first priority and
		// to be able to re-assign them later as manually mapped.
		int[] fragmentAtom = new int[mFragment.length];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			int fragment = mFragmentNo[atom];
			mFragment[fragment].setAtomMapNo(fragmentAtom[fragment], 0, false);
			if (mMol.getAtomMapNo(atom) != 0 && !mMol.isAutoMappedAtom(atom)) {
				int manualMapNo = mMol.getAtomMapNo(atom);

				// make sure that negative manual mapNo is in TreeMap
				Integer newMapNo = oldToNewMapNo.get(-manualMapNo);
				if (newMapNo == null)
					oldToNewMapNo.put(-manualMapNo, newMapNo = new Integer(nextMapNo++));

				mFragment[fragment].setAtomMass(fragmentAtom[fragment], fakeAtomMassBase + newMapNo);
			}
			fragmentAtom[fragment]++;
		}

		rxn = sMapper.mapReaction(rxn, sss);

		if (rxn != null) {
			// assign new mapping numbers: manually mapped atom starting from 1. Automapped atoms follow.

			// write manually mapped atoms' mapping number from rxn (i.e. fragments) into display molecule
			fragmentAtom = new int[mFragment.length];
			for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
				int fragment = mFragmentNo[atom];
				boolean hasFakeAtomMass = (mFragment[fragment].getAtomMass(fragmentAtom[fragment]) > fakeAtomMassBase);
				if (hasFakeAtomMass) {
					// rescue new mapNo
					int newMapNo = mFragment[fragment].getAtomMass(fragmentAtom[fragment]) - fakeAtomMassBase;

					// repair fake atom mass
					mFragment[fragment].setAtomMass(fragmentAtom[fragment], mMol.getAtomMass(atom));

					mMol.setAtomMapNo(atom, newMapNo, false);
					mFragment[fragment].setAtomMapNo(fragmentAtom[fragment], newMapNo, false);
					}
				else {
					// take generated mapNo from reaction
					int generatedMapNo = mFragment[fragment].getAtomMapNo(fragmentAtom[fragment]);

					Integer newMapNo = 0;
					if (generatedMapNo != 0) {
						newMapNo = oldToNewMapNo.get(generatedMapNo);
						if (newMapNo == null)
							oldToNewMapNo.put(generatedMapNo, newMapNo = new Integer(nextMapNo++));
					}

					mMol.setAtomMapNo(atom, newMapNo, true);
					mFragment[fragment].setAtomMapNo(fragmentAtom[fragment], newMapNo, true);
				}
				fragmentAtom[fragment]++;
			}
		}
		else {
			// restore original atom masses in fragments and copy molecule's mapping number into fragments
			fragmentAtom = new int[mFragment.length];
			for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
				int fragment = mFragmentNo[atom];
				mFragment[fragment].setAtomMass(fragmentAtom[fragment], mMol.getAtomMass(atom));
				mFragment[fragment].setAtomMapNo(fragmentAtom[fragment], mMol.getAtomMapNo(atom), mMol.isAutoMappedAtom(atom));
				fragmentAtom[fragment]++;
			}
		}




//		mMapper.resetFragments(mappedReaction);

//		AStarReactionMapper m = new AStarReactionMapper();
//		Reaction rxn = getReaction();
//		List<AStarReactionMapper.SlimMapping> l = m.map(rxn);
//		if (l != null && l.size() > 0) {
//			AStarReactionMapper.SlimMapping sm = l.get(0);
////			int[] mps = sm.getMapping();
////			for (int i = 0; i < mps.length; i++) {
////				System.out.printf("Maps %d -> %d\n", i, mps[i]);
////			}
//			m.activateMapping(sm);
////			for (int i = 0; i < mFragment.length; i++) {
////				StereoMolecule mol = mFragment[i];
////				for (int a = 0; a < mol.getAllAtoms(); a++) {
////					System.out.printf("T Map %d = %d\n", a, mol.getAtomMapNo(a));
////				}
////			}
//
//			int[] fragmentAtom = new int[mFragment.length];
//			for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
//				int fragment = mFragmentNo[atom];
//				if (mMol.getAtomMapNo(atom) == 0)
//					mMol.setAtomMapNo(atom, mFragment[fragment].getAtomMapNo(fragmentAtom[fragment]), true);
//				fragmentAtom[fragment]++;
//			}
//		}
	}

	/**
	 * Checks whether this bond is a stereo bond and whether it refers to a
	 * stereo center or BINAP bond, making it eligible for ESR information.
	 *
	 * @param stereoBond the up/down stereo bond
	 * @return
	 */
	private boolean qualifiesForESR(int stereoBond)
	{
		return mMol.isStereoBond(stereoBond) && (getESRAtom(stereoBond) != -1 || getESRBond(stereoBond) != -1);
	}

	/**
	 * Locates the stereo center with parity 1 or 2 that is defined by the stereo bond.
	 *
	 * @param stereoBond
	 * @return stereo center atom or -1 if no stereo center found
	 */
	private int getESRAtom(int stereoBond)
	{
		int atom = mMol.getBondAtom(0, stereoBond);
		if (mMol.getAtomParity(atom) != Molecule.cAtomParityNone) {
			return (mMol.isAtomParityPseudo(atom)
				|| (mMol.getAtomParity(atom) != Molecule.cAtomParity1
				&& mMol.getAtomParity(atom) != Molecule.cAtomParity2)) ? -1 : atom;
		}
		if (mMol.getAtomPi(atom) == 1) {
			for (int i = 0; i < mMol.getConnAtoms(atom); i++) {
				if (mMol.getConnBondOrder(atom, i) == 2) {
					int connAtom = mMol.getConnAtom(atom, i);
					if (mMol.getAtomPi(connAtom) == 2
						&& (mMol.getAtomParity(connAtom) == Molecule.cAtomParity1
						|| mMol.getAtomParity(connAtom) == Molecule.cAtomParity2)) {
						return connAtom;
					}
				}
			}
		}
		return -1;
	}

	private int getESRBond(int stereoBond)
	{
		int bond = mMol.findBINAPChiralityBond(mMol.getBondAtom(0, stereoBond));
		if (bond != -1
			&& mMol.getBondParity(bond) != Molecule.cBondParityEor1
			&& mMol.getBondParity(bond) != Molecule.cBondParityZor2) {
			bond = -1;
		}

		return bond;
	}

	/**
	 * Sets the ESR information of an atom or bond.
	 *
	 */
	private void setESRInfo(int stereoBond, int type)
	{
		int group = -1;

		int atom = getESRAtom(stereoBond);
		int bond = (atom == -1) ? getESRBond(stereoBond) : -1;

		// if type requires group information (type==And or type==Or)
		if (type != Molecule.cESRTypeAbs) {
			int maxGroup = -1;
			for (int i = 0; i < mMol.getAtoms(); i++) {
				if (i != atom
					&& mMol.getAtomESRType(i) == type
					&& (!mMol.isSelectedBond(stereoBond) || !mMol.isSelectedAtom(i))) {
					int grp = mMol.getAtomESRGroup(i);
					if (maxGroup < grp) {
						maxGroup = grp;
					}
				}
			}
			for (int i = 0; i < mMol.getBonds(); i++) {
				if (i != bond
					&& mMol.getBondESRType(i) == type
					&& (!mMol.isSelectedBond(stereoBond) || !mMol.isSelectedBond(i))) {
					int grp = mMol.getBondESRGroup(i);
					if (maxGroup < grp) {
						maxGroup = grp;
					}
				}
			}

			if ((atom == -1 ? mMol.getBondESRType(bond) : mMol.getAtomESRType(atom)) != type) {
				group = Math.min(maxGroup + 1, Molecule.cESRMaxGroups - 1);
			} else {
				group = (atom == -1) ? mMol.getBondESRGroup(bond) : mMol.getAtomESRGroup(atom);
				if (mMol.isSelectedBond(stereoBond)) {
					boolean selectedShareOneGroup = true;
					for (int i = 0; i < mMol.getAtoms(); i++) {
						if (i != atom && mMol.isSelectedAtom(i)
							&& mMol.getAtomESRType(i) == type
							&& mMol.getAtomESRGroup(i) != group) {
							selectedShareOneGroup = false;
							break;
						}
					}
					for (int i = 0; i < mMol.getBonds(); i++) {
						if (i != bond && mMol.isSelectedBond(i)
							&& mMol.getBondESRType(i) == type
							&& mMol.getBondESRGroup(i) != group) {
							selectedShareOneGroup = false;
							break;
						}
					}
					if (selectedShareOneGroup) {
						if (group <= maxGroup) {
							group++;
							if (group == Molecule.cESRMaxGroups) {
								group = 0;
							}
						} else {
							group = 0;
						}
					}
				} else {
					if (group <= maxGroup) {
						group++;
						if (group == Molecule.cESRMaxGroups) {
							group = 0;
						}
					} else {
						group = 0;
					}
				}
			}
		}

		if (mMol.isSelectedBond(stereoBond)) {
			for (int i = 0; i < mMol.getBonds(); i++) {
				if (mMol.isSelectedBond(i) && mMol.isStereoBond(i)) {
					int a = getESRAtom(i);
					int b = getESRBond(i);
					if (a != -1) {
						mMol.setAtomESR(a, type, group);
					} else if (b != -1) {
						mMol.setBondESR(b, type, group);
					}
				}
			}
		} else {
			if (atom != -1) {
				mMol.setAtomESR(atom, type, group);
			} else if (bond != -1) {
				mMol.setBondESR(bond, type, group);
			}
		}
	}

	private int findFragment(double x, double y)
	{
		int fragment = -1;
		double minDistance = Double.MAX_VALUE;
		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
			double dx = mX1 - mMol.getAtomX(atom);
			double dy = mY1 - mMol.getAtomY(atom);
			double distance = Math.sqrt(dx * dx + dy * dy);
			if (distance < FRAGMENT_MAX_CLICK_DISTANCE
				&& minDistance > distance) {
				minDistance = distance;
				fragment = mFragmentNo[atom];
			}
		}

		return fragment;
	}

	/**
	 * requires cHelperNeighbours
	 *
	 * @param atom
	 */
	private void suggestNewX2AndY2(int atom)
	{
		double newAngle = Math.PI * 2 / 3;
		if (atom != -1) {
			double angle[] = new double[MAX_CONNATOMS + 1];
			for (int i = 0; i < mMol.getAllConnAtomsPlusMetalBonds(atom); i++) {
				angle[i] = mMol.getBondAngle(atom, mMol.getConnAtom(atom, i));
			}

			if (mMol.getAllConnAtomsPlusMetalBonds(atom) == 1) {
				if (angle[0] < -Math.PI * 5 / 6) {
					newAngle = Math.PI / 3;
				} else if (angle[0] < -Math.PI / 2) {
					newAngle = Math.PI * 2 / 3;
				} else if (angle[0] < -Math.PI / 6) {
					newAngle = Math.PI / 3;
				} else if (angle[0] < 0.0) {
					newAngle = Math.PI * 2 / 3;
				} else if (angle[0] < Math.PI / 6) {
					newAngle = -Math.PI * 2 / 3;
				} else if (angle[0] < Math.PI / 2) {
					newAngle = -Math.PI / 3;
				} else if (angle[0] < Math.PI * 5 / 6) {
					newAngle = -Math.PI * 2 / 3;
				} else {
					newAngle = -Math.PI / 3;
				}
			} else {
				for (int i = mMol.getAllConnAtomsPlusMetalBonds(atom) - 1; i > 0; i--) {	// bubble sort
					for (int j = 0; j < i; j++) {
						if (angle[j] > angle[j + 1]) {
							double temp = angle[j];
							angle[j] = angle[j + 1];
							angle[j + 1] = temp;
						}
					}
				}
				angle[mMol.getAllConnAtomsPlusMetalBonds(atom)] = angle[0] + Math.PI * 2;

				int largestNo = 0;
				double largestDiff = 0.0;
				for (int i = 0; i < mMol.getAllConnAtomsPlusMetalBonds(atom); i++) {
					double angleDiff = angle[i + 1] - angle[i];
					if (largestDiff < angleDiff) {
						largestDiff = angleDiff;
						largestNo = i;
					}
				}
				newAngle = (angle[largestNo] + angle[largestNo + 1]) / 2;
			}
		}
		double avbl = mMol.getAverageBondLength();
		mX2 = ((atom == -1) ? mX1 : mMol.getAtomX(atom)) + avbl * (float) Math.sin(newAngle);
		mY2 = ((atom == -1) ? mY1 : mMol.getAtomY(atom)) + avbl * (float) Math.cos(newAngle);
	}

	private boolean areAtomsMappingCompatible(int atom1, int atom2)
	{
		if (mMol.isFragment()) {
			if ((mMol.getAtomQueryFeatures(atom1) & Molecule.cAtomQFExcludeGroup) != 0
			 || (mMol.getAtomQueryFeatures(atom1) & Molecule.cAtomQFExcludeGroup) != 0)
				return false;

			int[] atomList1 = mMol.getAtomList(atom1);
			int[] atomList2 = mMol.getAtomList(atom2);
			if (atomList1 == null ^ atomList2 == null) {
				return false;
			}

			if (atomList1 != null) {
				if (atomList1.length != atomList2.length) {
					return false;
				}
				for (int i = 0; i < atomList1.length; i++) {
					if (atomList1[i] != atomList2[i]) {
						return false;
					}
				}
			}

			boolean isAny1 = ((mMol.getAtomQueryFeatures(atom1) & Molecule.cAtomQFAny) != 0);
			boolean isAny2 = ((mMol.getAtomQueryFeatures(atom2) & Molecule.cAtomQFAny) != 0);
			if (isAny1 != isAny2) {
				return false;
			}
		}

		return mMol.getAtomicNo(atom1) == mMol.getAtomicNo(atom2);
	}

	private boolean trackHiliting(double x, double y, boolean isDragging)
	{
		int theAtom = mMol.findAtom(x, y);
		int theBond = -1;

		if (isDragging
			&& mPendingRequest == cRequestMapAtoms
			&& theAtom != -1
			&& (!areAtomsMappingCompatible(mAtom1, theAtom)
			|| (mMol.getAtomMapNo(mAtom1) != 0 && mMol.getAtomMapNo(mAtom1) == mMol.getAtomMapNo(theAtom) && !mMol.isAutoMappedAtom(mAtom1))
			|| shareSameReactionSide(mAtom1, theAtom))) {
			theAtom = -1;
		}

		if (theAtom != -1) {
			if (mCurrentTool == JDrawToolbar.cToolESRAbs
				|| mCurrentTool == JDrawToolbar.cToolESRAnd
				|| mCurrentTool == JDrawToolbar.cToolESROr) {
				theBond = mMol.getStereoBond(theAtom);
				theAtom = -1;
			} else if (mCurrentTool == JDrawToolbar.cToolMapper
				&& theAtom >= mMol.getAtoms()) {
				theAtom = -1;
			}
		}

		if (theBond == -1
			&& theAtom == -1
			&& mCurrentTool != JDrawToolbar.cToolChain
			&& mCurrentTool != JDrawToolbar.cToolMapper
			&& mCurrentTool != JDrawToolbar.cToolUnknownParity
			&& mCurrentTool != JDrawToolbar.cToolPosCharge
			&& mCurrentTool != JDrawToolbar.cToolNegCharge
			&& mCurrentTool != JDrawToolbar.cToolAtomH
			&& mCurrentTool != JDrawToolbar.cToolAtomC
			&& mCurrentTool != JDrawToolbar.cToolAtomN
			&& mCurrentTool != JDrawToolbar.cToolAtomO
			&& mCurrentTool != JDrawToolbar.cToolAtomSi
			&& mCurrentTool != JDrawToolbar.cToolAtomP
			&& mCurrentTool != JDrawToolbar.cToolAtomS
			&& mCurrentTool != JDrawToolbar.cToolAtomF
			&& mCurrentTool != JDrawToolbar.cToolAtomCl
			&& mCurrentTool != JDrawToolbar.cToolAtomBr
			&& mCurrentTool != JDrawToolbar.cToolAtomI
			&& mCurrentTool != JDrawToolbar.cToolAtomOther) {
			theBond = mMol.findBond(x, y);
		}

		if (theBond != -1
			&& (mCurrentTool == JDrawToolbar.cToolESRAbs
			|| mCurrentTool == JDrawToolbar.cToolESRAnd
			|| mCurrentTool == JDrawToolbar.cToolESROr)
			&& !qualifiesForESR(theBond)) {
			theBond = -1;
		}

		// don't change object hiliting while dragging
		AbstractDrawingObject hiliteObject = mCurrentHiliteObject;
		if (!isDragging && mDrawingObjectList != null) {
			hiliteObject = null;
			if (theAtom == -1 && theBond == -1
				&& (mCurrentTool == JDrawToolbar.cToolLassoPointer
				|| mCurrentTool == JDrawToolbar.cToolDelete
				|| mCurrentTool == JDrawToolbar.cToolText)) {
				for (AbstractDrawingObject theObject : mDrawingObjectList) {
					if (mCurrentTool == JDrawToolbar.cToolLassoPointer
						|| (mCurrentTool == JDrawToolbar.cToolDelete && !(theObject instanceof ReactionArrow))
						|| (mCurrentTool == JDrawToolbar.cToolText && theObject instanceof TextDrawingObject)) {
						if (theObject.checkHiliting(x, y)) {
							hiliteObject = theObject;
							if (mCurrentHiliteObject != null && mCurrentHiliteObject != theObject) {
								mCurrentHiliteObject.clearHiliting();
							}
							break;
						}
					}
				}
			}
		}

		boolean repaintNeeded = (mCurrentHiliteAtom != theAtom
			|| mCurrentHiliteBond != theBond
			|| mCurrentHiliteObject != hiliteObject
			|| hiliteObject != null);

		if (mCurrentHiliteAtom != theAtom) {
			if (mCurrentHiliteAtom != -1 && mAtomKeyStrokeBuffer.length() != 0)
				expandAtomKeyStrokes(mAtomKeyStrokeBuffer.toString());

			mCurrentHiliteAtom = theAtom;
			mAtomKeyStrokeBuffer.setLength(0);
			fireEvent(new EditorEvent(this, EditorEvent.WHAT_HILITE_ATOM_CHANGED, true));
		}
		if (mCurrentHiliteBond != theBond) {
			mCurrentHiliteBond = theBond;
			fireEvent(new EditorEvent(this, EditorEvent.WHAT_HILITE_BOND_CHANGED, true));
		}
		mCurrentHiliteObject = hiliteObject;

		return repaintNeeded;
	}

	private int getAtomKeyStrokeValidity(String s)
	{
		if (Molecule.getAtomicNoFromLabel(s) != 0)
			return KEY_IS_ATOM_LABEL;
		if (NamedSubstituents.getSubstituentIDCode(s) != null)
			return KEY_IS_SUBSTITUENT;
		if (isValidAtomKeyStrokeStart(s))
			return KEY_IS_VALID_START;
		return KEY_IS_INVALID;
	}

	/**
	 * @param s
	 * @return true if s is either a valid atom symbol or a valid substituent name
	 */
	private boolean isValidAtomKeyStroke(String s)
	{
		return Molecule.getAtomicNoFromLabel(s) != 0
			|| NamedSubstituents.getSubstituentIDCode(s) != null;
	}

	/**
	 * @param s
	 * @return true if adding one or more chars may still create a valid key stroke sequence
	 */
	private boolean isValidAtomKeyStrokeStart(String s)
	{
		if (s.length() < 3)
			for (int i = 1; i < Molecule.cAtomLabel.length; i++)
				if (Molecule.cAtomLabel[i].startsWith(s))
					return true;

		return NamedSubstituents.isValidSubstituentNameStart(s);
	}

	private void expandAtomKeyStrokes(String keyStrokes)
	{
		mAtomKeyStrokeBuffer.setLength(0);

		int atomicNo = Molecule.getAtomicNoFromLabel(keyStrokes);
		if (atomicNo != 0) {
			storeState();
			if (mMol.changeAtom(mCurrentHiliteAtom, atomicNo, 0, -1, 0)) {
				fireMoleculeChanged();
				update(UPDATE_CHECK_COORDS);
				return;
			}
		}

		StereoMolecule substituent = NamedSubstituents.getSubstituent(keyStrokes);
		if (substituent != null) {
			storeState();

			// Copy the the fragment containing the attachment point into a new molecule.
			// Then attach the substituent, create new atom coordinates for the substituent,
			// while retaining coordinates of the fragment.
			StereoMolecule fragment = new StereoMolecule();
			fragment.addFragment(mMol, mCurrentHiliteAtom, null);
			double sourceAVBL = fragment.getAverageBondLength();
			int firstAtomInFragment = fragment.getAllAtoms();
			for (int atom = 0; atom < fragment.getAllAtoms(); atom++)
				fragment.setAtomMarker(atom, true);
			fragment.addSubstituent(substituent, 0);
			new CoordinateInventor(CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS).invent(fragment);

			double dx = mMol.getAtomX(mCurrentHiliteAtom) - sourceAVBL * fragment.getAtomX(0);
			double dy = mMol.getAtomY(mCurrentHiliteAtom) - sourceAVBL * fragment.getAtomY(0);

			// Attach the substituent to the complete molecule and take coodinates from the
			// previously created fragment-substituent species.
			int firstAtomInMol = mMol.getAllAtoms();
			mMol.addSubstituent(substituent, mCurrentHiliteAtom);
			int substituentAtoms = mMol.getAllAtoms() - firstAtomInMol;
			for (int i = 0; i < substituentAtoms; i++) {
				mMol.setAtomX(firstAtomInMol + i, sourceAVBL * fragment.getAtomX(firstAtomInFragment + i) + dx);
				mMol.setAtomY(firstAtomInMol + i, sourceAVBL * fragment.getAtomY(firstAtomInFragment + i) + dy);
			}
			mMol.setStereoBondsFromParity();

			fireMoleculeChanged();
			update(UPDATE_CHECK_COORDS);
		}
	}

	private AbstractDrawingObject findDrawingObject(double x, double y, String type, boolean forDeletion)
	{
		if (mDrawingObjectList != null) {
			for (AbstractDrawingObject drawingObject : mDrawingObjectList) {
				if ((type == null || type.equals(drawingObject.getTypeString())
					&& !forDeletion || drawingObject.isDeletable())
					&& drawingObject.contains(x, y)) {
					return drawingObject;
				}
			}
		}

		return null;
	}

	private void editTextObject(TextDrawingObject object)
	{
		Component c = this;
		while (c.getParent() != null) {
			c = c.getParent();
		}
		new JTextDrawingObjectDialog((Frame) c, object);

		boolean nonWhiteSpaceFound = false;
		for (int i = 0; i < object.getText().length(); i++) {
			if (!Character.isWhitespace(object.getText().charAt(i))) {
				nonWhiteSpaceFound = true;
				break;
			}
		}
		if (!nonWhiteSpaceFound) {
			mDrawingObjectList.remove(object);
		}

		repaint();
	}

	private boolean shareSameReactionSide(int atom1, int atom2)
	{
		ReactionArrow arrow = (ReactionArrow) mDrawingObjectList.get(0);
		return !(arrow.isOnProductSide(mMol.getAtomX(atom1), mMol.getAtomY(atom1))
			^ arrow.isOnProductSide(mMol.getAtomX(atom2), mMol.getAtomY(atom2)));
	}

	protected void restoreState()
	{
		if (mUndoMol == null) {
			return;
		}
		mUndoMol.copyMolecule(mMol);
		mDrawingObjectList = (mUndoDrawingObjectList == null) ?
			null : new DrawingObjectList(mUndoDrawingObjectList);
	}

	public void storeState()
	{
		if (mUndoMol == null) {
			mUndoMol = new Molecule();
		}
		mMol.copyMolecule(mUndoMol);

		mUndoDrawingObjectList = (mDrawingObjectList == null) ?
			null : new DrawingObjectList(mDrawingObjectList);
	}

	private boolean deleteHilited()
	{
		if (mCurrentHiliteAtom != -1) {
			mMol.deleteAtom(mCurrentHiliteAtom);
			mCurrentHiliteAtom = -1;
			fireMoleculeChanged();
			update(UPDATE_REDRAW);
			return true;
		}

		if (mCurrentHiliteBond != -1) {
			mMol.deleteBondAndSurrounding(mCurrentHiliteBond);
			mCurrentHiliteBond = -1;
			fireMoleculeChanged();
			update(UPDATE_REDRAW);
			return true;
		}

		if (mCurrentHiliteObject != null && mCurrentHiliteObject.isDeletable()) {
			mDrawingObjectList.remove(mCurrentHiliteObject);
			mCurrentHiliteObject = null;
			update(UPDATE_REDRAW);
			return true;
		}
		return false;
	}

	private boolean deleteAt(double x, double y)
	{
		if (mMol.deleteAtomOrBond(x, y)) {
			fireMoleculeChanged();
			update(UPDATE_REDRAW);
			return true;
		}
		AbstractDrawingObject drawingObject = findDrawingObject(x, y, null, true);
		if (drawingObject != null) {
			mDrawingObjectList.remove(drawingObject);
			mCurrentHiliteObject = null;
			update(UPDATE_REDRAW);
			return true;
		}
		return false;
	}

	private void duplicateSelected()
	{
		int atomCount = 0;
		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
			if (mMol.isSelectedAtom(atom)) {
				atomCount++;
			}
		}
//		int bondCount = 0;
//		for (int bond=0; bond<mMol.getAllBonds(); bond++)
//			if (mMol.isSelectedBond(bond))
//				bondCount++;

		int originalAtoms = mMol.getAllAtoms();
		int originalBonds = mMol.getAllBonds();

		mX = Arrays.copyOf(mX, mX.length + atomCount);
		mY = Arrays.copyOf(mY, mY.length + atomCount);
		int[] atomMap = new int[mMol.getAllAtoms()];
		int esrGroupCountAND = mMol.renumberESRGroups(Molecule.cESRTypeAnd);
		int esrGroupCountOR = mMol.renumberESRGroups(Molecule.cESRTypeOr);
		for (int atom = 0; atom < originalAtoms; atom++) {
			if (mMol.isSelectedAtom(atom)) {
				int newAtom = mMol.getAllAtoms();
				mX[newAtom] = mX[atom];
				mY[newAtom] = mY[atom];
				atomMap[atom] = newAtom;
				mMol.copyAtom(mMol, atom, esrGroupCountAND, esrGroupCountOR);
			}
		}
		for (int bond = 0; bond < originalBonds; bond++) {
			if (mMol.isSelectedBond(bond)) {
				mMol.copyBond(mMol, bond, esrGroupCountAND, esrGroupCountOR, atomMap, false);
			}
		}
		for (int atom = 0; atom < originalAtoms; atom++) {
			mMol.setAtomSelection(atom, false);
		}
		for (int atom = originalAtoms; atom < mMol.getAllAtoms(); atom++) {
			mMol.setAtomMapNo(atom, 0, false);
		}

		if (mDrawingObjectList != null) {
			for (int i = mDrawingObjectList.size() - 1; i >= 0; i--) {
				AbstractDrawingObject object = mDrawingObjectList.get(i);
				if (object.isSelected() && !(object instanceof ReactionArrow)) {
					mDrawingObjectList.add(object.clone());
				}
			}
		}
	}

	private void fireMoleculeChanged()
	{
		fireEvent(new EditorEvent(this, EditorEvent.WHAT_MOLECULE_CHANGED, true));
	}

	private void fireEvent(EditorEvent e)
	{
		for (GenericEventListener<EditorEvent> l : mListeners) {
			l.eventHappened(e);
		}
	}

	/**
	 * Use this to inform the JDrawArea after changing its molecule from outside.
	 */
	public void moleculeChanged()
	{
		moleculeChanged(false);
	}

	/**
	 * Ideally don't use this from outside JDrawArea. Use moleculeChanged() instead.
	 *
	 * @param userChange is true if the change was done within the editor
	 */
	public void moleculeChanged(boolean userChange)
	{
		fireEvent(new EditorEvent(this, EditorEvent.WHAT_MOLECULE_CHANGED, userChange));
		update(UPDATE_SCALE_COORDS);
	}

	public StereoMolecule getMolecule()
	{
		return mMol;
	}

	public void setMolecule(StereoMolecule theMolecule)
	{
		if (mMol == theMolecule) {
			return;
		}
		mMol = theMolecule;
		storeState();
		moleculeChanged(false);
	}

	public StereoMolecule[] getFragments()
	{
		return mFragment;
	}

	public void setFragments(StereoMolecule[] fragment)
	{
		mMol.clear();
		mFragment = fragment;
		for (int i = 0; i < fragment.length; i++) {
			mMol.addMolecule(mFragment[i]);
		}
		storeState();

		mFragmentNo = new int[mMol.getAllAtoms()];
		for (int atom = 0, f = 0; f < mFragment.length; f++) {
			for (int j = 0; j < mFragment[f].getAllAtoms(); j++) {
				mFragmentNo[atom++] = f;
			}
		}

		fireEvent(new EditorEvent(this, EditorEvent.WHAT_MOLECULE_CHANGED, false));

		mMode = MODE_MULTIPLE_FRAGMENTS;
		update(UPDATE_SCALE_COORDS_USE_FRAGMENTS);
	}

	/**
	 * @return mapped reaction with absolute coordinates, but without drawing objects
	 */
	public Reaction getReaction()
	{
		if ((mMode & MODE_REACTION) == 0 || mFragment == null) {
			return null;
		}

		Reaction rxn = new Reaction();
		for (int i = 0; i < mFragment.length; i++) {
			if (i < mReactantCount) {
				rxn.addReactant(mFragment[i]);
			} else {
				rxn.addProduct(mFragment[i]);
			}
		}
		return rxn;
	}

	/**
	 * @return mapped reaction with absolute coordinates and drawing objects
	 */
	public Reaction getReactionAndDrawings()
	{
		Reaction rxn = getReaction();
		if (rxn != null) {
			rxn.setDrawingObjects(getDrawingObjects());
		}
		return rxn;
	}

	public void setReaction(Reaction rxn)
	{
		mMol.clear();
		mFragment = new StereoMolecule[rxn.getMolecules()];
		mReactantCount = rxn.getReactants();
		for (int i = 0; i < rxn.getMolecules(); i++) {
			mFragment[i] = rxn.getMolecule(i);
			mMol.addMolecule(mFragment[i]);
		}
		mMol.setFragment(rxn.isFragment());
		storeState();

		mFragmentNo = new int[mMol.getAllAtoms()];
		for (int atom = 0, f = 0; f < mFragment.length; f++) {
			for (int j = 0; j < mFragment[f].getAllAtoms(); j++) {
				mFragmentNo[atom++] = f;
			}
		}

		fireEvent(new EditorEvent(this, EditorEvent.WHAT_MOLECULE_CHANGED, false));

		mMode = MODE_MULTIPLE_FRAGMENTS | MODE_REACTION;
		update(UPDATE_SCALE_COORDS_USE_FRAGMENTS);
	}

	public MarkushStructure getMarkushStructure()
	{
		if ((mMode & MODE_MARKUSH_STRUCTURE) == 0) {
			return null;
		}

		MarkushStructure markush = new MarkushStructure();
		for (int i = 0; i < mFragment.length; i++) {
			if (i < mReactantCount) {
				markush.addCore(mFragment[i]);
			} else {
				markush.addRGroup(mFragment[i]);
			}
		}
		return markush;
	}

	public void setMarkushStructure(MarkushStructure markush)
	{
		mMol.clear();
		mFragment = new StereoMolecule[markush.getCoreCount() + markush.getRGroupCount()];
		mReactantCount = markush.getCoreCount();
		boolean isFragment = false;
		for (int i = 0; i < markush.getCoreCount() + markush.getRGroupCount(); i++) {
			mFragment[i] = (i < markush.getCoreCount()) ? markush.getCoreStructure(i)
				: markush.getRGroup(i - markush.getCoreCount());
			isFragment |= mFragment[i].isFragment();
			mMol.addMolecule(mFragment[i]);
		}
		mMol.setFragment(isFragment);
		storeState();

		mFragmentNo = new int[mMol.getAllAtoms()];
		for (int atom = 0, f = 0; f < mFragment.length; f++) {
			for (int j = 0; j < mFragment[f].getAllAtoms(); j++) {
				mFragmentNo[atom++] = f;
			}
		}

		fireEvent(new EditorEvent(this, EditorEvent.WHAT_MOLECULE_CHANGED, false));

		mMode = MODE_MULTIPLE_FRAGMENTS | MODE_MARKUSH_STRUCTURE;
		update(UPDATE_SCALE_COORDS_USE_FRAGMENTS);
	}

	public int getDisplayMode()
	{
		return mDisplayMode;
	}

	public void setDisplayMode(int dMode)
	{
		mDisplayMode = dMode;
		update(UPDATE_REDRAW);
	}

	/**
	 * If set to false then any query features will be removed from the molecule
	 * and any functionality that allows to define atom- or bond-query features
	 * won't be available. This feature is only relevant if the molecule is a fragment.
	 * @param allow
	 */
	public void setAllowQueryFeatures(boolean allow) {
		if (mAllowQueryFeatures != allow) {
			mAllowQueryFeatures = allow;
			if (!allow)
				mMol.removeQueryFeatures();
			}
		}

	/**
	 * Defines additional atom text to be displayed in top right
	 * position of some/all atom label. If the atom is charged, then
	 * the atom text follows the charge information.
	 * If using atom text make sure to update it accordingly, if atom
	 * indexes change due to molecule changes.
	 * Atom text is not supported for MODE_REACTION, MODE_MULTIPLE_FRAGMENTS or MODE_MARKUSH_STRUCTURE.
	 *
	 * @param atomText String[] matching atom indexes (may contain null entries)
	 */
	public void setAtomText(String[] atomText)
	{
		mAtomText = atomText;
	}

	public DrawingObjectList getDrawingObjects()
	{
		return mDrawingObjectList;
	}

	public void setDrawingObjects(DrawingObjectList drawingObjectList)
	{
		mDrawingObjectList = drawingObjectList;
		storeState();
		update(UPDATE_SCALE_COORDS);
	}

	public int getMode()
	{
		return mMode;
	}

	public int getHiliteAtom()
	{
		return mCurrentHiliteAtom;
	}

	public int getHiliteBond()
	{
		return mCurrentHiliteBond;
	}

	public void setHiliteBondSet(int[] bondSet)
	{
		mHiliteBondSet = bondSet;
		update(UPDATE_REDRAW);
	}

	public void setReactionMode(boolean rxn)
	{
		if (rxn) {
			Molecule m[] = this.getFragments();
			if (m == null) {
				setReaction(new Reaction(new StereoMolecule[]{new StereoMolecule(this.getMolecule())}, 1));
			} else {
				mMode = MODE_MULTIPLE_FRAGMENTS | MODE_REACTION;
				Reaction r = getReaction();
				setReaction(r);
			}
		} else {
			mMode &= ~MODE_REACTION;
		}
	}

	// CXR Need this method for derived components
	protected void setUpdateMode(int mode)
	{
		mUpdateMode = mode;
	}

	public boolean isAtomColorSupported()
	{
		return mAtomColorSupported;
	}

	public void setAtomColorSupported(boolean acs)
	{
		mAtomColorSupported = acs;
	}

	private void cleanupCoordinates(SwingDrawContext context, Graphics2D g)
	{
		int selectedAtomCount = 0;
		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
			if (mMol.isSelectedAtom(atom)) {
				selectedAtomCount++;
			}
		}
		boolean selectedOnly = (selectedAtomCount != 0 && selectedAtomCount != mMol.getAllAtoms());

		if ((mMode & MODE_MULTIPLE_FRAGMENTS) != 0)
			cleanupMultiFragmentCoordinates(g, selectedOnly);
		else
			cleanupMoleculeCoordinates(context, selectedOnly);
	}

	private void cleanupMoleculeCoordinates(SwingDrawContext context, boolean selectedOnly)
	{
		if (mUpdateMode == UPDATE_INVENT_COORDS) {
			if (selectedOnly)
				for (int atom = 0; atom < mMol.getAllAtoms(); atom++)
					mMol.setAtomMarker(atom, !mMol.isSelectedAtom(atom));

			new CoordinateInventor(selectedOnly ? CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS : 0).invent(mMol);

			if (selectedOnly)
				mMol.removeAtomMarkers();
		}

		mDepictor.updateCoords(context, new GenericRectangle(0, 0, getWidth(), getHeight()), maxUpdateMode());
	}

	private void cleanupMultiFragmentCoordinates(Graphics2D g, boolean selectedOnly)
	{
		if (selectedOnly && mUpdateMode == UPDATE_INVENT_COORDS) {
			int[] fragmentAtom = new int[mFragment.length];
			for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
				int fragment = mFragmentNo[atom];
				mFragment[fragment].setAtomMarker(fragmentAtom[fragment], !mMol.isSelectedAtom(atom));
				fragmentAtom[fragment]++;
			}
		}

		GenericRectangle[] boundingRect = new GenericRectangle[mFragment.length];
//		float fragmentWidth = 0.0f;
		for (int fragment = 0; fragment < mFragment.length; fragment++) {
			if (mUpdateMode == UPDATE_INVENT_COORDS) {
				new CoordinateInventor(selectedOnly ? CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS : 0).invent(mFragment[fragment]);
				mFragment[fragment].setStereoBondsFromParity();
			}
			Depictor2D d = new Depictor2D(mFragment[fragment]);
			d.updateCoords(g, null, AbstractDepictor.cModeInflateToMaxAVBL);
			boundingRect[fragment] = d.getBoundingRect();
//			fragmentWidth += boundingRect[fragment].width;
		}

		double spacing = FRAGMENT_CLEANUP_DISTANCE * AbstractDepictor.cOptAvBondLen;
		double avbl = mMol.getAverageBondLength();
		double arrowWidth = ((mMode & MODE_REACTION) == 0) ?
			0f
			: (mUpdateMode == UPDATE_SCALE_COORDS_USE_FRAGMENTS) ?
			DEFAULT_ARROW_LENGTH * getWidth()
			: ((ReactionArrow) mDrawingObjectList.get(0)).getLength() * AbstractDepictor.cOptAvBondLen / avbl;

		double rawX = 0.5 * spacing;
		for (int fragment = 0; fragment <= mFragment.length; fragment++) {
			if ((mMode & MODE_REACTION) != 0 && fragment == mReactantCount) {
				((ReactionArrow) mDrawingObjectList.get(0)).setCoordinates(
					rawX - spacing / 2, getHeight() / 2, rawX - spacing / 2 + arrowWidth, getHeight() / 2);
				rawX += arrowWidth;
			}

			if (fragment == mFragment.length) {
				break;
			}

			double dx = rawX - boundingRect[fragment].x;
			double dy = 0.5 * (getHeight() - boundingRect[fragment].height)
				- boundingRect[fragment].y;
			mFragment[fragment].translateCoords(dx, dy);

			rawX += spacing + boundingRect[fragment].width;
		}

		mDepictor.updateCoords(new SwingDrawContext(g), new GenericRectangle(0, 0, getWidth(), getHeight()), maxUpdateMode());

		int[] fragmentAtom = new int[mFragment.length];
		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
			int fragment = mFragmentNo[atom];
			mMol.setAtomX(atom, mFragment[fragment].getAtomX(fragmentAtom[fragment]));
			mMol.setAtomY(atom, mFragment[fragment].getAtomY(fragmentAtom[fragment]));

			fragmentAtom[fragment]++;
		}

		mMol.setStereoBondsFromParity();
	}

	private void analyzeFragmentMembership()
	{
		mMol.ensureHelperArrays(Molecule.cHelperParities);

		int[] fragmentNo = new int[mMol.getAllAtoms()];
		int fragments = mMol.getFragmentNumbers(fragmentNo, false, true);

		fragments = joinCloseFragments(fragmentNo, fragments);
		sortFragmentsByPosition(fragmentNo, fragments);
		mFragmentNo = fragmentNo;

		mFragment = mMol.getFragments(fragmentNo, fragments);
	}

	private int joinCloseFragments(int[] fragmentNo, int fragments)
	{
		if (fragments < 2) {
			return fragments;
		}

		boolean[][] mergeFragments = new boolean[fragments][];
		for (int i = 1; i < fragments; i++) {
			mergeFragments[i] = new boolean[i];
		}

		double avbl = mMol.getAverageBondLength();
		for (int atom1 = 1; atom1 < mMol.getAllAtoms(); atom1++) {
			for (int atom2 = 0; atom2 < atom1; atom2++) {
				double dx = mMol.getAtomX(atom2) - mMol.getAtomX(atom1);
				double dy = mMol.getAtomY(atom2) - mMol.getAtomY(atom1);
				double distance = Math.sqrt(dx * dx + dy * dy);
				if (distance < FRAGMENT_GROUPING_DISTANCE * avbl) {
					int fragment1 = fragmentNo[atom1];
					int fragment2 = fragmentNo[atom2];
					if (fragment1 != fragment2) {
						if (fragment1 > fragment2) {
							mergeFragments[fragment1][fragment2] = true;
						} else {
							mergeFragments[fragment2][fragment1] = true;
						}
					}
				}
			}
		}

		int[] newFragmentIndex = new int[fragments];
		for (int fragment = 0; fragment < fragments; fragment++) {
			newFragmentIndex[fragment] = fragment;
		}

		int mergeCount = 0;
		for (int i = 1; i < fragments; i++) {
			for (int j = 0; j < i; j++) {
				if (mergeFragments[i][j]) {
					int index1 = newFragmentIndex[i];
					int index2 = newFragmentIndex[j];
					if (index1 != index2) {
						mergeCount++;
						int minIndex = Math.min(index1, index2);
						int maxIndex = Math.max(index1, index2);
						for (int k = 0; k < fragments; k++) {
							if (newFragmentIndex[k] == maxIndex) {
								newFragmentIndex[k] = minIndex;
							} else if (newFragmentIndex[k] > maxIndex) {
								newFragmentIndex[k]--;
							}
						}
					}
				}
			}
		}

		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
			fragmentNo[atom] = newFragmentIndex[fragmentNo[atom]];
		}

		return fragments - mergeCount;
	}

	private void sortFragmentsByPosition(int[] fragmentNo, int fragments)
	{
		int[][] fragmentDescriptor = new int[fragments][((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE)) != 0) ? 2 : 1];
		for (int fragment = 0; fragment < fragments; fragment++) {
			fragmentDescriptor[fragment][0] = fragment;
		}

		Point[] fragmentCOG = calculateFragmentCenterOfGravity(fragmentNo, fragments);

		if ((mMode & MODE_REACTION) != 0) {
			mReactantCount = 0;
			ReactionArrow arrow = ((mMode & MODE_REACTION) != 0) ? (ReactionArrow) mDrawingObjectList.get(0) : null;
			for (int fragment = 0; fragment < fragments; fragment++) {
				fragmentDescriptor[fragment][1] = (arrow.isOnProductSide(fragmentCOG[fragment].x,
					fragmentCOG[fragment].y)) ? 1 : 0;
				if (fragmentDescriptor[fragment][1] == 0) {
					mReactantCount++;
				}
			}
		} else if ((mMode & MODE_MARKUSH_STRUCTURE) != 0) {
			mReactantCount = fragments;
			for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
				if (mMol.getAtomicNo(atom) == 0 && fragmentDescriptor[fragmentNo[atom]][1] == 0) {
					fragmentDescriptor[fragmentNo[atom]][1] = 1;
					mReactantCount--;
				}
			}
		}

		final Point[] cog = fragmentCOG;
		Arrays.sort(fragmentDescriptor, new Comparator<int[]>()
		{
			public int compare(int[] fragmentDescriptor1, int[] fragmentDescriptor2)
			{
				if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE)) != 0) {
					if (fragmentDescriptor1[1] != fragmentDescriptor2[1]) {
						return (fragmentDescriptor1[1] == 0) ? -1 : 1;
					}
				}

				return (cog[fragmentDescriptor1[0]].x
					+ cog[fragmentDescriptor1[0]].y
					< cog[fragmentDescriptor2[0]].x
					+ cog[fragmentDescriptor2[0]].y) ? -1 : 1;
			}
		});

		int[] newFragmentIndex = new int[fragments];
		Point[] centerOfGravity = new Point[fragments];
		for (int fragment = 0; fragment < fragments; fragment++) {
			int oldIndex = ((int[]) fragmentDescriptor[fragment])[0];
			newFragmentIndex[oldIndex] = fragment;
			centerOfGravity[fragment] = fragmentCOG[oldIndex];
		}

		fragmentCOG = centerOfGravity;
		for (int atom1 = 0; atom1 < mMol.getAllAtoms(); atom1++) {
			fragmentNo[atom1] = newFragmentIndex[fragmentNo[atom1]];
		}
	}

	private Point[] calculateFragmentCenterOfGravity(int[] fragmentNo, int fragments)
	{
		Point[] fragmentCOG = new Point[fragments];
		int[] fragmentAtoms = new int[fragments];
		for (int fragment = 0; fragment < fragments; fragment++) {
			fragmentCOG[fragment] = new Point(0, 0);
		}
		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
			fragmentCOG[fragmentNo[atom]].x += mMol.getAtomX(atom);
			fragmentCOG[fragmentNo[atom]].y += mMol.getAtomY(atom);
			fragmentAtoms[fragmentNo[atom]]++;
		}
		for (int fragment = 0; fragment < fragments; fragment++) {
			fragmentCOG[fragment].x /= fragmentAtoms[fragment];
			fragmentCOG[fragment].y /= fragmentAtoms[fragment];
		}
		return fragmentCOG;
	}

	private void updateCursor()
	{
		int cursor = -1;
		switch (mCurrentTool) {
			case JDrawToolbar.cToolZoom:
				cursor = SwingCursorHelper.cZoomCursor;
				break;
			case JDrawToolbar.cToolLassoPointer:
				if ((mCurrentHiliteAtom != -1 && mMol.isSelectedAtom(mCurrentHiliteAtom))
					|| (mCurrentHiliteBond != -1 && mMol.isSelectedBond(mCurrentHiliteBond))) {
					cursor = mMouseIsDown ? SwingCursorHelper.cFistCursor
						: mShiftIsDown ? SwingCursorHelper.cHandPlusCursor
						: SwingCursorHelper.cHandCursor;
				} else if (mCurrentHiliteAtom != -1
					|| mCurrentHiliteBond != -1) {
					cursor = SwingCursorHelper.cPointerCursor;
				} else if (mCurrentHiliteObject != null) {
					cursor = mMouseIsDown ? SwingCursorHelper.cFistCursor
						: (mShiftIsDown
						&& !(mCurrentHiliteObject instanceof ReactionArrow)) ?
						SwingCursorHelper.cHandPlusCursor : SwingCursorHelper.cHandCursor;
				} else {
					cursor = mShiftIsDown ?
						(mAltIsDown ? SwingCursorHelper.cSelectRectPlusCursor : SwingCursorHelper.cLassoPlusCursor)
						: (mAltIsDown ? SwingCursorHelper.cSelectRectCursor : SwingCursorHelper.cLassoCursor);
				}
				break;
			case JDrawToolbar.cToolDelete:
				cursor = SwingCursorHelper.cDeleteCursor;
				break;
			case JDrawToolbar.cToolChain:
				cursor = SwingCursorHelper.cChainCursor;
				break;
			case JDrawToolbar.cToolText:
				cursor = SwingCursorHelper.cTextCursor;
				break;
			default:
				cursor = SwingCursorHelper.cPointerCursor;
				break;
		}

		if (mCurrentCursor != cursor) {
			mCurrentCursor = cursor;
			setCursor(SwingCursorHelper.getCursor(cursor));
		}
	}

	private void initializeDragAndDrop(int dropAction)
	{
		if (dropAction != DnDConstants.ACTION_NONE) {
			MoleculeDropAdapter d = new MoleculeDropAdapter()
			{
				public void onDropMolecule(StereoMolecule m, Point pt)
				{
					if (m != null && m.getAllAtoms() != 0) {
						for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
							mMol.setAtomSelection(atom, false);				// deselect all current atoms
						}

						double avbl = mMol.getAverageBondLength();
						double dropAVBL = m.getAverageBondLength();
						m.scaleCoords(avbl / dropAVBL);
						Point cog = new Point();
						for (int atom = 0; atom < m.getAllAtoms(); atom++) {
							m.setAtomSelection(atom, true);					// select all new atoms
							cog.x += m.getAtomX(atom);
							cog.y += m.getAtomY(atom);
						}
						cog.x /= m.getAllAtoms();
						cog.y /= m.getAllAtoms();
						m.translateCoords(pt.x - cog.x, pt.y - cog.y);
						m.removeAtomColors();
						m.removeBondHiliting();
						mMol.addMolecule(m);
						fireMoleculeChanged();
						update(UPDATE_CHECK_COORDS);
					}
				}
			};

			new DropTarget(this, dropAction, d, true, new OurFlavorMap());
		}
	}

// This class is needed for inter-jvm drag&drop. Although not neccessary for standard environments, it prevents
// nasty "no native data was transfered" errors. It still might create ClassNotFoundException in the first place by
// the SystemFlavorMap, but as I found it does not hurt, since the context classloader will be installed after
// the first call. I know, that this depends heavely on a specific behaviour of the systemflavormap, but for now
// there's nothing I can do about it.
	static class OurFlavorMap implements java.awt.datatransfer.FlavorMap
	{
		@Override
		public java.util.Map<DataFlavor, String> getNativesForFlavors(DataFlavor[] dfs)
		{
			java.awt.datatransfer.FlavorMap m = java.awt.datatransfer.SystemFlavorMap.getDefaultFlavorMap();
			return m.getNativesForFlavors(dfs);
		}

		@Override
		public java.util.Map<String, DataFlavor> getFlavorsForNatives(String[] natives)
		{
			java.awt.datatransfer.FlavorMap m = java.awt.datatransfer.SystemFlavorMap.getDefaultFlavorMap();
			return m.getFlavorsForNatives(natives);
		}
	}


/*	private void dumpBytesOfGif(String imageFileName) {
		byte[] data = new byte[10000];
		try {
			BufferedInputStream in=new BufferedInputStream(getClass().getResourceAsStream(imageFileName));
			int count = in.read(data);
			Image image = Toolkit.getDefaultToolkit().createImage(data);
			PixelGrabber pg = new PixelGrabber(image, 0, 0, 16, 16, false);
			try {
				pg.grabPixels();
				byte[] pixels = (byte[])pg.getPixels();
				int pixelRow = 0;
				for (int i=0; i<pixels.length; i++) {
					if ((i%16) == 0)
						pixelRow = 0;
					pixelRow = (pixelRow<<2) + pixels[i];
					if ((i%16) == 15)
						System.out.print("0x"+Integer.toHexString(pixelRow)+", ");
					if ((i%80) == 79)
						System.out.println();
					}
				System.out.println("},{");
				}
			catch (InterruptedException e) {
				System.out.println("Pixel Grabbing interrupted: "+e);
				}
			}
		catch (Exception e) {
			System.out.println("Error reading data: "+e);
			}
		}*/

	private int maxUpdateMode() {
		return AbstractDepictor.cModeInflateToMaxAVBL + HiDPIHelper.scale(AbstractDepictor.cOptAvBondLen);
	}

	private Point2D calculateCenterOfGravity(boolean selectedOnly)
	{
		int atoms = 0;
		double sumx = 0;
		double sumy = 0;
		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
			if (!selectedOnly || mMol.isSelectedAtom(atom)) {
				sumx += mMol.getAtomX(atom);
				sumy += mMol.getAtomY(atom);
				atoms++;
			}
		}
		return atoms > 1 ? new Point2D.Double(sumx / atoms, sumy / atoms) : null;
	}

	private void flip(boolean horiz)
	{
		boolean selectedOnly = false;
		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
			if (mMol.isSelectedAtom(atom)) {
				selectedOnly = true;
				break;
			}
		}

		Point2D cog = calculateCenterOfGravity(selectedOnly);
		if (cog != null) {
			storeState();

			for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
				if (!selectedOnly || mMol.isSelectedAtom(atom)) {
					if (horiz) {
						mMol.setAtomX(atom, 2 * cog.getX() - mMol.getAtomX(atom));
					}
					else {
						mMol.setAtomY(atom, 2 * cog.getY() - mMol.getAtomY(atom));
					}
				}
			}

			// invert stereo bonds
			for (int bond=0; bond<mMol.getAllBonds(); bond++) {
				if (!selectedOnly || mMol.isSelectedBond(bond)) {
					if (mMol.getBondType(bond) == Molecule.cBondTypeUp)
						mMol.setBondType(bond, Molecule.cBondTypeDown);
					else if (mMol.getBondType(bond) == Molecule.cBondTypeDown)
						mMol.setBondType(bond, Molecule.cBondTypeUp);
				}
			}

			update(UPDATE_REDRAW);
		}
	}
}

