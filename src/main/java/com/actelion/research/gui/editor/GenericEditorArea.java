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

package com.actelion.research.gui.editor;

import com.actelion.research.chem.*;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.io.RDFileParser;
import com.actelion.research.chem.io.RXNFileParser;
import com.actelion.research.chem.name.StructureNameResolver;
import com.actelion.research.chem.reaction.*;
import com.actelion.research.gui.FileHelper;
import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.gui.clipboard.IClipboardHandler;
import com.actelion.research.gui.generic.*;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.gui.swing.SwingCursorHelper;
import com.actelion.research.util.ColorHelper;

import java.awt.*;
import java.awt.geom.Point2D;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.TreeMap;

public class GenericEditorArea implements GenericEventListener {
	public static final int MODE_MULTIPLE_FRAGMENTS = 1;
	public static final int MODE_MARKUSH_STRUCTURE = 2;
	public static final int MODE_REACTION = 4;
	public static final int MODE_DRAWING_OBJECTS = 8;

	public static final String TEMPLATE_TYPE_KEY = "TEMPLATE";
	public static final String TEMPLATE_TYPE_REACTION_QUERIES = "REACTION_QUERIES";
	public static final String TEMPLATE_SECTION_KEY = "SECTION";

	public static final int DEFAULT_ALLOWED_PSEUDO_ATOMS
			= Molecule.cPseudoAtomsHydrogenIsotops
			| Molecule.cPseudoAtomsAminoAcids
			| Molecule.cPseudoAtomR
			| Molecule.cPseudoAtomsRGroups
			| Molecule.cPseudoAtomAttachmentPoint;

	private static final int MAX_CONNATOMS = 8;
	private static final int MIN_BOND_LENGTH_SQUARE = 100;

	private static final int KEY_IS_ATOM_LABEL = 1;
	private static final int KEY_IS_SUBSTITUENT = 2;
	private static final int KEY_IS_VALID_START = 3;
	private static final int KEY_IS_INVALID = 4;

	private static final int RGB_BLUE = 0xFF0000FF;
	private static final int RGB_RED = 0xFFFF0000;
	private static final int RGB_DARK_RED = 0xFF800000;
	private static final int RGB_BLACK = 0xFF000000;
	private static final int RGB_GRAY = 0xFF808080;

	private static final String ITEM_COPY_STRUCTURE = "Copy Structure";
	private static final String ITEM_COPY_REACTION = "Copy Reaction";
	private static final String ITEM_PASTE_STRUCTURE = "Paste Structure";
	private static final String ITEM_PASTE_REACTION = "Paste Reaction";
	private static final String ITEM_USE_TEMPLATE = "Use Reaction Template...";
	private static final String ITEM_PASTE_WITH_NAME = ITEM_PASTE_STRUCTURE + " or Name";
	private static final String ITEM_LOAD_REACTION = "Open Reaction File...";
	private static final String ITEM_ADD_AUTO_MAPPING = "Auto-Map Reaction";
	private static final String ITEM_REMOVE_MAPPING = "Remove Manual Atom Mapping";
	private static final String ITEM_FLIP_HORIZONTALLY = "Flip Horizontally";
	private static final String ITEM_FLIP_VERTICALLY = "Flip Vertically";
	private static final String ITEM_FLIP_ROTATE180 = "Rotate 180°";
	private static final String ITEM_SHOW_HELP = "Help Me";

	// development items
	private static final String ITEM_SHOW_ATOM_BOND_NUMBERS = "Show Atom & Bond Numbers";
	private static final String ITEM_SHOW_SYMMETRY = "Show Symmetry";
	private static final String ITEM_SHOW_NORMAL = "Show Normal";

	private static final long WARNING_MILLIS = 1200;

	private static final float FRAGMENT_MAX_CLICK_DISTANCE = 24.0f;
	private static final float FRAGMENT_GROUPING_DISTANCE = 1.4f;    // in average bond lengths
	private static final float FRAGMENT_CLEANUP_DISTANCE = 1.5f;    // in average bond lengths
	private static final float DEFAULT_ARROW_LENGTH = 0.08f;        // relative to panel width

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

	private static final int DEFAULT_SELECTION_BACKGROUND = 0xFF80A4C0;

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
	private static String[][] sReactionQueryTemplates;
	private int mMode, mChainAtoms, mCurrentTool, mCustomAtomicNo, mCustomAtomMass, mCustomAtomValence, mCustomAtomRadical,
			mCurrentHiliteAtom, mCurrentHiliteBond, mPendingRequest, mEventsScheduled, mFirstAtomKey,
			mCurrentCursor, mReactantCount, mUpdateMode, mDisplayMode, mAtom1, mAtom2, mMaxAVBL, mAllowedPseudoAtoms;
	private int[] mChainAtom, mFragmentNo, mHiliteBondSet;
	private double mX1, mY1, mX2, mY2, mWidth, mHeight, mUIScaling, mTextSizeFactor;
	private double[] mX, mY, mChainAtomX, mChainAtomY;
	private boolean mAltIsDown, mShiftIsDown, mMouseIsDown, mIsAddingToSelection, mAtomColorSupported, mAllowQueryFeatures,
			mAllowFragmentChangeOnPasteOrDrop;
	private boolean[] mIsSelectedAtom, mIsSelectedObject;
	private String mCustomAtomLabel,mWarningMessage,mAtomKeyStrokeSuggestion;
	private String[] mAtomText;
	private ExtendedDepictor mDepictor;
	private StereoMolecule mMol;        // molecule being modified directly by the drawing editor
	private Molecule mUndoMol;          // molecule in undo buffer
	private StereoMolecule[] mFragment;    // in case of MODE_MULTIPLE_FRAGMENTS contains valid stereo fragments
	// for internal and external read-only-access (reconstructed at any change)
	private DrawingObjectList mDrawingObjectList, mUndoDrawingObjectList;
	private AbstractDrawingObject mCurrentHiliteObject;
	private GenericPolygon mLassoRegion;
	private ArrayList<GenericEventListener> mListeners;
	private IClipboardHandler mClipboardHandler;
	private StringBuilder mAtomKeyStrokeBuffer;
	private GenericUIHelper mUIHelper;
	private GenericCanvas mCanvas;

	/**
	 * @param mol  an empty or valid stereo molecule
	 * @param mode 0 or a meaningful combination of the mode flags, e.g. MODE_REACTION | MODE_DRAWING_OBJECTS
	 */
	public GenericEditorArea(StereoMolecule mol, int mode, GenericUIHelper helper, GenericCanvas canvas) {
		mMol = mol;
		mMode = mode;
		mUIHelper = helper;
		mCanvas = canvas;

		mListeners = new ArrayList<>();

		mCurrentTool = GenericEditorToolbar.cToolStdBond;
		mCurrentHiliteAtom = -1;
		mCurrentHiliteBond = -1;
		mCurrentHiliteObject = null;
		mAtom1 = -1;
		mCustomAtomicNo = 6;
		mCustomAtomMass = 0;
		mCustomAtomValence = -1;
		mCustomAtomRadical = 0;
		mCustomAtomLabel = null;
		mAllowQueryFeatures = true;
		mAllowFragmentChangeOnPasteOrDrop = false;
		mPendingRequest = cRequestNone;
		mCurrentCursor = SwingCursorHelper.cPointerCursor;
		mAtomKeyStrokeBuffer = new StringBuilder();

		mAllowedPseudoAtoms = DEFAULT_ALLOWED_PSEUDO_ATOMS;

		mTextSizeFactor = 1.0;

		mUIScaling = HiDPIHelper.getUIScaleFactor();
		mMaxAVBL = HiDPIHelper.scale(AbstractDepictor.cOptAvBondLen);

		if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE)) != 0) {
			mMode |= (MODE_MULTIPLE_FRAGMENTS);
		}

		if ((mMode & (MODE_DRAWING_OBJECTS | MODE_REACTION)) != 0) {
			mDrawingObjectList = new DrawingObjectList();
		}

		mUpdateMode = UPDATE_SCALE_COORDS;

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
	 * @return null or String[n][2] with pairs of reaction name and rxn-idcode
	 */
	public static String[][] getReactionQueryTemplates() {
		return sReactionQueryTemplates;
	}

	/**
	 * @param templates null or String[n][2] with pairs of reaction name and rxn-idcode
	 */
	public static void setReactionQueryTemplates(String[][] templates) {
		sReactionQueryTemplates = templates;
	}

	public GenericUIHelper getUIHelper() {
		return mUIHelper;
	}

	/**
	 * Call this after initialization to get clipboard support
	 *
	 * @param h
	 */
	public void setClipboardHandler(IClipboardHandler h) {
		mClipboardHandler = h;
	}

	private void update(int mode) {
		mUpdateMode = Math.max(mUpdateMode, mode);
		mCanvas.repaint();
	}

	public static void setReactionMapper(IReactionMapper mapper) {
		sMapper = mapper;
	}

	public void paintContent(GenericDrawContext context) {
		if (mWidth != mCanvas.getCanvasWidth() || mHeight != mCanvas.getCanvasHeight()) {
			mWidth = mCanvas.getCanvasWidth();
			mHeight = mCanvas.getCanvasHeight();
			if (mUpdateMode<UPDATE_CHECK_COORDS) {
				mUpdateMode = UPDATE_CHECK_COORDS;
			}
		}

		int background = context.getBackgroundRGB();
		int foreground = context.getForegroundRGB();

		context.setRGB(background);
		context.fillRectangle(0, 0, mWidth, mHeight);

		if ((mMode & MODE_REACTION) != 0 && mDrawingObjectList.size() == 0) {
			float mx = 0.5f * (float)mWidth;
			float my = 0.5f * (float)mHeight;
			float dx = 0.5f * DEFAULT_ARROW_LENGTH * (float)mWidth;
			ReactionArrow arrow = new ReactionArrow();
			arrow.setCoordinates(mx - dx, my, mx + dx, my);
			arrow.setDeletable(false);
			mDrawingObjectList.add(arrow);
		}

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
			mDepictor.setDefaultAVBL(mMaxAVBL);

			mDepictor.setFragmentNoColor(((mMode & MODE_MULTIPLE_FRAGMENTS) == 0) ? 0
					: LookAndFeelHelper.isDarkLookAndFeel() ? ColorHelper.brighter(background, 0.85f)
					: ColorHelper.darker(background, 0.85f));

			mDepictor.setFactorTextSize(mTextSizeFactor);

			mDepictor.setDisplayMode(mDisplayMode
					| AbstractDepictor.cDModeHiliteAllQueryFeatures
					| ((mCurrentTool == GenericEditorToolbar.cToolMapper) ?
					AbstractDepictor.cDModeShowMapping
							| AbstractDepictor.cDModeSuppressCIPParity : 0));

			if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE | MODE_MULTIPLE_FRAGMENTS)) == 0) {
				mDepictor.getMoleculeDepictor(0).setAtomText(mAtomText);
			}

			switch (mUpdateMode) {
				case UPDATE_INVENT_COORDS:
				case UPDATE_SCALE_COORDS:
				case UPDATE_SCALE_COORDS_USE_FRAGMENTS:
					cleanupCoordinates(context);
					break;
				case UPDATE_CHECK_COORDS:
					DepictorTransformation t1 = mDepictor.updateCoords(context, new GenericRectangle(0, 0, mWidth, mHeight), 0);
					if (t1 != null && (mMode & MODE_MULTIPLE_FRAGMENTS) != 0) {
						// in fragment mode depictor transforms mFragment[] rather than mMol
						t1.applyTo(mMol);
					}
					break;
				case UPDATE_CHECK_VIEW:
					DepictorTransformation t2 = mDepictor.validateView(context, new GenericRectangle(0, 0, mWidth, mHeight), 0);
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
			drawHiliting(context);
		}

		if (mDepictor != null) {
			mDepictor.paintStructures(context);
			mDepictor.paintDrawingObjects(context);
		}

		if (mCurrentHiliteAtom != -1 && mAtomKeyStrokeBuffer.length() != 0) {
			int x = (int)mMol.getAtomX(mCurrentHiliteAtom);
			int y = (int)mMol.getAtomY(mCurrentHiliteAtom);
			String s = mAtomKeyStrokeBuffer.toString();
			int validity = getAtomKeyStrokeValidity(s);
			context.setRGB((validity == KEY_IS_ATOM_LABEL) ? foreground
						   : (validity == KEY_IS_SUBSTITUENT) ? RGB_BLUE : RGB_GRAY);
			if (validity == KEY_IS_ATOM_LABEL)
				s = Molecule.cAtomLabel[Molecule.getAtomicNoFromLabel(s, mAllowedPseudoAtoms)];
			else if (validity == KEY_IS_SUBSTITUENT)
				s = mAtomKeyStrokeSuggestion.substring(0, s.length());
			int fontSize = 3*context.getFontSize()/2;
			context.setFont(fontSize, false, false);
			context.drawString(x, y, s);
			if (validity == KEY_IS_INVALID) {
				context.setRGB(RGB_RED);
				context.drawCenteredString(x+context.getBounds(s).getWidth()/2, y+fontSize, "<unknown>");
				}
			if (validity == KEY_IS_SUBSTITUENT) {
				context.setRGB(RGB_GRAY);
				context.drawString(x+context.getBounds(s).getWidth(), y, mAtomKeyStrokeSuggestion.substring(s.length()));
			}
		}

		context.setRGB(foreground);
		switch (mPendingRequest) {
			case cRequestNewBond:
				int x1, y1, x2, y2, xdiff, ydiff;
				x1 = (int)mX1;
				y1 = (int)mY1;
				if (mCurrentHiliteAtom == -1 || mCurrentHiliteAtom == mAtom1) {
					x2 = (int)mX2;
					y2 = (int)mY2;
				} else {
					x2 = (int)mMol.getAtomX(mCurrentHiliteAtom);
					y2 = (int)mMol.getAtomY(mCurrentHiliteAtom);
				}
				switch (mCurrentTool) {
					case GenericEditorToolbar.cToolStdBond:
						context.drawLine(x1, y1, x2, y2);
						break;
					case GenericEditorToolbar.cToolUpBond:
						xdiff = (y1 - y2) / 9;
						ydiff = (x2 - x1) / 9;
						GenericPolygon p = new GenericPolygon(3);
						p.addPoint( x1, y1);
						p.addPoint(x2 + xdiff, y2 + ydiff);
						p.addPoint(x2 - xdiff, y2 - ydiff);
						context.fillPolygon(p);
						break;
					case GenericEditorToolbar.cToolDownBond:
						int xx1, xx2, yy1, yy2;
						xdiff = x2 - x1;
						ydiff = y2 - y1;
						for (int i = 2; i<17; i += 2) {
							xx1 = x1 + i * xdiff / 17 - i * ydiff / 128;
							yy1 = y1 + i * ydiff / 17 + i * xdiff / 128;
							xx2 = x1 + i * xdiff / 17 + i * ydiff / 128;
							yy2 = y1 + i * ydiff / 17 - i * xdiff / 128;
							context.drawLine(xx1, yy1, xx2, yy2);
						}
						break;
				}
				break;
			case cRequestNewChain:
				if (mChainAtoms>0) {
					context.drawLine((int)mX1, (int)mY1, (int)mChainAtomX[0], (int)mChainAtomY[0]);
				}
				if (mChainAtoms>1) {
					for (int i = 1; i<mChainAtoms; i++) {
						context.drawLine((int)mChainAtomX[i - 1], (int)mChainAtomY[i - 1],
								(int)mChainAtomX[i], (int)mChainAtomY[i]);
					}
				}
				break;
			case cRequestLassoSelect:
				context.setRGB(lassoColor(context));
				context.drawPolygon(mLassoRegion);
				context.setRGB(foreground);
				break;
			case cRequestSelectRect:
				int x = (mX1<mX2) ? (int)mX1 : (int)mX2;
				int y = (mY1<mY2) ? (int)mY1 : (int)mY2;
				int w = (int)Math.abs(mX2 - mX1);
				int h = (int)Math.abs(mY2 - mY1);
				context.setRGB(lassoColor(context));
				context.drawRectangle(x, y, w, h);
				context.setRGB(foreground);
				break;
			case cRequestMapAtoms:
				x1 = (int)mX1;
				y1 = (int)mY1;
				if (mCurrentHiliteAtom == -1 || mCurrentHiliteAtom == mAtom1) {
					x2 = (int)mX2;
					y2 = (int)mY2;
				} else {
					x2 = (int)mMol.getAtomX(mCurrentHiliteAtom);
					y2 = (int)mMol.getAtomY(mCurrentHiliteAtom);
				}
				context.setRGB(mapToolColor(context));
				context.drawLine(x1, y1, x2, y2);
				context.setRGB(foreground);
				break;
		}

		if (mWarningMessage != null) {
			context.setFont(12, true, false);
			int saveRGB = context.getRGB();
			context.setRGB(RGB_RED);
			context.drawCenteredString(mWidth / 2, context.getFontSize(), mWarningMessage);
			context.setRGB(saveRGB);
		}
	}

	private double getScaledAVBL() {
		return mMol.getAverageBondLength(Molecule.cDefaultAVBL * mUIScaling);
	}

	public static int lassoColor(GenericDrawContext context) {
		int selectionColor = selectionColor(context);
		return ColorHelper.createColor(selectionColor, LookAndFeelHelper.isDarkLookAndFeel() ? 0.65f : 0.35f);
	}

	public static int selectionColor(GenericDrawContext context) {
		int selectionColor = context.getSelectionBackgroundRGB();
		return (selectionColor != 0) ? selectionColor : DEFAULT_SELECTION_BACKGROUND;
	}

	public static int mapToolColor(GenericDrawContext context) {
		int background = context.getBackgroundRGB();
		return ColorHelper.getContrastColor(RGB_DARK_RED, background);
	}

	public static int chainHiliteColor(GenericDrawContext context) {
		int background = context.getBackgroundRGB();
		int selectionColor = selectionColor(context);
		return ColorHelper.intermediateColor(selectionColor, background, 0.5f);
	}

	private void drawHiliting(GenericDrawContext context) {
		if (mHiliteBondSet != null) {
			context.setRGB(chainHiliteColor(context));
			for (int i = 0; i<mHiliteBondSet.length; i++) {
				hiliteBond(context, mHiliteBondSet[i]);
			}
		}

		if (mCurrentHiliteAtom != -1) {
			context.setRGB(selectionColor(context));
			hiliteAtom(context, mCurrentHiliteAtom);
			if (mCurrentTool == GenericEditorToolbar.cToolMapper) {
				int mapNo = mMol.getAtomMapNo(mCurrentHiliteAtom);
				if (mapNo != 0) {
					for (int atom = 0; atom<mMol.getAtoms(); atom++) {
						if (atom != mCurrentHiliteAtom
								&& mMol.getAtomMapNo(atom) == mapNo) {
							hiliteAtom(context, atom);
						}
					}
				}
			}
		}

		if (mCurrentHiliteBond != -1) {
			context.setRGB(selectionColor(context));
			hiliteBond(context, mCurrentHiliteBond);
		}

		if (mCurrentHiliteObject != null) {
			mCurrentHiliteObject.hilite(context);
		}
	}

	private void hiliteAtom(GenericDrawContext context, int atom) {
		int radius = (int)(0.32f * mMol.getAverageBondLength());

		int x = (int)mMol.getAtomX(atom);
		int y = (int)mMol.getAtomY(atom);
		context.fillCircle(x - radius, y - radius, 2 * radius);
		}

	private void hiliteBond(GenericDrawContext context, int bond) {
		int width = (int)(0.32f * mMol.getAverageBondLength());

		int x1 = (int)mMol.getAtomX(mMol.getBondAtom(0, bond));
		int y1 = (int)mMol.getAtomY(mMol.getBondAtom(0, bond));
		int x2 = (int)mMol.getAtomX(mMol.getBondAtom(1, bond));
		int y2 = (int)mMol.getAtomY(mMol.getBondAtom(1, bond));

		float oldWidth = context.getLineWidth();
		context.setLineWidth(width);
		context.drawLine(x1, y1, x2, y2);
		context.setLineWidth(oldWidth);
		}

	public void addDrawAreaListener(GenericEventListener<EditorEvent> l) {
		mListeners.add(l);
	}

	protected void buttonPressed(int button) {
		switch (button) {
			case GenericEditorToolbar.cButtonClear:
				clearAll();
				return;
			case GenericEditorToolbar.cButtonCleanStructure:
				storeState();
				updateAndFireEvent(UPDATE_INVENT_COORDS);
				return;
			case GenericEditorToolbar.cButtonUndo:
				restoreState();
				updateAndFireEvent(UPDATE_CHECK_VIEW);
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
		if (mUndoMol.getAllAtoms() != 0)
			updateAndFireEvent(UPDATE_REDRAW);
		else
			update(UPDATE_REDRAW);
	}

	public void toolChanged(int newTool) {
		if (mCurrentTool != newTool) {
			if (mCurrentTool == GenericEditorToolbar.cToolMapper
					|| newTool == GenericEditorToolbar.cToolMapper) {
				update(UPDATE_REDRAW);
			}

			mCurrentTool = newTool;
		}
	}

	public void setCustomAtom(int atomicNo, int mass, int valence, int radical, String customLabel) {
		mCustomAtomicNo = atomicNo;
		mCustomAtomMass = mass;
		mCustomAtomValence = valence;
		mCustomAtomRadical = radical;
		mCustomAtomLabel = customLabel;
	}

	@Override
	public void eventHappened(GenericEvent e) {
		if (e instanceof GenericActionEvent)
			eventHappened((GenericActionEvent)e);
		else if (e instanceof GenericKeyEvent)
			eventHappened((GenericKeyEvent)e);
		else if (e instanceof GenericMouseEvent)
			eventHappened((GenericMouseEvent)e);
	}

	private void eventHappened(GenericActionEvent e) {
		String command = e.getMessage();
		if (command.equals(ITEM_COPY_STRUCTURE) || command.equals(ITEM_COPY_REACTION)) {
			copy();
		} else if (command.equals(ITEM_PASTE_REACTION)) {
			pasteReaction();
		} else if (command.startsWith(ITEM_USE_TEMPLATE)) {
			useTemplate(command.substring(ITEM_USE_TEMPLATE.length()));
		} else if (command.startsWith(ITEM_PASTE_STRUCTURE)) {
			pasteMolecule();
		} else if (command.equals(ITEM_LOAD_REACTION)) {
			openReaction();
		} else if (command.equals(ITEM_ADD_AUTO_MAPPING)) {
			autoMapReaction();
			updateAndFireEvent(Math.max(mUpdateMode, UPDATE_REDRAW));
		} else if (command.equals(ITEM_REMOVE_MAPPING)) {
			removeManualMapping();
		} else if (command.equals(ITEM_FLIP_HORIZONTALLY)) {
			flip(true);
		} else if (command.equals(ITEM_FLIP_VERTICALLY)) {
			flip(false);
		} else if (command.equals(ITEM_FLIP_ROTATE180)) {
			rotate180();
		} else if (command.equals(ITEM_SHOW_HELP)) {
			showHelpDialog();
		} else if (command.startsWith("atomColor")) {
			int index = command.indexOf(':');
			int atom = Integer.parseInt(command.substring(9, index));
			int color = Integer.parseInt(command.substring(index + 1));
			if (mMol.isSelectedAtom(atom)) {
				for (int i = 0; i<mMol.getAtoms(); i++) {
					if (mMol.isSelectedAtom(i)) {
						mMol.setAtomColor(i, color);
					}
				}
			} else {
				mMol.setAtomColor(atom, color);
			}
		} else if (command.equals(ITEM_SHOW_ATOM_BOND_NUMBERS)) {
			setDisplayMode(AbstractDepictor.cDModeAtomNo | AbstractDepictor.cDModeBondNo);
		} else if (command.equals(ITEM_SHOW_SYMMETRY)) {
			setDisplayMode(AbstractDepictor.cDModeShowSymmetrySimple);
		} else if (command.equals(ITEM_SHOW_NORMAL)) {
			setDisplayMode(mCurrentTool == GenericEditorToolbar.cToolMapper ? AbstractDepictor.cDModeShowMapping : 0);
		}
	}

	private void removeManualMapping() {
		boolean changed = false;
		for (int atom = 0; atom<mMol.getAtoms(); atom++) {
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
			updateAndFireEvent(Math.max(mUpdateMode, UPDATE_REDRAW));
		}
	}

	/**
	 * Checks, whether a copy operation would copy a molecule or reaction.
	 *
	 * @param doCopy if true, then the chemistry object is copied to the clipboard
	 * @return true, if reaction
	 */
	private boolean analyseCopy(boolean doCopy) {
		boolean isReaction = ((mMode & MODE_REACTION) != 0);
		boolean selectionFound = false;
		boolean isBothSideSelection = false;
		boolean isOnProductSide = false;
		ReactionArrow arrow = null;

		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
			if (mMol.isSelectedAtom(atom)) {
				if (!selectionFound) {
					selectionFound = true;
					if (!isReaction) {
						break;
					}

					arrow = (ReactionArrow)mDrawingObjectList.get(0);
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

	private void copy() {
		analyseCopy(true);
	}

	private boolean copyReaction(boolean selectionOnly) {
		Reaction rxn = selectionOnly ? getSelectedReaction() : getReaction();
		if (rxn != null && mClipboardHandler != null) {
			return mClipboardHandler.copyReaction(rxn);
		}

		return false;
	}

	private Reaction getSelectedReaction() {
		Reaction rxn = new Reaction();
		for (int i = 0; i<mFragment.length; i++) {
			StereoMolecule selectedMol = getSelectedCopy(mFragment[i]);
			if (selectedMol != null) {
				if (i<mReactantCount) {
					rxn.addReactant(selectedMol);
				} else {
					rxn.addProduct(selectedMol);
				}
			}
		}
		return rxn;
	}

	private StereoMolecule getSelectedCopy(StereoMolecule sourceMol) {
		int atomCount = 0;
		for (int atom = 0; atom<sourceMol.getAllAtoms(); atom++) {
			if (sourceMol.isSelectedAtom(atom)) {
				atomCount++;
			}
		}

		if (atomCount == 0) {
			return null;
		}

		int bondCount = 0;
		for (int bond = 0; bond<sourceMol.getAllBonds(); bond++) {
			if (sourceMol.isSelectedBond(bond)) {
				bondCount++;
			}
		}

		boolean[] includeAtom = new boolean[sourceMol.getAllAtoms()];
		for (int atom = 0; atom<sourceMol.getAllAtoms(); atom++) {
			includeAtom[atom] = sourceMol.isSelectedAtom(atom);
		}

		StereoMolecule destMol = new StereoMolecule(atomCount, bondCount);
		sourceMol.copyMoleculeByAtoms(destMol, includeAtom, false, null);
		return destMol;
	}

	private boolean copyMolecule(boolean selectionOnly) {
		if (mMol.getAllAtoms() != 0 && mClipboardHandler != null) {
			return mClipboardHandler.copyMolecule(selectionOnly ? getSelectedCopy(mMol) : mMol.getCompactCopy());
		}

		return false;
	}

	private void paste() {
		if ((mMode & MODE_REACTION) != 0) {
			if (pasteReaction()) {
				return;
			}
		}

		pasteMolecule();
	}

	private boolean pasteReaction() {
		boolean ret = false;
		if (mClipboardHandler != null) {
			Reaction rxn = mClipboardHandler.pasteReaction();
			if (rxn != null) {
				if (!mAllowFragmentChangeOnPasteOrDrop)
					for (int i = 0; i<rxn.getMolecules(); i++)
						rxn.getMolecule(i).setFragment(mMol.isFragment());

				storeState();
				setReaction(rxn);
				ret = true;
			} else {
				showWarningMessage("No reaction on clipboard!");
			}
		}
		return ret;
	}

	private boolean pasteMolecule() {
		if (mClipboardHandler != null) {
			StereoMolecule mol = mClipboardHandler.pasteMolecule();
			if (addPastedOrDropped(mol, null))
				return true;

			showWarningMessage("No molecule on clipboard!");
			}
		return false;
		}

	public boolean addPastedOrDropped(StereoMolecule mol, GenericPoint p) {
		if (mol == null || mol.getAllAtoms() == 0)
			return false;

		if (mol.getAllBonds() != 0) {
			double avbl = getScaledAVBL();
			new GenericDepictor(mol).updateCoords(mCanvas.getDrawContext(),
					new GenericRectangle(0, 0, mCanvas.getCanvasWidth(), mCanvas.getCanvasHeight()),
					AbstractDepictor.cModeInflateToMaxAVBL + (int)avbl);

			while (atomCoordinatesCollide(mol, 0.2*avbl))
				mol.translateCoords(0.5*avbl, 0.5*avbl);
			}

		storeState();

		mol.removeAtomColors();
		mol.removeBondHiliting();

		boolean editorIsFragment = mMol.isFragment();
		if (mMol.getAllAtoms() == 0) {
			mol.copyMolecule(mMol);
			if (!mAllowFragmentChangeOnPasteOrDrop)
				mMol.setFragment(editorIsFragment);
			updateAndFireEvent(UPDATE_SCALE_COORDS);
		} else {
			if (p != null)
				mol.translateCoords(p.x - mCanvas.getCanvasWidth()/2, p.y - mCanvas.getCanvasHeight()/2);
			int originalAtoms = mMol.getAllAtoms();
			mMol.addMolecule(mol);
			if (!mAllowFragmentChangeOnPasteOrDrop)
				mMol.setFragment(editorIsFragment);
			for (int atom = 0; atom<mMol.getAllAtoms(); atom++)
				mMol.setAtomSelection(atom, atom>=originalAtoms);

			updateAndFireEvent(UPDATE_CHECK_COORDS);
			}

		return true;
		}

	private void useTemplate(String rxncode) {
		storeState();
		setReaction(ReactionEncoder.decode(rxncode, true));
		}

	private boolean atomCoordinatesCollide(StereoMolecule mol, double tolerance) {
		int count = 0;
		tolerance *= tolerance;
		for (int i=0; i<mol.getAllAtoms(); i++) {
			double x1 = mol.getAtomX(i);
			double y1 = mol.getAtomY(i);
			boolean found = false;
			for (int j=0; j<mMol.getAllAtoms(); j++) {
				double x2 = mMol.getAtomX(j);
				double y2 = mMol.getAtomY(j);
				double dx = x2 - x1;
				double dy = y2 - y1;
				if (dx * dx + dy * dy < tolerance) {
					found = true;
					break;
					}
				}
			if (found)
				count++;
			}
		return count == mol.getAllAtoms();
		}

	private void openReaction() {
		File rxnFile = mUIHelper.openChemistryFile(true);
		if (rxnFile != null) {
			try {
				Reaction reaction = null;

				if (FileHelper.getFileType(rxnFile.getName()) == FileHelper.cFileTypeRXN) {
					reaction = new RXNFileParser().getReaction(rxnFile);
				} else {
					RDFileParser rdfParser = new RDFileParser(rxnFile);
					if (rdfParser.isReactionNext())
						reaction = rdfParser.getNextReaction();
				}

				if (reaction != null) {
					for (int i = 0; i<reaction.getMolecules(); i++) {
						reaction.getMolecule(i).setFragment(mMol.isFragment());
					}
					storeState();
					setReaction(reaction);
				}
			} catch (Exception ex) {
			}
		}
	}

	private void showWarningMessage(String msg) {
		mWarningMessage = msg;
		mCanvas.repaint();
		new Thread(() -> {
			try {
				Thread.sleep(WARNING_MILLIS);
			} catch (InterruptedException ie) {
			}
			mWarningMessage = null;
			mCanvas.repaint();
		}).start();
	}

	private void eventHappened(GenericMouseEvent e) {
		if (e.getWhat() == GenericMouseEvent.MOUSE_PRESSED) {
			if (mCurrentHiliteAtom != -1 && mAtomKeyStrokeBuffer.length() != 0)
				expandAtomKeyStrokes(mAtomKeyStrokeBuffer.toString());

			mAtomKeyStrokeBuffer.setLength(0);

			if (e.isPopupTrigger()) {
				handlePopupTrigger(e.getX(), e.getY());
				return;
			}

			if (e.getButton() == 1) {
				if (e.getClickCount() == 2) {
					return;
				}

				mMouseIsDown = true;
				updateCursor();
				mousePressedButton1(e);
			}
		}

		if (e.getWhat() == GenericMouseEvent.MOUSE_RELEASED) {
			if (e.isPopupTrigger()) {
				handlePopupTrigger(e.getX(), e.getY());
				return;
			}

			if (e.getButton() == 1) {
				if (e.getClickCount() == 2) {
					handleDoubleClick(e.getX(), e.getY());
					return;
				}

				mMouseIsDown = false;
				updateCursor();
				mouseReleasedButton1();
			}
		}

		if (e.getWhat() == GenericMouseEvent.MOUSE_ENTERED) {
			mUIHelper.grabFocus();
			updateCursor();
		}

		if (e.getWhat() == GenericMouseEvent.MOUSE_MOVED) {
			mMouseIsDown = false;
			int x = e.getX();
			int y = e.getY();

			if (trackHiliting(x, y, false)) {
				mCanvas.repaint();
			}

			updateCursor();
		}

		if (e.getWhat() == GenericMouseEvent.MOUSE_DRAGGED) {
			mMouseIsDown = true;
			mX2 = e.getX();
			mY2 = e.getY();

			boolean repaintNeeded = trackHiliting(mX2, mY2, true);

			switch (mPendingRequest) {
				case cRequestNewChain:
					double lastX, lastY;
					if (mChainAtoms>0) {
						lastX = mChainAtomX[mChainAtoms - 1];
						lastY = mChainAtomY[mChainAtoms - 1];
					} else {
						lastX = 0;
						lastY = 0;
					}
					double avbl = getScaledAVBL();
					double s0 = (int)avbl;
					double s1 = (int)(0.866 * avbl);
					double s2 = (int)(0.5 * avbl);
					double dx = mX2 - mX1;
					double dy = mY2 - mY1;
					if (Math.abs(dy)>Math.abs(dx)) {
						mChainAtoms = (int)(2 * Math.abs(dy) / (s0 + s2));
						if (Math.abs(dy) % (s0 + s2)>s0) {
							mChainAtoms++;
						}
						mChainAtomX = new double[mChainAtoms];
						mChainAtomY = new double[mChainAtoms];
						if (mX2<mX1) {
							s1 = -s1;
						}
						if (mY2<mY1) {
							s0 = -s0;
							s2 = -s2;
						}
						for (int i = 0; i<mChainAtoms; i++) {
							mChainAtomX[i] = mX1 + ((i + 1) / 2) * s1;
							mChainAtomY[i] = mY1 + ((i + 1) / 2) * (s0 + s2);
							if ((i & 1) == 0) {
								mChainAtomY[i] += s0;
							}
						}
					} else {
						mChainAtoms = (int)(Math.abs(dx) / s1);
						mChainAtomX = new double[mChainAtoms];
						mChainAtomY = new double[mChainAtoms];
						if (mX2<mX1) {
							s1 = -s1;
						}
						if (mY2<mY1) {
							s2 = -s2;
						}
						for (int i = 0; i<mChainAtoms; i++) {
							mChainAtomX[i] = mX1 + (i + 1) * s1;
							mChainAtomY[i] = mY1;
							if ((i & 1) == 0) {
								mChainAtomY[i] += s2;
							}
						}
					}
					if (mChainAtoms>0) {
						mChainAtom = new int[mChainAtoms];
						for (int i = 0; i<mChainAtoms; i++) {
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
					if ((mX2 - mX1) * (mX2 - mX1) + (mY2 - mY1) * (mY2 - mY1)<MIN_BOND_LENGTH_SQUARE) {
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
					for (int i = 0; i<mMol.getAllAtoms(); i++) {
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
					for (int atom = 0; atom<mMol.getAllAtoms() && !selectedAtomsFound; atom++) {
						selectedAtomsFound = mMol.isSelectedAtom(atom);
					}
					boolean selectedObjectsFound = false;
					if (mDrawingObjectList != null) {
						for (int i = 0; i<mDrawingObjectList.size() && !selectedObjectsFound; i++) {
							selectedObjectsFound = mDrawingObjectList.get(i).isSelected();
						}
					}
					double magnification = (Math.abs(mY2 - mY1)<20) ? 1.0 : Math.exp((mY2 - mY1) / 100);
					double angleChange = (Math.abs(mX2 - mX1)<20) ? 0.0f : (mX2 - mX1) / 50;
					boolean selectedOnly = (selectedAtomsFound || selectedObjectsFound);
					if (mDrawingObjectList != null && (!selectedOnly || selectedObjectsFound)) {
						for (int i = 0; i<mDrawingObjectList.size(); i++) {
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
				mCanvas.repaint();
			}
		}
	}

	public void showHelpDialog() {
		mUIHelper.showHelpDialog("/html/editor/editor.html", "Structure Editor Help");
	}

	private void eventHappened(GenericKeyEvent e) {
		if (e.getWhat() == GenericKeyEvent.KEY_PRESSED) {
			if (e.getKey() == GenericKeyEvent.KEY_SHIFT) {
				mShiftIsDown = true;
				updateCursor();
			}
			if (e.getKey() == GenericKeyEvent.KEY_ALT) {
				mAltIsDown = true;
				updateCursor();
			}
			if (e.getKey() == GenericKeyEvent.KEY_CTRL) {
				updateCursor();
			}

			if (e.isMenuShortcut()) {
				if (e.getKey() == 'z') {
					restoreState();
					updateAndFireEvent(UPDATE_CHECK_VIEW);
				}
				else if (e.getKey() == 'c') {
					copy();
				}
				else if (e.getKey() == 'v') {
					paste();
				}
			} else if (e.getKey() == GenericKeyEvent.KEY_DELETE) {
				storeState();
				if (mCurrentTool == GenericEditorToolbar.cToolMapper) {
					boolean found = false;
					for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
						if (mMol.getAtomMapNo(atom) != 0) {
							mMol.setAtomMapNo(atom, 0, false);
							found = true;
						}
					}
					if (found) {
						updateAndFireEvent(UPDATE_REDRAW);
					}
				} else if (!deleteHilited()) {
					if (mMol.deleteSelectedAtoms()) {
						updateAndFireEvent(UPDATE_REDRAW);
					}
				}
			} else if (e.getKey() == GenericKeyEvent.KEY_HELP || (mCurrentHiliteAtom == -1 && e.getKey() == '?')) {
				showHelpDialog();
				return;
			} else if (mCurrentHiliteBond != -1) {
				int ch = e.getKey();
				if (ch == 'q' && mMol.isFragment()) {
					showBondQFDialog(mCurrentHiliteBond);
				} else if (ch == 'v') { // ChemDraw uses the same key
					if (mMol.addRingToBond(mCurrentHiliteBond, 3, false, getScaledAVBL()))
						updateAndFireEvent(UPDATE_CHECK_COORDS);
				} else if (ch>='4' && ch<='7') {
					if (mMol.addRingToBond(mCurrentHiliteBond, ch - '0', false, getScaledAVBL()))
						updateAndFireEvent(UPDATE_CHECK_COORDS);
				} else if (ch == 'a' || ch == 'b') {    // ChemDraw uses 'a', we use 'b' since a long time
					if (mMol.addRingToBond(mCurrentHiliteBond, 6, true, getScaledAVBL()))
						updateAndFireEvent(UPDATE_CHECK_COORDS);
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
					if (bondChanged)
						updateAndFireEvent(UPDATE_REDRAW);
				}
			} else if (mCurrentHiliteAtom != -1) {
				int ch = e.getKey();
				boolean isFirst = (mAtomKeyStrokeBuffer.length() == 0);
				if (isFirst)
					mFirstAtomKey = ch;
				else {
					if (mFirstAtomKey == 'l') { // if we don't want first 'l' to be a chlorine
						mAtomKeyStrokeBuffer.setLength(0);
						mAtomKeyStrokeBuffer.append('L');
						}
					mFirstAtomKey = -1;
					}

				if (isFirst && ch == 'l') { // if no chars are following, we interpret 'l' as chlorine analog to ChemDraw
					mAtomKeyStrokeBuffer.append("Cl");
					update(UPDATE_REDRAW);
				} else if (isFirst && (ch == '+' || ch == '-')) {
					storeState();
					if (mMol.changeAtomCharge(mCurrentHiliteAtom, ch == '+'))
						updateAndFireEvent(UPDATE_CHECK_COORDS);
				} else if (isFirst && ch == '.') {
					storeState();
					int newRadical = (mMol.getAtomRadical(mCurrentHiliteAtom) == Molecule.cAtomRadicalStateD) ?
							0 : Molecule.cAtomRadicalStateD;
					mMol.setAtomRadical(mCurrentHiliteAtom, newRadical);
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				} else if (isFirst && ch == ':') {
					storeState();
					int newRadical = (mMol.getAtomRadical(mCurrentHiliteAtom) == Molecule.cAtomRadicalStateT) ? Molecule.cAtomRadicalStateS
							: (mMol.getAtomRadical(mCurrentHiliteAtom) == Molecule.cAtomRadicalStateS) ? 0 : Molecule.cAtomRadicalStateT;
					mMol.setAtomRadical(mCurrentHiliteAtom, newRadical);
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				} else if (isFirst && ch == 'q' && mMol.isFragment()) {
					showAtomQFDialog(mCurrentHiliteAtom);

				} else if (isFirst && mMol.isFragment() && (ch == 'x' || ch == 'X')) {
					int[] list = { 9, 17, 35, 53 };
					mMol.setAtomList(mCurrentHiliteAtom, list);
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				} else if (isFirst && ch == '?') {
					storeState();
					if (mMol.changeAtom(mCurrentHiliteAtom, 0, 0, -1, 0)) {
						updateAndFireEvent(UPDATE_CHECK_COORDS);
					}
				} else if (isFirst && ch>48 && ch<=57) {
					if (mMol.getFreeValence(mCurrentHiliteAtom)>0) {
						storeState();
						int chainAtoms = ch - 47;
						int atom1 = mCurrentHiliteAtom;
						int hydrogenCount = mMol.getAllAtoms() - mMol.getAtoms();
						for (int i = 1; i<chainAtoms; i++) {
							suggestNewX2AndY2(atom1);
							int atom2 = mMol.addAtom(mX2, mY2);
							if (atom2 == -1) {
								break;
							}

							mMol.addBond(atom1, atom2);
							atom1 = atom2 - hydrogenCount;    // new atom was added behind all hydrogens and travels now to the front
							mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
						}
						updateAndFireEvent(UPDATE_CHECK_COORDS);
					}
				} else if (!isFirst && e.getKey() == GenericKeyEvent.KEY_ESCAPE) {
					mAtomKeyStrokeBuffer.setLength(0);
					update(UPDATE_REDRAW);
				} else if (!isFirst && e.getKey() == GenericKeyEvent.KEY_BACK_SPACE) {
					mAtomKeyStrokeBuffer.setLength(mAtomKeyStrokeBuffer.length() - 1);
					update(UPDATE_REDRAW);
				} else if ((ch>=65 && ch<=90)
						|| (ch>=97 && ch<=122)
						|| (ch>=48 && ch<=57)
						|| (ch == '-')) {
					mAtomKeyStrokeBuffer.append((char)ch);
					update(UPDATE_REDRAW);
				} else if (ch == '\n' || ch == '\r') {
					expandAtomKeyStrokes(mAtomKeyStrokeBuffer.toString());
				}
			} else if (mCurrentHiliteAtom == -1 && mCurrentHiliteBond == -1) {
				if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE | MODE_MULTIPLE_FRAGMENTS)) == 0) {
					int ch = e.getKey();
					if (ch == 'h')
						flip(true);
					if (ch == 'v')
						flip(false);
				}
			}
		}
		if (e.getWhat() == GenericKeyEvent.KEY_RELEASED) {
			if (e.getKey() == GenericKeyEvent.KEY_SHIFT) {
				mShiftIsDown = false;
				updateCursor();
			}
			if (e.getKey() == GenericKeyEvent.KEY_ALT) {
				mAltIsDown = false;
				updateCursor();
			}
			if (e.getKey() == GenericKeyEvent.KEY_CTRL) {
				updateCursor();
			}
		}
	}

	private boolean changeHighlightedBond(int type) {
		storeState();
		return mMol.changeBond(mCurrentHiliteBond, type);
		}

	private void handlePopupTrigger(int x, int y) {
		GenericPopupMenu popup = null;

		if (mClipboardHandler != null) {
			popup = mUIHelper.createPopupMenu(this);

			String item1 = analyseCopy(false) ? ITEM_COPY_REACTION : ITEM_COPY_STRUCTURE;
			popup.addItem(item1, null, mMol.getAllAtoms() != 0);

			String item3 = (StructureNameResolver.getInstance() == null) ? ITEM_PASTE_STRUCTURE : ITEM_PASTE_WITH_NAME;
			popup.addItem(item3, null, true);

			if ((mMode & MODE_REACTION) != 0) {
				popup.addItem(ITEM_PASTE_REACTION, null, true);
				popup.addSeparator();
				if (sReactionQueryTemplates != null && mMol.isFragment()) {
					boolean isSubMenu = false;
					for (String[] template: sReactionQueryTemplates) {
						if (TEMPLATE_SECTION_KEY.equals(template[0])) {
							if (isSubMenu)
								popup.endSubMenu();

							popup.startSubMenu("Use " + template[1] + " Template");
							isSubMenu = true;
							continue;
							}

						if (!isSubMenu) {
							popup.startSubMenu("Use Template");
							isSubMenu = true;
							}

						popup.addItem(template[0], ITEM_USE_TEMPLATE + template[1], true);
						}
					popup.endSubMenu();
					}
				}

			if ((mMode & MODE_REACTION) != 0)
				popup.addItem(ITEM_LOAD_REACTION, null, true);
			}

		if ((mMode & MODE_REACTION) != 0 && mCurrentTool == GenericEditorToolbar.cToolMapper) {
			if (popup == null)
				popup = mUIHelper.createPopupMenu(this);
			else
				popup.addSeparator();

			popup.addItem(ITEM_ADD_AUTO_MAPPING, null, true);
			popup.addItem(ITEM_REMOVE_MAPPING, null, true);
			}

		if (mCurrentTool == GenericEditorToolbar.cToolZoom) {
			if (popup == null)
				popup = mUIHelper.createPopupMenu(this);
			else
				popup.addSeparator();

			popup.addItem(ITEM_FLIP_HORIZONTALLY, null, true);
			popup.addItem(ITEM_FLIP_VERTICALLY, null, true);
			popup.addItem(ITEM_FLIP_ROTATE180, null, true);
			}

		if (mAtomColorSupported && mCurrentHiliteAtom != -1) {
			int atomColor = mMol.getAtomColor(mCurrentHiliteAtom);
			if (popup == null)
				popup = mUIHelper.createPopupMenu(this);
			else
				popup.addSeparator();

			popup.startSubMenu("Set Atom Color");
			popup.addRadioButtonItem("	  ", "atomColor" + mCurrentHiliteAtom + ":" + Molecule.cAtomColorNone, RGB_BLACK, atomColor == Molecule.cAtomColorNone);
			popup.addRadioButtonItem("	  ", "atomColor" + mCurrentHiliteAtom + ":" + Molecule.cAtomColorBlue, AbstractDepictor.COLOR_BLUE, atomColor == Molecule.cAtomColorBlue);
			popup.addRadioButtonItem("	  ", "atomColor" + mCurrentHiliteAtom + ":" + Molecule.cAtomColorDarkRed, AbstractDepictor.COLOR_DARK_RED, atomColor == Molecule.cAtomColorDarkRed);
			popup.addRadioButtonItem("	  ", "atomColor" + mCurrentHiliteAtom + ":" + Molecule.cAtomColorRed, AbstractDepictor.COLOR_RED, atomColor == Molecule.cAtomColorRed);
			popup.addRadioButtonItem("	  ", "atomColor" + mCurrentHiliteAtom + ":" + Molecule.cAtomColorDarkGreen, AbstractDepictor.COLOR_DARK_GREEN, atomColor == Molecule.cAtomColorDarkGreen);
			popup.addRadioButtonItem("	  ", "atomColor" + mCurrentHiliteAtom + ":" + Molecule.cAtomColorGreen, AbstractDepictor.COLOR_GREEN, atomColor == Molecule.cAtomColorGreen);
			popup.addRadioButtonItem("	  ", "atomColor" + mCurrentHiliteAtom + ":" + Molecule.cAtomColorMagenta, AbstractDepictor.COLOR_MAGENTA, atomColor == Molecule.cAtomColorMagenta);
			popup.addRadioButtonItem("	  ", "atomColor" + mCurrentHiliteAtom + ":" + Molecule.cAtomColorOrange, AbstractDepictor.COLOR_ORANGE, atomColor == Molecule.cAtomColorOrange);
			popup.endSubMenu();
			}

		if (popup == null)
			popup = mUIHelper.createPopupMenu(this);
		else
			popup.addSeparator();
		popup.addItem(ITEM_SHOW_HELP, null, true);

		if (System.getProperty("development") != null) {
			if (popup == null)
				popup = mUIHelper.createPopupMenu(this);
			else
				popup.addSeparator();

			popup.addItem(ITEM_SHOW_ATOM_BOND_NUMBERS, null, true);
			popup.addItem(ITEM_SHOW_SYMMETRY, null, true);
			popup.addItem(ITEM_SHOW_NORMAL, null, true);
			}

		if (popup != null)
			popup.show(x, y);
		}

	private void handleDoubleClick (int x, int y) {
		int atom = mMol.findAtom(x, y);
		int bond = mMol.findBond(x, y);

		if (mCurrentTool == GenericEditorToolbar.cToolLassoPointer) {
			if (mMol.isFragment()) {
				if (atom != -1) {
					showAtomQFDialog(atom);
				} else if (bond != -1) {
					showBondQFDialog(bond);
				} else if (mCurrentHiliteObject != null) {
					if (!mShiftIsDown) {
						for (int i = 0; i<mMol.getAllAtoms(); i++)
							mMol.setAtomSelection(i, false);
						for (AbstractDrawingObject ado:mDrawingObjectList)
							ado.setSelected(false);
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
						for (int i = 0; i<mMol.getAllAtoms(); i++) {
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
							for (int i = 0; i<mMol.getAllAtoms(); i++) {
								if (mFragmentNo[i] == fragment) {
									mMol.setAtomSelection(i, true);
								}
							}
						} else {
							int[] fragmentMember = mMol.getFragmentAtoms(rootAtom);
							for (int i = 0; i<fragmentMember.length; i++) {
								mMol.setAtomSelection(fragmentMember[i], true);
							}
						}
					} else {
						mCurrentHiliteObject.setSelected(true);
					}

					update(UPDATE_REDRAW);
				}
			}
		} else if (mCurrentTool == GenericEditorToolbar.cToolZoom) {
			int fragment = -2;
			if ((mMode & MODE_MULTIPLE_FRAGMENTS) != 0) {
				fragment = findFragment(x, y);
			}

			if (fragment != -1) {
				double minX = Integer.MAX_VALUE;
				double maxX = Integer.MIN_VALUE;
				for (int i = 0; i<mMol.getAllAtoms(); i++) {
					if (fragment == -2 || mFragmentNo[i] == fragment) {
						if (minX>mMol.getAtomX(i)) {
							minX = mMol.getAtomX(i);
						}
						if (maxX<mMol.getAtomX(i)) {
							maxX = mMol.getAtomX(i);
						}
					}
				}

				if (maxX>minX) {
					double centerX = (maxX + minX) / 2;
					for (int i = 0; i<mMol.getAllAtoms(); i++) {
						if (fragment == -2 || mFragmentNo[i] == fragment) {
							mMol.setAtomX(i, 2 * centerX - mMol.getAtomX(i));
						}
					}
					for (int i = 0; i<mMol.getAllBonds(); i++) {
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

				updateAndFireEvent(UPDATE_REDRAW);
			}
		} else if (mCurrentTool == GenericEditorToolbar.cToolCustomAtom) {
			mUIHelper.showMessage("To change current custom atom properties hold 'Ctrl'\nwhile clicking an atom with the left mouse button.");
		}
	}

	private void showAtomQFDialog(int atom) {
		if (mAllowQueryFeatures) {
			storeState();
			boolean showReactionHints = ((mMode & MODE_REACTION) != 0);
			AtomQueryFeatureDialogBuilder builder = new AtomQueryFeatureDialogBuilder(mUIHelper, mMol, atom, showReactionHints);
			if (builder.showDialog())
				updateAndFireEvent(UPDATE_REDRAW);
			}
		}

	private void showBondQFDialog(int bond) {
		if (mAllowQueryFeatures) {
			storeState();
			BondQueryFeatureDialogBuilder builder = new BondQueryFeatureDialogBuilder(mUIHelper, mMol, bond);
			if (builder.showDialog())
				updateAndFireEvent(UPDATE_REDRAW);
			}
		}

	public void showCustomAtomDialog(int atom) {
		storeState();
		CustomAtomDialogBuilder builder = (atom == -1) ?
				new CustomAtomDialogBuilder(mUIHelper, this, mCustomAtomicNo, mCustomAtomMass, mCustomAtomValence, mCustomAtomRadical, mCustomAtomLabel)
				: new CustomAtomDialogBuilder(mUIHelper, this, mMol, atom);
		if (builder.showDialog() && atom != -1)
			updateAndFireEvent(UPDATE_REDRAW);
		}

	private void mousePressedButton1(GenericMouseEvent gme) {
		mX1 = gme.getX();
		mY1 = gme.getY();

		switch (mCurrentTool) {
			case GenericEditorToolbar.cToolZoom:
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
			case GenericEditorToolbar.cToolLassoPointer:
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
					for (int i = 0; i<mMol.getAllAtoms(); i++) {
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

				mIsAddingToSelection = gme.isShiftDown();
				for (int i = 0; i<mMol.getAllAtoms(); i++) {
					mIsSelectedAtom[i] = mMol.isSelectedAtom(i);
				}
				if (mDrawingObjectList != null) {
					for (int i = 0; i<mDrawingObjectList.size(); i++) {
						mIsSelectedObject[i] = mDrawingObjectList.get(i).isSelected();
					}
				}

				if (gme.isAltDown()) {
					mPendingRequest = cRequestSelectRect;
				} else {
					mLassoRegion = new GenericPolygon();
					mLassoRegion.addPoint(mX1, mY1);
					mLassoRegion.addPoint(mX1, mY1);
					mPendingRequest = cRequestLassoSelect;
				}
				break;
			case GenericEditorToolbar.cToolDelete:
				storeState();
				deleteAt(mX1, mY1);
				break;
			case GenericEditorToolbar.cToolUnknownParity:
				int theAtom = mMol.findAtom(mX1, mY1);
				if (theAtom != -1) {
					storeState();
					mMol.setAtomConfigurationUnknown(theAtom, !mMol.isAtomConfigurationUnknown(theAtom));
					updateAndFireEvent(UPDATE_REDRAW);
				}
				break;
			case GenericEditorToolbar.cToolESRAbs:
			case GenericEditorToolbar.cToolESRAnd:
			case GenericEditorToolbar.cToolESROr:
				if (mCurrentHiliteBond != -1 && qualifiesForESR(mCurrentHiliteBond)) {
					storeState();
					setESRInfo(mCurrentHiliteBond, (mCurrentTool == GenericEditorToolbar.cToolESRAbs) ?
							Molecule.cESRTypeAbs :
							(mCurrentTool == GenericEditorToolbar.cToolESRAnd) ?
									Molecule.cESRTypeAnd : Molecule.cESRTypeOr);
					updateAndFireEvent(UPDATE_REDRAW);
				}
				break;
			case GenericEditorToolbar.cToolStdBond:
			case GenericEditorToolbar.cToolUpBond:
			case GenericEditorToolbar.cToolDownBond:
				mAtom1 = mMol.findAtom(mX1, mY1);
				if (mAtom1 == -1) {
					int bond = mMol.findBond(mX1, mY1);
					if (bond != -1) {
						storeState();
						int bondType = Molecule.cBondTypeIncreaseOrder;
						if (mCurrentTool == GenericEditorToolbar.cToolUpBond) {
							bondType = Molecule.cBondTypeUp;
						} else if (mCurrentTool == GenericEditorToolbar.cToolDownBond) {
							bondType = Molecule.cBondTypeDown;
						}
						if (mMol.changeBond(bond, bondType))
							updateAndFireEvent(UPDATE_REDRAW);
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
				mCanvas.repaint();
				break;
			case GenericEditorToolbar.cToolChain:
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
			case GenericEditorToolbar.cTool3Ring:
				storeState();
				if (mMol.addRing(mX1, mY1, 3, false, getScaledAVBL()))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cTool4Ring:
				storeState();
				if (mMol.addRing(mX1, mY1, 4, false, getScaledAVBL()))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cTool5Ring:
				storeState();
				if (mMol.addRing(mX1, mY1, 5, false, getScaledAVBL()))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cTool6Ring:
				storeState();
				if (mMol.addRing(mX1, mY1, 6, false, getScaledAVBL()))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cTool7Ring:
				storeState();
				if (mMol.addRing(mX1, mY1, 7, false, getScaledAVBL()))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAromRing:
				storeState();
				if (mMol.addRing(mX1, mY1, 6, true, getScaledAVBL()))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolPosCharge:
				storeState();
				if (mMol.changeAtomCharge(mX1, mY1, true))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolNegCharge:
				storeState();
				if (mMol.changeAtomCharge(mX1, mY1, false))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomH:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 1, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomC:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 6, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomN:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 7, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomO:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 8, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomSi:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 14, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomP:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 15, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomS:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 16, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomF:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 9, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomCl:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 17, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomBr:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 35, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolAtomI:
				storeState();
				if (mMol.addOrChangeAtom(mX1, mY1, 53, 0, -1, 0, null))
					updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case GenericEditorToolbar.cToolCustomAtom:
				if (gme.isControlDown()) {
					int atom = mMol.findAtom(mX1, mY1);
					if (atom != -1) {
						showCustomAtomDialog(atom);
					}
				} else {
					storeState();
					if (mMol.addOrChangeAtom(mX1, mY1, mCustomAtomicNo, mCustomAtomMass, mCustomAtomValence, mCustomAtomRadical, mCustomAtomLabel))
						updateAndFireEvent(UPDATE_CHECK_COORDS);
				}
				break;
			case GenericEditorToolbar.cToolMapper:
				mAtom1 = mMol.findAtom(mX1, mY1);
				if (mAtom1 != -1 && mAtom1<mMol.getAtoms()) {
					mX1 = mMol.getAtomX(mAtom1);
					mY1 = mMol.getAtomY(mAtom1);
					mPendingRequest = cRequestMapAtoms;
				}
				break;
			case GenericEditorToolbar.cToolText:
				TextDrawingObject object = null;
				if (mCurrentHiliteObject == null) {
					object = new TextDrawingObject();
					object.setCoordinates(mX1, mY1);
					mDrawingObjectList.add(object);
				} else if (mCurrentHiliteObject instanceof TextDrawingObject) {
					object = (TextDrawingObject)mCurrentHiliteObject;
				}
				editTextObject(object);
				storeState();
				update(UPDATE_CHECK_COORDS);
				break;
		}
	}

	private void mouseReleasedButton1() {
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
					if (mCurrentTool == GenericEditorToolbar.cToolUpBond) {
						bondType = Molecule.cBondTypeUp;
					}
					if (mCurrentTool == GenericEditorToolbar.cToolDownBond) {
						bondType = Molecule.cBondTypeDown;
					}
					mMol.addOrChangeBond(mAtom1, stopAtom, bondType);
				}

				updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case cRequestNewChain:
				storeState();
				if (mChainAtoms>0) {
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

				if (mChainAtoms>1) {
					for (int i = 1; i<mChainAtoms; i++) {
						if (mChainAtom[i] == -1) {
							mChainAtom[i] = mMol.addAtom(mChainAtomX[i],
									mChainAtomY[i]);
						}
						if (mChainAtom[i] != -1) {
							mMol.addBond(mChainAtom[i - 1], mChainAtom[i]);
						}
					}
				}
				updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case cRequestMoveSingle:
			case cRequestMoveSelected:
			case cRequestZoomAndRotate:
				updateAndFireEvent(UPDATE_CHECK_COORDS);
				break;
			case cRequestMoveObject:
				update(UPDATE_CHECK_COORDS);
				break;
			case cRequestLassoSelect:
			case cRequestSelectRect:
				boolean selectionChanged = false;
				for (int i = 0; i<mMol.getAllAtoms(); i++) {
					if (mIsSelectedAtom[i] != mMol.isSelectedAtom(i)) {
						selectionChanged = true;
						break;
					}
				}
				mCanvas.repaint();
				if (selectionChanged)
					fireEventLater(new EditorEvent(this, EditorEvent.WHAT_SELECTION_CHANGED, true));
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
						for (int atom = 0; atom<mMol.getAtoms(); atom++) {
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
						for (int atom = 0; atom<mMol.getAtoms(); atom++) {
							if (mMol.getAtomMapNo(atom) == mapNo) {
								mMol.setAtomMapNo(atom, 0, false);
							}
						}
					} else {
						// remove old mapping numbers of atom1 and atom2
						int mapNoAtom2 = mMol.getAtomMapNo(atom2);
						for (int atom = 0; atom<mMol.getAtoms(); atom++) {
							if (mMol.getAtomMapNo(atom) == mapNoAtom1
									|| mMol.getAtomMapNo(atom) == mapNoAtom2) {
								mMol.setAtomMapNo(atom, 0, false);
							}
						}
						int freeMapNo = 1;
						for (int atom = 0; atom<mMol.getAtoms(); atom++) {
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

				if (mapNoChanged)
					updateAndFireEvent(Math.max(mUpdateMode, UPDATE_REDRAW));

				mCanvas.repaint();
				break;
		}
	}

	/**
	 * Takes the manually mapped atom mapping numbers from the display molecule,
	 * copies them into the current fragments, creates a reaction from these,
	 * uses the MCS-mapper to map the reaction and returns the rxn's mapping
	 * into the display molecule.
	 */
	private void autoMapReaction () {
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
		TreeMap<Integer, Integer> oldToNewMapNo = new TreeMap<>();
		int nextMapNo = 1;

		final int fakeAtomMassBase = 512;

		// Mark the manually mapped atoms such that the mapper uses them first priority and
		// to be able to re-assign them later as manually mapped.
		int[] fragmentAtom = new int[mFragment.length];
		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
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
			for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
				int fragment = mFragmentNo[atom];
				boolean hasFakeAtomMass = (mFragment[fragment].getAtomMass(fragmentAtom[fragment])>fakeAtomMassBase);
				if (hasFakeAtomMass) {
					// rescue new mapNo
					int newMapNo = mFragment[fragment].getAtomMass(fragmentAtom[fragment]) - fakeAtomMassBase;

					// repair fake atom mass
					mFragment[fragment].setAtomMass(fragmentAtom[fragment], mMol.getAtomMass(atom));

					mMol.setAtomMapNo(atom, newMapNo, false);
					mFragment[fragment].setAtomMapNo(fragmentAtom[fragment], newMapNo, false);
				} else {
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
		} else {
			// restore original atom masses in fragments and copy molecule's mapping number into fragments
			fragmentAtom = new int[mFragment.length];
			for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
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
	private boolean qualifiesForESR ( int stereoBond){
		return mMol.isStereoBond(stereoBond) && (getESRAtom(stereoBond) != -1 || getESRBond(stereoBond) != -1);
	}

	/**
	 * Locates the stereo center with parity 1 or 2 that is defined by the stereo bond.
	 *
	 * @param stereoBond
	 * @return stereo center atom or -1 if no stereo center found
	 */
	private int getESRAtom ( int stereoBond){
		int atom = mMol.getBondAtom(0, stereoBond);
		if (mMol.getAtomParity(atom) != Molecule.cAtomParityNone) {
			return (mMol.isAtomParityPseudo(atom)
					|| (mMol.getAtomParity(atom) != Molecule.cAtomParity1
					&& mMol.getAtomParity(atom) != Molecule.cAtomParity2)) ? -1 : atom;
		}
		if (mMol.getAtomPi(atom) == 1) {
			for (int i = 0; i<mMol.getConnAtoms(atom); i++) {
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

	private int getESRBond ( int stereoBond)
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
	private void setESRInfo ( int stereoBond, int type)
	{
		int group = -1;

		int atom = getESRAtom(stereoBond);
		int bond = (atom == -1) ? getESRBond(stereoBond) : -1;

		// if type requires group information (type==And or type==Or)
		if (type != Molecule.cESRTypeAbs) {
			int maxGroup = -1;
			for (int i = 0; i<mMol.getAtoms(); i++) {
				if (i != atom
						&& mMol.getAtomESRType(i) == type
						&& (!mMol.isSelectedBond(stereoBond) || !mMol.isSelectedAtom(i))) {
					int grp = mMol.getAtomESRGroup(i);
					if (maxGroup<grp) {
						maxGroup = grp;
					}
				}
			}
			for (int i = 0; i<mMol.getBonds(); i++) {
				if (i != bond
						&& mMol.getBondESRType(i) == type
						&& (!mMol.isSelectedBond(stereoBond) || !mMol.isSelectedBond(i))) {
					int grp = mMol.getBondESRGroup(i);
					if (maxGroup<grp) {
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
					for (int i = 0; i<mMol.getAtoms(); i++) {
						if (i != atom && mMol.isSelectedAtom(i)
								&& mMol.getAtomESRType(i) == type
								&& mMol.getAtomESRGroup(i) != group) {
							selectedShareOneGroup = false;
							break;
						}
					}
					for (int i = 0; i<mMol.getBonds(); i++) {
						if (i != bond && mMol.isSelectedBond(i)
								&& mMol.getBondESRType(i) == type
								&& mMol.getBondESRGroup(i) != group) {
							selectedShareOneGroup = false;
							break;
						}
					}
					if (selectedShareOneGroup) {
						if (group<=maxGroup) {
							group++;
							if (group == Molecule.cESRMaxGroups) {
								group = 0;
							}
						} else {
							group = 0;
						}
					}
				} else {
					if (group<=maxGroup) {
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
			for (int i = 0; i<mMol.getBonds(); i++) {
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

	private int findFragment ( double x, double y)
	{
		int fragment = -1;
		double minDistance = Double.MAX_VALUE;
		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
			double dx = mX1 - mMol.getAtomX(atom);
			double dy = mY1 - mMol.getAtomY(atom);
			double distance = Math.sqrt(dx * dx + dy * dy);
			if (distance<FRAGMENT_MAX_CLICK_DISTANCE
					&& minDistance>distance) {
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
	private void suggestNewX2AndY2 ( int atom)
	{
		double newAngle = Math.PI * 2 / 3;
		if (atom != -1) {
			double angle[] = new double[MAX_CONNATOMS + 1];
			for (int i = 0; i<mMol.getAllConnAtomsPlusMetalBonds(atom); i++) {
				angle[i] = mMol.getBondAngle(atom, mMol.getConnAtom(atom, i));
			}

			if (mMol.getAllConnAtomsPlusMetalBonds(atom) == 1) {
				if (angle[0]<-Math.PI * 5 / 6) {
					newAngle = Math.PI / 3;
				} else if (angle[0]<-Math.PI / 2) {
					newAngle = Math.PI * 2 / 3;
				} else if (angle[0]<-Math.PI / 6) {
					newAngle = Math.PI / 3;
				} else if (angle[0]<0.0) {
					newAngle = Math.PI * 2 / 3;
				} else if (angle[0]<Math.PI / 6) {
					newAngle = -Math.PI * 2 / 3;
				} else if (angle[0]<Math.PI / 2) {
					newAngle = -Math.PI / 3;
				} else if (angle[0]<Math.PI * 5 / 6) {
					newAngle = -Math.PI * 2 / 3;
				} else {
					newAngle = -Math.PI / 3;
				}
			} else {
				for (int i = mMol.getAllConnAtomsPlusMetalBonds(atom) - 1; i>0; i--) {    // bubble sort
					for (int j = 0; j<i; j++) {
						if (angle[j]>angle[j + 1]) {
							double temp = angle[j];
							angle[j] = angle[j + 1];
							angle[j + 1] = temp;
						}
					}
				}
				angle[mMol.getAllConnAtomsPlusMetalBonds(atom)] = angle[0] + Math.PI * 2;

				int largestNo = 0;
				double largestDiff = 0.0;
				for (int i = 0; i<mMol.getAllConnAtomsPlusMetalBonds(atom); i++) {
					double angleDiff = angle[i + 1] - angle[i];
					if (largestDiff<angleDiff) {
						largestDiff = angleDiff;
						largestNo = i;
					}
				}
				newAngle = (angle[largestNo] + angle[largestNo + 1]) / 2;
			}
		}
		double avbl = getScaledAVBL();
		mX2 = ((atom == -1) ? mX1 : mMol.getAtomX(atom)) + avbl * (float)Math.sin(newAngle);
		mY2 = ((atom == -1) ? mY1 : mMol.getAtomY(atom)) + avbl * (float)Math.cos(newAngle);
	}

	private boolean areAtomsMappingCompatible ( int atom1, int atom2){
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
				for (int i = 0; i<atomList1.length; i++) {
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

	private boolean trackHiliting ( double x, double y, boolean isDragging){
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
			if (mCurrentTool == GenericEditorToolbar.cToolESRAbs
					|| mCurrentTool == GenericEditorToolbar.cToolESRAnd
					|| mCurrentTool == GenericEditorToolbar.cToolESROr) {
				theBond = mMol.getStereoBond(theAtom);
				theAtom = -1;
			} else if (mCurrentTool == GenericEditorToolbar.cToolMapper
					&& theAtom>=mMol.getAtoms()) {
				theAtom = -1;
			}
		}

		if (theBond == -1
				&& theAtom == -1
				&& mCurrentTool != GenericEditorToolbar.cToolChain
				&& mCurrentTool != GenericEditorToolbar.cToolMapper
				&& mCurrentTool != GenericEditorToolbar.cToolUnknownParity
				&& mCurrentTool != GenericEditorToolbar.cToolPosCharge
				&& mCurrentTool != GenericEditorToolbar.cToolNegCharge
				&& mCurrentTool != GenericEditorToolbar.cToolAtomH
				&& mCurrentTool != GenericEditorToolbar.cToolAtomC
				&& mCurrentTool != GenericEditorToolbar.cToolAtomN
				&& mCurrentTool != GenericEditorToolbar.cToolAtomO
				&& mCurrentTool != GenericEditorToolbar.cToolAtomSi
				&& mCurrentTool != GenericEditorToolbar.cToolAtomP
				&& mCurrentTool != GenericEditorToolbar.cToolAtomS
				&& mCurrentTool != GenericEditorToolbar.cToolAtomF
				&& mCurrentTool != GenericEditorToolbar.cToolAtomCl
				&& mCurrentTool != GenericEditorToolbar.cToolAtomBr
				&& mCurrentTool != GenericEditorToolbar.cToolAtomI
				&& mCurrentTool != GenericEditorToolbar.cToolCustomAtom) {
			theBond = mMol.findBond(x, y);
		}

		if (theBond != -1
				&& (mCurrentTool == GenericEditorToolbar.cToolESRAbs
				|| mCurrentTool == GenericEditorToolbar.cToolESRAnd
				|| mCurrentTool == GenericEditorToolbar.cToolESROr)
				&& !qualifiesForESR(theBond)) {
			theBond = -1;
		}

		// don't change object hiliting while dragging
		AbstractDrawingObject hiliteObject = mCurrentHiliteObject;
		if (!isDragging && mDrawingObjectList != null) {
			hiliteObject = null;
			if (theAtom == -1 && theBond == -1
					&& (mCurrentTool == GenericEditorToolbar.cToolLassoPointer
					|| mCurrentTool == GenericEditorToolbar.cToolDelete
					|| mCurrentTool == GenericEditorToolbar.cToolText)) {
				for (AbstractDrawingObject theObject : mDrawingObjectList) {
					if (mCurrentTool == GenericEditorToolbar.cToolLassoPointer
							|| (mCurrentTool == GenericEditorToolbar.cToolDelete && !(theObject instanceof ReactionArrow))
							|| (mCurrentTool == GenericEditorToolbar.cToolText && theObject instanceof TextDrawingObject)) {
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
			fireEventLater(new EditorEvent(this, EditorEvent.WHAT_HILITE_ATOM_CHANGED, true));
		}
		if (mCurrentHiliteBond != theBond) {
			mCurrentHiliteBond = theBond;
			fireEventLater(new EditorEvent(this, EditorEvent.WHAT_HILITE_BOND_CHANGED, true));
		}
		mCurrentHiliteObject = hiliteObject;

		return repaintNeeded;
	}

	private int getAtomKeyStrokeValidity(String s){
		if (Molecule.getAtomicNoFromLabel(s, mAllowedPseudoAtoms) != 0)
			return KEY_IS_ATOM_LABEL;
		mAtomKeyStrokeSuggestion = NamedSubstituents.identify(s);
		if (mAtomKeyStrokeSuggestion == null)
			return isValidAtomKeyStrokeStart(s) ? KEY_IS_VALID_START : KEY_IS_INVALID;
		if (mAtomKeyStrokeSuggestion.length() == 0)
			return KEY_IS_VALID_START;
		else
			return KEY_IS_SUBSTITUENT;
	}

	/**
	 * @param s
	 * @return true if adding one or more chars may still create a valid key stroke sequence
	 */
	private boolean isValidAtomKeyStrokeStart(String s){
		if (s.length()<3)
			for (int i=1; i<Molecule.cAtomLabel.length; i++)
				if (Molecule.cAtomLabel[i].startsWith(s))
					return true;

		return false;
	}

	private void expandAtomKeyStrokes(String keyStrokes){
		mAtomKeyStrokeBuffer.setLength(0);

		int atomicNo = Molecule.getAtomicNoFromLabel(keyStrokes, mAllowedPseudoAtoms);
		if (atomicNo != 0) {
			storeState();
			if (mMol.changeAtom(mCurrentHiliteAtom, atomicNo, 0, -1, 0)) {
				updateAndFireEvent(UPDATE_CHECK_COORDS);
				return;
			}
		}

		if (mAtomKeyStrokeSuggestion != null && mAtomKeyStrokeSuggestion.length() != 0)
			keyStrokes = mAtomKeyStrokeSuggestion;

		StereoMolecule substituent = NamedSubstituents.getSubstituent(keyStrokes);
		if (substituent != null) {
			storeState();

			// Copy the the fragment containing the attachment point into a new molecule.
			// Then attach the substituent, create new atom coordinates for the substituent,
			// while retaining coordinates of the fragment.
			StereoMolecule fragment = new StereoMolecule();
			fragment.addFragment(mMol, mCurrentHiliteAtom, null);
			double sourceAVBL = fragment.getAverageBondLength(mUIScaling * Molecule.cDefaultAVBL);
			int firstAtomInFragment = fragment.getAllAtoms();
			for (int atom = 0; atom<fragment.getAllAtoms(); atom++)
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
			for (int i = 0; i<substituentAtoms; i++) {
				mMol.setAtomX(firstAtomInMol + i, sourceAVBL * fragment.getAtomX(firstAtomInFragment + i) + dx);
				mMol.setAtomY(firstAtomInMol + i, sourceAVBL * fragment.getAtomY(firstAtomInFragment + i) + dy);
			}
			mMol.setStereoBondsFromParity();

			updateAndFireEvent(UPDATE_CHECK_COORDS);
		}
	}

	public int getAllowedPseudoAtoms() {
		return mAllowedPseudoAtoms;
	}

	public void setAllowedPseudoAtoms(int apa) {
		mAllowedPseudoAtoms = apa;
	}

	private AbstractDrawingObject findDrawingObject ( double x, double y, String type,boolean forDeletion){
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

	private void editTextObject (TextDrawingObject object){
		new TextDrawingObjectDialogBuilder(mUIHelper, object);

		boolean nonWhiteSpaceFound = false;
		for (int i = 0; i<object.getText().length(); i++) {
			if (!Character.isWhitespace(object.getText().charAt(i))) {
				nonWhiteSpaceFound = true;
				break;
			}
		}

		if (!nonWhiteSpaceFound)
			mDrawingObjectList.remove(object);

		mCanvas.repaint();
	}

	private boolean shareSameReactionSide ( int atom1, int atom2){
		ReactionArrow arrow = (ReactionArrow)mDrawingObjectList.get(0);
		return !(arrow.isOnProductSide(mMol.getAtomX(atom1), mMol.getAtomY(atom1))
				^ arrow.isOnProductSide(mMol.getAtomX(atom2), mMol.getAtomY(atom2)));
	}

	protected void restoreState () {
		if (mUndoMol == null) {
			return;
		}
		mUndoMol.copyMolecule(mMol);
		mDrawingObjectList = (mUndoDrawingObjectList == null) ?
				null : new DrawingObjectList(mUndoDrawingObjectList);
	}

	public void storeState ()
	{
		if (mUndoMol == null) {
			mUndoMol = new Molecule();
		}
		mMol.copyMolecule(mUndoMol);

		mUndoDrawingObjectList = (mDrawingObjectList == null) ?
				null : new DrawingObjectList(mDrawingObjectList);
	}

	private boolean deleteHilited () {
		if (mCurrentHiliteAtom != -1) {
			mMol.deleteAtom(mCurrentHiliteAtom);
			mCurrentHiliteAtom = -1;
			updateAndFireEvent(UPDATE_REDRAW);
			return true;
		}

		if (mCurrentHiliteBond != -1) {
			mMol.deleteBondAndSurrounding(mCurrentHiliteBond);
			mCurrentHiliteBond = -1;
			updateAndFireEvent(UPDATE_REDRAW);
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

	private boolean deleteAt ( double x, double y){
		if (mMol.deleteAtomOrBond(x, y)) {
			updateAndFireEvent(UPDATE_REDRAW);
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

	private void duplicateSelected () {
		int atomCount = 0;
		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
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
		for (int atom = 0; atom<originalAtoms; atom++) {
			if (mMol.isSelectedAtom(atom)) {
				int newAtom = mMol.getAllAtoms();
				mX[newAtom] = mX[atom];
				mY[newAtom] = mY[atom];
				atomMap[atom] = newAtom;
				mMol.copyAtom(mMol, atom, esrGroupCountAND, esrGroupCountOR);
			}
		}
		for (int bond = 0; bond<originalBonds; bond++) {
			if (mMol.isSelectedBond(bond)) {
				mMol.copyBond(mMol, bond, esrGroupCountAND, esrGroupCountOR, atomMap, false);
			}
		}
		for (int atom = 0; atom<originalAtoms; atom++) {
			mMol.setAtomSelection(atom, false);
		}
		for (int atom = originalAtoms; atom<mMol.getAllAtoms(); atom++) {
			mMol.setAtomMapNo(atom, 0, false);
		}

		if (mDrawingObjectList != null) {
			for (int i = mDrawingObjectList.size() - 1; i>=0; i--) {
				AbstractDrawingObject object = mDrawingObjectList.get(i);
				if (object.isSelected() && !(object instanceof ReactionArrow)) {
					mDrawingObjectList.add(object.clone());
				}
			}
		}
	}

	private void fireEventLater(EditorEvent e) {
		final int what = e.getWhat();
		if ((what & mEventsScheduled) == 0) {
			mUIHelper.runLater(() -> {
				mEventsScheduled &= ~what;
				for (GenericEventListener<EditorEvent> l : mListeners)
					l.eventHappened(e);
			} );
			mEventsScheduled |= what;
		}
	}

	/**
	 * Redraws the molecule(s) or the reaction after scaling coordinates.
	 * Then analyses fragment membership and recreate individual molecules, reaction, or markush structure
	 * Then, fires molecule change events with userChange=false, i.e. external change.
	 */
	public void moleculeChanged() {
		update(UPDATE_SCALE_COORDS);
		fireEventLater(new EditorEvent(this, EditorEvent.WHAT_MOLECULE_CHANGED, false));
	}

	/**
	 * Redraws the molecule(s) or the reaction after scaling coordinates.
	 * Then analyses fragment membership and recreate individual molecules, reaction, or markush structure
	 * Then, fires molecule change events with userChange=false, i.e. external change.
	 * @param updateMode
	 */
	private void updateAndFireEvent(int updateMode) {
		update(updateMode);
		fireEventLater(new EditorEvent(this, EditorEvent.WHAT_MOLECULE_CHANGED, true));
	}

	public StereoMolecule getMolecule ()
	{
		return mMol;
	}

	public void setMolecule (StereoMolecule theMolecule){
		if (mMol == theMolecule) {
			return;
		}
		storeState();
		mMol = theMolecule;
		mMode = 0;
		mDrawingObjectList = null;
		moleculeChanged();
	}

	public StereoMolecule[] getFragments () {
		return mFragment;
	}

	public void setFragments(StereoMolecule[]fragment) {
		mMol.clear();
		mFragment = fragment;
		for (int i = 0; i<fragment.length; i++) {
			mMol.addMolecule(mFragment[i]);
		}
		storeState();

		mFragmentNo = new int[mMol.getAllAtoms()];
		for (int atom = 0, f = 0; f<mFragment.length; f++) {
			for (int j = 0; j<mFragment[f].getAllAtoms(); j++) {
				mFragmentNo[atom++] = f;
			}
		}

		mMode = MODE_MULTIPLE_FRAGMENTS;
		mDrawingObjectList = null;
		update(UPDATE_SCALE_COORDS_USE_FRAGMENTS);
		fireEventLater(new EditorEvent(this, EditorEvent.WHAT_MOLECULE_CHANGED, false));
	}

	/**
	 * @return mapped reaction with absolute coordinates, but without drawing objects
	 */
	public Reaction getReaction () {
		if ((mMode & MODE_REACTION) == 0 || mFragment == null) {
			return null;
		}

		Reaction rxn = new Reaction();
		for (int i = 0; i<mFragment.length; i++) {
			if (i<mReactantCount) {
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
	public Reaction getReactionAndDrawings () {
		Reaction rxn = getReaction();
		if (rxn != null) {
			rxn.setDrawingObjects(getDrawingObjects());
		}
		return rxn;
	}

	public void setReaction (Reaction rxn) {
		mMol.clear();
		mFragment = new StereoMolecule[rxn.getMolecules()];
		mReactantCount = rxn.getReactants();
		for (int i = 0; i<rxn.getMolecules(); i++) {
			mFragment[i] = rxn.getMolecule(i);
			mMol.addMolecule(mFragment[i]);
		}
		mMol.setFragment(rxn.isFragment());
		storeState();

		mFragmentNo = new int[mMol.getAllAtoms()];
		for (int atom = 0, f = 0; f<mFragment.length; f++) {
			for (int j = 0; j<mFragment[f].getAllAtoms(); j++) {
				mFragmentNo[atom++] = f;
			}
		}

		mDrawingObjectList = new DrawingObjectList();

		mMode = MODE_MULTIPLE_FRAGMENTS | MODE_REACTION;
		update(UPDATE_SCALE_COORDS_USE_FRAGMENTS);
		fireEventLater(new EditorEvent(this, EditorEvent.WHAT_MOLECULE_CHANGED, false));
	}

	public MarkushStructure getMarkushStructure () {
		if ((mMode & MODE_MARKUSH_STRUCTURE) == 0) {
			return null;
		}

		MarkushStructure markush = new MarkushStructure();
		for (int i = 0; i<mFragment.length; i++) {
			if (i<mReactantCount) {
				markush.addCore(mFragment[i]);
			} else {
				markush.addRGroup(mFragment[i]);
			}
		}
		return markush;
	}

	public void setMarkushStructure (MarkushStructure markush){
		mMol.clear();
		mDrawingObjectList = null;
		mFragment = new StereoMolecule[markush.getCoreCount() + markush.getRGroupCount()];
		mReactantCount = markush.getCoreCount();
		boolean isFragment = false;
		for (int i = 0; i<markush.getCoreCount() + markush.getRGroupCount(); i++) {
			mFragment[i] = (i<markush.getCoreCount()) ? markush.getCoreStructure(i)
					: markush.getRGroup(i - markush.getCoreCount());
			isFragment |= mFragment[i].isFragment();
			mMol.addMolecule(mFragment[i]);
		}
		mMol.setFragment(isFragment);
		storeState();

		mFragmentNo = new int[mMol.getAllAtoms()];
		for (int atom = 0, f = 0; f<mFragment.length; f++) {
			for (int j = 0; j<mFragment[f].getAllAtoms(); j++) {
				mFragmentNo[atom++] = f;
			}
		}

		mMode = MODE_MULTIPLE_FRAGMENTS | MODE_MARKUSH_STRUCTURE;
		update(UPDATE_SCALE_COORDS_USE_FRAGMENTS);
		fireEventLater(new EditorEvent(this, EditorEvent.WHAT_MOLECULE_CHANGED, false));
	}

	public int getDisplayMode () {
		return mDisplayMode;
	}

	public void setDisplayMode ( int dMode){
		mDisplayMode = dMode;
		update(UPDATE_REDRAW);
	}

	public void setTextSizeFactor(double factor) {
		mTextSizeFactor = factor;
		update(UPDATE_REDRAW);
	}

	/**
	 * If set to false then any query features will be removed from the molecule
	 * and any functionality that allows to define atom- or bond-query features
	 * won't be available. This feature is only relevant if the molecule is a fragment.
	 * @param allow if false, then query feature editing is not allowed, even for molecules being a fragment
	 */
	public void setAllowQueryFeatures (boolean allow){
		if (mAllowQueryFeatures != allow) {
			mAllowQueryFeatures = allow;
			if (!allow)
				mMol.removeQueryFeatures();
		}
	}

	/**
	 * In case a molecule/reaction is pasted into or dropped onto the editor area,
	 * then the fragment state of the editor area stays unchanged even if it is empty.
	 * If this method has been called with parameter 'true' then the behaviour changes as follows:
	 * If a query molecule/reaction (fragment==true) is pasted/dropped, then the editor's
	 * fragment state is set to true. If a non-fragment molecule/reaction is pasted/dropped
	 * then the editor's fragment state stays unchanged unless the editor area is empty,
	 * in which case the fragment state of the dropped object will be retained.
	 * @param b
	 */
	public void setAllowFragmentChangeOnPasteOrDrop(boolean b) {
		mAllowFragmentChangeOnPasteOrDrop = b;
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
	public void setAtomText (String[]atomText){
		mAtomText = atomText;
	}

	public DrawingObjectList getDrawingObjects () {
		return mDrawingObjectList;
	}

	public void setDrawingObjects (DrawingObjectList drawingObjectList){
		mDrawingObjectList = drawingObjectList;
		storeState();
		update(UPDATE_SCALE_COORDS);
	}

	public int getMode () {
		return mMode;
	}

	public int getHiliteAtom () {
		return mCurrentHiliteAtom;
	}

	public int getHiliteBond () {
		return mCurrentHiliteBond;
	}

	public void setHiliteBondSet ( int[] bondSet)
	{
		mHiliteBondSet = bondSet;
		update(UPDATE_REDRAW);
	}

	public void setReactionMode ( boolean rxn)
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
	protected void setUpdateMode(int mode){
		mUpdateMode = mode;
	}

	public boolean isAtomColorSupported() {
		return mAtomColorSupported;
	}

	public void setAtomColorSupported(boolean acs){
		mAtomColorSupported = acs;
	}

	private void cleanupCoordinates(GenericDrawContext context) {
		int selectedAtomCount = 0;
		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
			if (mMol.isSelectedAtom(atom)) {
				selectedAtomCount++;
			}
		}
		boolean selectedOnly = (selectedAtomCount != 0 && selectedAtomCount != mMol.getAllAtoms());

		if ((mMode & MODE_MULTIPLE_FRAGMENTS) != 0)
			cleanupMultiFragmentCoordinates(context, selectedOnly);
		else
			cleanupMoleculeCoordinates(context, selectedOnly);
	}

	private void cleanupMoleculeCoordinates (GenericDrawContext context,boolean selectedOnly){
		if (mUpdateMode == UPDATE_INVENT_COORDS) {
			if (selectedOnly)
				for (int atom = 0; atom<mMol.getAllAtoms(); atom++)
					mMol.setAtomMarker(atom, !mMol.isSelectedAtom(atom));

			new CoordinateInventor(selectedOnly ? CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS : 0).invent(mMol);

			if (selectedOnly)
				mMol.removeAtomMarkers();
		}

		mDepictor.updateCoords(context, new GenericRectangle(0, 0, mCanvas.getCanvasWidth(),
				mCanvas.getCanvasHeight()), AbstractDepictor.cModeInflateToMaxAVBL | mMaxAVBL);
	}

	private void cleanupMultiFragmentCoordinates(GenericDrawContext context, boolean selectedOnly) {
		if (selectedOnly && mUpdateMode == UPDATE_INVENT_COORDS) {
			int[] fragmentAtom = new int[mFragment.length];
			for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
				int fragment = mFragmentNo[atom];
				mFragment[fragment].setAtomMarker(fragmentAtom[fragment], !mMol.isSelectedAtom(atom));
				fragmentAtom[fragment]++;
			}
		}

		GenericRectangle[] boundingRect = new GenericRectangle[mFragment.length];
//		float fragmentWidth = 0.0f;
		for (int fragment = 0; fragment<mFragment.length; fragment++) {
			if (mUpdateMode == UPDATE_INVENT_COORDS) {
				new CoordinateInventor(selectedOnly ? CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS : 0).invent(mFragment[fragment]);
				mFragment[fragment].setStereoBondsFromParity();
			}
			GenericDepictor d = new GenericDepictor(mFragment[fragment]);
			d.updateCoords(context, null, AbstractDepictor.cModeInflateToMaxAVBL | mMaxAVBL);
			boundingRect[fragment] = d.getBoundingRect();
//			fragmentWidth += boundingRect[fragment].width;
		}

		double spacing = FRAGMENT_CLEANUP_DISTANCE * mMaxAVBL;
		double avbl = getScaledAVBL();
		double arrowWidth = ((mMode & MODE_REACTION) == 0) ? 0f
				: (mUpdateMode == UPDATE_SCALE_COORDS_USE_FRAGMENTS) ?
				DEFAULT_ARROW_LENGTH * mCanvas.getCanvasWidth()
				: ((ReactionArrow)mDrawingObjectList.get(0)).getLength() * mMaxAVBL / avbl;

		double rawX = 0.5 * spacing;
		for (int fragment = 0; fragment<=mFragment.length; fragment++) {
			if ((mMode & MODE_REACTION) != 0 && fragment == mReactantCount) {
				((ReactionArrow)mDrawingObjectList.get(0)).setCoordinates(
						rawX - spacing / 2, mCanvas.getCanvasHeight() / 2, rawX - spacing / 2 + arrowWidth, mCanvas.getCanvasHeight() / 2);
				rawX += arrowWidth;
			}

			if (fragment == mFragment.length) {
				break;
			}

			double dx = rawX - boundingRect[fragment].x;
			double dy = 0.5 * (mCanvas.getCanvasHeight() - boundingRect[fragment].height)
					- boundingRect[fragment].y;
			mFragment[fragment].translateCoords(dx, dy);

			rawX += spacing + boundingRect[fragment].width;
		}

		mDepictor.updateCoords(context, new GenericRectangle(0, 0, mCanvas.getCanvasWidth(),
				mCanvas.getCanvasHeight()), AbstractDepictor.cModeInflateToMaxAVBL | mMaxAVBL);

		int[] fragmentAtom = new int[mFragment.length];
		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
			int fragment = mFragmentNo[atom];
			mMol.setAtomX(atom, mFragment[fragment].getAtomX(fragmentAtom[fragment]));
			mMol.setAtomY(atom, mFragment[fragment].getAtomY(fragmentAtom[fragment]));

			fragmentAtom[fragment]++;
		}

		mMol.setStereoBondsFromParity();
	}

	private void analyzeFragmentMembership() {
		mMol.ensureHelperArrays(Molecule.cHelperParities);

		int[] fragmentNo = new int[mMol.getAllAtoms()];
		int fragments = mMol.getFragmentNumbers(fragmentNo, false, true);

		fragments = joinCloseFragments(fragmentNo, fragments);
		sortFragmentsByPosition(fragmentNo, fragments);
		mFragmentNo = fragmentNo;

		mFragment = mMol.getFragments(fragmentNo, fragments);
	}

	private int joinCloseFragments ( int[] fragmentNo, int fragments){
		if (fragments<2) {
			return fragments;
		}

		boolean[][] mergeFragments = new boolean[fragments][];
		for (int i = 1; i<fragments; i++) {
			mergeFragments[i] = new boolean[i];
		}

		double avbl = getScaledAVBL();
		for (int atom1 = 1; atom1<mMol.getAllAtoms(); atom1++) {
			for (int atom2 = 0; atom2<atom1; atom2++) {
				double dx = mMol.getAtomX(atom2) - mMol.getAtomX(atom1);
				double dy = mMol.getAtomY(atom2) - mMol.getAtomY(atom1);
				double distance = Math.sqrt(dx * dx + dy * dy);
				if (distance<FRAGMENT_GROUPING_DISTANCE * avbl) {
					int fragment1 = fragmentNo[atom1];
					int fragment2 = fragmentNo[atom2];
					if (fragment1 != fragment2) {
						if (fragment1>fragment2) {
							mergeFragments[fragment1][fragment2] = true;
						} else {
							mergeFragments[fragment2][fragment1] = true;
						}
					}
				}
			}
		}

		int[] newFragmentIndex = new int[fragments];
		for (int fragment = 0; fragment<fragments; fragment++) {
			newFragmentIndex[fragment] = fragment;
		}

		int mergeCount = 0;
		for (int i = 1; i<fragments; i++) {
			for (int j = 0; j<i; j++) {
				if (mergeFragments[i][j]) {
					int index1 = newFragmentIndex[i];
					int index2 = newFragmentIndex[j];
					if (index1 != index2) {
						mergeCount++;
						int minIndex = Math.min(index1, index2);
						int maxIndex = Math.max(index1, index2);
						for (int k = 0; k<fragments; k++) {
							if (newFragmentIndex[k] == maxIndex) {
								newFragmentIndex[k] = minIndex;
							} else if (newFragmentIndex[k]>maxIndex) {
								newFragmentIndex[k]--;
							}
						}
					}
				}
			}
		}

		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
			fragmentNo[atom] = newFragmentIndex[fragmentNo[atom]];
		}

		return fragments - mergeCount;
	}

	private void sortFragmentsByPosition ( int[] fragmentNo, int fragments){
		int[][] fragmentDescriptor = new int[fragments][((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE)) != 0) ? 2 : 1];
		for (int fragment = 0; fragment<fragments; fragment++) {
			fragmentDescriptor[fragment][0] = fragment;
		}

		Point[] fragmentCOG = calculateFragmentCenterOfGravity(fragmentNo, fragments);

		if ((mMode & MODE_REACTION) != 0) {
			mReactantCount = 0;
			ReactionArrow arrow = ((mMode & MODE_REACTION) != 0) ? (ReactionArrow)mDrawingObjectList.get(0) : null;
			for (int fragment = 0; fragment<fragments; fragment++) {
				fragmentDescriptor[fragment][1] = (arrow.isOnProductSide(fragmentCOG[fragment].x,
						fragmentCOG[fragment].y)) ? 1 : 0;
				if (fragmentDescriptor[fragment][1] == 0) {
					mReactantCount++;
				}
			}
		} else if ((mMode & MODE_MARKUSH_STRUCTURE) != 0) {
			mReactantCount = fragments;
			for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
				if (mMol.getAtomicNo(atom) == 0 && fragmentDescriptor[fragmentNo[atom]][1] == 0) {
					fragmentDescriptor[fragmentNo[atom]][1] = 1;
					mReactantCount--;
				}
			}
		}

		final Point[] cog = fragmentCOG;
		Arrays.sort(fragmentDescriptor, new Comparator<int[]>() {
			public int compare(int[] fragmentDescriptor1, int[] fragmentDescriptor2) {
				if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE)) != 0) {
					if (fragmentDescriptor1[1] != fragmentDescriptor2[1]) {
						return (fragmentDescriptor1[1] == 0) ? -1 : 1;
					}
				}

				return (cog[fragmentDescriptor1[0]].x
						+ cog[fragmentDescriptor1[0]].y
						<cog[fragmentDescriptor2[0]].x
						+ cog[fragmentDescriptor2[0]].y) ? -1 : 1;
			}
		});

		int[] newFragmentIndex = new int[fragments];
		Point[] centerOfGravity = new Point[fragments];
		for (int fragment = 0; fragment<fragments; fragment++) {
			int oldIndex = ((int[])fragmentDescriptor[fragment])[0];
			newFragmentIndex[oldIndex] = fragment;
			centerOfGravity[fragment] = fragmentCOG[oldIndex];
		}

		fragmentCOG = centerOfGravity;
		for (int atom1 = 0; atom1<mMol.getAllAtoms(); atom1++) {
			fragmentNo[atom1] = newFragmentIndex[fragmentNo[atom1]];
		}
	}

	private Point[] calculateFragmentCenterOfGravity ( int[] fragmentNo, int fragments){
		Point[] fragmentCOG = new Point[fragments];
		int[] fragmentAtoms = new int[fragments];
		for (int fragment = 0; fragment<fragments; fragment++) {
			fragmentCOG[fragment] = new Point(0, 0);
		}
		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
			fragmentCOG[fragmentNo[atom]].x += mMol.getAtomX(atom);
			fragmentCOG[fragmentNo[atom]].y += mMol.getAtomY(atom);
			fragmentAtoms[fragmentNo[atom]]++;
		}
		for (int fragment = 0; fragment<fragments; fragment++) {
			fragmentCOG[fragment].x /= fragmentAtoms[fragment];
			fragmentCOG[fragment].y /= fragmentAtoms[fragment];
		}
		return fragmentCOG;
	}

	private void updateCursor () {
		int cursor = -1;
		switch (mCurrentTool) {
			case GenericEditorToolbar.cToolZoom:
				cursor = SwingCursorHelper.cZoomCursor;
				break;
			case GenericEditorToolbar.cToolLassoPointer:
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
			case GenericEditorToolbar.cToolDelete:
				cursor = SwingCursorHelper.cDeleteCursor;
				break;
			case GenericEditorToolbar.cToolChain:
				cursor = SwingCursorHelper.cChainCursor;
				break;
			case GenericEditorToolbar.cToolText:
				cursor = SwingCursorHelper.cTextCursor;
				break;
			default:
				cursor = SwingCursorHelper.cPointerCursor;
				break;
		}

		if (mCurrentCursor != cursor) {
			mCurrentCursor = cursor;
			mUIHelper.setCursor(cursor);
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

	private Point2D calculateCenterOfGravity ( boolean selectedOnly){
		int atoms = 0;
		double sumx = 0;
		double sumy = 0;
		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
			if (!selectedOnly || mMol.isSelectedAtom(atom)) {
				sumx += mMol.getAtomX(atom);
				sumy += mMol.getAtomY(atom);
				atoms++;
			}
		}
		return atoms>1 ? new Point2D.Double(sumx / atoms, sumy / atoms) : null;
	}

	private void rotate180() {
		boolean selectedOnly = false;
		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
			if (mMol.isSelectedAtom(atom)) {
				selectedOnly = true;
				break;
			}
		}

		Point2D cog = calculateCenterOfGravity(selectedOnly);
		if (cog != null) {
			storeState();

			for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
				if (!selectedOnly || mMol.isSelectedAtom(atom)) {
					mMol.setAtomX(atom, 2 * cog.getX() - mMol.getAtomX(atom));
					mMol.setAtomY(atom, 2 * cog.getY() - mMol.getAtomY(atom));
					}
				}
			update(UPDATE_REDRAW);
			}
		}

	private void flip(boolean horiz) {
		boolean selectedOnly = false;
		for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
			if (mMol.isSelectedAtom(atom)) {
				selectedOnly = true;
				break;
				}
			}

		Point2D cog = calculateCenterOfGravity(selectedOnly);
		if (cog != null) {
			storeState();

			for (int atom = 0; atom<mMol.getAllAtoms(); atom++) {
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
			for (int bond = 0; bond<mMol.getAllBonds(); bond++) {
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


