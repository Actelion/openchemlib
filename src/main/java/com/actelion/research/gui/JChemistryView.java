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
import com.actelion.research.chem.io.RXNFileParser;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.clipboard.ClipboardHandler;
import com.actelion.research.gui.dnd.MoleculeDropAdapter;
import com.actelion.research.gui.dnd.MoleculeTransferable;
import com.actelion.research.gui.dnd.ReactionDropAdapter;
import com.actelion.research.gui.dnd.ReactionTransferable;
import com.actelion.research.util.ColorHelper;
import com.actelion.research.util.CursorHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.util.ArrayList;

public class JChemistryView extends JComponent
				 implements ActionListener, DragGestureListener, DragSourceListener, MouseListener, MouseMotionListener {
	public static final int PASTE_AND_DROP_OPTION_ALLOW_FRAGMENT_STATE_CHANGE = 1;
	public static final int PASTE_AND_DROP_OPTION_KEEP_ATOM_COLORS = 2;
	public static final int PASTE_AND_DROP_OPTION_KEEP_BOND_HIGHLIGHTING = 4;
	public static final int PASTE_AND_DROP_OPTION_REMOVE_CATALYSTS = 8;
	public static final int PASTE_AND_DROP_OPTION_REMOVE_DRAWING_OBJECTS = 16;
	public static final int PASTE_AND_DROP_OPTION_LAYOUT_REACTION = 32;

	private static final String ITEM_COPY_RXN = "Copy Reaction";
	private static final String ITEM_PASTE_RXN = "Paste Reaction";
	private static final String ITEM_OPEN_RXN = "Open Reaction File";
	private static final String ITEM_SAVE_RXN = "Save Reaction File";
	private static final String ITEM_COPY_MOLS = "Copy Molecules";
	private static final String ITEM_PASTE_MOLS = "Paste Molecules";

	private static final long serialVersionUID = 20150204L;
	private static final int UPDATE_REDRAW_ONLY = 0;
	private static final int UPDATE_CHECK_COORDS = 1;
	private static final int UPDATE_SCALE_COORDS = 2;

	private static final int PASTE_AND_DROP_OPTIONS_DEFAULT = 0;

	private static final int ALLOWED_DRAG_ACTIONS = DnDConstants.ACTION_COPY_OR_MOVE;
	private static final int ALLOWED_DROP_ACTIONS = DnDConstants.ACTION_COPY_OR_MOVE;

	private ExtendedDepictor mDepictor;
	private ArrayList<StructureListener> mListener;
	private Dimension mSize;
	private int					mChemistryType,mUpdateMode,mDisplayMode,mDragMolecule,mCopyOrDragActions,mPasteOrDropActions,mPasteAndDropOptions;
	private boolean				mShowBorder,mIsDraggingThis,mAllowDropOrPasteWhenDisabled,mIsEditable;
	private Color				mFragmentNoColor;
	private MoleculeDropAdapter mMoleculeDropAdapter = null;
	private ReactionDropAdapter mReactionDropAdapter = null;


	/**
	 * Creates a new JChemistryView for showing a reaction or molecules.
	 * A JChemistryView uses an ExtendedDepictor to handle multiple molecules or a reaction.
	 * For showing one molecule use a JStructureView.
	 * This default implementation will support copy/paste and drag&drop.
	 * @param chemistryType one of ExtendedDepictor.TYPE_MOLECULES and ExtendedDepictor.TYPE_REACTION
	 */
	public JChemistryView(int chemistryType) {
		this(chemistryType, ALLOWED_DRAG_ACTIONS, ALLOWED_DROP_ACTIONS);
		}

	/**
	 * Creates a new JChemistryView for showing a reaction or molecules.
	 * A JChemistryView uses an ExtendedDepictor to handle multiple molecules or a reaction.
	 * For showing one molecule use a JStructureView.
	 * @param chemistryType one of the ExtendedDepictor.TYPE_... options
	 * @param allowedCopyOrDragActions DnDConstants.ACTION_xxx
	 * @param allowedPasteOrDropActions DnDConstants.ACTION_xxx
	 */
	public JChemistryView(int chemistryType, int allowedCopyOrDragActions, int allowedPasteOrDropActions) {
		mChemistryType = chemistryType;
		mCopyOrDragActions = allowedCopyOrDragActions;
		mPasteOrDropActions = allowedPasteOrDropActions;
		mPasteAndDropOptions = PASTE_AND_DROP_OPTIONS_DEFAULT;
		initializeDragAndDrop();
	    addMouseListener(this);
	    addMouseMotionListener(this);
	    mDragMolecule = -1;
		mIsEditable = true;
	    }

    public int getChemistryType() {
		return mChemistryType;
		}

	public void setContent(StereoMolecule mol) {
		setContent(mol, null);
		}

	public void setContent(StereoMolecule mol[]) {
		setContent(mol, null);
		}

	public void setContent(Reaction rxn) {
		setContent(rxn, null);
		}

	public void setContent(StereoMolecule mol, DrawingObjectList drawingObjectList) {
		mDepictor = new ExtendedDepictor(mol, drawingObjectList, true);
		mDepictor.setDisplayMode(mDisplayMode);
		mDepictor.setFragmentNoColor(mFragmentNoColor);
		mUpdateMode = UPDATE_SCALE_COORDS;
	    mDragMolecule = -1;
		repaint();
		}

	public void setContent(StereoMolecule[] mol, DrawingObjectList drawingObjectList) {
		mDepictor = new ExtendedDepictor(mol, drawingObjectList, true);
		mDepictor.setDisplayMode(mDisplayMode);
		mDepictor.setFragmentNoColor(mFragmentNoColor);
		mUpdateMode = UPDATE_SCALE_COORDS;
	    mDragMolecule = -1;
		repaint();
		}

	public void setContent(Reaction rxn, DrawingObjectList drawingObjectList) {
		mDepictor = new ExtendedDepictor(rxn, drawingObjectList, rxn == null ? false : rxn.isReactionLayoutRequired(), true);
		mDepictor.setDisplayMode(mDisplayMode);
		mDepictor.setFragmentNoColor(mFragmentNoColor);
		mUpdateMode = UPDATE_SCALE_COORDS;
	    mDragMolecule = -1;
		repaint();
		}

	/**
	 * fragment status change on drop is allowed then dropping a fragment (molecule)
	 * on a molecule (fragment) inverts the status of the view's chemical object.
	 * As default status changes are prohibited.
	 * @param options flag list of PASTE_AND_DROP_OPTION...
	 */
	public void setPasteAndDropOptions(int options) {
		mPasteAndDropOptions = options;
	}

	public void setAllowDropOrPasteWhenDisabled(boolean b) {
		mAllowDropOrPasteWhenDisabled = b;
		}

	@Override
	public void setEnabled(boolean enable) {
		boolean changed = (enable != isEnabled());
		super.setEnabled(enable);
		if (changed) {
			if (mMoleculeDropAdapter != null)
				mMoleculeDropAdapter.setActive(enable);
			if (mReactionDropAdapter != null)
				mReactionDropAdapter.setActive(enable);
			if (enable)
				mDepictor.setOverruleColor(null, null);
			else
				mDepictor.setOverruleColor(ColorHelper.getContrastColor(Color.GRAY, getBackground()), getBackground());
			repaint();
			}
		}

	public boolean isEditable() {
		return mIsEditable;
		}

	public void setEditable(boolean b) {
		if (mIsEditable != b)
			mIsEditable = b;
		}

	public void setFragmentNoColor(Color c) {
	        // use setFragmentNoColor(null) if you don't want fragment numbers to be shown
		mFragmentNoColor = c;
		if (mDepictor != null)
			mDepictor.setFragmentNoColor(c);
		}

	public void setDisplayMode(int displayMode) {
		mDisplayMode = displayMode;
		if (mDepictor != null)
			mDepictor.setDisplayMode(displayMode);
		}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);

		if (mDepictor == null)
			return;

		Dimension theSize = getSize();
		if (theSize.width == 0 || theSize.height == 0)
			return;

		Insets insets = getInsets();
		theSize.width -= insets.left + insets.right;
		theSize.height -= insets.top + insets.bottom;

		if (mSize == null
		 || mSize.width != theSize.width
		 || mSize.height != theSize.height
		 || mUpdateMode == UPDATE_SCALE_COORDS) {
			mDepictor.validateView(g, new Rectangle2D.Double(insets.left, insets.top, theSize.width, theSize.height),
									AbstractDepictor.cModeInflateToMaxAVBL);
			}
		else if (mUpdateMode == UPDATE_CHECK_COORDS) {
			mDepictor.validateView(g, new Rectangle2D.Double(insets.left, insets.top, theSize.width, theSize.height), 0);
			}

		mSize = theSize;

//		g.setColor(getBackground());
//		g.fillRect(0, 0, theSize.width, theSize.height);

		Graphics2D g2 = (Graphics2D)g;
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g2.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
		g2.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);

		mDepictor.setForegroundColor(getForeground(), getBackground());

		mDepictor.paint(g);

		if (mShowBorder) {
			g.setColor(Color.gray);
			g.drawRect(insets.left,insets.top,theSize.width - 1,theSize.height - 1);
			g.drawRect(insets.left + 1,insets.top + 1,theSize.width - 3,theSize.height - 3);
			}

		mUpdateMode = UPDATE_REDRAW_ONLY;
		}

	private void updateBorder(boolean showBorder) {
		if(mShowBorder != showBorder){
			mShowBorder = showBorder;
			repaint();
			}
		}

	@Override public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals(ITEM_COPY_MOLS)) {
			ClipboardHandler ch = new ClipboardHandler();
			if (mDepictor.getMoleculeCount() == 1) {
				ch.copyMolecule(mDepictor.getMolecule(0));
				}
			else {
				StereoMolecule mol = new StereoMolecule();
				for (int i=0; i<mDepictor.getMoleculeCount(); i++) {
					AbstractDepictor depictor = mDepictor.getMoleculeDepictor(i);
					int offset = mol.getAllAtoms();
					StereoMolecule m = mDepictor.getMolecule(i);
					m.copyMolecule(mol);
					DepictorTransformation dt = depictor.getTransformation();
					for (int atom=0; atom<m.getAllAtoms(); atom++) {
						mol.setAtomX(atom+offset, dt.transformX(m.getAtomX(atom)));
						mol.setAtomY(atom+offset, dt.transformY(m.getAtomY(atom)));
						}
					}
				ch.copyMolecule(mol);
				}
			}

		if (e.getActionCommand().equals(ITEM_COPY_RXN)) {
			ClipboardHandler ch = new ClipboardHandler();
			Reaction rxn = mDepictor.getReaction();
			ch.copyReaction(rxn);
			}

		if (e.getActionCommand().equals(ITEM_PASTE_MOLS) && mIsEditable) {
			ClipboardHandler ch = new ClipboardHandler();
			StereoMolecule mol = ch.pasteMolecule();
			if (mol != null)
				pasteOrDropMolecule(mol);
			}

		if (e.getActionCommand().equals(ITEM_PASTE_RXN) && mIsEditable) {
			ClipboardHandler ch = new ClipboardHandler();
			Reaction rxn = ch.pasteReaction();
			if (rxn != null)
				pasteOrDropReaction(rxn);
			}

		if (e.getActionCommand().equals(ITEM_OPEN_RXN) && mIsEditable) {
			File rxnFile = FileHelper.getFile(this, "Please select a reaction file", FileHelper.cFileTypeRXN);
			if (rxnFile != null) {
				try {
					Reaction reaction = new RXNFileParser().getReaction(rxnFile);
					pasteOrDropReaction(reaction);
					}
				catch (Exception ex) {}
				}
			}

		if (e.getActionCommand().equals(ITEM_SAVE_RXN)) {
			Reaction rxn = mDepictor.getReaction();
			new FileHelper(this).saveRXNFile(rxn);
			}
		}

	@Override
	public void mousePressed(MouseEvent e) {
		handlePopupTrigger(e);

		setCursor(CursorHelper.getCursor((mDragMolecule == -1) ? CursorHelper.cPointerCursor : CursorHelper.cFistCursor));
		}

	@Override
	public void mouseReleased(MouseEvent e) {
		handlePopupTrigger(e);
		}

	private void handlePopupTrigger(MouseEvent e) {
		if (e.isPopupTrigger()) {
			JPopupMenu popup = new JPopupMenu();
			if (mChemistryType == ExtendedDepictor.TYPE_MOLECULES) {
				if ((mCopyOrDragActions & DnDConstants.ACTION_COPY) != 0 && mDepictor.getMoleculeCount() != 0) {
					JMenuItem item = new JMenuItem(ITEM_COPY_MOLS);
					item.addActionListener(this);
					popup.add(item);
					}
				if ((isEnabled() || mAllowDropOrPasteWhenDisabled) && mIsEditable) {
					if (mCopyOrDragActions != DnDConstants.ACTION_NONE) {
						JMenuItem item = new JMenuItem(ITEM_PASTE_MOLS);
						item.addActionListener(this);
						popup.add(item);
						}
					}
				}

			if (mChemistryType == ExtendedDepictor.TYPE_REACTION) {
				if ((mCopyOrDragActions & DnDConstants.ACTION_COPY) != 0 && mDepictor.getReaction() != null && !mDepictor.getReaction().isEmpty()) {
					JMenuItem item = new JMenuItem(ITEM_COPY_RXN);
					item.addActionListener(this);
					popup.add(item);
					}
				if ((isEnabled() || mAllowDropOrPasteWhenDisabled) && mIsEditable) {
					if (mCopyOrDragActions != DnDConstants.ACTION_NONE) {
						JMenuItem item = new JMenuItem(ITEM_PASTE_RXN);
						item.addActionListener(this);
						popup.add(item);
						}

					if (popup.getComponentCount() != 0)
						popup.addSeparator();

					JMenuItem itemOpen = new JMenuItem(ITEM_OPEN_RXN);
					itemOpen.addActionListener(this);
					popup.add(itemOpen);
					}

				if (mDepictor.getReaction() != null && !mDepictor.getReaction().isEmpty()) {
					if (popup.getComponentCount() != 0 || isEnabled() || mAllowDropOrPasteWhenDisabled)
						popup.addSeparator();

					JMenuItem item = new JMenuItem(ITEM_SAVE_RXN);
					item.addActionListener(this);
					popup.add(item);
					}
				}

			if (popup.getComponentCount() != 0)
				popup.show(this, e.getX(), e.getY());

			return;
			}
		}

	@Override public void mouseClicked(MouseEvent e) {}
	@Override public void mouseEntered(MouseEvent e) {}
	@Override public void mouseExited(MouseEvent e) {}

	@Override
	public void mouseDragged(MouseEvent e) {
		if (mDragMolecule != -1)
			setCursor(CursorHelper.getCursor(CursorHelper.cFistCursor));
		}

	@Override
	public void mouseMoved(MouseEvent e) {
		int x = e.getX();
		int y = e.getY();
		if (mDepictor != null) {
			int dragMolecule = -1;
			for (int i=0; i<mDepictor.getMoleculeCount(); i++) {
				Rectangle2D.Double bounds = mDepictor.getMoleculeDepictor(i).getBoundingRect();
				if (bounds != null && bounds.contains(x, y)) {
					dragMolecule = i;
					break;
					}
				}
			if (mDragMolecule != dragMolecule) {
				mDragMolecule = dragMolecule;
				repaint();
				}
			}

		setCursor(CursorHelper.getCursor((mDragMolecule == -1) ? CursorHelper.cPointerCursor
															   : CursorHelper.cHandCursor));
		}

	public void addStructureListener(StructureListener l) {
		if(mListener == null)
			mListener = new ArrayList<StructureListener>();

		mListener.add(l);
		}

	public void removeStructureListener(StructureListener l) {
		if(mListener != null)
			mListener.remove(l);
		}

	public void informListeners() {
		if (mListener != null)
			for (StructureListener l:mListener)
				l.structureChanged(null);
		}

	public boolean canDrop() {
		return isEditable() && (isEnabled() || mAllowDropOrPasteWhenDisabled) && !mIsDraggingThis;
		}

	private void initializeDragAndDrop() {
 		if (mCopyOrDragActions != DnDConstants.ACTION_NONE)
			DragSource.getDefaultDragSource().createDefaultDragGestureRecognizer(this, mCopyOrDragActions, this);

		if(mPasteOrDropActions != DnDConstants.ACTION_NONE) {
			if (mChemistryType == ExtendedDepictor.TYPE_MOLECULES) {
				mMoleculeDropAdapter = new MoleculeDropAdapter() {
					public void onDropMolecule(StereoMolecule m,Point pt) {
						if (m != null && canDrop()) {
							pasteOrDropMolecule(m);
							onDrop();
							}
						updateBorder(false);
						}

					public void dragEnter(DropTargetDragEvent e) {
						boolean drop = canDrop() && isDropOK(e) ;
						if (!drop)
							e.rejectDrag();
						updateBorder(drop);
						}

					public void dragExit(DropTargetEvent e) {
						updateBorder(false);
						}
					};

				new DropTarget(this, mPasteOrDropActions, mMoleculeDropAdapter, true);
//			new DropTarget(this,mAllowedDropAction,mDropAdapter,true, getSystemFlavorMap());
				}

			if (mChemistryType == ExtendedDepictor.TYPE_REACTION) {
				mReactionDropAdapter = new ReactionDropAdapter() {
					public void onDropReaction(Reaction r, Point pt) {
						if (r != null && canDrop()) {
							pasteOrDropReaction(r);
							onDrop();
							}
						updateBorder(false);
						}

					public void dragEnter(DropTargetDragEvent e) {
						boolean drop = canDrop() && isDropOK(e) ;
						if (!drop)
							e.rejectDrag();
						updateBorder(drop);
						}

					public void dragExit(DropTargetEvent e) {
						updateBorder(false);
						}
					};

				new DropTarget(this, mPasteOrDropActions, mReactionDropAdapter, true);
//			new DropTarget(this,mAllowedDropAction,mDropAdapter,true, getSystemFlavorMap());
				}
			}
		}

	private void pasteOrDropMolecule(StereoMolecule m) {
		boolean isFragment = mDepictor.isFragment();
		StereoMolecule mol = new StereoMolecule(m);
		if ((mPasteAndDropOptions & PASTE_AND_DROP_OPTION_KEEP_ATOM_COLORS) == 0)
			mol.removeAtomColors();
		if ((mPasteAndDropOptions & PASTE_AND_DROP_OPTION_KEEP_BOND_HIGHLIGHTING) == 0)
			mol.removeBondHiliting();
		if ((mPasteAndDropOptions & PASTE_AND_DROP_OPTION_ALLOW_FRAGMENT_STATE_CHANGE) == 0)
			mol.setFragment(isFragment);
		setContent(mol);
		repaint();
		informListeners();
		}

	private void pasteOrDropReaction(Reaction r) {
		boolean isFragment = mDepictor.isFragment();
		Reaction rxn = new Reaction(r);
		if ((mPasteAndDropOptions & PASTE_AND_DROP_OPTION_KEEP_ATOM_COLORS) == 0)
			for (int m=0; m<rxn.getMolecules(); m++)
				rxn.getMolecule(m).removeAtomColors();
		if ((mPasteAndDropOptions & PASTE_AND_DROP_OPTION_KEEP_BOND_HIGHLIGHTING) == 0)
			for (int m=0; m<rxn.getMolecules(); m++)
				rxn.getMolecule(m).removeBondHiliting();
		if ((mPasteAndDropOptions & PASTE_AND_DROP_OPTION_ALLOW_FRAGMENT_STATE_CHANGE) == 0)
			rxn.setFragment(isFragment);
		if ((mPasteAndDropOptions & PASTE_AND_DROP_OPTION_REMOVE_CATALYSTS) != 0)
			rxn.removeCatalysts();
		if ((mPasteAndDropOptions & PASTE_AND_DROP_OPTION_REMOVE_DRAWING_OBJECTS) != 0)
			rxn.removeDrawingObjects();
		setContent(rxn);
		repaint();
		informListeners();
		}

	public Transferable getMoleculeTransferable(Point pt) {
		StereoMolecule mol = (mDragMolecule == -1) ? null
							 : mDepictor.getMolecule(mDragMolecule).getCompactCopy();
		return new MoleculeTransferable(mol);
		}

	public Transferable getReactionTransferable() {
		return new ReactionTransferable(mDepictor.getReaction());
		}

	public void dragEnter(DragSourceDragEvent e) {
		DragSourceContext context = e.getDragSourceContext();
		int dropAction = e.getDropAction();
		if ((dropAction & ALLOWED_DRAG_ACTIONS) != 0) {
			context.setCursor(DragSource.DefaultCopyDrop);
			}
		else {
            context.setCursor(DragSource.DefaultCopyNoDrop);
			}
		}

	public void dragOver(DragSourceDragEvent e) {}
	public void dragExit(DragSourceEvent e) {}
	public void dragDropEnd(DragSourceDropEvent e) {
		System.out.println("dragDropEnd()");
		mIsDraggingThis = false;
		}

	public void dropActionChanged(DragSourceDragEvent e) {
		DragSourceContext context = e.getDragSourceContext();
		int dropAction = e.getDropAction();
		if ((dropAction & ALLOWED_DRAG_ACTIONS) != 0) {
			context.setCursor(DragSource.DefaultCopyDrop);
			}
		else {
			context.setCursor(DragSource.DefaultCopyNoDrop);
			}
		}

	public void dragGestureRecognized(DragGestureEvent e) {
		if((e.getDragAction() & ALLOWED_DRAG_ACTIONS) != 0) {
			Transferable transferable = (mChemistryType == ExtendedDepictor.TYPE_REACTION) ?
					getReactionTransferable() : getMoleculeTransferable(e.getDragOrigin());
			if (transferable != null) {
				try {
					e.startDrag(CursorHelper.getCursor(CursorHelper.cFistCursor), transferable, this);
					mIsDraggingThis = true;
//					e.startDrag(DragSource.DefaultCopyNoDrop, transferable, this);
					}
				catch (InvalidDnDOperationException idoe) {
					System.err.println(idoe);
					}
				}
			}
		}

	// Drag notifications if needed by subclasses
	protected void onDragEnter() {}
	protected void onDragExit() {}
	protected void onDragOver() {}
	protected void onDrop() {}
	}

