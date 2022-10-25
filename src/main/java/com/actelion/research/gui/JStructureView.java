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
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import com.actelion.research.chem.*;
import com.actelion.research.chem.name.StructureNameResolver;
import com.actelion.research.gui.clipboard.IClipboardHandler;
import com.actelion.research.gui.dnd.MoleculeDragAdapter;
import com.actelion.research.gui.dnd.MoleculeDropAdapter;
import com.actelion.research.gui.dnd.MoleculeTransferable;
import com.actelion.research.gui.generic.GenericDepictor;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.gui.swing.SwingDrawContext;
import com.actelion.research.util.ColorHelper;
import com.actelion.research.gui.swing.SwingCursorHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;

public class JStructureView extends JComponent implements ActionListener,MouseListener,MouseMotionListener,StructureListener {
    static final long serialVersionUID = 0x20061113;

    private static final String ITEM_COPY = "Copy Structure";
	private static final String ITEM_COPY_SMILES = "Copy Structure As SMILES-String";
	private static final String ITEM_PASTE= "Paste Structure";
	private static final String ITEM_PASTE_WITH_NAME = ITEM_PASTE+" or Name";
	private static final String ITEM_CLEAR = "Clear Structure";

	private static final long WARNING_MILLIS = 1200;

	private static final int DRAG_MARGIN = 12;

	private ArrayList<StructureListener> mListener;
	private String mIDCode;
	private StereoMolecule mMol,mDisplayMol;
    private GenericDepictor mDepictor;
	private boolean mShowBorder,mAllowFragmentStatusChangeOnPasteOrDrop,mIsDraggingThis,mIsEditable,mDisableBorder;
	private int mChiralTextPosition,mDisplayMode;
	private String[] mAtomText;
	private IClipboardHandler mClipboardHandler;
	protected MoleculeDropAdapter mDropAdapter = null;
	protected int mAllowedDragAction;
	protected int mAllowedDropAction;
	private String mWarningMessage;
	private int[] mAtomHiliteColor;
	private float[] mAtomHiliteRadius;
	private double mTextSizeFactor;

	public JStructureView() {
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
	public JStructureView(StereoMolecule mol) {
        this(mol, DnDConstants.ACTION_COPY_OR_MOVE, DnDConstants.ACTION_COPY_OR_MOVE);
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
	 */
	public JStructureView(StereoMolecule mol, StereoMolecule displayMol) {
        this(mol, displayMol, DnDConstants.ACTION_COPY_OR_MOVE, DnDConstants.ACTION_COPY_OR_MOVE);
	    }

	public JStructureView(int dragAction, int dropAction) {
        this(null, dragAction, dropAction);
	    }

	/**
	 * This creates a standard structure view where the displayed molecule is
	 * used for D&D and clipboard transfer after removing atom colors and bond highlights.
	 * The default will support copy/paste and drag&drop from this view only,
	 * but dropping anything onto this view doesn't have an effect.
	 * Call setEditable(true) to allow changes through drag&drop and pasting.
	 * @param mol used for display, clipboard copy and d&d
	 * @param dragAction
	 * @param dropAction
	 */
	public JStructureView(StereoMolecule mol, int dragAction, int dropAction) {
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
	 * @param dragAction
	 * @param dropAction
	 */
	public JStructureView(StereoMolecule mol, StereoMolecule displayMol, int dragAction, int dropAction) {
		mMol = (mol == null) ? new StereoMolecule() : new StereoMolecule(mol);
		mDisplayMol = (displayMol == null) ? mMol : displayMol;
		mDisplayMode = AbstractDepictor.cDModeHiliteAllQueryFeatures;
		mTextSizeFactor = 1.0;
		mIsEditable = false;
		updateBackground();
		addMouseListener(this);
		addMouseMotionListener(this);
		initializeDragAndDrop(dragAction, dropAction);
	    }

	@Override
	public void updateUI() {
		super.updateUI();
		updateBackground();
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
	 * Atom text is not supported for MODE_REACTION, MODE_MULTIPLE_FRAGMENTS or MODE_MARKUSH_STRUCTURE.
	 * @param atomText null or String array matching atom indexes (may contain null entries)
	 */
	public void setAtomText(String[] atomText) {
		mAtomText = atomText;
		}

	public void setEnabled(boolean enable) {
		if (enable != isEnabled()) {
			updateBackground();
			repaint();
			if (mDropAdapter != null)
				mDropAdapter.setActive(enable);
			}
		super.setEnabled(enable);
		}


	public boolean isEditable() {
		return mIsEditable;
	}

	public void setEditable(boolean b) {
		if (mIsEditable != b)
			mIsEditable = b;
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
		return mIsEditable && isEnabled() && !mIsDraggingThis;
	    }

	@Override
	public synchronized void paintComponent(Graphics g) {
        super.paintComponent(g);

        Dimension theSize = getSize();
		Insets insets = getInsets();
		theSize.width -= insets.left + insets.right;
		theSize.height -= insets.top + insets.bottom;

        if (theSize.width <= 0 || theSize.height <= 0)
            return;

        Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
        g2.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);

		Color fg = g2.getColor();
		Color bg = getBackground();
		g2.setColor(bg);
		g2.fill(new Rectangle(insets.left, insets.top, theSize.width, theSize.height));

		if (mShowBorder && !mDisableBorder) {
			GenericRectangle rect = mDepictor.getBoundingRect();
			if (rect != null) {
				g.setColor(ColorHelper.perceivedBrightness(bg) < 0.5f ? ColorHelper.brighter(bg, 0.85f) : ColorHelper.darker(bg, 0.85f));
				int arc = (int)Math.min(rect.height/4, Math.min(rect.width/4, HiDPIHelper.scale(10)));
				g.fillRoundRect((int)rect.x, (int)rect.y, (int)rect.width, (int)rect.height, arc, arc);
				}
			}

		g2.setColor(fg);

		if (mDisplayMol != null && mDisplayMol.getAllAtoms() != 0) {
			mDepictor = new GenericDepictor(mDisplayMol);
            mDepictor.setDisplayMode(mDisplayMode);
            mDepictor.setFactorTextSize(mTextSizeFactor);
            mDepictor.setAtomText(mAtomText);
            mDepictor.setAtomHighlightColors(mAtomHiliteColor, mAtomHiliteRadius);

            int bgRGB = bg.getRGB();
			if (!isEnabled())
                mDepictor.setOverruleColor(ColorHelper.getContrastColor(0x808080, bgRGB), bgRGB);
			else
				mDepictor.setForegroundColor(getForeground().getRGB(), bgRGB);

			int avbl = HiDPIHelper.scale(AbstractDepictor.cOptAvBondLen);
			SwingDrawContext context = new SwingDrawContext((Graphics2D)g);
			mDepictor.validateView(context, new GenericRectangle(insets.left, insets.top, theSize.width,theSize.height),
								   AbstractDepictor.cModeInflateToMaxAVBL | mChiralTextPosition | avbl);
            mDepictor.paint(context);
			}

		if (mWarningMessage != null) {
			int fontSize = HiDPIHelper.scale(12);
			g.setFont(getFont().deriveFont(Font.BOLD, (float)fontSize));
			Color original = g.getColor();
			g.setColor(Color.RED);
			FontMetrics metrics = g.getFontMetrics();
			Rectangle2D bounds = metrics.getStringBounds(mWarningMessage, g);
			g.drawString(mWarningMessage, insets.left+(int)(theSize.width-bounds.getWidth())/2,
					insets.top+metrics.getHeight());
			g.setColor(original);
			}
		}

	public void setIDCode(String idcode) {
		setIDCode(idcode, null);
	    }

	public synchronized void setIDCode(String idcode, String coordinates) {
		if (idcode != null && idcode.length() == 0)
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
	 * @param mol new molecule used for display; may be null
	 * @param displayMol new molecule used for clipboard copy and d&d, may be null
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
	 * Should only be called if JStructureView's internal Molecule is changed
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
		if(mListener == null)
			mListener = new ArrayList<>();

		mListener.add(l);
		}

    public void removeStructureListener(StructureListener l) {
        if(mListener != null)
            mListener.remove(l);
        }

	public void setChiralDrawPosition(int p) {
		mChiralTextPosition = p;
		}

	@Override public void mouseClicked(MouseEvent e) {}
	@Override public void mouseEntered(MouseEvent e) {}
	@Override public void mouseExited(MouseEvent e) {}

	@Override
	public void mousePressed(MouseEvent e) {
		handlePopupTrigger(e);
		}

	@Override
	public void mouseReleased(MouseEvent e) {
		handlePopupTrigger(e);
		}

	@Override
	public void mouseMoved(MouseEvent e) {
		int x = e.getX();
		int y = e.getY();
		boolean isInRect = false;
		if (mDepictor != null && (mAllowedDragAction & DnDConstants.ACTION_COPY) != 0) {
			GenericRectangle bounds = shrink(mDepictor.getBoundingRect());
			if (bounds != null && bounds.contains(x, y))
				isInRect = true;
			}

		updateBorder(isInRect);
		setCursor(SwingCursorHelper.getCursor(isInRect ? SwingCursorHelper.cHandCursor : SwingCursorHelper.cPointerCursor));
		}

	private GenericRectangle shrink(GenericRectangle rect) {
		int margin = HiDPIHelper.scale(DRAG_MARGIN);
		int marginX = Math.min(margin, (int)rect.width / 6);
		int marginY = Math.min(margin, (int)rect.height / 6);
		return new GenericRectangle((int)rect.x+marginX, (int)rect.y+marginY, (int)rect.width-2*marginX, (int)rect.height-2*marginY);
	}

	@Override public void mouseDragged(MouseEvent e) {}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals(ITEM_COPY)) {
			mClipboardHandler.copyMolecule(mMol);
			}
		if (e.getActionCommand().equals(ITEM_COPY_SMILES)) {
			final String smiles = new IsomericSmilesCreator(mMol).getSmiles();
			final StringSelection data = new StringSelection(smiles);
			Toolkit.getDefaultToolkit().getSystemClipboard().setContents(data, data);
			}
		if (e.getActionCommand().startsWith(ITEM_PASTE) && mIsEditable) {
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
		if (e.getActionCommand().equals(ITEM_CLEAR) && mIsEditable) {
			mMol.clear();
			mDisplayMol = mMol;
			structureChanged();
			}
		}

	protected void showWarningMessage(String msg) {
		mWarningMessage = msg;
		repaint();
		new Thread(() -> {
			try { Thread.sleep(WARNING_MILLIS); } catch (InterruptedException ie) {}
			mWarningMessage = null;
			repaint();
			} ).start();
		}

	private void updateBackground() {
		Color bg = UIManager.getColor(isEditable() && isEnabled() ? "TextField.background" : "TextField.inactiveBackground");
		if (bg != null)
			setBackground(bg);
		}

	private void handlePopupTrigger(MouseEvent e) {
		if (mMol != null && e.isPopupTrigger() && mClipboardHandler != null) {
			JPopupMenu popup = new JPopupMenu();

			JMenuItem item1 = new JMenuItem(ITEM_COPY);
			item1.addActionListener(this);
			item1.setEnabled(mMol.getAllAtoms() != 0);
			popup.add(item1);

			JMenuItem itemCopySmiles = new JMenuItem(ITEM_COPY_SMILES);
			itemCopySmiles.addActionListener(this);
			itemCopySmiles.setEnabled(mMol.getAllAtoms() != 0);
			popup.add(itemCopySmiles);

			if (mIsEditable) {
				String itemText = StructureNameResolver.getInstance() == null ? ITEM_PASTE : ITEM_PASTE_WITH_NAME;
				JMenuItem item2 = new JMenuItem(itemText);
				item2.addActionListener(this);
				popup.add(item2);

				popup.addSeparator();

				JMenuItem item3 = new JMenuItem(ITEM_CLEAR);
				item3.addActionListener(this);
				item3.setEnabled(mMol.getAllAtoms() != 0);
				popup.add(item3);
				}

			popup.show(this, e.getX(), e.getY());
			}
		}

	private void informListeners() {
		if (mListener != null)
			for (int i = 0; i<mListener.size(); i++)
				mListener.get(i).structureChanged(mMol);
		}

	private void initializeDragAndDrop(int dragAction, int dropAction) {
		final JStructureView outer = this;
		mAllowedDragAction = dragAction;
		mAllowedDropAction = dropAction;
		mAllowFragmentStatusChangeOnPasteOrDrop = false;

		if(dragAction != DnDConstants.ACTION_NONE){
			new MoleculeDragAdapter(this) {
				public Transferable getTransferable(Point origin) {
					return getMoleculeTransferable(origin);
				}

				public void onDragEnter() {
					outer.onDragEnter();
				}

				public void dragIsValidAndStarts() {
					mIsDraggingThis = true;
					}

				/*	public void onDragOver() {
					 outer.onDragOver();
					 }
				 */
				public void onDragExit() {
					outer.onDragExit();
				}

				public void dragDropEnd(DragSourceDropEvent e) {
					mIsDraggingThis = false;
				}
			};
		}

		if(dropAction != DnDConstants.ACTION_NONE) {
			mDropAdapter = new MoleculeDropAdapter() {
				public void onDropMolecule(StereoMolecule m,Point pt) {
					if (m != null && canDrop()){
						boolean isFragment = mMol.isFragment();
						mMol = new StereoMolecule(m);
				        mMol.removeAtomColors();
				        mMol.removeBondHiliting();
				        if (!mAllowFragmentStatusChangeOnPasteOrDrop)
				        	mMol.setFragment(isFragment);
				        mDisplayMol = mMol;
						repaint();
						informListeners();
						onDrop();
					}
					updateBorder(false);
				}

				public void dragEnter(DropTargetDragEvent e) {
					boolean drop = canDrop() && isDropOK(e) ;
					if (!drop)
						e.rejectDrag();
//					updateBorder(drop);
				}

				public void dragExit(DropTargetEvent e) {
//					updateBorder(false);
				}
			};

			new DropTarget(this, mAllowedDropAction, mDropAdapter, true);
//			new DropTarget(this,mAllowedDropAction,mDropAdapter,true, getSystemFlavorMap());
		}
	}


	protected Transferable getMoleculeTransferable(Point pt) {
		return new MoleculeTransferable(mMol);
	}

	// Drag notifications if needed by subclasses
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

/*	public java.awt.datatransfer.FlavorMap getSystemFlavorMap() {
	    return new OurFlavorMap();
	    }

    // This class is needed for inter-jvm drag&drop. Although not neccessary for standard environments, it prevents
    // nasty "no native data was transfered" errors. It still might create ClassNotFoundException in the first place by
    // the SystemFlavorMap, but as I found it does not hurt, since the context classloader will be installed after
    // the first call. I know, that this depends heavely on a specific behaviour of the systemflavormap, but for now
    // there's nothing I can do about it.
    static class OurFlavorMap implements FlavorMap, FlavorTable {
    	public java.util.Map<DataFlavor,String> getNativesForFlavors(DataFlavor[] dfs) {
    //	    System.out.println("getNativesForFlavors " + dfs.length);
    //	    for (int i = 0; i < dfs.length; i++)
    //		    System.out.println(" -> " + dfs[i]);
    //
    	    return SystemFlavorMap.getDefaultFlavorMap().getNativesForFlavors(dfs);
    	    }
    
    	public java.util.Map<String,DataFlavor> getFlavorsForNatives(String[] natives) {
    //	    System.out.println("getFlavorsForNatives " + natives.length);
    //	    for (int i = 0; i < natives.length; i++)
    //	        System.out.println(" -> " + natives[i]);
    //
    	    return SystemFlavorMap.getDefaultFlavorMap().getFlavorsForNatives(natives);
    	    }
    
    	public synchronized java.util.List<DataFlavor> getFlavorsForNative(String nat) {
    //	    System.out.println("getFlavorsForNative " + nat);
    	    
    	    return ((SystemFlavorMap)SystemFlavorMap.getDefaultFlavorMap()).getFlavorsForNative(nat);
    	    }
    
    	public synchronized java.util.List<String> getNativesForFlavor(DataFlavor flav) {
    //	    System.out.println("getNativesForFlavor " + flav);
    	    return ((SystemFlavorMap)SystemFlavorMap.getDefaultFlavorMap()).getNativesForFlavor(flav);
    	    }
        }*/
    }
