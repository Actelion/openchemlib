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
import com.actelion.research.gui.clipboard.IClipboardHandler;
import com.actelion.research.gui.dnd.MoleculeDragAdapter;
import com.actelion.research.gui.dnd.MoleculeDropAdapter;
import com.actelion.research.gui.dnd.MoleculeTransferable;
import com.actelion.research.gui.editor.SwingEditorDialog;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.util.ColorHelper;
import com.actelion.research.gui.swing.SwingCursorHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DropTarget;
import java.awt.dnd.DropTargetDragEvent;
import java.awt.dnd.DropTargetEvent;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;

public class CompoundCollectionPane<T> extends JScrollPane
			implements ActionListener,CompoundCollectionListener,MouseListener,MouseMotionListener,StructureListener {
	private static final long serialVersionUID = 0x20060904;

	private static final String[] MESSAGE = { "<to add compounds use popup menu,", "drag&drop, or paste structure(s) or name(s)>" };

	private static final String ADD = "Add...";
	private static final String EDIT = "Edit...";
	private static final String REMOVE = "Remove";
	private static final String REMOVE_ALL = "Remove All";
	private static final String COPY = "Copy";
	private static final String PASTE = "Paste";
	private static final String OPEN = "Add From File...";
	private static final String SAVE_DWAR = "Save DataWarrior-File...";
	private static final String SAVE_SDF2 = "Save SD-File V2...";
	private static final String SAVE_SDF3 = "Save SD-File V3...";

	public static final int FILE_SUPPORT_NONE = 0;
	public static final int FILE_SUPPORT_OPEN_FILES = 1;
	public static final int FILE_SUPPORT_SAVE_FILES = 2;
	public static final int FILE_SUPPORT_OPEN_AND_SAVE_FILES = 3;

	private static final int ALLOWED_DRAG_ACTIONS = DnDConstants.ACTION_COPY_OR_MOVE;
	private static final int ALLOWED_DROP_ACTIONS = DnDConstants.ACTION_COPY_OR_MOVE;

	private final static int cWhiteSpace = 4;

	private CompoundCollectionModel<T> mModel;
	private IClipboardHandler   mClipboardHandler;
	private MoleculeFilter		mCompoundFilter;
	private int					mDisplayMode,mSelectedIndex,mHighlightedIndex,
								mEditedIndex,mFileSupport,mStructureSize,mDragIndex,mDropIndex;
	private Dimension		    mContentSize,mCellSize;
	private JPanel			    mContentPanel;
	private boolean			    mIsVertical,mIsEditable,mIsSelectable,mCreateFragments,
								mIsEnabled,mShowValidationError,mInternalDragAndDropIsMove;
	private String[]            mMessage;
	private ArrayList<JMenuItem> mCustomPopupItemList;
	private ScrollPaneAutoScrollerWhenDragging mScroller;


	/**
	 * This is a visual component to display and edit a compound collection maintained
	 * by a CompoundCollectionModel. Three variations of DefaultCompoundCollectionModel
	 * (.Native, .Molecule, and .IDCode) are available, which internally keep molecule
	 * instances as Object, StereoMolecule or String, respectively. If one of these
	 * default model is used, than the CompoundCollectionPane's T must match this class.
	 * @param model
	 * @param isVertical
	 */
	public CompoundCollectionPane(CompoundCollectionModel<T> model, boolean isVertical) {
		this(model, isVertical, 0, ALLOWED_DRAG_ACTIONS, ALLOWED_DROP_ACTIONS);
		}

	public CompoundCollectionPane(CompoundCollectionModel<T> model, boolean isVertical, int displayMode,
								  int dragAction, int dropAction) {
		mModel = model;
		mModel.addCompoundCollectionListener(this);
		mIsEnabled = true;
		mIsVertical = isVertical;
		mDisplayMode = displayMode;
		mFileSupport = FILE_SUPPORT_OPEN_AND_SAVE_FILES;
		mStructureSize = 0;
		mSelectedIndex = -1;
		mHighlightedIndex = -1;
		mDragIndex = -1;
		mDropIndex = -1;
		init();
		initializeDragAndDrop(dragAction, dropAction);
		mScroller = new ScrollPaneAutoScrollerWhenDragging(this, isVertical);
		}

	public CompoundCollectionModel<T> getModel() {
		return mModel;
		}

	public void setEnabled(boolean b) {
		super.setEnabled(b);
		if (mIsVertical)
			getVerticalScrollBar().setEnabled(b);
		else
			getHorizontalScrollBar().setEnabled(b);
		mIsEnabled = b;
		repaint();
		}

	/**
	 * Defines the behaviour for internal drag&drop. The default is false,
	 * which means that a dragged structure is copied to the new internal position.
	 * @param b
	 */
	public void setInternalDragAndDropIsMove(boolean b) {
		mInternalDragAndDropIsMove = b;
		}

	/**
	 * Defines the width or height of individual structure cells,
	 * depending on whether the the CompoundCollectionPane is horizontal
	 * or vertical, respectively. Setting size to 0 (default) causes
	 * an automatic behavior with a square cell areas and width and height
	 * being implicitly defined by the CompoundCollectionPane component size.
	 */
	public void setStructureSize(int size) {
		mStructureSize = size;
		validateSize();
		repaint();
		}

	/**
	 * Defines, whether the list and individual structures can be edited.
	 * @param editable
	 */
	public void setEditable(boolean editable) {
		mIsEditable = editable;
		updateMouseListening();
		}

	/**
	 * Defines, whether the popup menu contains 'Open' and/or 'Save' items.
	 * As default both, OPEN and SAVE options are active.
	 * @param fileSupport one of the FILE_SUPPORT_... options
	 */
	public void setFileSupport(int fileSupport) {
		mFileSupport = fileSupport;
		}

	public void setCompoundFilter(MoleculeFilter filter) {
		if (mCompoundFilter != filter) {
			mCompoundFilter = filter;
			if (mCompoundFilter != null)
				for (int i=mModel.getSize()-1; i>=0; i--)
					if (!mCompoundFilter.moleculeQualifies(mModel.getMolecule(i)))
						mModel.remove(i);
			}
		}

	public void setSelectable(boolean selectable) {
		mIsSelectable = selectable;
		updateMouseListening();
		}

	/**
	 * Defines whether new created structures are fragments of molecules.
	 * @param createFragments
	 */
	public void setCreateFragments(boolean createFragments) {
		mCreateFragments = createFragments;
		}

	/**
	 * Defines whether a large red question mark is shown
	 * in case of a structure validation error. 
	 * @param showError
	 */
	public void setShowValidationError(boolean showError) {
		mShowValidationError = showError;
		}

	/**
	 *  call this in order to get clipboard support
	 */
	public void setClipboardHandler(IClipboardHandler h) {
		mClipboardHandler = h;
		}

	public IClipboardHandler getClipboardHandler() {
		return mClipboardHandler;
		}

	public void addCustomPopupItem(JMenuItem customItem) {
		if (mCustomPopupItemList == null)
			mCustomPopupItemList = new ArrayList<>();

		mCustomPopupItemList.add(customItem);
		}

	/**
	 * @param msg null for default message or custom message to show if this pane is empty
	 */
	public void setMessage(String msg) {
		if (msg == null) {
			mMessage = null;
			}
		else {
			mMessage = new String[1];
			mMessage[0] = msg;
			}
		if (mModel.getSize() == 0)
			repaint();
		}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals(COPY) && mHighlightedIndex != -1) {
			mClipboardHandler.copyMolecule(mModel.getMolecule(mHighlightedIndex));
			}
		else if (e.getActionCommand().equals(PASTE)) {
			int index = (mHighlightedIndex == -1) ? mModel.getSize() : mHighlightedIndex;
			ArrayList<StereoMolecule> molList = mClipboardHandler.pasteMolecules();
			if (molList != null) {
				int errorCount = 0;
				for (StereoMolecule mol:molList) {
					mol.setFragment(mCreateFragments);
					if (mCompoundFilter == null || mCompoundFilter.moleculeQualifies(mol))
						mModel.addMolecule(index, mol);
					else
						errorCount++;
					}
				if (errorCount != 0)
					JOptionPane.showMessageDialog(getParentFrame(), errorCount+" compound(s) could not be added, because they doesn't qualify.");
				}
			}
		else if (e.getActionCommand().equals(ADD)) {
			editStructure(-1);
			}
		else if (e.getActionCommand().equals(EDIT) && mHighlightedIndex != -1) {
			editStructure(mHighlightedIndex);
			}
		else if (e.getActionCommand().equals(REMOVE) && mHighlightedIndex != -1) {
			mModel.remove(mHighlightedIndex);
			mHighlightedIndex = -1;
			}
		else if (e.getActionCommand().equals(REMOVE_ALL)) {
			mModel.clear();
			mHighlightedIndex = -1;
			}
		else if (e.getActionCommand().equals(OPEN)) {
			ArrayList<StereoMolecule> compounds = new FileHelper(getParentFrame()).readStructuresFromFile(true);
			if (compounds != null) {
				for (StereoMolecule compound:compounds)
					compound.setFragment(mCreateFragments);
				if (mCompoundFilter != null) {
					int count = 0;
					for (int i = compounds.size() - 1; i >= 0; i--) {
						if (!mCompoundFilter.moleculeQualifies(compounds.get(i))) {
							compounds.remove(i);
							count++;
							}
						}
					if (count != 0) {
						JOptionPane.showMessageDialog(getParentFrame(),Integer.toString(count).concat(" compounds were removed, because they don't qualify."));
						}
					}
				mModel.addMoleculeList(compounds);
				}
			}
		else if (e.getActionCommand().equals(SAVE_DWAR)) {
			String filename = new FileHelper(getParentFrame()).selectFileToSave(
					"Save DataWarrior File", FileHelper.cFileTypeDataWarrior, "Untitled");
			if (filename != null) {
				try {
					String title = mCreateFragments ? "Fragment" : "Structure";
					BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename),"UTF-8"));
					writer.write("<datawarrior-fileinfo>");
					writer.newLine();
					writer.write("<version=\"3.1\">");
					writer.newLine();
					writer.write("<rowcount=\""+mModel.getSize()+"\">");
					writer.newLine();
					writer.write("</datawarrior-fileinfo>");
					writer.newLine();
					writer.write("<column properties>");
					writer.newLine();
					writer.write("<columnName=\""+title+"\">");
					writer.newLine();
					writer.write("<columnProperty=\"specialType\tidcode\">");
					writer.newLine();
					writer.write("<columnName=\"coords\">");
					writer.newLine();
					writer.write("<columnProperty=\"specialType\tidcoordinates2D\">");
					writer.newLine();
					writer.write("<columnProperty=\"parent\t"+title+"\">");
					writer.newLine();
					writer.write("</column properties>");
					writer.newLine();
					writer.write(title+"\tcoords");
					writer.newLine();
					for (int i=0; i<mModel.getSize(); i++) {
						if (mModel instanceof DefaultCompoundCollectionModel.IDCode) {
							String idcode = (String)mModel.getCompound(i);
							int index = idcode.indexOf(' ');
							if (index == -1) {
								writer.write(idcode.substring(0, index));
								writer.write('\t');
								writer.write(idcode.substring(index+1));
								}
							else {
								writer.write(idcode);
								writer.write('\t');
								}
							}
						else {
							Canonizer canonizer = new Canonizer(mModel.getMolecule(i));
							writer.write(canonizer.getIDCode());
							writer.write('\t');
							writer.write(canonizer.getEncodedCoordinates());
							}
						writer.newLine();
						}
					writer.close();
					}
				catch (IOException ioe) {
					JOptionPane.showMessageDialog(getParentFrame(), ioe.toString());
					}
				}
			}
		else if (e.getActionCommand().equals(SAVE_SDF2)
			  || e.getActionCommand().equals(SAVE_SDF3)) {
			String version = "Version " + (e.getActionCommand().equals(SAVE_SDF2) ? "2" : "3");
			String filename = new FileHelper(getParentFrame()).selectFileToSave(
					"Save SD-File "+version, FileHelper.cFileTypeSD, "Untitled");
			if (filename != null) {
				try {
					BufferedWriter theWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename),"UTF-8"));
	
					for (int i=0; i<mModel.getSize(); i++) {
						StereoMolecule mol = mModel.getMolecule(i);
	
						if (e.getActionCommand().equals(SAVE_SDF3))
							new MolfileV3Creator(mol).writeMolfile(theWriter);
						else
							new MolfileCreator(mol).writeMolfile(theWriter);
	
						theWriter.write("$$$$");
						theWriter.newLine();
						}
					theWriter.close();
					}
				catch (IOException ioe) {
					JOptionPane.showMessageDialog(getParentFrame(), ioe.toString());
					}
				}
			}
		}

	private void editStructure(int index) {
		mEditedIndex = index;
		StereoMolecule mol = null;
		if (index == -1) {
			mol = new StereoMolecule();
			mol.setFragment(mCreateFragments);
			}
		else {
			mol = mModel.getMolecule(mEditedIndex);
			}
		Component c = getParentFrame();
		SwingEditorDialog theDialog = (c instanceof Frame) ? new SwingEditorDialog((Frame)c, mol) : new SwingEditorDialog((Dialog)c, mol);
		theDialog.addStructureListener(this);
		theDialog.setVisible(true);
		}

	private void updateMouseListening() {
		if (mIsSelectable || mIsEditable) {
			addMouseListener(this);
			addMouseMotionListener(this);
			}
		else {
			removeMouseListener(this);
			removeMouseMotionListener(this);
			}
		}

	private void init() {
		mContentSize = new Dimension();
		mContentPanel = new JPanel() {
			private static final long serialVersionUID = 0x20060904;

			public void paintComponent(Graphics g) {
				super.paintComponent(g);
				((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
				((Graphics2D)g).setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);

				validateSize();

				Rectangle clipRect = g.getClipBounds();

				Color background = UIManager.getColor(mIsEnabled ? "TextArea.background" : "TextArea.inactiveBackground");
				Color foreground = UIManager.getColor(mIsEnabled ? "TextArea.foreground" : "TextArea.inactiveForeground");
				g.setColor(background);
				g.fillRect(clipRect.x, clipRect.y, clipRect.width, clipRect.height);

				int i1 = 0;
				int i2 = 0;

				if (mModel.getSize() != 0) {
					i1 = Math.max(0, mIsVertical ? clipRect.y / mCellSize.height
													 : clipRect.x / mCellSize.width);
					i2 = Math.min(mModel.getSize(), mIsVertical ? 1+(clipRect.y+clipRect.height) / mCellSize.height
																	: 1+(clipRect.x+clipRect.width) / mCellSize.width);
					Color warningColor = (ColorHelper.perceivedBrightness(background) < 0.5f) ?
							ColorHelper.brighter(background, 0.8f) : ColorHelper.darker(background, 0.8f);

					for (int i=i1; i<i2; i++) {
						Rectangle bounds = getMoleculeBounds(i);

						StereoMolecule compound = mModel.getMoleculeForDisplay(i);
						if (mShowValidationError) {
							try {
								compound.validate();
								}
							catch (Exception e) {
								int size = Math.min(bounds.width, bounds.height);
								g.setColor(warningColor);
								g.setFont(g.getFont().deriveFont(Font.BOLD, size));
								FontMetrics m = g.getFontMetrics();
								Rectangle2D b = m.getStringBounds("?", g);
								g.drawString("?", bounds.x+(bounds.width-(int)b.getWidth())/2, bounds.y+(bounds.height-(int)b.getHeight())/2+m.getAscent());
								}
							}

						Depictor2D d = new Depictor2D(compound, mDisplayMode);
						d.validateView((Graphics2D)g,
									   new GenericRectangle(bounds.x, bounds.y, bounds.width, bounds.height),
									   AbstractDepictor.cModeInflateToMaxAVBL);

						d.setForegroundColor(foreground, background);
						d.paint((Graphics2D)g);

						if (mSelectedIndex == i || mHighlightedIndex == i) {
							g.setColor(!mIsEnabled ? ColorHelper.getContrastColor(Color.GRAY, background)
									 : (mSelectedIndex != i) ? Color.BLUE
									 : (mHighlightedIndex != i) ? Color.RED : Color.MAGENTA);
							g.drawRect(bounds.x-2, bounds.y-2, bounds.width+3, bounds.height+3);
							g.drawRect(bounds.x-1, bounds.y-1, bounds.width+1, bounds.height+1);
							}
						}
					}
				else {
					Rectangle bounds = getViewportBorderBounds();
					g.setColor(foreground);
					g.setFont(g.getFont().deriveFont(Font.PLAIN, HiDPIHelper.scale(12)));
					FontMetrics m = g.getFontMetrics();
					String[] message = (mMessage == null) ? MESSAGE : mMessage;
					for (int i=0; i<message.length; i++) {
						Rectangle2D b = m.getStringBounds(message[i], g);
						g.drawString(message[i], bounds.x + (bounds.width - (int)b.getWidth()) / 2,
								bounds.y + i*m.getHeight() + (bounds.height - message.length * m.getHeight()) / 2 + m.getAscent());
						}
					}

				if (mIsEnabled && mIsEditable && mDropIndex>=i1 && mDropIndex<=i2) {
					Rectangle bounds = getMoleculeBounds(mDropIndex);
					g.setColor(ColorHelper.getContrastColor(Color.GRAY, background));
					if (mIsVertical)
						g.fillRect(bounds.x-2, bounds.y-4, bounds.width+4, 5);
					else
						g.fillRect(bounds.x-4, bounds.y-2, 5, bounds.height+4);
					}
				}
			};
		setHorizontalScrollBarPolicy(mIsVertical ? JScrollPane.HORIZONTAL_SCROLLBAR_NEVER : JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		setVerticalScrollBarPolicy(mIsVertical ? JScrollPane.VERTICAL_SCROLLBAR_ALWAYS : JScrollPane.VERTICAL_SCROLLBAR_NEVER);
		setViewportView(mContentPanel);
		}

	public void collectionUpdated(int fromIndex, int toIndex) {
		if (mSelectedIndex >= fromIndex && mSelectedIndex <= toIndex)
			mSelectedIndex = -1;
		if (mHighlightedIndex >= fromIndex && mHighlightedIndex <= toIndex)
			mHighlightedIndex = -1;

		repaint();
		}

	private Rectangle getMoleculeBounds(int molIndex) {
		int x = cWhiteSpace/2;
		int y = cWhiteSpace/2;

		if (mIsVertical)
			y += molIndex * mCellSize.height;
		else
			x += molIndex * mCellSize.width;

		return new Rectangle(x, y, mCellSize.width-cWhiteSpace, mCellSize.height-cWhiteSpace);
		}

	private int getMoleculeIndex(int x, int y) {
		if (mModel.getSize() == 0 || mCellSize.width == 0 || mCellSize.height == 0)
			return -1;

		Point p = getViewport().getViewPosition();
		int index = (mIsVertical) ? (y+p.y) / mCellSize.height
								  : (x+p.x) / mCellSize.width;
		return (index < mModel.getSize()) ? index : -1;
		}

	public void mouseClicked(MouseEvent e) {
		if (mIsEnabled && mIsEditable && e.getClickCount() == 2 && mHighlightedIndex != -1) {
			editStructure(mHighlightedIndex);
			}
		}

	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}

	public void mousePressed(MouseEvent e) {
		if (mIsEnabled) {
			if (e.isPopupTrigger()) {
				handlePopupTrigger(e);
				}
			else if (mIsSelectable) {
				int index = getMoleculeIndex(e.getX(), e.getY());
				if (mSelectedIndex != index) {
					mSelectedIndex = index;
					setSelection(index);
					repaint();
					}
				}
			}
		}

	public void mouseReleased(MouseEvent e) {
		if (mIsEnabled && e.isPopupTrigger())
			handlePopupTrigger(e);
		}

	public void mouseDragged(MouseEvent e) {}

	public void mouseMoved(MouseEvent e) {
		if (mIsEnabled) {
			int index = getMoleculeIndex(e.getX(), e.getY());
			if (mHighlightedIndex != index) {
				mHighlightedIndex = index;
				setCursor(SwingCursorHelper.getCursor(index == -1 ? SwingCursorHelper.cPointerCursor : SwingCursorHelper.cHandCursor));
				repaint();
				}
			}
		mDragIndex = -1;
		}

	private void handlePopupTrigger(MouseEvent e) {
		JPopupMenu popup = new JPopupMenu();
		JMenuItem item = new JMenuItem(ADD);
		item.addActionListener(this);
		popup.add(item);
		if (mHighlightedIndex != -1) {
			item = new JMenuItem(EDIT);
			item.addActionListener(this);
			popup.add(item);
			item = new JMenuItem(REMOVE);
			item.addActionListener(this);
			popup.add(item);
			}

		if (mModel.getSize() != 0) {
			item = new JMenuItem(REMOVE_ALL);
			item.addActionListener(this);
			popup.add(item);
			}

		if (mClipboardHandler != null) {
			popup.addSeparator();
			if (mHighlightedIndex != -1) {
				item = new JMenuItem(COPY);
				item.addActionListener(this);
				popup.add(item);
				}
			item = new JMenuItem(PASTE);
			item.addActionListener(this);
			popup.add(item);
			}

		if (mFileSupport != 0) {
			popup.addSeparator();
			if ((mFileSupport & FILE_SUPPORT_OPEN_FILES) != 0) {
				item = new JMenuItem(OPEN);
				item.addActionListener(this);
				popup.add(item);
				}
			if ((mFileSupport & FILE_SUPPORT_SAVE_FILES) != 0 && mModel.getSize() != 0) {
				item = new JMenuItem(SAVE_DWAR);
				item.addActionListener(this);
				popup.add(item);
				item = new JMenuItem(SAVE_SDF2);
				item.addActionListener(this);
				popup.add(item);
				item = new JMenuItem(SAVE_SDF3);
				item.addActionListener(this);
				popup.add(item);
				}
			}

		if (mCustomPopupItemList != null) {
			popup.addSeparator();
			for (JMenuItem customItem:mCustomPopupItemList)
				popup.add(customItem);
			}

		popup.show(this, e.getX(), e.getY());
		}

	public void structureChanged(StereoMolecule mol) {
		String reason = (mCompoundFilter == null) ? null
				: (mCompoundFilter instanceof SubstructureFilter) ? "match the substructure" : "qualify";
		if (mEditedIndex == -1) {	// new structure
			if (mol.getAllAtoms() != 0) {
				if (mCompoundFilter == null || mCompoundFilter.moleculeQualifies(mol))
					mModel.addMolecule(mModel.getSize(), mol);
				else
					JOptionPane.showMessageDialog(getParentFrame(),"The compound could not be added, because it doesn't "+reason+".");
				}
			}
		else {
			if (mol.getAllAtoms() == 0)
				mModel.remove(mEditedIndex);
			else {
				if (mCompoundFilter == null || mCompoundFilter.moleculeQualifies(mol))
					mModel.setMolecule(mEditedIndex, mol);
				else
					JOptionPane.showMessageDialog(getParentFrame(),"The compound could not be changed, because the changed structure doesn't "+reason+".");
				}
			}
		}

	/**
	 * May be overridden to act on selection changes
	 * @param molIndex
	 */
	public void setSelection(int molIndex) {}

	private void validateSize() {
		Rectangle viewportBounds = getViewportBorderBounds();

		int width = mIsVertical ? viewportBounds.width : mStructureSize == 0 ? viewportBounds.height : mStructureSize;
		int height = !mIsVertical ? viewportBounds.height : mStructureSize == 0 ? viewportBounds.width : mStructureSize;

		mCellSize = new Dimension(width, height);

		if (mIsVertical) {
			height *= mModel.getSize();
			if (height < viewportBounds.height)
				height = viewportBounds.height;
			}
		else {
			width *= mModel.getSize();
			if (width < viewportBounds.width)
				width = viewportBounds.width;
			}

		if (mContentSize.width != width
		 || mContentSize.height != height) {
			mContentSize.width = width;
			mContentSize.height = height;
			mContentPanel.setPreferredSize(mContentSize);
			mContentPanel.revalidate();
			}
		}

	private void initializeDragAndDrop(int dragAction, int dropAction) {
		if (dragAction != DnDConstants.ACTION_NONE) {
			new MoleculeDragAdapter(this) {
				public Transferable getTransferable(Point p) {
					if (mHighlightedIndex == -1)
						return null;
					setCursor(SwingCursorHelper.getCursor(SwingCursorHelper.cFistCursor));
					mDragIndex = mHighlightedIndex;
					return new MoleculeTransferable(mModel.getMolecule(mHighlightedIndex));
					}
				};
			}

		if (dropAction != DnDConstants.ACTION_NONE) {
			MoleculeDropAdapter d = new MoleculeDropAdapter() {
				@Override
				public void onDropMolecule(StereoMolecule mol, Point pt) {
					if (mIsEnabled && mIsEditable && mol != null && mol.getAllAtoms() != 0 && mDropIndex != -1) {
						for (int atom=0; atom<mol.getAllAtoms(); atom++)
							mol.setAtomColor(atom, Molecule.cAtomColorNone); // don't copy atom coloring

						mol.setFragment(mCreateFragments);
						if (mCompoundFilter == null || mCompoundFilter.moleculeQualifies(mol)) {
							if (mDragIndex != -1 && mInternalDragAndDropIsMove) {
								mModel.remove(mDragIndex);
								if (mDropIndex > mDragIndex)
									mDropIndex--;
								}

							mModel.addMolecule(mDropIndex, mol);
							}
						else {
							String reason = (mCompoundFilter instanceof SubstructureFilter) ? "match the substructure" : "qualify";
							JOptionPane.showMessageDialog(getParentFrame(),"The compound could not be added, because it doesn't "+reason+".");
							}
						}
					updateDropPosition(-1);
					}

				@Override
				public void dragEnter(DropTargetDragEvent e) {
					boolean drop = mIsEnabled && mIsEditable && isDropOK(e) ;
					if (!drop) {
						e.rejectDrag();
						}
					else {
						updateDropPosition(getDropIndex(e));
						}
					}

				@Override
				public void dragOver(DropTargetDragEvent e) {
					mScroller.autoScroll();
					updateDropPosition(getDropIndex(e));
					}

				@Override
				public void dragExit(DropTargetEvent e) {
					updateDropPosition(-1);
					}

				private int getDropIndex(DropTargetDragEvent e) {
					int x = e.getLocation().x + (mIsVertical ? 0 : mCellSize.width / 2);
					int y = e.getLocation().y + (mIsVertical ? mCellSize.height / 2 : 0);
					int dropIndex = getMoleculeIndex(x, y);
					if (dropIndex == -1)
						dropIndex = mModel.getSize();

					// if we move internally onto the same position, don't indicate it
					if (mInternalDragAndDropIsMove
					 && (dropIndex == mDragIndex || dropIndex == mDragIndex+1))
						dropIndex = -1;

					return dropIndex;
					}
				};

			new DropTarget(this, dropAction, d, true, new OurFlavorMap());
			}
		}

	private void updateDropPosition(int dropIndex) {
		if (mIsEnabled && mIsEditable && mDropIndex != dropIndex) {
			mDropIndex = dropIndex;
			repaint();
			}
		}

	private Component getParentFrame() {
		Component c = this;
		while (c != null && !(c instanceof Frame) && !(c instanceof Dialog))
			c = c.getParent();
		return c;
		}

		// This class is needed for inter-jvm drag&drop. Although not neccessary for standard environments, it prevents
		// nasty "no native data was transfered" errors. It still might create ClassNotFoundException in the first place by
		// the SystemFlavorMap, but as I found it does not hurt, since the context classloader will be installed after
		// the first call. I know, that this depends heavely on a specific behaviour of the systemflavormap, but for now
		// there's nothing I can do about it.
	static class OurFlavorMap implements java.awt.datatransfer.FlavorMap {
		public java.util.Map<DataFlavor,String> getNativesForFlavors(DataFlavor[] dfs) {
			java.awt.datatransfer.FlavorMap m = java.awt.datatransfer.SystemFlavorMap.getDefaultFlavorMap();
			return m.getNativesForFlavors(dfs);
			}

		public java.util.Map<String,DataFlavor> getFlavorsForNatives(String[] natives) {
			java.awt.datatransfer.FlavorMap m = java.awt.datatransfer.SystemFlavorMap.getDefaultFlavorMap();
			return m.getFlavorsForNatives(natives);
			}
		}
	}
