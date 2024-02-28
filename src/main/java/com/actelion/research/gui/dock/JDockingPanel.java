package com.actelion.research.gui.dock;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;
import java.util.*;

public class JDockingPanel extends JPanel implements ActionListener {
	private static final long serialVersionUID = 0x20070720;

	public  static final int DOCK_CENTER = 0;
	public  static final int DOCK_TOP = 1;
	public  static final int DOCK_LEFT = 2;
	public  static final int DOCK_BOTTOM = 3;
	public  static final int DOCK_RIGHT = 4;

	public static final int ALLOWED_DRAG_DROP_ACTIONS = DnDConstants.ACTION_MOVE;   // currently no copy

	private TreeMap<String,Dockable> mDockableMap;
	private TreeMap<String,TreeLeaf> mLeafMap;
	private Dockable mPreviousTargetDockable;
	private TreeRoot mTreeRoot;
	private TreeLeaf mTargetLeaf;
	private int mTargetPosition,mPreviousTargetPosition;
	private GhostPreview mPreview;
	private Dockable mMaximizedView;
	private Vector<DividerChangeListener> mDividerChangeListeners;

	/**
	 * Creates a docking panel to which any Dockables may be added by
	 * the respective dock methods.
	 */
	public JDockingPanel() {
		mDockableMap = new TreeMap<>();
		mLeafMap = new TreeMap<>();
		mTargetPosition = -1;
		mPreview = new GhostPreview();

		// Unfortunately, this cannot be used to catch drag events of default JTable with active DropTarget
//		Toolkit.getDefaultToolkit().addAWTEventListener(this, AWTEvent.MOUSE_MOTION_EVENT_MASK);

		setLayout(new OverlayLayout(this));

		new DropTarget(this, ALLOWED_DRAG_DROP_ACTIONS, new DropTargetAdapter() {
			@Override
			public void dragOver(DropTargetDragEvent dtde) {
				Dockable dockable = getDraggedDockable(dtde);
				if (dockable != null) {
					Point p = dtde.getLocation();
					Rectangle b = dockable.getHeader().getBounds();
					p.translate(b.x, b.y);  // make p relative to Dockable
					updatePreview(p, dockable);
					}
				}

			@Override
			public void drop(DropTargetDropEvent dtde) {
				Transferable transferable = dtde.getTransferable();
				if (mTargetPosition == -1
				 || !transferable.isDataFlavorSupported(TransferableDockable.DF_DOCKABLE_DEF)) {
// TODO do we need this?				fireDockableSelected(mDraggedDockable);
					dtde.rejectDrop();
					}
				else {
					dtde.acceptDrop(ALLOWED_DRAG_DROP_ACTIONS);
					try {
						String draggedTitle = (String)transferable.getTransferData(TransferableDockable.DF_DOCKABLE_DEF);
						for (int i=0; i<mTargetLeaf.getDockableCount(); i++) {
							if (!draggedTitle.equals(mTargetLeaf.getDockable(i).getTitle())) {
								String targetTitle = mTargetLeaf.getDockable(i).getTitle();
								relocateView(draggedTitle, targetTitle, mTargetPosition, 0.5f);
								break;
								}
							}
						dtde.dropComplete(true);
						}
					catch (Exception e) {}
					}

				mTargetPosition = -1;
				}
			}, true);
		}

	public void addDividerChangeLister(DividerChangeListener l) {
		if (mDividerChangeListeners == null)
			mDividerChangeListeners = new Vector<>();

		mDividerChangeListeners.add(l);
		if (mTreeRoot != null)
			mTreeRoot.setDividerChangeListeners(mDividerChangeListeners);
		}

	public void removeDividerChangeLister(DividerChangeListener l) {
		if (mDividerChangeListeners != null) {
			mDividerChangeListeners.remove(l);
			if (mTreeRoot != null)
				mTreeRoot.setDividerChangeListeners(mDividerChangeListeners);
			}
		}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().startsWith("close_")) {
			String title = e.getActionCommand().substring(6);
			undock(title, false);
			}
		else if (e.getActionCommand().startsWith("max_")) {
			String title = e.getActionCommand().substring(4);
			maximize(title, null);
			}
		}

	/**
	 * Toggle the maximization state of the dockable in the docking panel.
	 * If it is maximizing, then a toolbar may be provided that is shown in the maximized state.
	 * Note that a null toolbar causes the maximized state to use the default toolbar with
	 * maximize and close buttons only.
	 * @param maximizeToolBar null for default or toolbar that shall be shown in maximized state
	 */
	public void maximize(String title, JToolBar maximizeToolBar) {
		Dockable dockable = mDockableMap.get(title);
		if (mMaximizedView != null) {
			remove(mMaximizedView);
			dockable.endBorrowContent();
			mTreeRoot.getChild().getComponent().setVisible(true);
			mMaximizedView = null;
			validate();
			}
		else if (!(mTreeRoot.getComponent() instanceof Dockable)) {
			mTreeRoot.getChild().getComponent().setVisible(false);
			selectDockable(dockable);
			JComponent content = dockable.borrowContent();

			// create dummy dockable to hold the maximized view
			mMaximizedView = new Dockable(this, content, title, maximizeToolBar);
			mMaximizedView.setPopupProvider(dockable.getPopupProvider());
			mMaximizedView.setVisible(true);
			add(mMaximizedView);

			validate();
			}
		}

	/**
	 * Sends ActionEvent to listeners when the user interactively selected another dockable
	 * @param dockable
	 */
	protected void fireDockableSelected(Dockable dockable) {
		actionPerformed(new ActionEvent(this, ActionEvent.ACTION_PERFORMED, "selected_"+dockable.getTitle()));
		}

	private void updatePreview(Point p, Dockable draggedDockable) {
		mTargetPosition = -1;
		Dockable targetDockable = null;

		for (Dockable d:mDockableMap.values()) {
			if (d.isVisible()) {
				Rectangle bounds = getAbsoluteBounds(d.getContent());
				if (bounds.contains(p)) {
					targetDockable = d;
					mTargetPosition = getPosition(p, bounds);
					break;
					}
				}
			}

		updatePreview(draggedDockable, targetDockable);
		}

	/**
	 * Calculates bounds of the dockable relative to this JDockingPanel
	 * @param jc
	 * @return
	 */
	public Rectangle getAbsoluteBounds(Component jc) {
		Rectangle b = jc.getBounds();
		Component c = jc.getParent();
		while (c != this) {
			Point p = c.getLocation();
			b.translate(p.x, p.y);
			c = c.getParent();
			}
		return b;
		}

	private int getPosition(Point p, Rectangle bounds) {
		int x = p.x - bounds.x;
		int y = p.y - bounds.y;
		if (x > bounds.width/4  && x < bounds.width*3/4
		 && y > bounds.height/4 && y < bounds.height*3/4)
			return DOCK_CENTER;

		boolean isTopOrRight = (((double)x / (double)bounds.width) > ((double)y / (double)bounds.height));
		boolean isTopOrLeft =  (((double)x / (double)bounds.width) + ((double)y / (double)bounds.height) < 1.0);

		if (isTopOrRight)
			return isTopOrLeft ? DOCK_TOP : DOCK_RIGHT;

		return isTopOrLeft ? DOCK_LEFT : DOCK_BOTTOM;
		}

	private void updatePreview(Dockable draggedDockable, Dockable targetDockable) {
		if (targetDockable == null) {
			mTargetPosition = -1;
			}
		else if (targetDockable == draggedDockable) {
			TreeLeaf sourceLeaf = mLeafMap.get(draggedDockable.getTitle());
			if (mTargetPosition == DOCK_CENTER
			 || sourceLeaf.getDockableCount() == 1)
				mTargetPosition = -1;
			else
				mTargetLeaf = sourceLeaf;
			}
		else {
			mTargetLeaf = mLeafMap.get(targetDockable.getTitle());
			}

		if (mPreviousTargetDockable != targetDockable
		 || mPreviousTargetPosition != mTargetPosition) {
			if (mTargetPosition != -1)
				mPreview.createPreview(draggedDockable, targetDockable, mTargetPosition, this);
			mPreviousTargetDockable = targetDockable;
			mPreviousTargetPosition = mTargetPosition;
			repaint();
			}
		}

	@Override
	public void paint(Graphics g) {
		super.paint(g);
		if (mTargetPosition != -1)
			mPreview.drawPreview((Graphics2D)g);
		}

/*	private void performDrag() {
		undock(mDraggedDockable.getTitle(), true);
		dock(mDraggedDockable, mTargetPosition, mTargetLeaf, 0.5, true);
		} */

	/**
	 * May be overriden to catch user initiated drag&drop
	 * @param movedDockableName
	 * @param targetDockableName
	 * @param targetPosition
	 */
	public void relocateView(String movedDockableName, String targetDockableName, int targetPosition, float dividerPosition) {
		Dockable movedDockable = getDockable(movedDockableName);
		undock(movedDockableName, true);
		dock(movedDockable, targetPosition, mLeafMap.get(targetDockableName), dividerPosition, true);
		}

	/**
	 * Dock the dockable to the center of the selected view.
	 * @param dockable
	 */
	public void dock(Dockable dockable) {
		dock(dockable, DOCK_CENTER);
		}

	/**
	 * Dock the Dockable at the given position to the selected view.
	 * @param dockable
	 * @param position
	 */
	public void dock(Dockable dockable, int position) {
		dock(dockable, position, getSelectedLeaf(), 0.5, false);
		}

	/**
	 * Dock the dockable either as a new root dockable or
	 * to the treeElement that contains dockable specified as part dockInfo.
	 * The dockInfo parameter is either "root" or has the following format:
	 * "title[\tposition[\tdividerlocation]]" with position is one of center,top,left,bottom,right
	 * and dividerlocation is between 0.0 and 1.0 giving the percentage of space given to the left/top
	 * element.
	 * @param dockable
	 * @param dockInfo encoded position and parent obtainable by getDockInfo(title)
	 */
	public void dock(Dockable dockable, String dockInfo) {
		if (dockInfo == null || dockInfo.equals("root")) {
			dock(dockable, DOCK_CENTER, getSelectedLeaf(), Double.NaN, false);
			return;
			}

		int index1 = dockInfo.indexOf('\t');
		int index2 = dockInfo.indexOf('\t', index1+1);
		String title = dockInfo.substring(0, index1);
		String p = (index2 == -1) ? dockInfo.substring(index1+1) : dockInfo.substring(index1+1, index2);
		String dp = (index2 == -1) ? null : dockInfo.substring(index2);

		int position = "center".equals(p) ? DOCK_CENTER
					 : "top".equals(p)	? DOCK_TOP
					 : "left".equals(p)   ? DOCK_LEFT
					 : "bottom".equals(p) ? DOCK_BOTTOM
					 : "right".equals(p)  ? DOCK_RIGHT
					 : -1;
		double dividerPosition = (dp != null) ? Double.parseDouble(dp)
				: (position == DOCK_CENTER) ? Double.NaN : 0.5;

		dock(dockable, position, mLeafMap.get(title), dividerPosition, false);
		}

	/**
	 * Dock the dockable at the given position to the treeElement's component.
	 * If component is the first to dock then treeElement may be null.
	 * @param dockable
	 * @param position DOCK_CENTER,DOCK_LEFT,...
	 * @param treeLeaf may be null, if the dockable is the first in the JDockingPanel
	 * @param dividerPosition position of divider of new JSplitPane if position!=DOCK_CENTER
	 */
	private void dock(Dockable dockable, int position, TreeLeaf treeLeaf, double dividerPosition, boolean isDragging) {
		mDockableMap.put(dockable.getTitle(), dockable);

		if (treeLeaf == null) {
			TreeLeaf leaf = new TreeLeaf(dockable, this, isDragging);
			mTreeRoot = new TreeRoot(this, leaf);
			mLeafMap.put(dockable.getTitle(), leaf);
			}
		else if (position == DOCK_CENTER) {
			treeLeaf.addContent(dockable, isDragging);
			mLeafMap.put(dockable.getTitle(), treeLeaf);
			}
		else {
			TreeLeaf newLeaf = new TreeLeaf(dockable, this, isDragging);
			TreeContainer parent = treeLeaf.getParent();
			TreeFork treeFork = new TreeFork(treeLeaf, newLeaf, position, dividerPosition, mDividerChangeListeners);
			parent.replaceChildElement(treeLeaf, treeFork);
			mLeafMap.put(dockable.getTitle(), newLeaf);
			}

		selectDockable(dockable);
		validate();
//System.out.println("*****Status after dock()");
//mTreeRoot.printStatus();
		}

	public void undock(String title) {
		undock(title, false);
		}

	public void undock(String title, boolean isDragging) {
		Dockable dockable = mDockableMap.get(title);
		if (dockable != null) {
			if(mDockableMap.size()==1) {
				undockAll();
				return;				
			}
			boolean isSelected = dockable.isSelected();
			TreeLeaf leaf = mLeafMap.get(title);
			mLeafMap.remove(title);
			if (leaf.removeContent(dockable, isDragging)) {
				if (isSelected) {
					TreeLeaf newSelected = mLeafMap.get(mLeafMap.firstKey());
					if (newSelected != null)
						newSelected.setSelected(true);
					}
				}
			else {
				if (isSelected)
					leaf.setSelected(true);
				}

			mDockableMap.remove(title);

			validate();
			repaint();
			}
//System.out.println("*****Status after undock()");
//mTreeRoot.printStatus();
		}

	public void undockAll() {
		mTreeRoot = null;
		mDockableMap.clear();
		mLeafMap.clear();
		removeAll();
		validate();
		repaint();
		}

	public Set<String> getDockableTitles() {
		return mDockableMap.keySet();
		}

	public Collection<Dockable> getDockables() {
		return mDockableMap.values();
		}

	public Dockable getDockable(String title) {
		return mDockableMap.get(title);
		}

	public String getTitle(Dockable dockable) {
		for (String title:mDockableMap.keySet())
			if (mDockableMap.get(title) == dockable)
				return title;
		return null;
		}

	public int getDockableCount() {
		return mDockableMap.size();
		}

	/**
	 * returns the persistence information necessary to recreate the current docking state.
	 * To recreate a docking state create an empty JDockingPanel and call dock(dockable, dockInfo)
	 * once for every dockable keeping the order of the title[] array and pass the respective
	 * dockInfo from the array returned by this function.
	 * @return dockInfo[] array with a docking instruction for every dockable
	 */
	public String[] getDockInfoSequence() {
		if (mTreeRoot == null)
			return null;

		return mTreeRoot.createStateInfo().toArray(new String[0]);
		}

	/**
	 * @param title
	 * @return true if a dockable with the given title exists and if it is the selected component of a tabbed pane
	 */
	public boolean isInFrontInTabbedPane(String title) {
		Dockable dockable = mDockableMap.get(title);
		if (dockable == null)
			return false;

		Component parent = dockable.getParent();
		if (parent instanceof JTabbedPane)
			return (dockable == ((JTabbedPane)parent).getSelectedComponent());

		return false;
		}

	/**
	 * @param title
	 */
	public void setToFrontInTabbedPane(String title) {
		Dockable dockable = mDockableMap.get(title);
		if (dockable != null) {
			Component parent = dockable.getParent();
			if (parent instanceof JTabbedPane)
				((JTabbedPane)parent).setSelectedComponent(dockable);

			}
		}

	/**
	 * Checks, whether a dockable with the given title exists, and whether it is visible (selected if in a JTabbedPane)
	 * and whether it is in the left or right branch (depending on parameter left) of a JSplitPane.
	 * @param title
	 * @param left if true, then
	 * @return true the respective dockable is visible and in the left part of a JSplitPane
	 */
	public boolean isVisibleInSplitPane(String title, boolean left) {
		Component view = mDockableMap.get(title);
		if (view == null)
			return false;

		Component pane = view.getParent();
		if (pane instanceof JTabbedPane) {
			if (view != ((JTabbedPane)pane).getSelectedComponent())
				return false;

			view = pane;
			pane = pane.getParent();
			}

		if (pane instanceof JSplitPane)
			return (left && view == ((JSplitPane)pane).getLeftComponent())
				|| (!left && view == ((JSplitPane)pane).getRightComponent());

		return false;
		}

	/**
	 * changes a title of a docked Dockable
	 * @param oldTitle existing title
	 * @param newTitle must be a unique title and must not contain a TAB
	 * @return true in case of success
	 */
	public boolean changeTitle(String oldTitle, String newTitle) {
		if (mDockableMap.containsKey(newTitle) || newTitle.indexOf('\t') != -1)
			return false;

		Dockable dockable = mDockableMap.get(oldTitle);
		mDockableMap.remove(oldTitle);
		mDockableMap.put(newTitle, dockable);
		dockable.setTitle(newTitle);

		TreeLeaf leaf = mLeafMap.get(oldTitle);
		leaf.changeTitle(oldTitle, newTitle);
		mLeafMap.remove(oldTitle);
		mLeafMap.put(newTitle, leaf);

		if (mMaximizedView != null
		 && mMaximizedView.getTitle().equals(oldTitle))
			mMaximizedView.setTitle(newTitle);

		return true;
		}

	/**
	 * Returns the Dockable that contains the given content.
	 * @param content the content to be looked for
	 * @return content owning Dockable or null if content is not found
	 */
	public Dockable getDockable(Component content) {
		if (content != null)
			for (Dockable d:mDockableMap.values())
				if (d.getContent() == content)
					return d;

		return null;
		}

	/**
	 * Returns the currently active Dockable, i.e. the one with the
	 * highlighted header area.
	 * @return active Dockable or null if there are no Dockables
	 */
	public Dockable getSelectedDockable() {
		for (Dockable d:mDockableMap.values())
			if (d.isSelected())
				return d;

		return null;
		}

	public Dockable getDraggedDockable(DropTargetDragEvent dtde) {
		Transferable transferable = dtde.getTransferable();
		if (transferable.isDataFlavorSupported(TransferableDockable.DF_DOCKABLE_DEF)) {
			try {
				return getDockable((String)transferable.getTransferData(TransferableDockable.DF_DOCKABLE_DEF));
				}
			catch (Exception e) {}
			}
		return null;
		}

	/**
	 * This is called when a Dockable's visibility changes, i.e.<br>
	 * - when a new Dockable is docked<br>
	 * - when a new Dockable is undocked<br>
	 * - when the user drags & docks a Dockable<br>
	 * - when the user actively switches tabs in a tabbed pane<br>
	 * It may be overridden in order to react on a visibility change due
	 * to the user selecting a different tab of dragging & docking a
	 * Dockable.
	 * @param dockable
	 * @param isVisible
	 */
	public void visibilityChanged(Dockable dockable, boolean isVisible) {
		dockable.notifyVisibility(isVisible);
		}

	public void selectDockable(Dockable dockable) {
		if (dockable != null && !isMaximized()) {
			for (Dockable d:mDockableMap.values())
				d.setSelected(d == dockable);
	
			TreeLeaf leaf = mLeafMap.get(dockable.getTitle());
			if (leaf != null)
				leaf.setSelectedDockable(dockable);
			}
		}

	private TreeLeaf getSelectedLeaf() {
		for (TreeLeaf leaf:mLeafMap.values())
			if (leaf.isSelected())
				return leaf;

		return null;
		}

	/**
	 * @return whether one of the dockable views is currently maximized
	 */
	public boolean isMaximized() {
		return (mMaximizedView != null);
		}

	public void redistribute() {
		List<Dockable> dockables = new ArrayList<Dockable>(getDockables());
		Collections.sort(dockables, new Comparator<Dockable>() {
			public int compare(Dockable o1, Dockable o2) {
				return o1.getTitle().compareTo(o2.getTitle());
				}
			});

		if(dockables.size()<=1)
			return;

		mTreeRoot = null;
		mLeafMap.clear();
		mDockableMap.clear();
		removeAll();
		validate();
		
		for (Dockable dockable : dockables) {
			dockable.borrowContent();
			}
		
		if(dockables.size()>0) {
			Dockable dockable = dockables.get(0);
			TreeLeaf leaf = new TreeLeaf(dockable, this, false);
			mTreeRoot = new TreeRoot(this, leaf);
			mLeafMap.put(dockable.getTitle(), leaf);
			mDockableMap.put(dockable.getTitle(), dockable);
			redistribute(dockables, 0, dockables.size(), leaf);
			}

		for (Dockable dockable : dockables) {
			dockable.endBorrowContent();
			}

		validate();
		repaint();
		}
	
	private void redistribute(List<Dockable> dockables, int startIndex, int endIndex, TreeLeaf topLeft) {

		//int size = (endIndex-startIndex + 3)/4;
		//if(size>1 && size%2==1) size++;
		//size = (int)(Math.sqrt(size)+.99);
		//size*= size;
		//System.out.println("JDockingPanel.redistribute() "+startIndex+" "+endIndex+" size="+size);

		int size1 = startIndex + (endIndex-startIndex+3)/4;
		int size2 = size1 + (endIndex-size1+2)/3;
		int size3 = size2 + (endIndex-size2+1)/2;
		System.out.println(startIndex+" "+size1+" "+size2+" "+size3+" "+endIndex );

		TreeLeaf topRight = null;
		TreeLeaf bottomLeft = null;
		TreeLeaf bottomRight = null;

		if(size2>size1 && size2<dockables.size()) {
			Dockable dockable = dockables.get(size2);
			bottomLeft = new TreeLeaf(dockable, this, false);
			TreeContainer parent = topLeft.getParent();
			TreeFork treeFork = new TreeFork(topLeft, bottomLeft, DOCK_BOTTOM, .5, mDividerChangeListeners);
			parent.replaceChildElement(topLeft, treeFork);
			mLeafMap.put(dockable.getTitle(), bottomLeft);
			mDockableMap.put(dockable.getTitle(), dockable);
			}

		if(size1>startIndex && size1<dockables.size()) {
			Dockable dockable = dockables.get(size1);
			topRight = new TreeLeaf(dockable, this, false);
			TreeContainer parent = topLeft.getParent();
			TreeFork treeFork = new TreeFork(topLeft, topRight, DOCK_RIGHT, .5, mDividerChangeListeners);
			parent.replaceChildElement(topLeft, treeFork);
			mLeafMap.put(dockable.getTitle(), topRight);
			mDockableMap.put(dockable.getTitle(), dockable);
			}

		if(size3>size2 && size3<dockables.size()) {
			Dockable dockable = dockables.get(size3);
			bottomRight = new TreeLeaf(dockable, this, false);
			TreeContainer parent = bottomLeft.getParent();
			TreeFork treeFork = new TreeFork(bottomLeft, bottomRight, DOCK_RIGHT, .5, mDividerChangeListeners);
			parent.replaceChildElement(bottomLeft, treeFork);
			mLeafMap.put(dockable.getTitle(), bottomRight);
			mDockableMap.put(dockable.getTitle(), dockable);
			}

		if (size1>startIndex+1)
			redistribute(dockables, startIndex, size1, topLeft);
		if (size2>size1+1)
			redistribute(dockables, size1, size2, topRight);
		if (size3>size2+1)
			redistribute(dockables, size2, size3, bottomLeft);
		if (endIndex>size3+1)
			redistribute(dockables, size3, endIndex, bottomRight);
		}
	}
