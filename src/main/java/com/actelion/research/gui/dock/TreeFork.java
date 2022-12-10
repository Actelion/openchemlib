package com.actelion.research.gui.dock;

import javax.swing.*;
import javax.swing.plaf.SplitPaneUI;
import javax.swing.plaf.basic.BasicSplitPaneUI;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Vector;

public class TreeFork extends TreeContainer {
	private TreeElement mLeftChildElement;
	private TreeElement mRightChildElement;
	private Vector<DividerChangeListener> mDividerChangeListeners;

	/**
	 * Constructor to create a fork element that is inserted between the specified
	 * parent and its former oldLeaf to carry newLeaf at the specified position.
	 * @param oldLeaf
	 * @param newLeaf
	 * @param newLeafPosition
	 * @param dividerLocation
	 */
	public TreeFork(TreeLeaf oldLeaf, TreeLeaf newLeaf, int newLeafPosition, double dividerLocation, Vector<DividerChangeListener> listeners) {
		JSplitPane splitPane = null;
		switch (newLeafPosition) {
		case JDockingPanel.DOCK_TOP:
			splitPane = new MySplitPane(JSplitPane.VERTICAL_SPLIT,
					newLeaf.getComponent(), oldLeaf.getComponent(), dividerLocation);
			mLeftChildElement = newLeaf;
			mRightChildElement = oldLeaf;
			break;
		case JDockingPanel.DOCK_LEFT:
			splitPane = new MySplitPane(JSplitPane.HORIZONTAL_SPLIT,
					newLeaf.getComponent(), oldLeaf.getComponent(), dividerLocation);
			mLeftChildElement = newLeaf;
			mRightChildElement = oldLeaf;
			break;
		case JDockingPanel.DOCK_BOTTOM:
			splitPane = new MySplitPane(JSplitPane.VERTICAL_SPLIT,
					oldLeaf.getComponent(), newLeaf.getComponent(), dividerLocation);
			mLeftChildElement = oldLeaf;
			mRightChildElement = newLeaf;
			break;
		case JDockingPanel.DOCK_RIGHT:
			splitPane = new MySplitPane(JSplitPane.HORIZONTAL_SPLIT,
					oldLeaf.getComponent(), newLeaf.getComponent(), dividerLocation);
			mLeftChildElement = oldLeaf;
			mRightChildElement = newLeaf;
			break;
			}
		splitPane.addPropertyChangeListener(JSplitPane.DIVIDER_LOCATION_PROPERTY, e -> {
			if (mDividerChangeListeners != null && ((MySplitPane)mComponent).isDragged())
				for (DividerChangeListener l:mDividerChangeListeners)
					l.dividerLocationChanged(this);
			} );
		mDividerChangeListeners = listeners;
		mComponent = splitPane;
		oldLeaf.setParent(this);
		newLeaf.setParent(this);
		}

	public void updateDividerChangeListeners(Vector<DividerChangeListener> listeners) {
		mDividerChangeListeners = listeners;
		if (mLeftChildElement instanceof TreeFork)
			((TreeFork)mLeftChildElement).updateDividerChangeListeners(listeners);
		if (mLeftChildElement instanceof TreeFork)
			((TreeFork)mRightChildElement).updateDividerChangeListeners(listeners);
		}

	public void removeWithLeaf(TreeLeaf leaf) {
		assert(leaf == mLeftChildElement || leaf == mRightChildElement);
		if (leaf == mLeftChildElement)
			mParent.replaceChildElement(this, mRightChildElement);
		else
			mParent.replaceChildElement(this, mLeftChildElement);
		}

	public void replaceChildElement(TreeElement oldElement, TreeElement newElement) {
		JSplitPane splitPane = (JSplitPane)mComponent;
		int dl = splitPane.getDividerLocation();
		if (oldElement == mLeftChildElement) {
			splitPane.setLeftComponent(newElement.getComponent());
			mLeftChildElement = newElement;
			}
		else {
			splitPane.setRightComponent(newElement.getComponent());
			mRightChildElement = newElement;
			}
		splitPane.setDividerLocation(dl);
		newElement.setParent(this);
		}

	public TreeElement getLeftChild() {
		return mLeftChildElement;
		}

	public TreeElement getRightChild() {
		return mRightChildElement;
		}

	public boolean isVertical() {
		return (((JSplitPane)mComponent).getOrientation() == JSplitPane.VERTICAL_SPLIT);
		}

	public double getDividerLocation() {
		JSplitPane sp = (JSplitPane)mComponent;
		int size = (isVertical() ? sp.getHeight() : sp.getWidth()) - sp.getDividerSize();
		return (double)sp.getDividerLocation() / (double)size;
		}

	public void updateChildElement(Component oldContent, TreeElement childElement) {
		JSplitPane splitPane = (JSplitPane)mComponent;
		int dl = splitPane.getDividerLocation();
		if (childElement == mLeftChildElement)
			splitPane.setLeftComponent(childElement.getComponent());
		else
			splitPane.setRightComponent(childElement.getComponent());
		splitPane.setDividerLocation(dl);
		}

	protected void addStateInfo(ArrayList<String> stateInfo, String firstDockInfo) {
		String leftTitle = addOneLeafToStateInfo(mLeftChildElement, stateInfo, firstDockInfo);
		addOneLeafToStateInfo(mRightChildElement, stateInfo, getRightDockInfo(leftTitle));
		if (mLeftChildElement instanceof TreeFork)
			((TreeFork)mLeftChildElement).addStateInfo(stateInfo, firstDockInfo);
		if (mRightChildElement instanceof TreeFork)
			((TreeFork)mRightChildElement).addStateInfo(stateInfo, firstDockInfo);
		}

	private String addOneLeafToStateInfo(ArrayList<String> stateInfoList, String firstDockInfo) {
		return addOneLeafToStateInfo(mLeftChildElement, stateInfoList, firstDockInfo);
		}

	private String addOneLeafToStateInfo(TreeElement child, ArrayList<String> stateInfoList, String firstDockInfo) {
		if (child instanceof TreeLeaf)
			return ((TreeLeaf)child).addStateInfo(stateInfoList, firstDockInfo);
		else
			return ((TreeFork)child).addOneLeafToStateInfo(stateInfoList, firstDockInfo);
		}

	protected void clearStateInfo() {
		mLeftChildElement.clearStateInfo();
		mRightChildElement.clearStateInfo();
		}

	private String getRightDockInfo(String leftTitle) {
		return leftTitle
			 + "\t" + (isVertical()?"bottom":"right")
			 + "\t" + ((double)((int)(1000*getDividerLocation()))/1000);
		}

	public void printStatus() {
		JSplitPane splitPane = (JSplitPane)mComponent;
		String status1 = (splitPane.getLeftComponent() == mLeftChildElement.getComponent()) ? "OK" : "failed";
		String status2 = (splitPane.getRightComponent() == mRightChildElement.getComponent()) ? "OK" : "failed";
		System.out.println("Fork leftChildStatus:"+status1+" rightChildStatus:"+status2);
		mLeftChildElement.printStatus();
		mRightChildElement.printStatus();
		}
	}

/**
 * JSplitPane.setDividerLocation() only works reliably after the layout is finished
 * and first painting was done.
 * MySplitPane allows to set the divider location before the splitpane is visible.
 * JSplitPane does its one funny stuff to determine a divider location during layout
 * and first painting, which sometimes results in locations substantially different
 * from the defined proportional value. We re-enforce setDividerLocation() until the
 * wanted location is reached.
 * @author sandert
 */
class MySplitPane extends JSplitPane {
	private static final long serialVersionUID = 0x20070726;

	private double mProportionalLocation;
	private boolean mMouseDown;

	public MySplitPane(int newOrientation, Component newLeftComponent, Component newRightComponent, double proportionalLocation) {
		super(newOrientation, true, newLeftComponent, newRightComponent);
		proportionalLocation = Math.min(1, Math.max(0, proportionalLocation));
		setDividerLocation(proportionalLocation);
		setResizeWeight(proportionalLocation);
		setBorder(null);
		mProportionalLocation = proportionalLocation;

		SplitPaneUI spui = getUI();
		if (spui instanceof BasicSplitPaneUI) {
			((BasicSplitPaneUI) spui).getDivider().addMouseListener(new MouseAdapter() {
				public void mousePressed(MouseEvent e) {
					mMouseDown = true;
					}
				public void mouseReleased(MouseEvent e) {
					mMouseDown = false;
					}
				} );
			}
		}

	public boolean isDragged() {
		return mMouseDown;
		}

	@Override
	public void paintComponent(Graphics g) {
		if (mProportionalLocation != -1) {
			int currentSize = (getOrientation() == JSplitPane.VERTICAL_SPLIT ? getHeight() : getWidth()) - getDividerSize();
			int wantedDeviderLocation = (int)(mProportionalLocation * currentSize);
			if (Math.abs(wantedDeviderLocation - getDividerLocation()) > 1)
				super.setDividerLocation(wantedDeviderLocation);
			else
				mProportionalLocation = -1;	// setting was successful
			}
		super.paintComponent(g);
		}
	}
