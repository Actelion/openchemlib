package com.actelion.research.gui.dock;

import com.actelion.research.gui.hidpi.HiDPIHelper;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;

public class TreeLeaf extends TreeElement implements ChangeListener {
	private JDockingPanel   mDockingPanel;
	private Dockable		mVisibleDockable;
	private String		  mStateTitle;
	private boolean		 mIsAutomatedStateChange;

	/**
	 * Constructor to create a leaf element with a visible component on the screen.
	 * @param dockable the first Dockable of this TreeLeaf
	 * @param dockingPanel
	 * @param isDragging
	 */
	public TreeLeaf(Dockable dockable, JDockingPanel dockingPanel, boolean isDragging) {
		mComponent = dockable;
		mVisibleDockable = dockable;
		mDockingPanel = dockingPanel;

		if (!isDragging)
			mDockingPanel.visibilityChanged(dockable, true);
		}
 
	public void stateChanged(ChangeEvent e) {
		if (mComponent instanceof JTabbedPane) {
			JTabbedPane tabbedPane = (JTabbedPane)e.getSource();
			if (mVisibleDockable != tabbedPane.getSelectedComponent()) {
				Dockable deselected = mVisibleDockable;
				mVisibleDockable = (Dockable)tabbedPane.getSelectedComponent();

				if (!mIsAutomatedStateChange) {
					mDockingPanel.visibilityChanged(deselected, false);
					mDockingPanel.visibilityChanged(mVisibleDockable, true);
					mDockingPanel.fireDockableSelected(mVisibleDockable);
					}

				mDockingPanel.selectDockable(mVisibleDockable);
				}
			}
		// if after removal of second last tab component the JTabbedPane was removed
		else if (mComponent instanceof Dockable) {
			Dockable deselected = mVisibleDockable;
			mVisibleDockable = (Dockable)mComponent;

			if (!mIsAutomatedStateChange) {
				mDockingPanel.visibilityChanged(deselected, false);
				mDockingPanel.visibilityChanged(mVisibleDockable, true);
				}
			mDockingPanel.selectDockable(mVisibleDockable);
			}
		}

	public boolean isSelected() {
		if (mComponent instanceof JTabbedPane) {
			JTabbedPane tabbedPane = (JTabbedPane)mComponent;
			for (int i=0; i<tabbedPane.getTabCount(); i++)
				if (((Dockable)tabbedPane.getComponentAt(i)).isSelected())
					return true;
			return false;
			}

		if (mComponent instanceof Dockable)
			return ((Dockable)mComponent).isSelected();

		return false;   // should not occur
		}

	public void setSelected(boolean b) {
		mVisibleDockable.setSelected(b);
		}

	public void addContent(Dockable dockable, boolean isDragging) {
		if (mComponent instanceof Dockable) {
			Dockable existingDockable = (Dockable)mComponent;
			JTabbedPane tabbedPane = new JTabbedPane(JTabbedPane.BOTTOM) {
				@Override public Dimension getMinimumSize() {
					// of not modified, the minimum size == preferred size.
					// Then width is at least the width of the longest tab name.
					return new Dimension(HiDPIHelper.scale(100),HiDPIHelper.scale(100));
					}
				};
			tabbedPane.putClientProperty("Quaqua.TabbedPane.contentBorderPainted", Boolean.FALSE);
			tabbedPane.add(existingDockable, existingDockable.getTitle());
			tabbedPane.add(dockable, dockable.getTitle());
			tabbedPane.addChangeListener(this);
			tabbedPane.addMouseListener(new MouseAdapter() {
				public void mousePressed(MouseEvent e) {
					handlePopupTrigger(e);
					}

				public void mouseReleased(MouseEvent e) {
					handlePopupTrigger(e);
					}

				private void handlePopupTrigger(MouseEvent e) {
					if (e.isPopupTrigger()) {
						JTabbedPane tp = (JTabbedPane)e.getSource();
						int index = tp.indexAtLocation(e.getX(), e.getY());
						if (index != -1) {
							Dockable dockable = (Dockable)tp.getComponentAt(index);
							PopupProvider pp = dockable.getPopupProvider();
							if (pp != null) {
								JPopupMenu popup = pp.createPopupMenu(dockable.getTitle(), dockable.isMaximized());
								popup.show(tp, e.getX(), e.getY());
								}
							}
						}
					}
				});

			mComponent = tabbedPane;
			mParent.updateChildElement(existingDockable, this);
			}
		if (mComponent instanceof JTabbedPane) {
			JTabbedPane tabbedPane = (JTabbedPane)mComponent;
			mDockingPanel.visibilityChanged(mVisibleDockable, false);

			// we need to recognize stateChanged() calls caused by tabbedPane.add()
			mIsAutomatedStateChange = true;
			tabbedPane.add(dockable, dockable.getTitle());
			tabbedPane.setSelectedComponent(dockable);
			mIsAutomatedStateChange = false;
			}

		if (!isDragging)
			mDockingPanel.visibilityChanged(dockable, true);

		mVisibleDockable = dockable;
		}

	/**
	 * undockes and removes a Dockable from the docking tree and from the component tree
	 * @param dockable
	 * @param isDragging true if the removed dockable is immediately docked somewere else as result of a drag&dock operation
	 * @return true if this TreeLeaf is empty after the removal of the Dockable
	 */
	public boolean removeContent(Dockable dockable, boolean isDragging) {
		if (!isDragging)
			mDockingPanel.visibilityChanged(dockable, false);

		if (mComponent instanceof JTabbedPane) {
			JTabbedPane tabbedPane = (JTabbedPane)mComponent;

			// we need to recognize stateChanged() calls caused by tabbedPane.remove()
			mIsAutomatedStateChange = true;

			try {   // JFXPanel.getInputMethodRequests() may cause a NullPointerException, e.g. when:
					// Open simple file with table view, create explanation view on top of table view,
					// drag explanation down to show both views at the same time: the JTabbedPane.remove()
					// method cascade contains JFXPanel.getInputMethodRequests() causing a NullPointerException
					// (unfixed in Bellsoft OpenJDK 8 232)
				tabbedPane.remove(dockable);
				}
			catch (Exception e) {}

			mIsAutomatedStateChange = false;

			if (tabbedPane.getTabCount() == 1) {
				mComponent = (JComponent)tabbedPane.getComponentAt(0);
				mVisibleDockable = (Dockable)mComponent;
				mParent.updateChildElement(tabbedPane, this);
				}
			else {
				mVisibleDockable = (Dockable)tabbedPane.getSelectedComponent();
				}
			mDockingPanel.visibilityChanged(mVisibleDockable, true);
			return false;
			}

		mParent.removeWithLeaf(this);
		return true;
		}

	public int getDockableCount() {
		return (mComponent instanceof JTabbedPane) ? ((JTabbedPane)mComponent).getTabCount() : 1;
		}

	/**
	 * adds docking state entries of this leaf's Dockables
	 * @param stateInfoList
	 * @param firstDockableState
	 * @return title of last Dockable to serve as reference for further stateInfo entries
	 */
	protected String addStateInfo(ArrayList<String> stateInfoList, String firstDockableState) {
		if (mStateTitle == null) {
			if (mComponent instanceof JTabbedPane) {
				JTabbedPane tabbedPane = (JTabbedPane)mComponent;
				for (int i=0; i<tabbedPane.getTabCount(); i++) {
					String title = tabbedPane.getTitleAt(i);
					if (i == 0)
						stateInfoList.add(title+"\t"+firstDockableState);
					else
						stateInfoList.add(title+"\t"+mStateTitle+"\tcenter");
					mStateTitle = title;
					}
				
				}
			else if (mComponent instanceof Dockable) {
				mStateTitle = ((Dockable)mComponent).getTitle();
				stateInfoList.add(mStateTitle+"\t"+firstDockableState);
				}
			}

		return mStateTitle;
		}

	public Dockable getDockable(int index) {
		if (mComponent instanceof JTabbedPane) {
			JTabbedPane tabbedPane = (JTabbedPane)mComponent;
			return (Dockable)tabbedPane.getComponentAt(index);
			}

		if (mComponent instanceof Dockable) {
			return (index == 0) ? (Dockable)mComponent : null;
			}

		return null;
		}

	public Dockable getDragable(Point p) {
		if (mComponent instanceof JTabbedPane) {
			JTabbedPane tabbedPane = (JTabbedPane)mComponent;
			for (int i=0; i<tabbedPane.getTabCount(); i++) {
				Dockable dockable = (Dockable)tabbedPane.getComponentAt(i);
				if (dockable.getDragHandle().getBounds().contains(p))
					return dockable;
				}
			}

		if (mComponent instanceof Dockable) {
			Dockable dockable = (Dockable)mComponent;
			if (dockable.getDragHandle().getBounds().contains(p))
				return dockable;
			}

		return null;
		}

	public Rectangle getBounds() {
		return mComponent.getBounds();
		}

	public void setSelectedDockable(Dockable dockable) {
		if (mComponent instanceof JTabbedPane
		 && mVisibleDockable != dockable) {
			mDockingPanel.visibilityChanged(mVisibleDockable, false);

			JTabbedPane tabbedPane = (JTabbedPane)mComponent;
			mIsAutomatedStateChange = true;
			tabbedPane.setSelectedComponent(dockable);	// updates mVisibleDockable through stateChanged()
			mIsAutomatedStateChange = false;

			mDockingPanel.visibilityChanged(dockable, true);
			}
		}

	public void changeTitle(String oldTitle, String newTitle) {
		if (mComponent instanceof JTabbedPane) {
			JTabbedPane tabbedPane = (JTabbedPane)mComponent;
			for (int i=0; i<tabbedPane.getTabCount(); i++) {
				if (tabbedPane.getTitleAt(i).equals(oldTitle)) {
					tabbedPane.setTitleAt(i, newTitle);
					break;
					}
				}
			}
		}

	protected void clearStateInfo() {
		mStateTitle = null;
		}

	public void printStatus() {
		System.out.print("Leaf ");
		if (mComponent instanceof JTabbedPane) {
			JTabbedPane tabbedPane = (JTabbedPane)mComponent;
			for (int i=0; i<tabbedPane.getTabCount(); i++) {
				Dockable dockable = (Dockable)tabbedPane.getComponentAt(i);
				System.out.print(dockable.getTitle()+",");
				}
			}
		else {
			Dockable dockable = (Dockable)mComponent;
			System.out.print(dockable.getTitle());
			}
		System.out.println();
		}
	}
