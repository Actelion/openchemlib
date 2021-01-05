package com.actelion.research.gui.dock;

import com.actelion.research.gui.HeaderPaintHelper;
import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.gui.hidpi.HiDPIIconButton;
import info.clearthought.layout.TableLayout;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

public class DockableHeader extends JPanel implements ActionListener {
	private static final long serialVersionUID = 0x20070723;

	private static final double[][] SIZE = {{TableLayout.FILL,TableLayout.PREFERRED},{TableLayout.FILL,TableLayout.PREFERRED,TableLayout.FILL}};

	private Dockable mDockable;
	private JLabel mTitleLabel;
	private ActionListener mActionListener;
	private HeaderMouseAdapter mMouseAdapter;
	private boolean mIsSelected;

	public DockableHeader(Dockable dockable, String title, ActionListener al, boolean isClosable, boolean hasPopupButton) {
		super(new TableLayout(SIZE));

		mDockable = dockable;
		mTitleLabel = new JLabel(title, SwingConstants.LEADING) /*{
			private static final long serialVersionUID = 0x20080423;
			public Dimension getPreferredSize() {
				return new Dimension(super.getPreferredSize().width, 16);
				}
			}*/;
		mTitleLabel.setBorder(BorderFactory.createEmptyBorder(0, 4, 0, 0));
		mTitleLabel.setOpaque(false);
		add(mTitleLabel, "0,0,0,2");

		add(createToolBar(isClosable, hasPopupButton),"1,1");

		setOpaque(true);
		mActionListener = al;
		mMouseAdapter = new HeaderMouseAdapter(this);
		addMouseListener(mMouseAdapter);
		addMouseMotionListener(mMouseAdapter);
		}

	@Override
	public void updateUI() {
		super.updateUI();
		updateToolbarSeparators();
		}

	private JToolBar createToolBar(boolean isClosable, boolean hasPopupButton) {
		JToolBar toolbar = new JToolBar();
		if (LookAndFeelHelper.isSubstance())
			toolbar.addSeparator();

		if (hasPopupButton) {
			JButton popupButton = new HiDPIIconButton("axisButton.png", null, "popup", 0, null);
			popupButton.addActionListener(this);
			toolbar.add(popupButton);
			if (LookAndFeelHelper.isSubstance())
				toolbar.addSeparator();
		}

		JButton maxButton = new HiDPIIconButton("maxButton.png", "Maximize view", "max", 0, null);
		maxButton.addActionListener(this);
		toolbar.add(maxButton);
		if (LookAndFeelHelper.isSubstance())
			toolbar.addSeparator();

		if (isClosable) {
			JButton closeButton = new HiDPIIconButton("closeButton.png", null, "close", 0, null);
			closeButton.addActionListener(this);
			toolbar.add(closeButton);
			if (LookAndFeelHelper.isSubstance())
				toolbar.addSeparator();
		}

		toolbar.setFloatable(false);
		toolbar.setRollover(true);

		toolbar.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
		toolbar.setOpaque(false);
		return toolbar;
		}

	/**
	 * Depending on the look&feel, this method puts separators between
	 * the toolbar buttons or removes them. For 'substance' we use use
	 * separators, for 'quaqua' we don't.
	 */
	private void updateToolbarSeparators() {
		for (Component c:getComponents()) {
			if (c instanceof JPanel && ((JPanel) c).getComponentCount() == 1) {
				Component t = ((JPanel)c).getComponent(0);
				if (t instanceof JToolBar) {
					JToolBar toolbar = (JToolBar)t;
					ArrayList<JButton> buttonList = new ArrayList<JButton>();
					for (Component b:toolbar.getComponents())
						if (b instanceof JButton)
							buttonList.add((JButton) b);
					toolbar.removeAll();
					if (LookAndFeelHelper.isSubstance())
						toolbar.addSeparator();
					for (JButton b:buttonList) {
						toolbar.add(b);
						if (LookAndFeelHelper.isSubstance())
							toolbar.addSeparator();
						else if (LookAndFeelHelper.isQuaQua())
							// after switch from substance to quaqua the toolbar button style is not re-established
							// (instead if only showing aqua design when rolling over, it uses steady round style)
							// "toolBarTab" seems the best what we can do.
							b.putClientProperty("Quaqua.Button.style", "toolBarTab");
						}
					}
				}
			}
		}

	public String getTitle() {
		return mTitleLabel.getText();
		}

	public void setTitle(String title) {
		mTitleLabel.setText(title);
		}

	public Dockable getDockable() {
		return mDockable;
		}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);

		int width = getWidth();
		int height = getHeight();

		Graphics2D g2 = (Graphics2D) g;
		Paint storedPaint = g2.getPaint();

		g2.setPaint(HeaderPaintHelper.getHeaderPaint(mIsSelected, height));
		g2.fillRect(0, 0, width, height);

		g2.setPaint(storedPaint);
		}

	public void update(boolean isSelected) {
		mIsSelected = isSelected;
		repaint();
		}

	public PopupProvider getPopupProvider() {
		return mMouseAdapter.getPopupProvider();
		}

	public void setPopupProvider(PopupProvider p) {
		mMouseAdapter.setPopupProvider(p);
		}

	public void actionPerformed(ActionEvent e) {
		String title = mTitleLabel.getText();
		if (e.getActionCommand().equals("popup")) {
			mActionListener.actionPerformed(new ActionEvent(this, ActionEvent.ACTION_PERFORMED, "popup_"+title));
			return;
			}

		if (e.getActionCommand().equals("max")) {
			mActionListener.actionPerformed(new ActionEvent(this, ActionEvent.ACTION_PERFORMED, "max_"+title));
			return;
			}

		if (e.getActionCommand().equals("close")) {
			mActionListener.actionPerformed(new ActionEvent(this, ActionEvent.ACTION_PERFORMED, "close_"+title));
			return;
			}
		}
	}
