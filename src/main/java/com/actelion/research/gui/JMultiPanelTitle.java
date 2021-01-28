/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import com.actelion.research.gui.hidpi.HiDPIHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

public class JMultiPanelTitle extends JComponent implements MouseListener,MouseMotionListener {
    private static final long serialVersionUID = 0x20100813;

    private static final int HEIGHT = 10;

	private MultiPanelDragListener	mDragListener;
	private String					mTitle;
	private boolean					mDragEnabled;

	/**
	 * @return default title panel height potentially adapted for HiDPI devices
	 */
	public static int height() {
		return HiDPIHelper.scale(HEIGHT);
		}

	public JMultiPanelTitle(MultiPanelDragListener parent, String title) {
		mDragListener = parent;
		mTitle = title;
		mDragEnabled = true;
		addMouseListener(this);
		addMouseMotionListener(this);
		}

	public void setDragEnabled(boolean b) {
		mDragEnabled = b;
		}

	public void setTitle(String title) {
		mTitle = title;
		repaint();
		}

	public void paintComponent(Graphics g) {
		((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

		Dimension size = getSize();

        Graphics2D g2 = (Graphics2D) g;
        Paint storedPaint = g2.getPaint();

        g2.setPaint(HeaderPaintHelper.getHeaderPaint(true, size.height));
        g2.fillRect(0, 0, size.width, size.height);

        g2.setPaint(storedPaint);

		g.setColor(UIManager.getColor("Label.foreground"));
		g.setFont(UIManager.getFont("Label.font").deriveFont(Font.PLAIN, height()-1));
		int stringWidth = (int)g.getFontMetrics().getStringBounds(mTitle, g).getWidth();
		g.drawString(mTitle, (size.width-stringWidth)/2, height()-2);
		}

	public void mouseClicked(MouseEvent e) {}

	public void mousePressed(MouseEvent e) {
		if (mDragEnabled)
			mDragListener.dragStarted(e.getY(), this);
		}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {
		if (mDragEnabled)
			setCursor(Cursor.getPredefinedCursor(Cursor.S_RESIZE_CURSOR));
		}

	public void mouseExited(MouseEvent e) {
		setCursor(Cursor.getDefaultCursor());
		}

	public void mouseDragged(MouseEvent e) {
		if (mDragEnabled)
			mDragListener.dragContinued(e.getY());
		}

	public void mouseMoved(MouseEvent e) {}
	}
