/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * 3. Neither the name of the copyright holder nor the
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

import com.actelion.research.gui.hidpi.HiDPIHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

public class JMultiPanelTitle extends JComponent implements MouseListener,MouseMotionListener {
    private static final long serialVersionUID = 0x20100813;

    public static final int HEIGHT = HiDPIHelper.scale(10);

	private MultiPanelDragListener	mDragListener;
	private String					mTitle;
	private boolean					mDragEnabled;

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
		g.setFont(UIManager.getFont("Label.font").deriveFont(Font.PLAIN, HEIGHT-1));
		int stringWidth = (int)g.getFontMetrics().getStringBounds(mTitle, g).getWidth();
		g.drawString(mTitle, (size.width-stringWidth)/2, HEIGHT-2);
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
