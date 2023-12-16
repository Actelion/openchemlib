/*
 * Copyright (c) 1997 - 2017
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
import com.actelion.research.util.ColorHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;

/**
 * Lightweight Slider for maintaining two values within given limits
 */

public class JPruningBar extends JPanel implements MouseListener, MouseMotionListener {
	static final long serialVersionUID = 0x20070307;

	private static final int cRedWaterRGB = Color.RED.getRGB();
	private static final int cBlueWaterRGB = 0x156be4;

	private static final int cThumbSize = 2 * HiDPIHelper.scale(7);	// even integer
	private static final int cBarWidth = cThumbSize - 2 * HiDPIHelper.scale(3);  // even integer
	private static final int cBorder = HiDPIHelper.scale(2);

	private float[] mPosition;
	private float	mX1,mY1,mLowValue,mMinValue,mHighValue,mMaxValue, mValuePerPixel;
	private boolean	mIsHorizontal,mUpdateNeeded,mUseRedColor,mAllowDoubleClickChange,mWasDragged;
	private int		mID,mMousePosition,mClickedArea,mActiveThumb;
	private ArrayList<PruningBarListener> mListener;

	public JPruningBar() {
		this(0, 100, true, 0, false);
		}

	public JPruningBar(boolean isHorizontal) {
		this(0, 100, isHorizontal, 0, false);
		}

	public JPruningBar(boolean isHorizontal, int id) {
		this(0, 100, isHorizontal, id, false);
		}

	public JPruningBar(float min, float max, boolean isHorizontal, int id) {
		this(min, max, isHorizontal, id, false);
		}

	public JPruningBar(float min, float max, boolean isHorizontal, int id, boolean allowDoubleClick) {
		init();
		mMinValue = min;
		mLowValue = min;
		mHighValue = max;
		mMaxValue = max;
		mIsHorizontal = isHorizontal;
		mID = id;
		mAllowDoubleClickChange = allowDoubleClick;
		mPosition = new float[2];
		mActiveThumb = -1;
		}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		Graphics2D g2D = (Graphics2D)g;
		Dimension theSize = getSize();
		float lineWidth = 2.5f*HiDPIHelper.getUIScaleFactor();
		g2D.setStroke(new BasicStroke(lineWidth, BasicStroke.CAP_BUTT,BasicStroke.JOIN_MITER));
		Color shadow = LookAndFeelHelper.isDarkLookAndFeel() ?
				ColorHelper.brighter(getBackground(), 0.7f) : ColorHelper.darker(getBackground(), 0.85f);
		Color grayThumb = LookAndFeelHelper.isDarkLookAndFeel() ?
				ColorHelper.brighter(getBackground(), 0.5f) : ColorHelper.darker(getBackground(), 0.75f);

		if (mUpdateNeeded) {
			float x2,y2;  // display area without border
			if (mIsHorizontal) {
				mX1 = cBorder;
				x2 = theSize.width - cBorder;
				mY1 = (theSize.height - cThumbSize) / 2;
				y2 = mY1 + cThumbSize;
				}
			else {
				mX1 = (theSize.width - cThumbSize) / 2;
				x2 = mX1 + cThumbSize;
				mY1 = cBorder;
				y2 = theSize.height - cBorder;
				}

			setColor(g, shadow);
			g2D.fill(new Ellipse2D.Float(mX1, mY1, cThumbSize, cThumbSize));
			g2D.fill(new Ellipse2D.Float(x2-cThumbSize, y2-cThumbSize, cThumbSize, cThumbSize));
			if (mIsHorizontal)
				g2D.fill(new Rectangle2D.Float(mX1+cThumbSize-1, mY1+(cThumbSize-cBarWidth)/2, x2-mX1-2*cThumbSize+2, cBarWidth));
			else
				g2D.fill(new Rectangle2D.Float(mX1+(cThumbSize-cBarWidth)/2, mY1+cThumbSize-1, cBarWidth, y2-mY1-2*cThumbSize+2));

			float zoomSpace = (mIsHorizontal ? x2-mX1 : y2-mY1) - 2 * cThumbSize;
			float valueSpace = mMaxValue - mMinValue;
			mValuePerPixel = valueSpace / zoomSpace;
			if (mIsHorizontal) {
				mPosition[0] = (mLowValue - mMinValue) / mValuePerPixel;
				mPosition[1] = (mHighValue - mMinValue) / mValuePerPixel;
				}
			else {   // tribute to inverted Y-scale in java
				mPosition[0] = (mMaxValue - mHighValue) / mValuePerPixel;
				mPosition[1] = (mMaxValue - mLowValue) / mValuePerPixel;
				}

			int rgb = mUseRedColor ? cRedWaterRGB : HeaderPaintHelper.getThemeColors() == null ? cBlueWaterRGB : HeaderPaintHelper.getThemeColors()[0];
			Color c1 = new Color(ColorHelper.createColor(rgb, 0.7f));
			Color c2 = new Color(ColorHelper.createColor(rgb, 0.5f));
			Color c3 = new Color(ColorHelper.createColor(rgb, 0.3f));

			if (mPosition[1] > mPosition[0]) {
				float waterStart = (mIsHorizontal ? mX1 : mY1) + mPosition[0] + cThumbSize;
				float waterStop = (mIsHorizontal ? mX1 : mY1) + mPosition[1] + cThumbSize;

				Paint storedPaint = g2D.getPaint();
				if (mIsHorizontal) {
					float yy1 = mY1 + (cThumbSize - cBarWidth) / 2;
					float yy2 = mY1 + (cThumbSize + cBarWidth) / 2;
					g2D.setPaint(new GradientPaint(0, yy1, c1, 0, yy2, c3));
					g2D.fill(new Rectangle2D.Float(waterStart, yy1, waterStop - waterStart, cBarWidth));
					}
				else {
					float xx1 = mX1 + (cThumbSize - cBarWidth) / 2;
					float xx2 = mY1 + (cThumbSize + cBarWidth) / 2;
					g2D.setPaint(new GradientPaint(xx1, 0, c1, xx2, 0, c3));
					g2D.fill(new Rectangle2D.Float(xx1, waterStart, cBarWidth, waterStop - waterStart));
					}

				g2D.setPaint(storedPaint);
				}

			for (int i=0; i<2; i++)
				drawThumb(g2D, i, i == mActiveThumb ? c2 : grayThumb, c3);
			}
		}

	private void drawThumb(Graphics2D g, int i, Color color1, Color color2) {
		Ellipse2D.Float elipse = getThumbElipse(i);
		setColor(g, color1);
		g.fill(elipse);
		setColor(g, color2);
		g.draw(elipse);
		}

	private Ellipse2D.Float getThumbElipse(int i) {
		float shift = mPosition[i] + (i == 1 ? cThumbSize : 0);
		float x = mX1 + (mIsHorizontal ? shift : 0);
		float y = mY1 + (mIsHorizontal ? 0 : shift);
		return new Ellipse2D.Float(x, y, cThumbSize, cThumbSize);
		}

	private void setColor(Graphics g, Color color) {
		g.setColor(isEnabled() ? color : ColorHelper.intermediateColor(color, Color.LIGHT_GRAY, 0.7f));
		}

	public void update(Graphics g) {
		paint(g);
		}

	public float getLowValue() {
		return mLowValue;
		}

	public float getHighValue() {
		return mHighValue;
		}

	public float getMaximumValue() {
		return mMaxValue;
		}

	public float getMinimumValue() {
		return mMinValue;
		}

	/**
	 * Sets high and low values to max and min, respectively.
	 * Sends PruningBarEvents in case of a change.
	 */
	public void reset() {
		boolean highChanged = setHigh(mMaxValue);
		boolean lowChanged = setLow(mMinValue);

		if (highChanged || lowChanged) {
			informListeners(false);
			mUpdateNeeded = true;
			repaint();
			}
		}

	/**
	 * Updates the low and high values simultaniously, provided the low <= high,
	 * low >= current min, and high <= current max.
	 * @param low
	 * @param high
	 * @param silent if true, PruningBarEvents are suppressed
	 */
	public void setLowAndHigh(float low, float high, boolean silent) {
		if (low < mMinValue)
			low = mMinValue;
		if (high > mMaxValue)
			high = mMaxValue;
		if ((low != mLowValue || high != mHighValue)
		 && low <= high) {
			mLowValue = low;
			mHighValue = high;
			if (!silent)
				informListeners(false);
			mUpdateNeeded = true;
			repaint();
			}
		}

	/**
	 * Updates the low value provided the new value is >= current min and <= current high.
	 * Sends PruningBarEvents in case of a successful change.
	 * @param low
	 */
	public void setLowValue(float low) {
		if (setLow(low)) {
			informListeners(false);
			mUpdateNeeded = true;
			repaint();
			}
		}

	public void setID(int id) {
		mID = id;
		}

	/**
	 * Updates the high value provided the new value is >= current low and <= current max.
	 * Sends PruningBarEvents in case of a successful change.
	 * @param high
	 */
	public void setHighValue(float high) {
		if (setHigh(high)) {
			informListeners(false);
			mUpdateNeeded = true;
			repaint();
			}
		}

	/**
	 * Updates the maximum value; may update low and high to stay within limits
	 * Sends PruningBarEvents in case of a successful change.
	 * @param max
	 */
	public void setMaximumValue(float max) {
		if (max < mMinValue)
			max = mMinValue;

		if (mMaxValue != max) {
			mMaxValue = max;
			if (mHighValue > mMaxValue)
				mHighValue = mMaxValue;
			if (mLowValue > mHighValue)
				mLowValue = mHighValue;

			informListeners(false);
			mUpdateNeeded = true;
			repaint();
			}
		}

	/**
	 * Updates the minimum value; may update low and high to stay within limits
	 * Sends PruningBarEvents in case of a successful change.
	 * @param min
	 */
	public void setMinimumValue(float min) {
			// changes the allowed min value; may update low and high to stay within limits
		if (min > mMaxValue)
			min = mMaxValue;

		if (mMinValue != min) {
			mMinValue = min;
			if (mLowValue < mMinValue)
				mLowValue = mMinValue;
			if (mHighValue < mLowValue)
				mHighValue = mLowValue;

			informListeners(false);
			mUpdateNeeded = true;
			repaint();
			}
		}

	/**
	 * Initializes the bar by setting min and low/max and high to the values given.
	 * Does not(!) send PruningBarEvents.
	 * @param min
	 * @param max
	 */
	public void setMinAndMax(float min, float max) {
		mLowValue = mMinValue = min;
		mHighValue = mMaxValue = max;
		mUpdateNeeded = true;
		repaint();
		}

	public void setUseRedColor(boolean useRed) {
		if (mUseRedColor != useRed) {
			mUseRedColor = useRed;
			mUpdateNeeded = true;
			repaint();
			}
		}
	
	public Dimension getMinimumSize() {
		if (mIsHorizontal)
			return new Dimension(4* cThumbSize, cThumbSize+2*cBorder);
		else
			return new Dimension(cThumbSize+2*cBorder, 4* cThumbSize);
		}

	public Dimension getPreferredSize() {
		if (mIsHorizontal)
			return new Dimension(100, cThumbSize+2*cBorder);
		else
			return new Dimension(cThumbSize+2*cBorder, 100);
		}

	public Dimension getMaximumSize() {
		if (mIsHorizontal)
			return new Dimension(Short.MAX_VALUE, cThumbSize+2*cBorder);
		else
			return new Dimension(cThumbSize+2*cBorder, Short.MAX_VALUE);
		}

	public void mousePressed(MouseEvent e) {
		if (!isEnabled())
			return;

		if (mIsHorizontal)
			mMousePosition = e.getX();
		else
			mMousePosition = e.getY();

		mWasDragged = false;

		mClickedArea = 0;
		if (mMousePosition < mPosition[0])
			mClickedArea = 0;
		else if (mMousePosition <= mPosition[0] + cThumbSize)
			mClickedArea = 1;
		else if (mMousePosition <= mPosition[1]+ cThumbSize)
			mClickedArea = 2;
		else if (mMousePosition <= mPosition[1] + 2 * cThumbSize)
			mClickedArea = 3;
		else
			mClickedArea = 0;

		if (!mIsHorizontal  // tribute to inverted Y-scale in java
		 && (mClickedArea == 1 || mClickedArea == 3))
			mClickedArea = 4 - mClickedArea;
		}

	public void mouseReleased(MouseEvent e) {
		if (mWasDragged)
			informListeners(false);
		}

	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}

	/**
	 * If this pruning bar allows changes via double click and text input, then such an action will
	 * cause a PruningBarEvent with type=TYPE_TYPED. Its low or high value will contain the number that
	 * the user has typed in. The other value will be Float.NaN. The JPruningBar will <b>not</b> update
	 * internal values, because JPruningBar scale and user visible scale may be different. The receiver
	 * of the event is responsible to translate the user number into JPuningBar scale and call the
	 * pruning bar's setLowValue() or setHighValue().
	 * @param e
	 */
	public void mouseClicked(MouseEvent e) {
		if (mAllowDoubleClickChange && e.getClickCount() == 2 && (mClickedArea == 1 || mClickedArea == 3)) {
			Component c = this;
			while (c != null && !(c instanceof Window || c instanceof Frame))
				c = c.getParent();
			if (c != null) {
				try {
					String s = JOptionPane.showInputDialog(c, "Please type in a value!", "Set Value", JOptionPane.QUESTION_MESSAGE);
					if (s != null) {
						float d = Float.parseFloat(s);
						if (mClickedArea == 1) {
							informListeners(d, Float.NaN);
							}
						else if (mClickedArea == 3) {
							informListeners(Float.NaN, d);
							}
						}
					}
				catch (NumberFormatException nfe) {}
				}
			}
		}

	public synchronized void mouseDragged(MouseEvent e) {
		if (!isEnabled())
			return;

		int position;
		if (mIsHorizontal)
			position = e.getX();
		else
			position = e.getY();

		if (position == mMousePosition)
			return;

		mWasDragged = true;

		float change = mValuePerPixel * (float)(position - mMousePosition);
		if (!mIsHorizontal) // tribute to inverted Y-scale in java
			change = -change;
		if (e.isControlDown())
			change /= 10.0;

		boolean valuesChanged = false;
		switch (mClickedArea) {
		case 1:
			if ((change < 0.0 && mLowValue > mMinValue)
			 || (change > 0.0 && mLowValue < mHighValue))
				valuesChanged = true;

			if (mLowValue + change < mMinValue)
				mLowValue = mMinValue;
			else if (mLowValue + change > mHighValue)
				mLowValue = mHighValue;
			else
				mLowValue += change;
			break;
		case 2:
			if ((change < 0.0 && mLowValue > mMinValue)
			 || (change > 0.0 && mHighValue < mMaxValue))
				valuesChanged = true;

			if (mLowValue + change < mMinValue) {
				mHighValue -= mLowValue - mMinValue;
				mLowValue = mMinValue;
				}
			else if (mHighValue + change > mMaxValue) {
				mLowValue += mMaxValue - mHighValue;
				mHighValue = mMaxValue;
				}
			else {
				mLowValue += change;
				mHighValue += change;
				}
			break;
		case 3:
			if ((change < 0.0 && mHighValue > mLowValue)
			 || (change > 0.0 && mHighValue < mMaxValue))
				valuesChanged = true;

			if (mHighValue + change > mMaxValue)
				mHighValue = mMaxValue;
			else if (mHighValue + change < mLowValue)
				mHighValue = mLowValue;
			else
				mHighValue += change;
			break;
			}

		if (valuesChanged) {
			informListeners(true);
			mMousePosition = position;
			mUpdateNeeded = true;
			repaint();
			}
		}

	public void mouseMoved(MouseEvent e) {
		int old = mActiveThumb;
		mActiveThumb = -1;
		for (int i=0; i<2; i++) {
			if (getThumbElipse(i).contains(e.getX(), e.getY())) {
				mActiveThumb = i;
				break;
				}
			}
		if (old != mActiveThumb)
			repaint();
		}

	public void addPruningBarListener(PruningBarListener listener) {
		mListener.add(listener);
		}

	public void removePruningBarListener(PruningBarListener listener) {
		mListener.remove(listener);
		}

	public void firePruningBarChanged() {
		informListeners(false);
		}

	private void init() {
		this.setOpaque(false);
		mUpdateNeeded = true;
		addMouseListener(this);
		addMouseMotionListener(this);
		mListener = new ArrayList<>();
		}

	private boolean setHigh(float value) {
		if (value < mLowValue)
			value = mLowValue;
		else if (value > mMaxValue)
			value = mMaxValue;

		if (value == mHighValue)
			return false;

		mHighValue = value;
		return true;
		}

	private boolean setLow(float value) {
		if (value < mMinValue)
			value = mMinValue;
		else if (value > mHighValue)
			value = mHighValue;

		if (value == mLowValue)
			return false;

		mLowValue = value;
		return true;
		}

	private void informListeners(float low, float high) {
		for (int i=0; i<mListener.size(); i++)
			mListener.get(i).pruningBarChanged(
				new PruningBarEvent(this, low, high, false, mID, PruningBarEvent.TYPE_TYPED));
		}

	private void informListeners(boolean isAdjusting) {
		for (int i=0; i<mListener.size(); i++)
			mListener.get(i).pruningBarChanged(
				new PruningBarEvent(this, mLowValue, mHighValue, isAdjusting, mID, PruningBarEvent.TYPE_DRAGGED));
		}
	}
