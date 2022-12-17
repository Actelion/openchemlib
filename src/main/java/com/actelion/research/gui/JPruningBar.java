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

import com.actelion.research.util.ColorHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;

/**
 * Lightweight Slider for maintaining two values within given limits
 */

public class JPruningBar extends JPanel implements MouseListener, MouseMotionListener {
	static final long serialVersionUID = 0x20070307;

	private static final Color cShadowColor		 = new Color(153, 153, 153);
	private static final Color cDarkShadowColor	 = new Color(102, 102, 102);
	private static final Color cThumbColor		  = new Color(153, 153, 204);
	private static final Color cThumbShadowColor	= new Color(102, 102, 153);
	private static final Color cThumbHighlightColor = new Color(204, 204, 255);
	private static final Color[] cWaterBlueColor = {
							new Color(188, 216, 248), new Color(184, 215, 247),
							new Color(170, 205, 245), new Color(150, 194, 242),
							new Color(116, 172, 235), new Color(138, 190, 243),
							new Color(146, 196, 247), new Color(155, 203, 252),
							new Color(164, 208, 255), new Color(170, 215, 255),
							new Color(177, 223, 255), new Color(193, 235, 255),
							new Color(207, 245, 255), new Color(207, 255, 255) };
	private static final Color[] cWaterRedColor = {
							new Color(255, 189, 219), new Color(247, 189, 215),
							new Color(247, 173, 207), new Color(247, 148, 195),
							new Color(239, 115, 174), new Color(247, 140, 190),
							new Color(247, 148, 199), new Color(255, 156, 203),
							new Color(255, 165, 211), new Color(255, 173, 215),
							new Color(255, 181, 223), new Color(255, 198, 235),
							new Color(255, 206, 247), new Color(255, 206, 255) };
	private static final int cOverallWidth = 16;
	private static final int cThumbHeight = 14;
	private static final int cBorder = 2;

	private float	mLowValue,mMinValue,mHighValue,mMaxValue,mSegmentSize;
	private boolean	mIsHorizontal,mUpdateNeeded,mUseRedColor,mAllowDoubleClickChange,mWasDragged;
	private int		mID,mMousePosition,mClickedArea,mPosition1,mPosition2;
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
		}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);

		Dimension theSize = getSize();

		if (mUpdateNeeded) {
			int x1,y1,x2,y2;
			if (mIsHorizontal) {
				x1 = cBorder;
				x2 = theSize.width - cBorder - 1;
				y1 = (theSize.height - cOverallWidth) / 2;
				y2 = y1 + cOverallWidth - 1;
				}
			else {
				x1 = (theSize.width - cOverallWidth) / 2;
				x2 = x1 + cOverallWidth - 1;
				y1 = cBorder;
				y2 = theSize.height - cBorder - 1;
				}

			setColor(g, cShadowColor);
			g.drawLine(x1+1, y1+1, x1+1, y2-1);
			g.drawLine(x1+1, y1+1, x2-1, y1+1);

			int freePixel = (mIsHorizontal) ? x2-x1 : y2-y1;
			freePixel -= 2 * cThumbHeight + 2;
			float allSize = mMaxValue - mMinValue;
			mSegmentSize = allSize / (float)freePixel;
			if (mIsHorizontal) {
				mPosition1 = (int)((mLowValue - mMinValue + mSegmentSize/2) / mSegmentSize) + 1;
				mPosition2 = (int)((mHighValue - mMinValue + mSegmentSize/2) / mSegmentSize) + 1 + cThumbHeight + 1;
				}
			else {   // tribute to inverted Y-scale in java
				mPosition1 = (int)((mMaxValue - mHighValue + mSegmentSize/2) / mSegmentSize) + 1;
				mPosition2 = (int)((mMaxValue - mLowValue + mSegmentSize/2) / mSegmentSize) + 1 + cThumbHeight + 1;
				}

			drawThumb(g, x1, y1, mPosition1);
			drawThumb(g, x1, y1, mPosition2);

			if (mPosition2 > mPosition1+cThumbHeight+2) {
				int waterStart = mPosition1 + cThumbHeight + 1;
				int waterStop = mPosition2 - 2;
				for (int i=0; i<14; i++) {
					setColor(g, mUseRedColor ? cWaterRedColor[i] : cWaterBlueColor[i]);
					if (mIsHorizontal)
						g.drawLine(x1+waterStart, y1+i+1, x1+waterStop, y1+i+1);
					else
						g.drawLine(x1+i+1, y1+waterStart, x1+i+1, y1+waterStop);
					}
				}

			setColor(g, cDarkShadowColor);
			g.drawRect(x1, y1, x2-x1, y2-y1);
			}
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
			return new Dimension(4*cThumbHeight, cOverallWidth+2*cBorder);
		else
			return new Dimension(cOverallWidth+2*cBorder, 4*cThumbHeight);
		}

	public Dimension getPreferredSize() {
		if (mIsHorizontal)
			return new Dimension(100, cOverallWidth+2*cBorder);
		else
			return new Dimension(cOverallWidth+2*cBorder, 100);
		}

	public Dimension getMaximumSize() {
		if (mIsHorizontal)
			return new Dimension(Short.MAX_VALUE, cOverallWidth+2*cBorder);
		else
			return new Dimension(cOverallWidth+2*cBorder, Short.MAX_VALUE);
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
		if (mMousePosition < mPosition1)
			mClickedArea = 0;
		else if (mMousePosition <= mPosition1 + cThumbHeight)
			mClickedArea = 1;
		else if (mMousePosition <= mPosition2)
			mClickedArea = 2;
		else if (mMousePosition <= mPosition2 + cThumbHeight)
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

		float change = mSegmentSize * (float)(position - mMousePosition);
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

	public void mouseMoved(MouseEvent e) {}

	public void addPruningBarListener(PruningBarListener listener) {
		mListener.add(listener);
		}

	public void removePruningBarListener(PruningBarListener listener) {
		mListener.remove(listener);
		}

	public void firePruningBarChanged() {
		informListeners(false);
		}

	private void drawThumb(Graphics g, int x1, int y1, int position) {
		int x,y;

		setColor(g, cThumbColor);
		if (mIsHorizontal) {
			g.fillRect(x1+position+1, y1+2, cThumbHeight-1, cOverallWidth-3);
			x = x1+position+2;
			y = y1+3;
			}
		else {
			g.fillRect(x1+2, y1+position+1, cOverallWidth-3, cThumbHeight-1);
			x = x1+3;
			y = y1+position+2;
			}

		setColor(g, cThumbHighlightColor);
		if (mIsHorizontal) {
			g.drawLine(x1+position, y1+1, x1+position+cThumbHeight-2, y1+1);
			g.drawLine(x1+position, y1+1, x1+position, y1+cOverallWidth-2);
			}
		else {
			g.drawLine(x1+1, y1+position, x1+1, y1+position+cThumbHeight-2);
			g.drawLine(x1+1, y1+position, x1+cOverallWidth-2, y1+position);
			}
		for (int i=0; i<5; i++)
			for (int j=0; j<5; j++)
				if (((i + j) & 1) == 0)
					g.drawLine(x+i*2, y+j*2, x+i*2, y+j*2);

		setColor(g, cThumbShadowColor);
		if (mIsHorizontal) {
			g.drawLine(x1+position-1, y1+1, x1+position-1, y1+cOverallWidth-2);
			g.drawLine(x1+position+cThumbHeight, y1+1, x1+position+cThumbHeight, y1+cOverallWidth-2);
			}
		else {
			g.drawLine(x1+1, y1+position-1, x1+cOverallWidth-2, y1+position-1);
			g.drawLine(x1+1, y1+position+cThumbHeight, x1+cOverallWidth-2, y1+position+cThumbHeight);
			}
		for (int i=0; i<5; i++)
			for (int j=0; j<5; j++)
				if (((i + j) & 1) == 0)
					g.drawLine(x+i*2+1, y+j*2+1, x+i*2+1, y+j*2+1);
		}

	private void setColor(Graphics g, Color color) {
		g.setColor(isEnabled() ? color : ColorHelper.intermediateColor(color, Color.LIGHT_GRAY, 0.7f));
		}

	private void init() {
		this.setOpaque(false);
		mUpdateNeeded = true;
		addMouseListener(this);
		addMouseMotionListener(this);
		mListener = new ArrayList<PruningBarListener>();
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
