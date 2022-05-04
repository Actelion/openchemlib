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
*/

package com.actelion.research.chem;


import com.actelion.research.gui.generic.GenericPolygon;

import java.awt.*;
import java.awt.font.GlyphVector;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.ArrayList;

public class Depictor2D extends AbstractDepictor<Graphics2D> {
    private static boolean isMac = (System.getProperty("os.name").toLowerCase().indexOf("mac") >= 0);

    private int			mpTextSize;
	private double      mCurrentStringWidth;
	private float		mLineWidth;
	private ArrayList<Font>	mFonts;
	private String      mCurrentString;
	private Font        mCurrentFont;
	private GlyphVector mCurrentGlyphVector;

	public Depictor2D(StereoMolecule mol) {
		super(mol);
		}

	public Depictor2D(StereoMolecule mol, int displayMode) {
		super(mol, displayMode);
		}

	@Override
	protected void init() {
		super.init();
		mFonts = new ArrayList<Font>();
		mLineWidth = 1.0f;
		}

	@Override
	protected void drawBlackLine(DepictorLine theLine) {
		// Lines on OSX are shifted left and down when drawn in Graphics2D by 0.5. We must compensate.
		if (isMac)
			mContext.draw(new Line2D.Double(theLine.x1-0.5, theLine.y1-0.5, theLine.x2-0.5, theLine.y2-0.5));
		else
			mContext.draw(new Line2D.Double(theLine.x1, theLine.y1, theLine.x2, theLine.y2));
		}

	@Override
    protected void drawDottedLine(DepictorLine theLine) {
        Stroke stroke = mContext.getStroke();
        mContext.setStroke(new BasicStroke(
        						mLineWidth,
                                BasicStroke.CAP_ROUND,
                                BasicStroke.JOIN_ROUND,
                                mLineWidth,
                                new float[] {3.0f*mLineWidth},
                                0f));

        drawBlackLine(theLine);

        mContext.setStroke(stroke);
        }

	@Override
	protected void drawString(String theString, double x, double y) {
		double strWidth = getStringWidth(theString);
		mContext.drawGlyphVector(mCurrentGlyphVector, (float)(x-strWidth/2.0),
										(float)(y+(float)mpTextSize/3.0));
		}

	@Override
	protected void drawPolygon(GenericPolygon p) {
		GeneralPath polygon = new GeneralPath(GeneralPath.WIND_NON_ZERO, p.getSize());
		polygon.moveTo(p.getX(0), p.getY(0));
		for (int i=1; i<p.getSize(); i++)
			polygon.lineTo(p.getX(i), p.getY(i));
		polygon.closePath();

		mContext.fill(polygon);

		if (isMac) {
			polygon = new GeneralPath(GeneralPath.WIND_NON_ZERO, p.getSize());
			polygon.moveTo(p.getX(0)-0.5f, p.getY(0)-0.5f);
			for (int i=1; i<p.getSize(); i++)
				polygon.lineTo(p.getX(i)-0.5f, p.getY(i)-0.5f);
			polygon.closePath();
			}
		
		mContext.draw(polygon);
		}

	@Override
	protected void fillCircle(double x, double y, double d) {
		if (isMac)
			mContext.fill(new Ellipse2D.Double(x-0.5f, y-0.5f, d, d));
		else
			mContext.fill(new Ellipse2D.Double(x, y, d, d));
		}

	@Override
	protected double getStringWidth(String theString) {
		if (mContext != null) {
			if (!theString.equals(mCurrentString)
				|| mCurrentFont != ((Graphics2D) mContext).getFont()) {
				mCurrentString = theString;
				mCurrentFont = ((Graphics2D) mContext).getFont();
				mCurrentGlyphVector = ((Graphics2D) mContext).getFont().createGlyphVector(((Graphics2D) mContext).getFontRenderContext(), theString);
				mCurrentStringWidth = mCurrentGlyphVector.getLogicalBounds().getWidth();
			}
		}
		return mCurrentStringWidth;
		}

	@Override
	protected void setTextSize(int theSize) {
		mpTextSize = (int)theSize;
		if (mContext != null) {
			if (mContext.getFont().getSize() != theSize) {
				for (int i = 0; i < mFonts.size(); i++) {
					if ((mFonts.get(i)).getSize() == theSize) {
						mContext.setFont(mFonts.get(i));
						return;
						}
					}
				Font newFont = mContext.getFont().deriveFont(0, theSize);
//				Font newFont = new Font("Helvetica", 0, (int) theSize);
				mFonts.add(newFont);
				mContext.setFont(newFont);
				}
			}
		}

	@Override
    public int getTextSize() {
        return mpTextSize;
        }

	@Override
	protected double getLineWidth() {
		return mLineWidth;
		}

	@Override
	protected void setLineWidth(double lineWidth) {
		mLineWidth = (float)lineWidth;
		mContext.setStroke(new BasicStroke((float)lineWidth, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
		}

	@Override
	protected void setRGB(int rgb) {
		Color color = new Color(rgb);
	    mContext.setColor(color);
		mContext.setPaint(color);
		}
	}

