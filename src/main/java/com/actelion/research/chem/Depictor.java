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
import java.util.ArrayList;

@Deprecated
public class Depictor extends AbstractDepictor<Graphics> {
	private static final int MAX_TEXTSIZE = 16;

	private int			mpTextSize,mMaxTextSize;
	private float		mLineWidth;
	private ArrayList<Font>	mFonts;
    private Font        currentFont;

	public Depictor(StereoMolecule mol) {
		super(mol);
		}

	public Depictor(StereoMolecule mol, int displayMode) {
		super(mol, displayMode);
		}

	public void setMaximumTextSize(int maxTextSize) {
		mMaxTextSize = maxTextSize;
		}

	@Override
	protected void init() {
		super.init();
		mFonts = new ArrayList<>();
		mMaxTextSize = MAX_TEXTSIZE;
		mLineWidth = 1.0f;
		}


	@Override
	protected void drawBlackLine(DepictorLine theLine) {
		mContext.drawLine((int)Math.round(theLine.x1),(int)Math.round(theLine.y1),
						  (int)Math.round(theLine.x2),(int)Math.round(theLine.y2));
		}


	@Override
    protected void drawDottedLine(DepictorLine theLine) {
        Stroke stroke = ((Graphics2D)mContext).getStroke();
        ((Graphics2D)mContext).setStroke(new BasicStroke(mLineWidth, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,
		                        mLineWidth, new float[] {3.0f*mLineWidth}, 0f));
        mContext.drawLine((int)theLine.x1,(int)theLine.y1, (int)theLine.x2,(int)theLine.y2);
        ((Graphics2D)mContext).setStroke(stroke);
        }


	@Override
    public void drawString(String theString,double x, double y) {
	    double strWidth = getStringWidth(theString);
		mContext.drawString(theString,(int)Math.round(x-strWidth/2),(int)Math.round(y+1+mpTextSize/3));
		}


	@Override
	protected void drawPolygon(GenericPolygon p) {
		int[] px = new int[p.getSize()];
		int[] py = new int[p.getSize()];
		for (int i=0; i<p.getSize(); i++) {
			px[i] = (int)Math.round(p.getX(i));
			py[i] = (int)Math.round(p.getY(i));
			}
		mContext.drawPolygon(px, py, p.getSize());
		mContext.fillPolygon(px, py, p.getSize());
		}


	@Override
	protected void fillCircle(double x, double y, double d) {
	    mContext.fillOval((int)Math.round(x), (int)Math.round(y), (int)Math.round(d), (int)Math.round(d));
		}


/*	protected void drawBlackArc(MyRect theArc,double arcAngle,boolean inverted) {
		double xdif = (theArc.x2 - theArc.x1);
		double ydif = (theArc.y2 - theArc.y1);
		double length = Math.sqrt(xdif*xdif + ydif*ydif);
		double theAngle = (inverted) ? (Math.PI - cArcAngle) / 2
									 : (cArcAngle - Math.PI) / 2;
		double radius = (length / 2) / Math.cos(theAngle);
		double centerX = theArc.x1 + radius * Math.sin(arcAngle + theAngle);
		double centerY = theArc.y1 + radius * Math.cos(arcAngle + theAngle);

		mContext.drawArc((int)(centerX - radius + 0.5),(int)(centerY - radius + 0.5),
				   (int)(2 * radius + 0.5),(int)(2 * radius + 0.5),
				   (int)((arcAngle - theAngle) * 180 / Math.PI) - 90,
				   (inverted) ? (int)(-cArcAngle * 180 / Math.PI + 0.5)
				              : (int)(cArcAngle * 180 / Math.PI + 0.5));
		}	*/


	@Override
	protected double getStringWidth(String theString) {
		return mContext.getFontMetrics().stringWidth(theString);
		}

	@Override
	public void setTextSize(int theSize) {
        mpTextSize = Math.min(theSize, mMaxTextSize);
        if (mContext != null) {
            if (currentFont == null || currentFont.getSize() != mpTextSize) {
                for (int i = 0; i < mFonts.size(); i++) {
                    if ((mFonts.get(i)).getSize() == mpTextSize) {
                        ((Graphics) mContext).setFont(mFonts.get(i));
                        return;
                    }
                }
				Font newFont = mContext.getFont().deriveFont(0, mpTextSize);
//                Font newFont = new Font("Helvetica", 0, mpTextSize);
                mFonts.add(newFont);
                currentFont = newFont;
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
		if (lineWidth <= 1.5f)
			lineWidth = 1.0f;
		if (mLineWidth != lineWidth) {
			mLineWidth = (float)lineWidth;
			((Graphics2D)mContext).setStroke(new BasicStroke((float)lineWidth, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
			}
		}

	@Override
	public void setRGB(int rgb) {
	    mContext.setColor(new Color(rgb));
		}
	}
