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

import java.awt.*;
import java.util.ArrayList;

public class Depictor extends AbstractDepictor {
	private static final int MAX_TEXTSIZE = 16;

	private int			mpTextSize,mMaxTextSize;
	private float		mLineWidth;
	private ArrayList<Font>	mFonts;
    private Font currentFont;

	public Depictor(StereoMolecule mol) {
		super(mol);
		}


	public Depictor(StereoMolecule mol, int displayMode) {
		super(mol, displayMode);
		}

	
	public void setMaximumTextSize(int maxTextSize) {
		mMaxTextSize = maxTextSize;
		}


	protected void init() {
		super.init();
		mFonts = new ArrayList<Font>();
		mMaxTextSize = MAX_TEXTSIZE;
		mLineWidth = 1.0f;
		}


	protected void drawBlackLine(DepictorLine theLine) {
		((Graphics)mG).drawLine((int)Math.round(theLine.x1),(int)Math.round(theLine.y1),
								(int)Math.round(theLine.x2),(int)Math.round(theLine.y2));
		}


    protected void drawDottedLine(DepictorLine theLine) {
        Stroke stroke = ((Graphics2D)mG).getStroke();
        ((Graphics2D)mG).setStroke(new BasicStroke(mLineWidth, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,
		                        mLineWidth, new float[] {3.0f*mLineWidth}, 0f));
        ((Graphics)mG).drawLine((int)theLine.x1,(int)theLine.y1,
                                (int)theLine.x2,(int)theLine.y2);
        ((Graphics2D)mG).setStroke(stroke);
        }


    public void drawString(String theString,double x, double y) {
	    double strWidth = getStringWidth(theString);
		((Graphics)mG).drawString(theString,(int)Math.round(x-strWidth/2),(int)Math.round(y+1+mpTextSize/3));
		}


	protected void drawPolygon(double[] x, double[] y, int count) {
		int[] px = new int[count];
		int[] py = new int[count];
		for (int i=0; i<count; i++) {
			px[i] = (int)Math.round(x[i]);
			py[i] = (int)Math.round(y[i]);
			}
		((Graphics)mG).drawPolygon(px, py, count);
		((Graphics)mG).fillPolygon(px, py, count);
		}


	protected void fillCircle(double x, double y, double r) {
	    ((Graphics)mG).fillOval((int)Math.round(x), (int)Math.round(y), (int)Math.round(r), (int)Math.round(r));
		}


	protected double getStringWidth(String theString) {
		return ((Graphics)mG).getFontMetrics().stringWidth(theString);
		}


	public void setTextSize(int theSize)
    {
        mpTextSize = Math.min(theSize, mMaxTextSize);
        if (mG != null) {
            if (currentFont == null || currentFont.getSize() != mpTextSize) {
                for (int i = 0; i < mFonts.size(); i++) {
                    if ((mFonts.get(i)).getSize() == mpTextSize) {
                        ((Graphics) mG).setFont(mFonts.get(i));
                        return;
                    }
                }
				Font newFont = ((Graphics)mG).getFont().deriveFont(0, mpTextSize);
//                Font newFont = new Font("Helvetica", 0, mpTextSize);
                mFonts.add(newFont);
                currentFont = newFont;
                ((Graphics) mG).setFont(newFont);
            }
        }
    }

	public int getTextSize() {
	    return mpTextSize;
	    }


	protected double getLineWidth() {
		return mLineWidth;
		}


	protected void setLineWidth(double lineWidth) {
		if (lineWidth <= 1.5f)
			lineWidth = 1.0f;
		if (mLineWidth != lineWidth) {
			mLineWidth = (float)lineWidth;
			((Graphics2D)mG).setStroke(new BasicStroke((float)lineWidth, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
			}
		}


	public void setColor(Color theColor) {
	    ((Graphics)mG).setColor(theColor);
		}
	}
