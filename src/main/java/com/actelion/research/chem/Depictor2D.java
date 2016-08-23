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
import java.awt.font.GlyphVector;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.ArrayList;

public class Depictor2D extends AbstractDepictor {
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

	
	protected void init() {
		super.init();
		mFonts = new ArrayList<Font>();
		mLineWidth = 1.0f;
		}


	protected void drawBlackLine(DepictorLine theLine) {
		// Lines on OSX are shifted left and down when drawn in Graphics2D by 0.5. We must compensate.
		if (isMacintosh())
			((Graphics2D)mG).draw(new Line2D.Double(theLine.x1-0.5, theLine.y1-0.5, theLine.x2-0.5, theLine.y2-0.5));
		else
			((Graphics2D)mG).draw(new Line2D.Double(theLine.x1, theLine.y1, theLine.x2, theLine.y2));
		}


    protected void drawDottedLine(DepictorLine theLine) {
        Stroke stroke = ((Graphics2D)mG).getStroke();
        ((Graphics2D)mG).setStroke(new BasicStroke(
        						mLineWidth,
                                BasicStroke.CAP_ROUND,
                                BasicStroke.JOIN_ROUND,
                                mLineWidth,
                                new float[] {3.0f*mLineWidth},
                                0f));

        drawBlackLine(theLine);

        ((Graphics2D)mG).setStroke(stroke);
        }


	protected void drawString(String theString, double x, double y) {
		double strWidth = getStringWidth(theString);
		((Graphics2D)mG).drawGlyphVector(mCurrentGlyphVector, (float)(x-strWidth/2.0),
										(float)(y+(float)mpTextSize/3.0));
		}


	protected void drawPolygon(double[] x, double[] y, int count) {
		GeneralPath polygon = new GeneralPath(GeneralPath.WIND_NON_ZERO, count);
		polygon.moveTo((float)x[0], (float)y[0]);
		for (int i=1; i<count; i++)
			polygon.lineTo((float)x[i], (float)y[i]);
		polygon.closePath();

		((Graphics2D)mG).fill(polygon);
		
		if (isMacintosh()) {
			polygon = new GeneralPath(GeneralPath.WIND_NON_ZERO, count);
			polygon.moveTo((float)x[0]-0.5f, (float)y[0]-0.5f);
			for (int i=1; i<count; i++)
				polygon.lineTo((float)x[i]-0.5f, (float)y[i]-0.5f);
			polygon.closePath();
			}
		
		((Graphics2D)mG).draw(polygon);
		}


	protected void fillCircle(double x, double y, double r) {
		if (isMacintosh())
			((Graphics2D)mG).fill(new Ellipse2D.Double(x-0.5f, y-0.5f, r, r));
		else
			((Graphics2D)mG).fill(new Ellipse2D.Double(x, y, r, r));
		}


	protected double getStringWidth(String theString) {
		if (mG != null) {
			if (!theString.equals(mCurrentString)
				|| mCurrentFont != ((Graphics2D) mG).getFont()) {
				mCurrentString = theString;
				mCurrentFont = ((Graphics2D) mG).getFont();
				mCurrentGlyphVector = ((Graphics2D) mG).getFont().createGlyphVector(((Graphics2D) mG).getFontRenderContext(), theString);
				mCurrentStringWidth = mCurrentGlyphVector.getLogicalBounds().getWidth();
			}
		}
		return mCurrentStringWidth;
		}


	protected void setTextSize(int theSize) {
		mpTextSize = (int)theSize;
		if (mG != null) {
			if (((Graphics2D) mG).getFont().getSize() != theSize) {
				for (int i = 0; i < mFonts.size(); i++) {
					if ((mFonts.get(i)).getSize() == theSize) {
						((Graphics2D) mG).setFont(mFonts.get(i));
						return;
						}
					}
				Font newFont = ((Graphics)mG).getFont().deriveFont(0, theSize);
//				Font newFont = new Font("Helvetica", 0, (int) theSize);
				mFonts.add(newFont);
				((Graphics2D) mG).setFont(newFont);
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
		mLineWidth = (float)lineWidth;
		((Graphics2D)mG).setStroke(new BasicStroke((float)lineWidth, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
		}


	protected void setColor(Color theColor) {
	    ((Graphics2D)mG).setColor(theColor);
		((Graphics2D)mG).setPaint(theColor);
		}
	
	private static boolean isMacintosh()
	{
		return isMac;
	}
	}

