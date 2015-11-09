/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
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
	private float      mCurrentStringWidth;
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
			((Graphics2D)mG).draw(new Line2D.Float(theLine.x1-0.5f, theLine.y1-0.5f, theLine.x2-0.5f, theLine.y2-0.5f));
		else
			((Graphics2D)mG).draw(new Line2D.Float(theLine.x1, theLine.y1, theLine.x2, theLine.y2));
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


	protected void drawString(String theString, float x, float y) {
		float strWidth = getStringWidth(theString);
		((Graphics2D)mG).drawGlyphVector(mCurrentGlyphVector,
									     (float)(x-strWidth/2.0),
									     (float)(y+(float)mpTextSize/3.0));
		}


	protected void drawPolygon(float[] x, float[] y, int count) {
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


	protected void fillCircle(float x, float y, float r) {
		if (isMacintosh())
			((Graphics2D)mG).fill(new Ellipse2D.Float(x-0.5f, y-0.5f, r, r));
		else
			((Graphics2D)mG).fill(new Ellipse2D.Float(x, y, r, r));
		}


	protected float getStringWidth(String theString) {
		if (!theString.equals(mCurrentString)
		 || mCurrentFont != ((Graphics2D)mG).getFont()) {
			mCurrentString = theString;
			mCurrentFont = ((Graphics2D)mG).getFont();
			mCurrentGlyphVector = ((Graphics2D)mG).getFont().createGlyphVector(((Graphics2D)mG).getFontRenderContext(), theString);
			mCurrentStringWidth = (float)mCurrentGlyphVector.getLogicalBounds().getWidth();
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
				Font newFont = new Font("Helvetica", 0, (int) theSize);
				mFonts.add(newFont);
				((Graphics2D) mG).setFont(newFont);
				}
			}
		}


    public int getTextSize() {
        return mpTextSize;
        }


	protected float getLineWidth() {
		return mLineWidth;
		}


	protected void setLineWidth(float lineWidth) {
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

