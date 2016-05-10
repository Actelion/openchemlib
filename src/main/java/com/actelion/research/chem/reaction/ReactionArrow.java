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

package com.actelion.research.chem.reaction;

import com.actelion.research.chem.AbstractDrawingObject;
import com.actelion.research.chem.DepictorTransformation;
import com.actelion.research.chem.Molecule;
import com.actelion.research.util.ColorHelper;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

public class ReactionArrow extends AbstractDrawingObject {
	public static final String TYPE_STRING = "arrow";
    
	private static final int PART_NONE = 0;
	private static final int PART_ARROW_START = 1;
	private static final int PART_ARROW_END = 2;
	private static final int PART_ARROW = 3;

	int		mHiliteStatus;

	public ReactionArrow() {
		mPoint = new Point2D.Double[2];
		mPoint[0] = new Point2D.Double();
		mPoint[1] = new Point2D.Double();
		mHiliteStatus = PART_NONE;
		}

    public ReactionArrow(String descriptorDetail) {
        this();

        int index1 = 0;
        while (index1 != -1) {
            // text="this is a test"
            // 012345678901234567890
            int index2 = descriptorDetail.indexOf("=\"", index1);
            if (index2 == -1)
                break;  // should never happen
            String key = descriptorDetail.substring(index1+1, index2);
            index1 = descriptorDetail.indexOf("\"", index2+2);
            String value = (index1 == -1) ? descriptorDetail.substring(index2+1)
                                : descriptorDetail.substring(index2+1, index1);
            if (key.equals("x1"))
                try { mPoint[0].x = Float.parseFloat(value); } catch (NumberFormatException nfe) {}
            else if (key.equals("y1"))
                try { mPoint[0].y = Float.parseFloat(value); } catch (NumberFormatException nfe) {}
            else if (key.equals("x2"))
                try { mPoint[1].x = Float.parseFloat(value); } catch (NumberFormatException nfe) {}
            else if (key.equals("y2"))
                try { mPoint[1].y = Float.parseFloat(value); } catch (NumberFormatException nfe) {}
            }
        }

    public String getTypeString() {
        return TYPE_STRING;   
    	}

    public String getDescriptorDetail() {
		StringBuilder detail = new StringBuilder();
        detail.append(" x1=\""+mPoint[0].x  + "\"");
        detail.append(" y1=\""+mPoint[0].y  + "\"");
        detail.append(" x2=\""+mPoint[1].x  + "\"");
        detail.append(" y2=\""+mPoint[1].y  + "\"");
        return detail.toString();
        }

	public AbstractDrawingObject clone() {
		ReactionArrow arrow = new ReactionArrow();
		arrow.mPoint[0].x = this.mPoint[0].x;
		arrow.mPoint[0].y = this.mPoint[0].y;
		arrow.mPoint[1].x = this.mPoint[1].x;
		arrow.mPoint[1].y = this.mPoint[1].y;
		arrow.mIsSelected = this.mIsSelected;
		return arrow;
		}

	public float getLength() {
		int dx = (int)(mPoint[1].x - mPoint[0].x);
		int dy = (int)(mPoint[1].y - mPoint[0].y);
		return (float)Math.sqrt(dx*dx+dy*dy);
		}

	public void setCoordinates(double x1, double y1, double x2, double y2) {
		mPoint[0].x = x1;
		mPoint[0].y = y1;
		mPoint[1].x = x2;
		mPoint[1].y = y2;
		}

	public void translate(double x, double y) {
		switch (mHiliteStatus) {
		case PART_ARROW_START:
			mPoint[0].x = mTransformationValue1[0] + x - mTransformationReferenceX;
			mPoint[0].y = mTransformationValue2[0] + y - mTransformationReferenceY;
			break;
		case PART_ARROW_END:
			mPoint[1].x = mTransformationValue1[1] + x - mTransformationReferenceX;
			mPoint[1].y = mTransformationValue2[1] + y - mTransformationReferenceY;
			break;
		default:
			super.translate(x, y);
			break;
			}
		}

	public void draw(Graphics g, DepictorTransformation t) {
		g.setColor(mIsSelected ? ColorHelper.getContrastColor(Color.red, UIManager.getColor("TextArea.background"))
							   : UIManager.getColor("TextArea.foreground"));
		double x1 = (t == null) ? mPoint[0].x : t.transformX(mPoint[0].x);
		double y1 = (t == null) ? mPoint[0].y : t.transformY(mPoint[0].y);
		double x2 = (t == null) ? mPoint[1].x : t.transformX(mPoint[1].x);
		double y2 = (t == null) ? mPoint[1].y : t.transformY(mPoint[1].y);
		g.drawLine((int)x1, (int)y1, (int)x2, (int)y2);
		int dx = (int)(x2 - x1);
		int dy = (int)(y2 - y1);
		int[] px = new int[4];
		int[] py = new int[4];
		px[0] = (int)x2;
		py[0] = (int)y2;
		px[1] = (int)(x2 - dx/5 + dy/10);
		py[1] = (int)(y2 - dy/5 - dx/10);
		px[2] = (int)(x2 - dx/15);
		py[2] = (int)(y2 - dy/15);
		px[3] = (int)(x2 - dx/5 - dy/10);
		py[3] = (int)(y2 - dy/5 + dx/10);
		g.fillPolygon(px, py, 4);
		}

	public void hilite(Graphics g) {
		g.setColor(SELECTION_COLOR);
		switch (mHiliteStatus) {
		case PART_ARROW_START:
			g.fillOval((int)mPoint[0].x-8, (int)mPoint[0].y-8, 16, 16);
			break;
		case PART_ARROW_END:
			g.fillOval((int)mPoint[1].x-8, (int)mPoint[1].y-8, 16, 16);
			break;
		case PART_ARROW:
			double length = getLength();
			double f = Math.max(length/8.0f, 3.0f);
			double angle = Molecule.getAngle(mPoint[0].x, mPoint[0].y, mPoint[1].x, mPoint[1].y);
			int dx = (int)(f * Math.cos(angle));
			int dy = -(int)(f * Math.sin(angle));
			int x[] = new int[4];
			int y[] = new int[4];
			x[0] = (int)(mPoint[0].x + dx);
			y[0] = (int)(mPoint[0].y + dy);
			x[1] = (int)(mPoint[1].x + dx);
			y[1] = (int)(mPoint[1].y + dy);
			x[2] = (int)(mPoint[1].x - dx);
			y[2] = (int)(mPoint[1].y - dy);
			x[3] = (int)(mPoint[0].x - dx);
			y[3] = (int)(mPoint[0].y - dy);
			g.fillPolygon(x, y, 4);
			break;
			}
		}

	public void draw2D(Graphics2D g, DepictorTransformation t) {
		draw(g, t);
		}

	public boolean contains(double x, double y) {
		return (findPart(x, y) != PART_NONE);
		}

	public boolean checkHiliting(double x, double y) {
		mHiliteStatus = findPart(x, y);
		return (mHiliteStatus != PART_NONE);
		}

	public void clearHiliting() {
		mHiliteStatus = PART_NONE;
		}

	private int findPart(double x, double y) {
		double distanceToStart = Math.sqrt((mPoint[0].x-x)*(mPoint[0].x-x)+(mPoint[0].y-y)*(mPoint[0].y-y));
		if (distanceToStart < 8.0)
			return PART_ARROW_START;

		double distanceToEnd = Math.sqrt((mPoint[1].x-x)*(mPoint[1].x-x)+(mPoint[1].y-y)*(mPoint[1].y-y));
		if (distanceToEnd < 8.0)
			return PART_ARROW_END;

		double arrowLength = Math.sqrt((mPoint[1].x-mPoint[0].x)*(mPoint[1].x-mPoint[0].x)+(mPoint[1].y-mPoint[0].y)*(mPoint[1].y-mPoint[0].y));
		if (distanceToStart + distanceToEnd < arrowLength + 5)
			return PART_ARROW;

		return PART_NONE;
		}

	public Rectangle2D.Double getBoundingRect() {
		double length = getLength();
		double f = Math.max(length/8.0, 3.0);
		double angle = Molecule.getAngle(mPoint[0].x, mPoint[0].y, mPoint[1].x, mPoint[1].y);
		double dx = Math.abs(f * Math.cos(angle));
		double dy = Math.abs(f * Math.sin(angle));

		Rectangle2D.Double bounds = new Rectangle2D.Double();
		if (mPoint[0].x < mPoint[1].x) {
			bounds.x = mPoint[0].x - dx;
			bounds.width = mPoint[1].x - mPoint[0].x + 2*dx;
			}
		else {
			bounds.x = mPoint[1].x - dx;
			bounds.width = mPoint[0].x - mPoint[1].x + 2*dx;
			}
		if (mPoint[0].y < mPoint[1].y) {
			bounds.y = mPoint[0].y - dy;
			bounds.height = mPoint[1].y - mPoint[0].y + 2*dy;
			}
		else {
			bounds.y = mPoint[1].y - dy;
			bounds.height = mPoint[0].y - mPoint[1].y + 2*dy;
			}

		return bounds;
		}

	public boolean isOnProductSide(double x, double y) {
		double dx = mPoint[1].x - mPoint[0].x;
		double dy = mPoint[1].y - mPoint[0].y;
		double mx = (mPoint[0].x + mPoint[1].x) / 2.0;
		double my = (mPoint[0].y + mPoint[1].y) / 2.0;

		if (dx == 0.0)
			return (dy < 0.0) ^ (y > my);

		if (dy == 0.0)
			return (dx < 0.0) ^ (x > mx);

		double m = -dx/dy;	// m of orthogonal line through S

		double sx = (mPoint[0].x + m*m*x - m*y + m*mPoint[0].y) / (1 + m*m);
		return (dx < 0.0) ^ (sx > mx);
		}
	}
