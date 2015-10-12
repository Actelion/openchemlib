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

package com.actelion.research.chem.reaction;

import java.awt.*;
import java.awt.geom.*;

import com.actelion.research.chem.AbstractDrawingObject;
import com.actelion.research.chem.DepictorTransformation;
import com.actelion.research.chem.Molecule;

public class ReactionArrow extends AbstractDrawingObject {
	public static final String TYPE_STRING = "arrow";
    
	private static final int PART_NONE = 0;
	private static final int PART_ARROW_START = 1;
	private static final int PART_ARROW_END = 2;
	private static final int PART_ARROW = 3;

	int		mHiliteStatus;

	public ReactionArrow() {
		mPoint = new Point2D.Float[2];
		mPoint[0] = new Point2D.Float();
		mPoint[1] = new Point2D.Float();
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

	public void setCoordinates(float x1, float y1, float x2, float y2) {
		mPoint[0].x = x1;
		mPoint[0].y = y1;
		mPoint[1].x = x2;
		mPoint[1].y = y2;
		}

	public void translate(float x, float y) {
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
		g.setColor(mIsSelected ? Color.red : Color.black);
		float x1 = (t == null) ? mPoint[0].x : t.transformX(mPoint[0].x);
		float y1 = (t == null) ? mPoint[0].y : t.transformY(mPoint[0].y);
		float x2 = (t == null) ? mPoint[1].x : t.transformX(mPoint[1].x);
		float y2 = (t == null) ? mPoint[1].y : t.transformY(mPoint[1].y);
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
			float length = getLength();
			float f = Math.max(length/8.0f, 3.0f);
			float angle = Molecule.getAngle(mPoint[0].x, mPoint[0].y, mPoint[1].x, mPoint[1].y);
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

	public boolean contains(float x, float y) {
		return (findPart(x, y) != PART_NONE);
		}

	public boolean checkHiliting(float x, float y) {
		mHiliteStatus = findPart(x, y);
		return (mHiliteStatus != PART_NONE);
		}

	public void clearHiliting() {
		mHiliteStatus = PART_NONE;
		}

	private int findPart(float x, float y) {
		float distanceToStart = (float)Math.sqrt((mPoint[0].x-x)*(mPoint[0].x-x)+(mPoint[0].y-y)*(mPoint[0].y-y));
		if (distanceToStart < 8.0)
			return PART_ARROW_START;

		float distanceToEnd = (float)Math.sqrt((mPoint[1].x-x)*(mPoint[1].x-x)+(mPoint[1].y-y)*(mPoint[1].y-y));
		if (distanceToEnd < 8.0)
			return PART_ARROW_END;

		float arrowLength = (float)Math.sqrt((mPoint[1].x-mPoint[0].x)*(mPoint[1].x-mPoint[0].x)+(mPoint[1].y-mPoint[0].y)*(mPoint[1].y-mPoint[0].y));
		if (distanceToStart + distanceToEnd < arrowLength + 5)
			return PART_ARROW;

		return PART_NONE;
		}

	public Rectangle2D.Float getBoundingRect() {
		float length = getLength();
		float f = Math.max(length/8.0f, 3.0f);
		float angle = Molecule.getAngle(mPoint[0].x, mPoint[0].y, mPoint[1].x, mPoint[1].y);
		float dx = Math.abs(f * (float)Math.cos(angle));
		float dy = Math.abs(f * (float)Math.sin(angle));

		Rectangle2D.Float bounds = new Rectangle2D.Float();
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

	public boolean isOnProductSide(float x, float y) {
		float dx = mPoint[1].x - mPoint[0].x;
		float dy = mPoint[1].y - mPoint[0].y;
		float mx = (mPoint[0].x + mPoint[1].x) / 2.0f;
		float my = (mPoint[0].y + mPoint[1].y) / 2.0f;

		if (dx == 0.0)
			return (dy < 0.0) ^ (y > my);

		if (dy == 0.0)
			return (dx < 0.0) ^ (x > mx);

		float m = -dx/dy;	// m of orthogonal line through S

		float sx = (mPoint[0].x + m*m*x - m*y + m*mPoint[0].y) / (1 + m*m);
		return (dx < 0.0) ^ (sx > mx);
		}
	}
