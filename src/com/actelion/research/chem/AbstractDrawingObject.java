/*
* Copyright (c) 1997 - 2015
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

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

import javax.swing.UIManager;


public abstract class AbstractDrawingObject {
	protected static final String DESCRIPTOR_START = "<DrawingObject";
	protected static final String DESCRIPTOR_END = "></DrawingObject>";
	protected static final String DESCRIPTOR_TYPE = " type=\"";

    protected static final Color SELECTION_COLOR = UIManager.getColor("TextArea.selectionBackground");

	protected Point2D.Float[]	mPoint;
	protected boolean			mIsSelected,mProtectedFromDeletion;

	protected float			mTransformationReferenceX,mTransformationReferenceY;
	protected float			mTransformationValue1[];
	protected float			mTransformationValue2[];

	abstract public void draw(Graphics g, DepictorTransformation t);
	abstract public void draw2D(Graphics2D g, DepictorTransformation t);
	abstract public void hilite(Graphics g);
	abstract public void clearHiliting();

	/**
	 * Checks, whether this drawing object contains the point at x,y
	 * @param x
	 * @param y
	 * @return
	 */
	abstract public boolean contains(float x, float y);

	abstract public boolean checkHiliting(float x, float y);
	abstract public AbstractDrawingObject clone();

	/**
	 * Creates a string encoding all properties specific to this drawing object.
	 * Its type is not part of the descriptor detail. The detail must start, but not
	 * end with a space. Example: ' size="12" x="50.0" y="80.0"'
	 * @return concatenated property list, each property with a leading space
	 */
	abstract public String getDescriptorDetail();
	abstract public String getTypeString();

        
	public static AbstractDrawingObject instantiate(String descriptor) {
		return DrawingObjectFactory.createObject(descriptor);
		}

	public void move(float dx, float dy) {
		if (mPoint != null) {
			for (int i=0; i<mPoint.length; i++) {
				mPoint[i].x += dx;
				mPoint[i].y += dy;
				}
			}
		}

	public void scale(float f) {
		if (mPoint != null) {
			for (int i=0; i<mPoint.length; i++) {
				mPoint[i].x *= f;
				mPoint[i].y *= f;
				}
			}
		}

	public boolean isDeletable() {
		return !mProtectedFromDeletion;
		}

	public void setDeletable(boolean d) {
		mProtectedFromDeletion = !d;
		}

	public boolean isSelected() {
		return mIsSelected;
		}

	public void setSelected(boolean s) {
		mIsSelected = s;
		}

	public Rectangle2D.Float getBoundingRect() {
		if (mPoint == null)
			return null;

		Rectangle2D.Float bounds = new Rectangle2D.Float();
		bounds.x = mPoint[0].x;
		bounds.y = mPoint[0].y;

		for (int i=1; i<mPoint.length; i++) {
			if (bounds.x > mPoint[i].x) {
				bounds.width += bounds.x - mPoint[i].x;
				bounds.x = mPoint[i].x;
				}
			else if (bounds.width < mPoint[i].x - bounds.x) {
				bounds.width = mPoint[i].x - bounds.x;
				}
			if (bounds.y > mPoint[i].y) {
				bounds.height += bounds.y - mPoint[i].y;
				bounds.y = mPoint[i].y;
				}
			else if (bounds.height < mPoint[i].y - bounds.y) {
				bounds.height = mPoint[i].y - bounds.y;
				}
			}

		return bounds;
		}

	public boolean isSurroundedBy(Shape shape) {
		if (mPoint == null)
			return false;

		for (int i=0; i<mPoint.length; i++)
			if (!shape.contains(mPoint[i]))
				return false;

		return true;
		}

	public void translateInit(float x, float y) {
		mTransformationReferenceX = x;
		mTransformationReferenceY = y;
		if (mPoint != null) {
			int pointCount = mPoint.length;
			mTransformationValue1 = new float[pointCount];
			mTransformationValue2 = new float[pointCount];
			for (int i=0; i<pointCount; i++) {
				mTransformationValue1[i] = mPoint[i].x;
				mTransformationValue2[i] = mPoint[i].y;
				}
			}
		}

	public void translate(float x, float y) {
			// overwrite this if only hilited parts of the object shall be moved
		if (mPoint != null) {
			for (int i=0; i<mPoint.length; i++) {
				mPoint[i].x = mTransformationValue1[i] + x - mTransformationReferenceX;
				mPoint[i].y = mTransformationValue2[i] + y - mTransformationReferenceY;
				}
			}
		}

	public void zoomAndRotateInit(float x, float y) {
		mTransformationReferenceX = x;
		mTransformationReferenceY = y;
		if (mPoint != null) {
			int pointCount = mPoint.length;
			mTransformationValue1 = new float[pointCount];
			mTransformationValue2 = new float[pointCount];
			for (int i=0; i<pointCount; i++) {
				float dx = x - mPoint[i].x;
				float dy = y - mPoint[i].y;
				mTransformationValue2[i] = (float)Math.sqrt(dx*dx+dy*dy);	// distance to center of gravity
				mTransformationValue1[i] = Molecule.getAngle(x,y,mPoint[i].x,mPoint[i].y);
				}
			}
		}

	public void zoomAndRotate(float zoom,float angle) {
		if (mPoint != null) {
			for (int i=0; i<mPoint.length; i++) {
				float newDistance = mTransformationValue2[i] * zoom;
				float newAngle = mTransformationValue1[i] - angle;
				mPoint[i].x = mTransformationReferenceX + newDistance*(float)Math.sin(newAngle);
				mPoint[i].y = mTransformationReferenceY + newDistance*(float)Math.cos(newAngle);
				}
			}
		}

	public String getDescriptor() {
		return DESCRIPTOR_START+" type=\""+getTypeString()+"\""+getDescriptorDetail()+DESCRIPTOR_END;
		}

	public String toString() {
		StringBuffer objectString = new StringBuffer();
        objectString.append(getDescriptor());
		return objectString.toString();
		}

}
