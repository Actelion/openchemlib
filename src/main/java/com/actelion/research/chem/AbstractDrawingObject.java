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

import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericPoint;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.gui.generic.GenericShape;


public abstract class AbstractDrawingObject {
	protected static final String DESCRIPTOR_START = "<DrawingObject";
	protected static final String DESCRIPTOR_END = "></DrawingObject>";
	protected static final String DESCRIPTOR_TYPE = " type=\"";

	protected GenericPoint[]	mPoint;
	protected boolean			mIsSelected,mProtectedFromDeletion;

	protected double			mTransformationReferenceX,mTransformationReferenceY;
	protected double			mTransformationValue1[];
	protected double			mTransformationValue2[];

	abstract public void draw(GenericDrawContext context, DepictorTransformation t);
	abstract public void hilite(GenericDrawContext context);
	abstract public void clearHiliting();

	/**
	 * Checks, whether this drawing object contains the point at x,y
	 * @param x
	 * @param y
	 * @return
	 */
	abstract public boolean contains(double x, double y);

	abstract public boolean checkHiliting(double x, double y);
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

	public void move(double dx, double dy) {
		if (mPoint != null) {
			for (int i=0; i<mPoint.length; i++) {
				mPoint[i].x += dx;
				mPoint[i].y += dy;
				}
			}
		}

	public void scale(double f) {
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

	public GenericRectangle getBoundingRect(GenericDrawContext context) {
		if (mPoint == null)
			return null;

		GenericRectangle bounds = new GenericRectangle(mPoint[0].x, mPoint[0].y, 0, 0);

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

	public boolean isSurroundedBy(GenericShape shape) {
		if (mPoint == null)
			return false;

		for (int i=0; i<mPoint.length; i++)
			if (!shape.contains(mPoint[i].x, mPoint[i].y))
				return false;

		return true;
		}

	public void translateInit(double x, double y) {
		mTransformationReferenceX = x;
		mTransformationReferenceY = y;
		if (mPoint != null) {
			int pointCount = mPoint.length;
			mTransformationValue1 = new double[pointCount];
			mTransformationValue2 = new double[pointCount];
			for (int i=0; i<pointCount; i++) {
				mTransformationValue1[i] = mPoint[i].x;
				mTransformationValue2[i] = mPoint[i].y;
				}
			}
		}

	public void translate(double x, double y) {
			// overwrite this if only hilited parts of the object shall be moved
		if (mPoint != null) {
			for (int i=0; i<mPoint.length; i++) {
				mPoint[i].x = mTransformationValue1[i] + x - mTransformationReferenceX;
				mPoint[i].y = mTransformationValue2[i] + y - mTransformationReferenceY;
				}
			}
		}

	public void zoomAndRotateInit(double x, double y) {
		mTransformationReferenceX = x;
		mTransformationReferenceY = y;
		if (mPoint != null) {
			int pointCount = mPoint.length;
			mTransformationValue1 = new double[pointCount];
			mTransformationValue2 = new double[pointCount];
			for (int i=0; i<pointCount; i++) {
				double dx = x - mPoint[i].x;
				double dy = y - mPoint[i].y;
				mTransformationValue2[i] = Math.sqrt(dx*dx+dy*dy);	// distance to center of gravity
				mTransformationValue1[i] = Molecule.getAngle(x,y,mPoint[i].x,mPoint[i].y);
				}
			}
		}

	public void zoomAndRotate(double zoom, double angle) {
		if (mPoint != null) {
			for (int i=0; i<mPoint.length; i++) {
				double newDistance = mTransformationValue2[i] * zoom;
				double newAngle = mTransformationValue1[i] - angle;
				mPoint[i].x = mTransformationReferenceX + newDistance*Math.sin(newAngle);
				mPoint[i].y = mTransformationReferenceY + newDistance*Math.cos(newAngle);
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
