/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
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

package com.actelion.research.gui;

import java.util.EventObject;

public class PruningBarEvent extends EventObject {
    private static final long serialVersionUID = 0x20101217;

    public static final int TYPE_DRAGGED = 1;
	public static final int TYPE_TYPED = 2;

    private float mLowValue,mHighValue;
	private int mType,mID;
	private boolean mAdjusting;

	public PruningBarEvent(Object source, float low, float high, boolean adjusting, int id) {
		this(source, low, high, adjusting, id, TYPE_DRAGGED);
		}

    public PruningBarEvent(Object source, float low, float high, boolean adjusting, int id, int type) {
		super(source);
		mLowValue = low;
		mHighValue = high;
		mAdjusting = adjusting;
		mID = id;
		mType = type;
	    }

	public float getLowValue() {
		return mLowValue;
		}

	public float getHighValue() {
		return mHighValue;
		}

	public boolean isAdjusting() {
		return mAdjusting;
		}

	public int getID() {
		return mID;
		}

	public int getType() {
		return mType;
		}
	}