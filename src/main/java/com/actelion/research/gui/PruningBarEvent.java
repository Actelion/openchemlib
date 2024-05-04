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
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import java.util.EventObject;

public class PruningBarEvent extends EventObject {
    private static final long serialVersionUID = 0x20101217;

    public static final int TYPE_DRAGGED = 1;
	public static final int TYPE_TYPED = 2;

    private double mLowValue,mHighValue;
	private int mType,mID;
	private boolean mAdjusting;

	public PruningBarEvent(Object source, double low, double high, boolean adjusting, int id) {
		this(source, low, high, adjusting, id, TYPE_DRAGGED);
		}

    public PruningBarEvent(Object source, double low, double high, boolean adjusting, int id, int type) {
		super(source);
		mLowValue = low;
		mHighValue = high;
		mAdjusting = adjusting;
		mID = id;
		mType = type;
	    }

	public double getLowValue() {
		return mLowValue;
		}

	public double getHighValue() {
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