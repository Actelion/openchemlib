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

package com.actelion.research.gui;

import java.util.EventObject;

public class DrawAreaEvent extends EventObject {
    private static final long serialVersionUID = 0x20090611;

    public static final int TYPE_MOLECULE_CHANGED = 1;
    public static final int TYPE_SELECTION_CHANGED = 2;
    public static final int TYPE_HILITE_ATOM_CHANGED = 3;
    public static final int TYPE_HILITE_BOND_CHANGED = 4;

	private boolean	mIsUserChange;
	private int mType;

    public DrawAreaEvent(Object source, int type, boolean isUserChange) {
		super(source);
		mType = type;
		mIsUserChange = isUserChange;
	    }

	public int getType() {
		return mType;
		}

	public boolean isUserChange() {
		return mIsUserChange;
		}
	}
