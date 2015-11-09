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

package com.actelion.research.chem.prediction;

import java.util.ArrayList;


public class ParameterizedStringList /* implements Serializable */ {
	public static final int cStringTypeIDCode = 1;
	public static final int cStringTypeText = 2;
	public static final int cStringTypeDouble = 3;

	private ArrayList<ParameterizedString> mList = new ArrayList<ParameterizedString>();

	public ParameterizedStringList() {
		}

	public void add(String s, int t) {
        mList.add(new ParameterizedString(s, t));
		}

	public int getSize(){
		return mList.size();
		}

	public String getStringAt(int i) {
		return mList.get(i).mString;
		}

	public int getStringTypeAt(int i) {
		return mList.get(i).mType;
		}
	}


class ParameterizedString /* implements Serializable */ {
	protected String	mString;
	protected int		mType;

    protected ParameterizedString() {
		mString = null;
		mType = 0;
    	}

	protected ParameterizedString(String s, int t) {
		mString = s;
		mType = t;
		}
	}
