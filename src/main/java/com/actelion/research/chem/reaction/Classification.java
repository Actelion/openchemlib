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

package com.actelion.research.chem.reaction;

public class Classification {
	public static final int MAXINDICES	= 8;		// max number of search strings and masks
	public static final int INDEXLEN	= 4;		// max length (longs) of entry in index files
	public static final int MAXFLASH	= 10;
	public static final int MAXUNITS	= 16;		// max number of unit rxns in query
	public static final int MAXENDUR	= 10;		// max# of handled unchanging functionalities

	public String mClassName;					// name of main class(es)
	public String[] mUnitName;					// names of recognized reaction types
	public int[][][]  mIndex;					// all permutations of all query unit rxns
	public int[] mEnduringFG;					// class no's of enduring FGs
	public int[] mClassResult;					// subClass of each unit rxn (only refu)
	public int[][] mChngGrps;					// changing groups in query
	public int[][] mFlashMol;					// hilite mol/atom list of
	public int[][] mFlashAtom;					// unit rxns for retrieve
	public int[] mEFlashMol;					// hilite mol/atom list of
	public int[] mEFlashAtom;					// unchanged functionalities
	public int[] mMainClass;					// mainClass of each unit rxn
	public int[] mStereoInfo;					// stereo info (if applicable)
	public int mUnitRxns;						// no of detected unit rxns
	public int mEnduringFGs;					// #of unchanged functionality
	public int[] mRearStrandLen;				// rarrangement strand length

	public Classification() {
		mUnitName = new String[MAXUNITS];
		mIndex = new int[MAXUNITS][MAXINDICES][INDEXLEN];
		mEnduringFG = new int[MAXENDUR];
		mClassResult = new int[MAXUNITS];
		mChngGrps = new int[MAXUNITS][4];
		mFlashMol = new int[MAXUNITS][MAXFLASH];
		mFlashAtom = new int[MAXUNITS][MAXFLASH];
		mEFlashMol = new int[MAXENDUR+1];
		mEFlashAtom = new int[MAXENDUR+1];
		mMainClass = new int[MAXUNITS];
		mStereoInfo = new int[MAXUNITS];
		mRearStrandLen = new int[2];
		}
	};
