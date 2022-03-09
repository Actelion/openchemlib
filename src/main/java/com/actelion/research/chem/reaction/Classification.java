/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * 3. Neither the name of the copyright holder nor the
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
