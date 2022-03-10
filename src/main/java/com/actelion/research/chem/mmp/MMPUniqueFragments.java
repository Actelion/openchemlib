/*
 * Copyright (c) 2017
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
 * @author Gregori Gerebtzoff
 */

package com.actelion.research.chem.mmp;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

public class MMPUniqueFragments {
	
	private StereoMolecule fragmentMolecule = new StereoMolecule();
	private StereoMolecule fragment = null;
	private int biggestFragment = 0;
	private LinkedHashMap<String, MMPUniqueFragment> uniqueFragments;
	private IDCodeParser idCodeParser = new IDCodeParser();
	private int fragmentIndexes = 0;
	
	public class MMPUniqueFragment {
		private int fragmentIndex;    // index of the whole fragment
		private int fragmentAtoms;    // number of heavy atoms (without considering the R-Group(s))
		private String[] fragmentFP;  // fragments of size 1-6 rooted at the R-Group(s)
		
		private MMPUniqueFragment(String fragmentID) {
			fragmentIndex = fragmentIndexes++;
			fragmentAtoms = idCodeParser.getAtomCount(fragmentID) - 1;
			fragmentFP = null;
		}
		
		private MMPUniqueFragment(int fragmentAtoms, String[] fragmentFP) {
			this.fragmentIndex = fragmentIndexes++;
			this.fragmentAtoms = fragmentAtoms;
			this.fragmentFP = fragmentFP;
		}
		
		/**
		 * Generates fragments of size 2-6 rooted at the R-Group(s)
		 * @param fragmentID idCode of the whole fragment
		 */
		private void addFP(String fragmentID) {
			if (fragmentFP == null) {
				fragmentFP = new String[5];
				idCodeParser.parse(fragmentMolecule, fragmentID);
				fragmentMolecule.ensureHelperArrays(StereoMolecule.cHelperRings);
				if (fragment == null) {
					fragment = new StereoMolecule(fragmentMolecule.getAllAtoms(), fragmentMolecule.getAllBonds());
				}
				int[] atomList = new int[fragmentMolecule.getAtoms()];
				boolean[] atomMask = new boolean[fragmentMolecule.getAtoms()];
				int[] atomMap = new int[fragmentMolecule.getAtoms()];
				Arrays.fill(atomMask, false);
				int rGroupCounter = 0;
				for (int atom=0; atom<fragmentMolecule.getAtoms(); atom++) {
					if (fragmentMolecule.getAtomicNo(atom) == 0 || fragmentMolecule.getAtomicNo(atom) >= 142) {
						atomList[rGroupCounter] = atom;
						atomMask[atom] = true;
						rGroupCounter++;
					}
				}
				int min = 0;
	            int max = rGroupCounter;
	            for (int sphere=0; sphere<5; sphere++) {
	                int newMax = max;
	                for (int i=min; i<max; i++) {
	                    int atom = atomList[i];
	                    for (int j=0; j<fragmentMolecule.getConnAtoms(atom); j++) {
	                        int connAtom = fragmentMolecule.getConnAtom(atom, j);
	                        if (!atomMask[connAtom]) {
	                            atomMask[connAtom] = true;
	                            atomList[newMax++] = connAtom;
	                        }
	                    }
	                }
	                min = max;
	                max = newMax;
	                fragmentMolecule.copyMoleculeByAtoms(fragment, atomMask, true, atomMap);
	                for (int i=0; i<atomMap.length; i++) {
	                	if (atomMap[i] != -1) {
	                		if (fragmentMolecule.isAromaticAtom(atomMap[i])) {
	                			fragment.setAtomQueryFeature(i, StereoMolecule.cAtomQFAromatic, true);
	                		}
	                	}
	                }
	                fragmentFP[sphere] = fragment.getIDCode();
	            }
			}
		}
		
		public int getFragmentIndex()  {
			return fragmentIndex;
		}
		
		public int getFragmentAtoms()  {
			return fragmentAtoms;
		}
		
		public String[] getFragmentFP() {
			return fragmentFP;
		}
	}
	
	public MMPUniqueFragments() {
		uniqueFragments = new LinkedHashMap<String, MMPUniqueFragment>();
	}
	
	/**
	 * Returns the number of heavy atoms of a fragment
	 * @param fragmentID idCode of the fragment
	 * @return number of heavy atoms
	 */
	public Integer getFragmentAtoms(String fragmentID) {
		if (uniqueFragments.containsKey(fragmentID)) {
			return uniqueFragments.get(fragmentID).fragmentAtoms;
		}
		return null;
	}
	
	/**
	 * Returns a MMPUniqueFragment object from a fragment idCode
	 * @param fragmentsID idCode of the fragment
	 * @return a MMPUniqueFragment object
	 */
	public MMPUniqueFragment fragmentIDToFragment(String[] fragmentsID) {
		String fragmentID = fragmentsID[0];
		if (fragmentsID.length == 2) {
			fragmentID = fragmentsID[0] + "\t" + fragmentsID[1];
		}
		return fragmentIDToFragment(fragmentID);
	}
	
	/**
	 * Returns a MMPUniqueFragment object from a fragment idCode
	 * @param fragmentID idCode of the fragment
	 * @return a MMPUniqueFragment object
	 */
	public MMPUniqueFragment fragmentIDToFragment(String fragmentID) {
		MMPUniqueFragment mmpUniqueFragment = null;
		if (uniqueFragments.containsKey(fragmentID)) {
			mmpUniqueFragment = uniqueFragments.get(fragmentID); 
		}
		if (mmpUniqueFragment == null) {
			addFragment(fragmentID);
		}
		if (mmpUniqueFragment != null && mmpUniqueFragment.fragmentFP == null) {
			mmpUniqueFragment.addFP(fragmentID);
		}
		return mmpUniqueFragment;
	}
	
	/**
	 * Adds a new unique fragment to the LinkedHashMap
	 * @param fragmentID idCode of the fragment
	 * @return index of the inserted unique fragment
	 */
	public int addFragment(String fragmentID) {
		if (!uniqueFragments.containsKey(fragmentID)) {
			MMPUniqueFragment mmpFragment = new MMPUniqueFragment(fragmentID);
			uniqueFragments.put(fragmentID, mmpFragment);
			if (mmpFragment.fragmentAtoms > biggestFragment) {
				biggestFragment = mmpFragment.fragmentAtoms; 
			}
			return mmpFragment.fragmentIndex;
		}
		else {
			return uniqueFragments.get(fragmentID).fragmentIndex;
		}
	}
	
	/**
	 * Adds a new unique fragment to the LinkedHashMap
	 * @param fragmentID idCode of the fragment
	 * @param mmpFragment MMPUniqueFragment object
	 */
	public void addFragment(String fragmentID, MMPUniqueFragment mmpFragment) {
		if (!uniqueFragments.containsKey(fragmentID)) {
			uniqueFragments.put(fragmentID, mmpFragment);
		}
	}
	
	/**
	 * Adds a new unique fragment to the LinkedHashMap
	 * @param fragmentID idCode of the whole fragment
	 * @param fragmentAtoms number of heavy atoms of the fragment
	 * @param fragmentFP fragments of size 1-6 rooted at the R-group(s)
	 * @return index of the inserted fragment
	 */
	public int addFragment(String fragmentID, int fragmentAtoms, String[] fragmentFP) {
		if (!uniqueFragments.containsKey(fragmentID)) {
			MMPUniqueFragment mmpFragment = new MMPUniqueFragment(fragmentAtoms, fragmentFP);
			uniqueFragments.put(fragmentID, mmpFragment);
			if (mmpFragment.fragmentAtoms > biggestFragment) {
				biggestFragment = mmpFragment.fragmentAtoms; 
			}
			return mmpFragment.fragmentIndex;
		}
		else {
			return uniqueFragments.get(fragmentID).fragmentIndex;
		}
	}
	
	/**
	 * Writes the Unique Fragments block
	 * @param printWriter
	 */
	public void writeUniqueFragments(PrintWriter printWriter) throws IOException {
		printWriter.println("<mmpUniqueFragments>");
		printWriter.println("<column properties>");
		printWriter.println("<columnName=\"fragmentID\">");
		printWriter.println("<columnProperty=\"specialType	idcode\">");
		printWriter.println("<columnName=\"fragmentAtoms\">");
		printWriter.println("<columnName=\"fragmentFP1\">");
		printWriter.println("<columnProperty=\"specialType	idcode\">");
		printWriter.println("<columnName=\"fragmentFP2\">");
		printWriter.println("<columnProperty=\"specialType	idcode\">");
		printWriter.println("<columnName=\"fragmentFP3\">");
		printWriter.println("<columnProperty=\"specialType	idcode\">");
		printWriter.println("<columnName=\"fragmentFP4\">");
		printWriter.println("<columnProperty=\"specialType	idcode\">");
		printWriter.println("<columnName=\"fragmentFP5\">");
		printWriter.println("<columnProperty=\"specialType	idcode\">");
		printWriter.println("</column properties>");
		printWriter.println("fragmentID\tfragmentAtoms\tfragmentFP1\tfragmentFP2\tfragmentFP3\tfragmentFP4\tfragmentFP5");
		if (biggestFragment != 0) {
			fragment = new StereoMolecule(biggestFragment+1, biggestFragment+1);
		}
		for (Map.Entry<String, MMPUniqueFragment> entry: uniqueFragments.entrySet()) {
			String fragmentID = entry.getKey();
			MMPUniqueFragment fragmentFP = entry.getValue();
			if (fragmentFP.fragmentFP == null) {
				fragmentFP.addFP(fragmentID);
				entry.setValue(fragmentFP);
			}
			String printLine = fragmentID + "\t" + Integer.toString(fragmentFP.fragmentAtoms) + "\t" + fragmentFP.fragmentFP[0] + "\t" + fragmentFP.fragmentFP[1] + "\t" + fragmentFP.fragmentFP[2] + "\t" + fragmentFP.fragmentFP[3] + "\t" + fragmentFP.fragmentFP[4];
			printWriter.println(printLine);
		}
		printWriter.println("</mmpUniqueFragments>");
	}
	
	public int getUniqueFragmentsCount() {
		return uniqueFragments.size();
	}
}
