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

package com.actelion.research.chem.conf;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.phesa.PheSAAlignment;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

public class BondRotationHelper {
	
	private StereoMolecule mMol;
	private int[] mRotatableBonds;
	private boolean[] mIsRotatableBond;
	private int[][] mSmallerSideAtomLists;
	private int[][] mBiggerSideAtomLists;
	private int[][] mTorsionAtoms;
	private int[][] mRearAtoms;
	private int[] mRotationCenters;
	private int[] mRotationCentersBig;
	private String[] mTorsionIDs;
	private boolean includeTerminalPolarH;
	private int[] terminalPolarHBond;
	
	public BondRotationHelper(StereoMolecule mol) {
		this(mol,false);
	}
	
	public BondRotationHelper(StereoMolecule mol, boolean includeTerminalPolarH) {
		mMol = mol;
		this.includeTerminalPolarH = includeTerminalPolarH;
		initialize();
	}
	
	public void initialize() {
		int[] disconnectedFragmentNo = new int[mMol.getAllAtoms()];
		int disconnectedFragmentCount = mMol.getFragmentNumbers(disconnectedFragmentNo, false, true);
		int disconnectedFragmentSize[] = new int[disconnectedFragmentCount];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			disconnectedFragmentSize[disconnectedFragmentNo[atom]]++;
		mIsRotatableBond = new boolean[mMol.getBonds()];
		TorsionDB.findRotatableBonds(mMol,true, mIsRotatableBond);
		if(includeTerminalPolarH)
			terminalPolarHBond = findTerminalBondsPolarHs(mIsRotatableBond);
		List<Integer> rotBonds = new ArrayList<Integer>();
		IntStream.range(0, mIsRotatableBond.length).forEach(e -> {
			if(mIsRotatableBond[e])
				rotBonds.add(e);
		});
		mRotatableBonds = rotBonds.stream().mapToInt(i->i).toArray();
		mTorsionAtoms = new int[mRotatableBonds.length][];
		mRearAtoms = new int[mRotatableBonds.length][];
		mSmallerSideAtomLists = new int[mRotatableBonds.length][];
		mBiggerSideAtomLists = new int[mRotatableBonds.length][];
		mRotationCenters = new int[mRotatableBonds.length];
		mRotationCentersBig = new int[mRotatableBonds.length];
		mTorsionIDs = new String[mRotatableBonds.length];
		
		for(int i=0;i<mRotatableBonds.length;i++) {
			int bond = mRotatableBonds[i];
			boolean isSpecialBond = false;
			if(terminalPolarHBond!=null) {
				for(int b : terminalPolarHBond) {
					if(b==bond)
						isSpecialBond=true;
				}
			}
			int[] torsionAtoms = new int[4];
			int[] rearAtoms = new int[2];
			TorsionDetail detail = new TorsionDetail();
			String torsionID = TorsionDB.getTorsionID(mMol, bond, torsionAtoms, detail);
			mTorsionIDs[i] = torsionID;
			if(isSpecialBond) {
				int a1 = mMol.getBondAtom(0, bond);
				int a2 = mMol.getBondAtom(1, bond);
				int a0 = mMol.getConnAtom(a1, 0)==a2 ? mMol.getConnAtom(a1, 1) : mMol.getConnAtom(a1, 0);
				int a3 = mMol.getConnAtom(a2, 0)==a1 ? mMol.getConnAtom(a2, 1) : mMol.getConnAtom(a2, 0);
				rearAtoms[0] = a2;
				rearAtoms[1] = a1;
				torsionAtoms[0] = a0;
				torsionAtoms[1] = a1;
				torsionAtoms[2] = a2;
				torsionAtoms[3] = a3;
			}
			else if (torsionID != null) {
				rearAtoms[0] = detail.getRearAtom(0);
				rearAtoms[1] = detail.getRearAtom(1);
				}
			else {
				predictAtomSequence(mMol,bond, torsionAtoms, rearAtoms);
				}
			mTorsionAtoms[i] = torsionAtoms;
			mRearAtoms[i] = rearAtoms;
			findSmallerSideAtomList(disconnectedFragmentSize[disconnectedFragmentNo[mMol.getBondAtom(0, bond)]],
					disconnectedFragmentNo, i);
		}
	}
	
	private int[] findTerminalBondsPolarHs(boolean[] isRotatableBond) {
		List<Integer> terminalHBonds = new ArrayList<>();
		for(int b=0;b<mMol.getBonds();b++) {
			int a1 = mMol.getBondAtom(0, b);
			int a2 = mMol.getBondAtom(1, b);
			if(mMol.getAtomicNo(a1)==7 || mMol.getAtomicNo(a1)==8) {
				if(mMol.getNonHydrogenNeighbourCount(a1)==1 && mMol.getAllHydrogens(a1)>0) {
					isRotatableBond[b]=true;
					terminalHBonds.add(b);
					continue;
				}
			}
			if(mMol.getAtomicNo(a2)==7 || mMol.getAtomicNo(a2)==8) {
				if(mMol.getNonHydrogenNeighbourCount(a2)==1 && mMol.getAllHydrogens(a2)>0) {
					isRotatableBond[b]=true;
					terminalHBonds.add(b);
				}
			}
		}
		return terminalHBonds.stream().mapToInt(Integer::intValue).toArray();
	}
	
	//populates smallerSideAtomList array and returns the rotation center
	private void findSmallerSideAtomList(int disconnectedFragmentSize, int[] disconnectedFragmentNo, int bondIndex) {
		boolean[] isMember = new boolean[mMol.getAllAtoms()];
		int memberCount = mMol.getSubstituent(mRearAtoms[bondIndex][0], mTorsionAtoms[bondIndex][1], isMember, null, null);
		
		int alkyneAtoms = 0;	// if we have an extended linear sp-atom strain
		if (mRearAtoms[bondIndex][0] != mTorsionAtoms[bondIndex][2])
			alkyneAtoms = mMol.getPathLength(mRearAtoms[bondIndex][0], mTorsionAtoms[bondIndex][2]);

		boolean invert = false;
		if (memberCount > disconnectedFragmentSize-alkyneAtoms-memberCount) {
			memberCount = disconnectedFragmentSize-alkyneAtoms-memberCount;
			invert = true;
			}

		// if invert, then flag all linear alkyne atoms to be avoided
		if (invert && alkyneAtoms != 0) {
			int spAtom = mRearAtoms[bondIndex][0];
			int backAtom = mTorsionAtoms[bondIndex][1];
        	while (mMol.getAtomPi(spAtom) == 2
           		&& mMol.getConnAtoms(spAtom) == 2
           		&& mMol.getAtomicNo(spAtom) < 10) {
        		isMember[spAtom] = true;
           		for (int j=0; j<2; j++) {
           			int connAtom = mMol.getConnAtom(spAtom, j);
           			if (connAtom != backAtom) {
           				backAtom = spAtom;
           				spAtom = connAtom;
           				break;
           				}
           			}
           		}
			}

		int memberNo = 0;
		int fragmentNo = disconnectedFragmentNo[mTorsionAtoms[bondIndex][1]];
		int [] smallerSideAtomList = new int[memberCount];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			if (disconnectedFragmentNo[atom] == fragmentNo && (isMember[atom] ^ invert))
				smallerSideAtomList[memberNo++] = atom;


		int rotationCenter = mTorsionAtoms[bondIndex][invert ? 2 : 1];
		mRotationCenters[bondIndex] = rotationCenter;
		mRotationCentersBig[bondIndex] = mTorsionAtoms[bondIndex][invert ? 1 : 2];
		mSmallerSideAtomLists[bondIndex] = smallerSideAtomList;
		Set<Integer> bigSideAtoms = new HashSet<Integer>();
		Set<Integer> smallSideAtoms = new HashSet<Integer>();
		for(int a : smallerSideAtomList)
			smallSideAtoms.add(a);
		for(int a=0;a<mMol.getAllAtoms();a++) {
			if(!smallSideAtoms.contains(a)) {
				bigSideAtoms.add(a);
			}
		}
		mBiggerSideAtomLists[bondIndex] = bigSideAtoms.stream().mapToInt(Integer::intValue).toArray();
		}
	
	
	public static void predictAtomSequence(StereoMolecule mol,int bond, int[] torsionAtoms, int[] rearAtoms) {
        for (int i=0; i<2; i++) {
    		int centralAtom = mol.getBondAtom(i, bond);
        	int rearAtom = mol.getBondAtom(1-i, bond);

        	// walk along sp-chains to first sp2 or sp3 atom
        	while (mol.getAtomPi(centralAtom) == 2
        		&& mol.getConnAtoms(centralAtom) == 2
        		&& mol.getAtomicNo(centralAtom) < 10) {
        		for (int j=0; j<2; j++) {
        			int connAtom = mol.getConnAtom(centralAtom, j);
        			if (connAtom != rearAtom) {
        				rearAtom = centralAtom;
        				centralAtom = connAtom;
        				break;
        				}
        			}
        		}

        	torsionAtoms[i+1] = centralAtom;
           	rearAtoms[i] = rearAtom;
        	}

    	// A TorsionPrediction does not distinguish hetero atoms from carbons a positions 0 and 3.
        // Therefore we can treat two sp2 neighbors as equivalent when predicting torsions.
        if (mol.getAtomPi(torsionAtoms[1]) == 0 && mol.getConnAtoms(torsionAtoms[1]) == 3) {
        	torsionAtoms[0] = -1;
        	}
        else {
			for (int i=0; i<mol.getConnAtoms(torsionAtoms[1]); i++) {
				int connAtom = mol.getConnAtom(torsionAtoms[1], i);
				if (connAtom != torsionAtoms[2]) {
					torsionAtoms[0] = connAtom;
					break;
					}
				}
        	}
        if (mol.getAtomPi(torsionAtoms[2]) == 0 && mol.getConnAtoms(torsionAtoms[2]) == 3) {
        	torsionAtoms[3] = -1;
        	}
        else {
			for (int i=0; i<mol.getConnAtoms(torsionAtoms[2]); i++) {
				int connAtom = mol.getConnAtom(torsionAtoms[2], i);
				if (connAtom != torsionAtoms[1]) {
					torsionAtoms[3] = connAtom;
					break;
					}
				}
        	}
		}
	
	public boolean isRotatableBond(int bond) {
		return mIsRotatableBond[bond];
	}

	public int[] getRotatableBonds() {
		return mRotatableBonds;
	}

	public void rotateSmallerSide(int bond, double alpha) {
		if(!isRotatableBond(bond))
			return;
		int bondIndex=-1;
		for(int i=0;i<mRotatableBonds.length;i++) {
			if(mRotatableBonds[i]==bond)
				bondIndex = i;
		}
		if(bondIndex==-1)
			return;
		Coordinates t2 = mMol.getCoordinates(mTorsionAtoms[bondIndex][2]);
		Coordinates unit = t2.subC(mMol.getCoordinates(mTorsionAtoms[bondIndex][1])).unit();
		double[][] m = new double[3][3];
		int rotationCenter = mRotationCenters[bondIndex];
		PheSAAlignment.getRotationMatrix((rotationCenter == mTorsionAtoms[bondIndex][1]) ? alpha : -alpha,unit,m);
		Coordinates t2Neg = t2.scaleC(-1.0);
		for (int atom:mSmallerSideAtomLists[bondIndex]) {
			if (atom != rotationCenter) {
				mMol.getCoordinates(atom).add(t2Neg);
				mMol.getCoordinates(atom).rotate(m);
				mMol.getCoordinates(atom).add(t2);

				}
			}
		}
	
	/**
	 * rotate torsion angle of a conformer
	 * @param bondIndex
	 * @param alpha
	 * @param conf
	 * @param biggerSide
	 */
	public void rotateAroundBond(int bondIndex, double alpha, Conformer conf, boolean biggerSide) {
		int[] atomList;
		Coordinates t2;
		Coordinates unit;
		Coordinates t2Neg;
		double[][] m;
		int rotationCenter;
		if(biggerSide) {
			atomList = mBiggerSideAtomLists[bondIndex];
			t2 = conf.getCoordinates(mTorsionAtoms[bondIndex][1]);
			unit = t2.subC(conf.getCoordinates(mTorsionAtoms[bondIndex][2])).unit();
			m = new double[3][3];
			rotationCenter = mRotationCentersBig[bondIndex];
			PheSAAlignment.getRotationMatrix((rotationCenter == mTorsionAtoms[bondIndex][2]) ? alpha : -alpha,unit,m);
			t2Neg = t2.scaleC(-1.0);
		}
		else {
			
			t2 = conf.getCoordinates(mTorsionAtoms[bondIndex][2]);
			unit = t2.subC(conf.getCoordinates(mTorsionAtoms[bondIndex][1])).unit();
			m = new double[3][3];
			rotationCenter = mRotationCenters[bondIndex];
			PheSAAlignment.getRotationMatrix((rotationCenter == mTorsionAtoms[bondIndex][1]) ? alpha : -alpha,unit,m);
			t2Neg = t2.scaleC(-1.0);
			atomList = mSmallerSideAtomLists[bondIndex];
		}
		for (int atom:atomList) {
			if (atom != rotationCenter) {
				Coordinates coords = conf.getCoordinates(atom);

				coords.add(t2Neg);
				coords.rotate(m);
				coords.add(t2);
				}
			}
		}
		
	public void setRotatableBonds(int[] rotatableBonds) {
		mRotatableBonds = rotatableBonds;
	}

	public int[][] getSmallerSideAtomLists() {
		return mSmallerSideAtomLists;
	}

	public void setSmallerSideAtomLists(int[][] smallerSideAtomLists) {
		mSmallerSideAtomLists = smallerSideAtomLists;
	}

	public int[][] getTorsionAtoms() {
		return mTorsionAtoms;
	}

	public void setTorsionAtoms(int[][] torsionAtoms) {
		mTorsionAtoms = torsionAtoms;
	}

	public int[][] getRearAtoms() {
		return mRearAtoms;
	}

	public void setRearAtoms(int[][] rearAtoms) {
		this.mRearAtoms = rearAtoms;
	}

	public int[] getRotationCenters() {
		return mRotationCenters;
	}

	public void setRotationCenters(int[] rotationCenters) {
		this.mRotationCenters = rotationCenters;
	}

	public String[] getTorsionIDs() {
		return mTorsionIDs;
	}

	public void setTorsionIDs(String[] torsionIDs) {
		this.mTorsionIDs = torsionIDs;
	}
}
