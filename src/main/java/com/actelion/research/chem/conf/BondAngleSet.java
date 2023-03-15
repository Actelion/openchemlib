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

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;

public class BondAngleSet {
	private static float TO_RADIAN = (float)Math.PI/180;

	private final StereoMolecule	mMol;
	private final BondLengthSet		mBondLengthSet;
	private final float[][][]		mBondAngle;
	private final int[]				mDefinedAngleCount,mTinyRingSizeSum;
	private final float[]			mDefinedAngleSum;

	/**
	 * Calculates and caches a list of bond angle estimates for any two neighbours
	 * of any atom of the molecule.
	 * Internally it requires a valid BondLengthSet of this molecule,
	 * which can be passed if there is already one available. Otherwise,
	 * one is created internally.
	 * @param mol
	 * @param set null or a valid BondLengthSet of the molecule
	 */
	public BondAngleSet(final StereoMolecule mol, final BondLengthSet set) {
		// initialize bond angle arrays
		mBondLengthSet = (set == null) ? new BondLengthSet(mol) : set;
		mMol = mol;

		mBondAngle = new float[mMol.getAtoms()][][];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			mBondAngle[atom] = new float[mMol.getAllConnAtoms(atom)][];
			for (int i=1; i<mMol.getAllConnAtoms(atom); i++)
				mBondAngle[atom][i] = new float[i];
			}
	
		mDefinedAngleCount = new int[mMol.getAtoms()];
		mDefinedAngleSum = new float[mMol.getAtoms()];
		mTinyRingSizeSum = new int[mMol.getBonds()];
	
			// define all angles between bonds of flat rings (aromatic or ringsize<=4)
		RingCollection ringSet = mMol.getRingSet();
		boolean[] isAromaticRing = new boolean[ringSet.getSize()];
		ringSet.determineAromaticity(isAromaticRing, new boolean[ringSet.getSize()], new int[ringSet.getSize()], true);
		for (int consideredRingSize=3; consideredRingSize<=7; consideredRingSize++) {
			for (int ring=0; ring<ringSet.getSize(); ring++) {
				int ringSize = ringSet.getRingSize(ring);
				if (ringSize == consideredRingSize) {
					if (isAromaticRing[ring])
						calculateBondAnglesOfAromaticRing(ring);
					else if (ringSize <= 4)
						calculateBondAnglesOfSmallRing(ring);
					}
				}
			}
		// We include tautomeric bond states to detect aromaticity leading to flatness of rings.
		// When later checking atoms for being aromatic, we need to apply the same concept.
		boolean[] isAromaticAtom = new boolean[mMol.getAtoms()];
		for (int ring=0; ring<ringSet.getSize(); ring++)
			if (isAromaticRing[ring])
				for (int ringAtom:ringSet.getRingAtoms(ring))
					isAromaticAtom[ringAtom] = true;

			// define remaining angles between yet undefined angles between vicinal bonds
		final int[] cTotalAngleCount = { 0, 0, 1, 3, 6, 10, 15, 21 };
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			int connAtoms = mMol.getAllConnAtoms(atom);
			if (connAtoms > 4) {
				for (int i=1; i<connAtoms; i++)
					for (int j=0; j<i; j++)
						mBondAngle[atom][i][j] = (float)(Math.PI / 2.0);
	
				mDefinedAngleCount[atom] = cTotalAngleCount[connAtoms];
				continue;
				}
	
			if (mDefinedAngleCount[atom] == cTotalAngleCount[connAtoms])
				continue;
	
				// exocyclic bonds of flat ring atoms (aromatic or sp2 in ringsize<=4)
			if (mMol.isSmallRingAtom(atom)
			 && (isAromaticAtom[atom]
			  || (mMol.getAtomRingSize(atom) <= 4
			   && mMol.getAtomPi(atom) > 0))) {
				if (connAtoms > 2) {
					float angle;
					if (mDefinedAngleCount[atom] == 1) {
						if (mMol.getAtomicNo(atom) <= 14)   // up to Si
							angle = (2f * (float)Math.PI - mDefinedAngleSum[atom]) / 2.0f;
						else
							angle = calculateRemainingTetrahedralAngle(mDefinedAngleSum[atom]);
						}
					else {
						// mDefinedAngleCount==2: we assume two annelated aromatic rings
						// (and not spiro, which theoretically might be possible with e.g. sulphur)
						angle = 2f * (float)Math.PI - mDefinedAngleSum[atom];
						if (connAtoms > 3) {
							// strange case with P,S,Se,etc sharing two/three annelated flat rings and additional substituent(s) sticking out of plane
							if (mDefinedAngleCount[atom] == 2) {  // we need to find the third angle in the plane
								boolean[] isMissingIndex = new boolean[mMol.getAllConnAtoms(atom)];
								for (int i=1; i<connAtoms; i++) {
									for (int j=0; j<i; j++) {
										if (mBondAngle[atom][i][j] != 0.0) {
											isMissingIndex[i] = !isMissingIndex[i];
											isMissingIndex[j] = !isMissingIndex[j];
											}
										}
									}
								for (int i=0; i<connAtoms; i++) {
									if (isMissingIndex[i]) {
										for (int j=i+1; j<connAtoms; j++) {
											if (isMissingIndex[j]) {
												mBondAngle[atom][j][i] = angle;
												break;
												}
											}
										break;
										}
									}
								}
							angle = (float)Math.PI / 2;
							}
						}
					for (int i=1; i<connAtoms; i++)
						for (int j=0; j<i; j++)
							if (mBondAngle[atom][i][j] == 0.0)
								mBondAngle[atom][i][j] = angle;
					}
				}
	
				// sp3 atoms within 3- and 4-membered rings
			else if (mMol.isSmallRingAtom(atom)
			      && mMol.getAtomRingSize(atom) <= 4) {
				switch (getRingStrainClass(atom)) {
				case 0x00000303:	// cyclopropane
					setRingStrainAngles(atom, 0, 3, 2.0654f);
					setRingStrainAngles(atom, 0, 0, 113.53f * TO_RADIAN);	// calculated from MOPAC dimethylcyclopropane
					break;
				case 0x00000404:	// cyclobutane
					setRingStrainAngles(atom, 0, 4, 1.9798f);
					setRingStrainAngles(atom, 0, 0, 111.55f * TO_RADIAN);	// calculated from MOPAC dimethylcyclobutane
					break;
				case 0x00060303:	// bicyclobutane
					setRingStrainAngles(atom, 0, 3, 2.2640f);
					setRingStrainAngles(atom, 0, 6, 2.2640f);
					setRingStrainAngles(atom, 3, 3, 98.715f * TO_RADIAN);	// calculated from MOPAC 1-methylbicyclobutane
					break;
				case 0x00070403:	// bicyclopentane
					setRingStrainAngles(atom, 0, 3, 2.1676f);
					setRingStrainAngles(atom, 0, 4, 2.1676f);
					setRingStrainAngles(atom, 0, 7, 2.1676f);
					setRingStrainAngles(atom, 3, 4, 110.71f * TO_RADIAN);	// calculated from MOPAC 1-methylbicyclopentane
					break;
				case 0x00080404:	// bicyclohexane
					setRingStrainAngles(atom, 0, 4, 2.0663f);
					setRingStrainAngles(atom, 0, 8, 2.0663f);
					setRingStrainAngles(atom, 4, 4, 114.56f * TO_RADIAN);	// calculated from MOPAC
				case 0x00060606:	// tricyclobutane
					setRingStrainAngles(atom, 0, 6, 2.5261f);
					break;
				case 0x00070706:	// tricyclopentane
					setRingStrainAngles(atom, 0, 6, 2.3562f);
					setRingStrainAngles(atom, 0, 7, 2.3562f);
					break;
				case 0x00080707:	// tricyclohexane
					setRingStrainAngles(atom, 0, 7, 2.2845f);
					setRingStrainAngles(atom, 0, 8, 2.2845f);
					break;
				case 0x00080808:	// tricycloheptane
					setRingStrainAngles(atom, 0, 8, 2.1863f);
					break;
				case 0x03030303:	// 3,3-spiro-pentan
					setRingStrainAngles(atom, 3, 3, 2.4189f);
					break;
				case 0x04040303:	// 4,3-spiro-hexan
					setRingStrainAngles(atom, 3, 4, 2.2299f);
					break;
				case 0x04040404:	// 4,4-spiro-heptan
					setRingStrainAngles(atom, 4, 4, 2.0944f);
					break;
				case 0x06060303:	// tricyclopentane	sKp@H~JYFjj@@
					setRingStrainAngles(atom, 3, 6, 105.42f * TO_RADIAN);	// calculated from MOPAC
					setRingStrainAngles(atom, 3, 3, 166.51f * TO_RADIAN);	// calculated from MOPAC
					break;
				case 0x07060403:	// tricyclohexane	pfH@BUSDckUUP@@
					setRingStrainAngles(atom, 3, 4, 161.13f * TO_RADIAN);	// calculated from MOPAC
					setRingStrainAngles(atom, 3, 7, 102.62f * TO_RADIAN);	// calculated from MOPAC
					setRingStrainAngles(atom, 4, 6, 121.61f * TO_RADIAN);	// calculated from MOPAC
					break;
				case 0x08070403:	// tricycloheptane	poH@BUSHjyuUU@@
					setRingStrainAngles(atom, 3, 4, 151.40f * TO_RADIAN);	// calculated from MOPAC
					setRingStrainAngles(atom, 3, 8, 116.18f * TO_RADIAN);	// calculated from MOPAC
					setRingStrainAngles(atom, 4, 7, 129.00f * TO_RADIAN);	// calculated from MOPAC
					break;
				case 0x07070303:	// tricyclohexane	pfH@BUSLbkUUP@@
					setRingStrainAngles(atom, 3, 7, 120.88f * TO_RADIAN);	// calculated from MOPAC
					setRingStrainAngles(atom, 3, 3, 177.92f * TO_RADIAN);	// calculated from MOPAC
					break;
				case 0x07070404:	// tricycloheptane	poH@BURqZeuUU@@
					setRingStrainAngles(atom, 4, 7, 119.79f * TO_RADIAN);	// calculated from MOPAC
					setRingStrainAngles(atom, 4, 4, 146.20f * TO_RADIAN);	// calculated from MOPAC
					break;
				case 0x08080404:	// tricyclooctane	HaT@@DjfTkFKjjjh@@
					setRingStrainAngles(atom, 4, 8, 122.57f * TO_RADIAN);	// calculated from MOPAC
					setRingStrainAngles(atom, 4, 4, 134.76f * TO_RADIAN);	// calculated from MOPAC
					break;
					}
				}
				// acyclic atoms and non-aromatic 6-membered or larger rings
			else {
				float angle = (mMol.getAtomicNo(atom) > 10) ? 109.47f * (float)Math.PI/180.0f
							: (mMol.getAtomPi(atom) == 2)   ? (float)Math.PI
				            :  mMol.isFlatNitrogen(atom)    ? 120.00f * (float)Math.PI/180.0f
							: (mMol.getAtomPi(atom) == 0)   ? 109.47f * (float)Math.PI/180.0f
															: 120.00f * (float)Math.PI/180.0f;
				for (int i=1; i<connAtoms; i++)
					for (int j=0; j<i; j++)
						mBondAngle[atom][i][j] = angle;
				}
			}
		}

	// TODO set query feature allylic in TorsionID for flat non-pi nitrogens to distinguish them
	// from the tetrahedral case !!!

	/**
	 * Returns the preferred angle between to three atoms in a row as positive value <= pi.
	 * @param atom central atom
	 * @param conn1 one neighbour atom
	 * @param conn2 other neighbour atom
	 * @return
	 */
	public float getAngle(int atom, int conn1, int conn2) {
		int connIndex1 = -1;
		int connIndex2 = -1;
		for (int i=0; i<mMol.getAllConnAtoms(atom); i++) {
		    int conn = mMol.getConnAtom(atom, i);
		    if (conn == conn1 || conn == conn2) {
		        if (connIndex2 == -1) {
		            connIndex2 = i;
		        	}
		        else {
		            connIndex1 = i;
		        	break;
		    		}
		    	}
			}
		return mBondAngle[atom][connIndex1][connIndex2];
		}

	/**
	 * Returns the preferred angle between to three atoms in a row as positive value <= pi.
	 * @param atom central atom
	 * @param i1 index of neighbour atom
	 * @param i2 index of neighbour atom
	 * @return
	 */
	public float getConnAngle(int atom, int i1, int i2) {
		return (i1 < i2) ? mBondAngle[atom][i2][i1] : mBondAngle[atom][i1][i2];
		}

	private void calculateBondAnglesOfSmallRing(int ring) {
		RingCollection ringSet = mMol.getRingSet();
	    int ringSize = ringSet.getRingSize(ring);
		int[] ringBond = ringSet.getRingBonds(ring);
		boolean untouchedBondFound = false;
		for (int i=0; i<ringSize; i++) {
			if (mTinyRingSizeSum[ringBond[i]] == 0) {
				untouchedBondFound = true;
				break;
				}
			}

		if (untouchedBondFound) {
			float angle = (180.0f * ringSize - 360.0f) / ringSize * TO_RADIAN;
			int[] ringAtom = ringSet.getRingAtoms(ring);
			for (int i=0; i<ringSize; i++) {
			    setBondAngle(ringAtom[i], ringBond[(i==0)?ringSize-1:i-1], ringBond[i], angle);
				}

			if (ringSize <= 4) {
				for (int i=0; i<ringSize; i++)
					mTinyRingSizeSum[ringBond[i]] += ringSize;
				}
			}
		}

	private void calculateBondAnglesOfAromaticRing(int ring) {
		RingCollection ringSet = mMol.getRingSet();
	    int ringSize = ringSet.getRingSize(ring);
		int[] ringAtom = ringSet.getRingAtoms(ring);
		int[] ringBond = ringSet.getRingBonds(ring);
	    boolean isRegularRing = true;
	    for (int i=1; i<ringSize; i++) {
	        if (mBondLengthSet.getLength(ringBond[i]) != mBondLengthSet.getLength(ringBond[0])) {
	            isRegularRing = false;
	            break;
	        	}
	    	}
	    if (isRegularRing) {
	        float bondAngle = (ringSize - 2.0f) * (float)Math.PI / ringSize;
			for (int i=0; i<ringSize; i++)
			    setBondAngle(ringAtom[i], ringBond[i], ringBond[(i==0)?ringSize-1:i-1], bondAngle);
		    return;
			}
	    
	    float[] optAngle = new float[ringSize];
	    float angleSum = 0.0f;
		for (int i=0; i<ringSize; i++) {
		    optAngle[i] = (mMol.getAtomPi(ringAtom[i]) == 0) ? 109.47f * TO_RADIAN
		            	: (mMol.getAtomPi(ringAtom[i]) == 1) ? 0.6667f * (float)Math.PI : (float)Math.PI;
		    angleSum += optAngle[i];
			}
		float angleInc = ((ringSize - 2.0f) * (float)Math.PI - angleSum) / ringSize;
		for (int i=0; i<ringSize; i++)
		    optAngle[i] += angleInc;
		float[] direction = new float[ringSize];
		for (int i=1; i<ringSize; i++)
		    direction[i] = direction[i-1] + (float)Math.PI - optAngle[i];

		float[] dError_ddir = new float[ringSize];
		final int cycles = 100;
		for (int cycle=0; cycle<cycles; cycle++) {
		    float sx = 0.0f;
		    float sy = 0.0f;
			for (int i=0; i<ringSize; i++) {
			    sx += mBondLengthSet.getLength(ringBond[i]) * Math.sin(direction[i]);
			    sy += mBondLengthSet.getLength(ringBond[i]) * Math.cos(direction[i]);
				}
			float gapDirection = (float)Molecule.getAngle(0.0f, 0.0f, sx, sy);
			float gapSize = (float)Math.sqrt(sx*sx+sy*sy);

			int maxForceIndex = -1;
			float maxForce = 0;
			for (int i=0; i<ringSize; i++) {
			    int im1 = (i==0) ? ringSize-1 : i-1;
			    int ip1 = (i+1==ringSize) ? 0 : i+1;
			    float dirDif1 = (float)Molecule.getAngleDif(direction[i], direction[im1]);
			    float dirDif2 = (float)Molecule.getAngleDif(direction[ip1], direction[i]);
			    float optDif = (float)Molecule.getAngleDif(optAngle[i], optAngle[ip1]);
			    dError_ddir[i] = 2*dirDif1-2*dirDif2+2*optDif;

			    float gapGain = (float)Math.cos(direction[i] - (float)Math.PI/2 - gapDirection);
						// a positive force calls for a increase of the direction
			    float force = gapSize * gapGain - 0.03f * dError_ddir[i];

				if (Math.abs(force) > Math.abs(maxForce)) {
				    maxForce = force;
					maxForceIndex = i;
					}
				}
			float factor = (float)Math.exp(-5*(float)cycle/cycles);
			direction[maxForceIndex] += factor * maxForce;
			}

		for (int i=0; i<ringSize; i++) {
	        int im1 = (i == 0) ? ringSize-1 : i-1;
	        float angle = direction[im1] + (float)Math.PI - direction[i];
	        if (angle > 2*Math.PI)
	            angle -= 2*Math.PI;
		    setBondAngle(ringAtom[i], ringBond[im1], ringBond[i], angle);
			}

/*	    double sx = 0.0;
	    double sy = 0.0;
	    System.out.println("------------------");
		for (int i=0; i<ringSize; i++) {
	        int im1 = (i == 0) ? ringSize-1 : i-1;
	        double angle = direction[im1] + Math.PI - direction[i];
	        if (angle > 2*Math.PI)
	            angle -= 2*Math.PI;

	        System.out.println("optAngle:"+optAngle[i]+"; angle("+i+"):"+com.actelion.research.util.DoubleFormat.toString(angle));
		    sx += mBondLengthSet.getLength(ringBond[i]) * Math.sin(direction[i]);
		    sy += mBondLengthSet.getLength(ringBond[i]) * Math.cos(direction[i]);
			}
		System.out.println("dsx:"+com.actelion.research.util.DoubleFormat.toString(sx));
		System.out.println("dsy:"+com.actelion.research.util.DoubleFormat.toString(sy));
*/		}

	private int getRingStrainClass(int atom) {
		// the ring strain class of an atom is the sorted list
		// of non-zero tinyRingSizeSums of all attached bonds
		boolean[] handled = new boolean[mMol.getConnAtoms(atom)];
		int strainClass = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int largestSum = 0;
			int largestIndex = -1;
			for (int j=0; j<mMol.getConnAtoms(atom); j++) {
				if (!handled[j]) {
					int bond = mMol.getConnBond(atom, j);
					if (largestSum < mTinyRingSizeSum[bond]) {
						largestSum = mTinyRingSizeSum[bond];
						largestIndex = j;
						}
					}
				}
			if (largestSum == 0)
				return strainClass;

			strainClass <<= 8;
			strainClass += largestSum;
			handled[largestIndex] = true;
			}

		return strainClass;
		}

	private void setRingStrainAngles(int atom, int tinyRingSizeSum1, int tinyRingSizeSum2, float angle) {
		int allConnAtoms = mMol.getAllConnAtoms(atom);
		int connAtoms = mMol.getConnAtoms(atom);
		for (int i=1; i<allConnAtoms; i++) {
			int bondSum1 = (i<connAtoms) ? mTinyRingSizeSum[mMol.getConnBond(atom, i)] : 0;
			for (int j=0; j<i; j++) {
				if (mBondAngle[atom][i][j] == 0.0) {
					int bondSum2 = (j<connAtoms) ? mTinyRingSizeSum[mMol.getConnBond(atom, j)] : 0;
					if ((bondSum1 == tinyRingSizeSum1 && bondSum2 == tinyRingSizeSum2)
					 || (bondSum1 == tinyRingSizeSum2 && bondSum2 == tinyRingSizeSum1))
						mBondAngle[atom][i][j] = angle;
					}
				}
			}
		}

	private void setBondAngle(int atom, int bond1, int bond2, float angle) {
	    int conn1 = -1;
	    int conn2 = -1;
	    for (int i=0; i<mMol.getAllConnAtoms(atom); i++) {
	        int connBond = mMol.getConnBond(atom, i);
	        if (connBond == bond1 || connBond == bond2) {
	            if (conn1 == -1) {
	                conn1 = i;
	            	}
	            else {
	                conn2 = i;
	                break;
	            	}
	        	}
	        }
        if (mBondAngle[atom][conn2][conn1] == 0.0) {
    	    mBondAngle[atom][conn2][conn1] = angle;
			mDefinedAngleSum[atom] += angle;
			mDefinedAngleCount[atom]++;
	    	}
		}

	private float calculateRemainingTetrahedralAngle(float firstAngle) {
		float a109 = 109.47f * (float)Math.PI/180.0f;
		return a109 + (a109 - firstAngle) * 0.18f;  // this is an aproximation
		}
	}
