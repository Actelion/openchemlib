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

package com.actelion.research.chem.reaction;

import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;

import java.util.TreeMap;


public class FunctionalGroupClassifier {
	private StereoMolecule mMol;
	private TreeMap<Integer,Integer> mFunctionalGroupCountMap;

	public FunctionalGroupClassifier(StereoMolecule mol) {
		mMol = mol;
		mMol.ensureHelperArrays(Molecule.cHelperParities);
		}

	/**
	 * Applying a predefined dictionary of logically arranged 1024 functional groups
	 * this method determines all functional groups present in the molecule, counts
	 * how often they occurr in the molecule. The first 512 function group IDs refer
	 * to metal-related functions, while the second half of the IDs refer to traditional
	 * organic functional group fragments. These fragments are logically arranged in a tree,
	 * such that functional groups, whose ID only differs in the lowest significant bits,
	 * are structurally related. Functional groups as usually perceived by the chemist are
	 * often described by multiple overlapping fragments.
	 * @return counts mapped to existing functional group IDs
	 */
	public TreeMap<Integer,Integer> getFunctionGroupCountMap() {
		if (mFunctionalGroupCountMap == null)
			classifyFunctionalGroups();

		return mFunctionalGroupCountMap;
		}

	/**
	 * @return int[n][2] n by ID sorted pairs of functional group ID and occurance counts
	 */
	public int[][] getOrganicFunctionalGroupCounts() {
		if (mFunctionalGroupCountMap == null)
			classifyFunctionalGroups();

		int[][] v = new int[mFunctionalGroupCountMap.size()][2];
		int index = 0;
		for (Integer key:mFunctionalGroupCountMap.keySet()) {
			v[index][0] = key;
			v[index][1] = mFunctionalGroupCountMap.get(key).byteValue();
			index++;
			}

		return v;
		}

	/**
	 * Determines the similarity of two functional groups based on their
	 * location in the binary functional group similarity tree.
	 * If the functional belong into the same class (e.g. carbonyl), then the number
	 * of steps upwards in the tree is determined where both groups find a common
	 * still meaningful node. If the functional groups are identical, then 0 is returned.
	 * If they belong into unrelated classes, -1 is returned.
	 * @return -1:unrelated; 0:equal; 1...7:higher value means lower similarity
	 */
	public static int getFunctionalGroupEquivalenceLevel(int fg1, int fg2) {
		if (fg1 == fg2)
			return 0;

		int mask = 1;
		for (int i=1; i<8; i++) {
			if (nodeExists(fg1, i) || nodeExists(fg2, i))
				return -1;

			fg1 |= mask;
			fg2 |= mask;
			if (fg1 == fg2)
				return i;

			mask *= 2;
			}

		return -1;
		}

	private static boolean nodeExists(int fgID, int level) {
		return ClassificationData.getInstance().getEFGName(fgID, 8-level) != null;
		}

	private void classifyFunctionalGroups() {
		mFunctionalGroupCountMap = new TreeMap<Integer,Integer>();

		// mark all atoms that are non-carbon, have pi electrons or are stereo centers
		for (int atm=0; atm<mMol.getAtoms(); atm++)
			mMol.setAtomMarker(atm, mMol.getAtomicNo( atm ) != 6
								 || getAtomPi(atm) != 0
								 || mMol.isAtomStereoCenter(atm));

		// mark all atoms in pure-carbon 3-membered rings
		for (int i=0; i<mMol.getRingSet().getSize(); i++) {
			if (mMol.getRingSet().getRingSize(i) == 3) {
				int[] ringAtom = mMol.getRingSet().getRingAtoms(i);
				boolean heteroFound = false;
				for (int j=0; j<3; j++)
					if (mMol.getAtomicNo(ringAtom[j]) != 6)
						heteroFound = true;
				if (!heteroFound);
					for (int j=0; j<3; j++)
						mMol.setAtomMarker(ringAtom[j], true);
				}
			}
			
		for (int atm=0; atm<mMol.getAtoms(); atm++) {				// metal
			if (mMol.isMarkedAtom(atm)
			 && mMol.isMetalAtom(atm))
				classMet(atm);
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {				// Boron
			if (mMol.isMarkedAtom(atm)
			 && mMol.getAtomicNo(atm) == 5)
				storeEClass(classB(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {				// Silicon
			if (mMol.isMarkedAtom(atm)
			 && mMol.getAtomicNo(atm) == 14)
				storeEClass(classSi(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {					// C=X
			if (mMol.isMarkedAtom(atm)
			 && mMol.isElectronegative(atm)) {
				if (mMol.getConnAtoms( atm ) == 1) {
					if (mMol.getConnBondOrder( atm, 0 ) != 2) continue;
					int connAtm = mMol.getConnAtom( atm, 0 );
					if (mMol.getAtomicNo( connAtm ) != 6) continue;
					storeEClass(classCX(atm,connAtm));
					}
				else {
					if (mMol.getAtomicNo( atm ) != 7) continue;
					for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
						if (mMol.getConnBondOrder( atm, i ) != 2) continue;
						int connAtm = mMol.getConnAtom( atm, i );
						if (mMol.getAtomicNo( connAtm ) != 6) continue;
						storeEClass(classCX(atm,connAtm));
						}
					}
				}
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 33) continue;				// Arsenium
			if (mMol.isMarkedAtom(atm))
				storeEClass(classAs(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 52) continue;				// Tellurium
			if (mMol.isMarkedAtom(atm))
				storeEClass(classTe(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 34) continue;				// Selenium
			if (mMol.isMarkedAtom(atm))
				storeEClass(classSe(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 15) continue;				// Phosphorous
			if (mMol.isMarkedAtom(atm))
				storeEClass(classP(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 16) continue;				// Sulfur
			if (mMol.isMarkedAtom(atm))
				storeEClass(classS(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 53) continue;				// Iodine
			if (mMol.isMarkedAtom(atm))
				storeEClass(classI(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 35) continue;				// Bromine
			if (mMol.isMarkedAtom(atm))
				storeEClass(classBr(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 17) continue;				// Chlorine
			if (mMol.isMarkedAtom(atm))
				storeEClass(classCl(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 9) continue;				// Fluorine
			if (mMol.isMarkedAtom(atm))
				storeEClass(classF(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 7) continue;				// Nitrogen
			if (mMol.isMarkedAtom(atm))
				storeEClass(classN(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 8) continue;				// Oxygen
			if (mMol.isMarkedAtom(atm))
				storeEClass(classO(atm));
			}
		
		for (int atm=0; atm<mMol.getAtoms(); atm++) {
			if (mMol.getAtomicNo( atm ) != 6) continue;				// Carbon
			if (mMol.isMarkedAtom(atm))
				storeEClass(classC(atm));
			}
		}
		
		
	private void classMet(int atm) {
		int[] carbon = new int[ExtendedMolecule.cMaxConnAtoms];
		
		int nrofCarbs = getCarbons(atm, carbon);
		
		mMol.setAtomMarker(atm, false);
		
		if (nrofCarbs > 1) {
			for (int i=0; i<nrofCarbs; i++)
				classCMet(carbon[i],atm);
			}
		else if (nrofCarbs == 1)
			{
			classCMet(carbon[0], atm);
			}
		}
		
	private void classCMet(int carbon, int metal) {
		final byte[] metClass = { -2,
		 -2,																 -2,
		  0,  8,										 -2, -2, -2, -2, -2, -2,
		  4,  8,										 16, -2, -2, -2, -2, -2,
		  4, 12, 40, 48, 52, 52, 56, 56, 56, 56, 32, 36, 20, 24, -2, -2, -2, -2,
		  4, 12, 40, 48, 52, 52, -2, 60, 60, 60, 32, 36, 20, 24, 28, -2, -2, -2,
		  4, 12, 40,	 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
					 48, 52, 52, 60, 60, 60, 60, 32, 36, 20, 24, 28, -2, -2, -2,
		 -2, -2, 40,	 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44 };
		
		int eClass = metClass[mMol.getAtomicNo(metal)];
		if (eClass == -2) { storeEClass(-2); return; }   // shouldn't happen
		if (getAtomZ( carbon ) != 0) eClass += 2;
		if (mMol.isAllylicAtom(carbon)) eClass += 1;
		storeEClass(eClass);
		}

	private int classB(int atm) {
		int[] hetAtm = new int[ExtendedMolecule.cMaxConnAtoms];
		
		mMol.setAtomMarker(atm, false);
		
		if (getAtomPi(atm) != 0) return -2;
		
		if (mMol.getAtomCharge( atm ) == -1) {
			int nrofHets = getHeteros(atm,hetAtm,-1);
			if (nrofHets != getAtomZ( atm )) return -2;
		
			switch (nrofHets) {
			case 4:
				return 0x0047;											// B(-)X4
			case 3:
				if (getAtomSigma( atm ) != 0)
					return 0x0046;									// C-B(-)X3
				else
					return 0x0040;									// BH(-)X3
			case 2:
				if (getAtomSigma( atm ) == 2)
					return 0x0045;								// (C-)2B(-)X2
				else
					return 0x0040;									// R-BH(-)X2
			case 1:
				if (getAtomSigma( atm ) == 3)
					return 0x0044;									// (C-)3B(-)X
				else
					return 0x0041;									// R2-BH(-)X
			default:
				if (getAtomSigma( atm ) == 4)
					return 0x0043;									// (C-)4B(-)
				else
					return 0x0042;									// R3BH(-)
				}
			}
		
		if (mMol.getAtomCharge( atm ) != 0) return -2;
		
		int nrofHets = getHeteros(atm,hetAtm,-1);
		if (nrofHets != getAtomZ( atm )) return 0x0057;					// B=X
		
		if (mMol.getConnAtoms( atm ) - getAtomZ( atm ) - getAtomSigma( atm ) != 0)
			return 0x004F;											// Met an B
		
		switch (nrofHets) {
		case 3:
			if (mMol.getAtomicNo( hetAtm[0] ) == mMol.getAtomicNo( hetAtm[1] )
			 && mMol.getAtomicNo( hetAtm[1] ) == mMol.getAtomicNo( hetAtm[2] )) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 7)
					return 0x0054;									// B(NR2)3
				if (mMol.getAtomicNo( hetAtm[0] ) == 8)
					return 0x0055;										// B(OR)3
				}
			return 0x0056;											// rare BX3
		case 2:
			if (mMol.getConnAtoms( atm ) < 3)
				{
				if (mMol.getAtomicNo( hetAtm[0] )
				 == mMol.getAtomicNo( hetAtm[1] ))
					{
					if (mMol.getAtomicNo( hetAtm[0] ) == 7)
						return 0x004C;								// BH(NR2)2
					if (mMol.getAtomicNo( hetAtm[0] ) == 8)
						return 0x004D;								// BH(OR)2
					}
				return 0x004E;										// rare BHX2
				}

			int connAtm = -1;
			for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
				connAtm = mMol.getConnAtom( atm, i );
				if (!mMol.isElectronegative( connAtm )) break;
				}
		
			if (getAtomSigma( atm ) != 0) {
				if (mMol.getAtomicNo( hetAtm[0] )
				 == mMol.getAtomicNo( hetAtm[1] )) {
					switch (mMol.getAtomicNo( hetAtm[0] )) {
					case  7: return 0x005C;							// C-B(NR2)2
					case  8:
						if (mMol.isAromaticAtom( connAtm ))
							return 0x0058;							// Ar-B(OR)2
						if (getAtomPi(connAtm) != 0)
							return 0x0059;						// C(pi>0)-B(OR)2
						if (mMol.isAllylicAtom( connAtm ))
							return 0x005A;						// C=C-C-B(OR)2
						else
							return 0x005B;					// C(aliph)-B(OR)2
						}
					}
				return 0x005D;										// rare C-BX2
				}
		
			return -2;  // Met-B-Het2
		case 1:
			if (mMol.getConnAtoms( atm ) == 1)
				return 0x004B;											// BH2X
			if (mMol.getConnAtoms( atm ) == 2)
				return 0x004A;											// C-BHX
			if (mMol.getConnAtoms( atm ) == 3)
				{
				switch (mMol.getAtomicNo( hetAtm[0] ))
					{
				case  7: return 0x0050;								// (C-)2B-N
				case  8: return 0x0051;								// (C-)2B-O
				case 16: return 0x0052;								// (C-)2B-S
				default: return 0x0053;						// (C-)2B-X  X!=N,O,S
					}
				}
		
		default:
			switch (getAtomSigma( atm ))
				{
			case  3: return 0x005E;									// (C-)3B
			case  2: return 0x0048;									// (C-)2BH
			default: return 0x0049;										// C-BH2
				}
			}
		}

	private int classSi(int atm) {
		int[] hetAtm = new int[ExtendedMolecule.cMaxConnAtoms];

		mMol.setAtomMarker(atm, false);
		
		if (getAtomPi(atm) != 0) return -2;

		boolean alkinyl = false;
		boolean aryl = false;
		boolean vinyl = false;
		boolean acyl = false;
		boolean allyl = false;
		
		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
			int connAtm = mMol.getConnAtom( atm, i );
			if (mMol.getAtomicNo( connAtm ) != 6) continue;
			if (getAtomPi(connAtm) == 2)	{ alkinyl = true; continue; }
			if (mMol.isAromaticAtom( connAtm ))		{ aryl = true; continue; }
			if (getAtomPi(connAtm) == 1)	{ vinyl = true; continue; }
			if (getAtomZ( connAtm ) > 1) {
				for (int j=0; j<mMol.getConnAtoms( connAtm ); j++) {
					int nextConn = mMol.getConnAtom( connAtm, j );
					if (mMol.getAtomicNo( nextConn ) == 8)
						if (mMol.getConnBondOrder( connAtm, j ) == 2) acyl = true;
					}
				}
			if (acyl) continue;
			if (mMol.isAllylicAtom( connAtm )) {   allyl = true; continue; }
			}
		
		if (getAtomZ( atm ) == 0) {
			switch (getAtomSigma( atm )) {
			case 4:
				if (acyl)	return 0x00B4;						// acylsilane
				if (alkinyl) return 0x00B0;						// alkinylsilane
				if (vinyl)   return 0x00B1;						// alkenylsilane
				if (allyl)   return 0x00B3;						// allylsilane
				if (aryl)	return 0x00B2;						// arylsilane
				return 0x00B5;								// tetraalkylsilane
			case 3:
				return 0x00A0;										// (C-)3SiH
			case 2:
				return 0x00A1;										// (C-)2SiH2
			default:
				return 0x00A2;											// C-SiH3
				}
			}
		
		if (getAtomZ( atm ) == 1) {
			int betaCs = 0;			// count beta carbons as measure of steric hindrance
			for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
				int connAtm = mMol.getConnAtom( atm, i );
				if (mMol.isElectronegative( connAtm )) {
					hetAtm[0] = connAtm;
					continue;
					}
				betaCs += getAtomSigma( connAtm );
				}
		
			switch (getAtomSigma( atm )) {
			case 3:
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {					// O-SiR3
					switch (betaCs) {
					case 0: return 0x0088;								// O-TMS
					case 1: return 0x0089;							// O-SiEtMe2
					case 2: return 0x008A;							// O-SiiPrMe2
					case 3: return 0x008B;							// O-SitBuMe2
					case 4: return 0x008C;							// O-SiMeiPr2
					case 5: return 0x008D;							// O-SiEtiPr2
					case 6: return 0x008E;							// O-SiMetBu2
					default: return 0x008F;						// O-SitBuPh2
						}
					}
				if (mMol.getAtomicNo( hetAtm[0] ) == 7) {					// N-SiR3
					switch (betaCs) {
					case 0: return 0x0080;								// N-TMS
					case 1: return 0x0081;							// N-SiEtMe2
					case 2: return 0x0082;							// N-SiiPrMe2
					case 3: return 0x0083;							// N-SitBuMe2
					case 4: return 0x0084;							// N-SiMeiPr2
					case 5: return 0x0085;							// N-SiEtiPr2
					case 6: return 0x0086;							// N-SiMetBu2
					default: return 0x0087;						// N-SitBuPh2
						}
					}
				switch (mMol.getAtomicNo( hetAtm[0] )) {
				case 15: return 0x0098;								// (C-)3Si-P
				case 16: return 0x0099;								// (C-)3Si-S
				default: return 0x009A;					// (C-)3Si-X  X!=O,N,P,S
					}
			case 2:
				if (mMol.getAtomicNo( hetAtm[0] ) == 8)
					return 0x00A4;									// (C-)2SiH-O
				else
					return 0x00A5;							// (C-)2SiH-X  X!=O
			default:
				if (mMol.getAtomicNo( hetAtm[0] ) == 8)
					return 0x00A6;									// C-SiH2-O
				else
					return 0x00A7;								// C-SiH2-X  X!=O
				}
			}
		
		if (getAtomZ( atm ) == 2) {
			if (getHeteros(atm,hetAtm,-1) != 2) return -2;
			
			switch (getAtomSigma( atm )) {
			case 2:
				if (mMol.getAtomicNo( hetAtm[0] )
				 == mMol.getAtomicNo( hetAtm[1] )) {
					switch (mMol.getAtomicNo( hetAtm[0] )) {
					case  7: return 0x0091;						// (C-)2Si(-N)2
					case  8: return 0x0090;						// (C-)2Si(-O)2
					case 16: return 0x0092;						// (C-)2Si(-S)2
						}
					}
				return 0x0093;									// rare (C-)2SiX2
			default:
				if (mMol.getAtomicNo( hetAtm[0] )
				 == mMol.getAtomicNo( hetAtm[1] )) {
					switch (mMol.getAtomicNo( hetAtm[0] )) {
					case  7: return 0x00A9;							// RSiH(-N)2
					case  8: return 0x00A8;							// RSiH(-O)2
					case 16: return 0x00AA;							// RSiH(-S)2
						}
					}
				return 0x00AB;									// rare RSiHX2
				}
			}
		
		if (getAtomZ( atm ) == 3) {
			if (getHeteros(atm,hetAtm,-1) != 3) return -2;
			
			switch (getAtomSigma( atm )) {
			case 1:
				if (mMol.getAtomicNo( hetAtm[0] ) == 8
				 && mMol.getAtomicNo( hetAtm[1] ) == 8
				 && mMol.getAtomicNo( hetAtm[2] ) == 8)
					return 0x0094;									// C-Si(-O)3
				else
					return 0x0095;								// rare C-SiX3
			default:
				if (mMol.getAtomicNo( hetAtm[0] ) == 8
				 && mMol.getAtomicNo( hetAtm[1] ) == 8
				 && mMol.getAtomicNo( hetAtm[2] ) == 8)
					return 0x00AC;									// SiH(-O)3
				else
					return 0x00AD;									// rare SiHX3
				}
			}
		
		if (getAtomZ( atm ) == 4) {
			if (getHeteros(atm,hetAtm,-1) != 4) return -2;
			if (mMol.getAtomicNo( hetAtm[0] ) == 8
			 && mMol.getAtomicNo( hetAtm[1] ) == 8
			 && mMol.getAtomicNo( hetAtm[2] ) == 8
			 && mMol.getAtomicNo( hetAtm[3] ) == 8)
				return 0x0096;										// Si(-O)4
			else
				return 0x0097;										// rare SiX4
			}
		
		return -1;
		}

	private int classCX(int atm, int carbonylC) {
		int[] carbon = new int[ExtendedMolecule.cMaxConnAtoms];
		int[] hetAtm = new int[ExtendedMolecule.cMaxConnAtoms];
		
		if (bondToMet(carbonylC)) return -1;
		
		if (mMol.getAtomicNo( atm ) == 7) {				// classify N-part in C=N
			if (mMol.isAromaticAtom( atm ))
				return -1;
		
			if (!bondToMet(atm)) {
				int nrofHets = getHeteros(atm,hetAtm,-1);
			
				if (mMol.getAtomCharge( atm ) == 1) {
					if (mMol.getConnAtoms( atm ) == 1) {
						storeEClass(0x0214);						// C=NH2(+)
						}
					else {
						switch (getAtomZ( atm )) {
						case 2:
							if (nrofHets == 1) {
								storeEClass(0x0204);					// C=N(+)=N(-) type
								}
							else if (mMol.getAtomicNo( hetAtm[0] ) == 8
							 && mMol.getAtomicNo( hetAtm[1] ) == 8) {
								for (int i=0; i<2; i++)
									if (mMol.getConnAtoms( hetAtm[i] ) == 1)
										mMol.setAtomMarker(hetAtm[i], false);
								storeEClass(0x0216);					// C=N(+)(-O)2 type
								}
							else {
								storeEClass(0x0217);						// C=N(+)X2  X!=O
								}
							break;
						case 1:
							if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
								if (mMol.getConnAtoms( hetAtm[0] ) == 1)
									mMol.setAtomMarker(hetAtm[0], false);
								storeEClass(0x020C);					// C=N(+)(-C)-O type
								}
							else {
								storeEClass(0x020D);						// rare C=N(+)(-R)-X
								}
							break;
						default:
							storeEClass(0x0215);						// C=N(+)R-C
							break;
							}
						}
					}
				else {
					if (mMol.getConnAtoms( atm ) == 1) {
						storeEClass(0x0250);						// C=NH
						}
					else if (nrofHets != 0) {
						if (hasDBondToHetero( hetAtm[0] )) {
							storeEClass(0x0259);					// C=N-SO2-? type
							}
						else if (mMol.getAtomicNo( hetAtm[0] ) == 7) {
							if (getAtomPi(hetAtm[0]) != 0)
								storeEClass(0x025C);				// C=N-N=C type
							else
								storeEClass(0x025D);				// C=N-NR2 type
							}
						else if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
							if (mMol.getConnAtoms( hetAtm[0] ) == 1) {
								mMol.setAtomMarker(hetAtm[0], false);
								storeEClass(0x025E);				// C=N-OH
								}
							else {
								storeEClass(0x025F);				// C=N-O-?  ?!=H
								}
							}
						else {
							storeEClass(0x025A);						// rare C=N-X
							}
						}
					else {
						int nextCarbon = getNextConn(atm,carbonylC);
						if (mMol.isAromaticAtom( nextCarbon ))
							storeEClass(0x0257);					// C=N-Ar
						else if (hasDBondToHetero( nextCarbon ))
							storeEClass(0x0258);					// C=N-C=X
						else if (getAtomPi(nextCarbon) != 0) {
							if (getAtomSigma( nextCarbon ) == 2)
								storeEClass(0x0254);				// C=N-C(-C)=C
							else
								storeEClass(0x0255);				// C=N-CH=C type
							}
						else {
							if (getAtomZ( nextCarbon ) > 1)
								storeEClass(0x0256);				// C=N-CXR2 type
							else if (getAtomSigma( nextCarbon ) > 1)
								storeEClass(0x0252);				// C=N-C(-C)2-R
							else
								storeEClass(0x0253);				// C=N-CH2-R
							}
						}
					}
				}
			}
		
		mMol.setAtomMarker(atm, false);
		
		if (getAtomPi(carbonylC) != 0) {			// C=C=X
			int nextCarbon = getNextConn(carbonylC,atm);
			mMol.setAtomMarker(carbonylC, false);
			mMol.setAtomMarker(nextCarbon, false);
			
			switch (mMol.getAtomicNo( atm ))
				{
			case 7:
				return 0x0218;									// C=C=NR type
			case 8:
				return 0x0219;									// C=C=O
			case 16:
				return 0x021A;									// C=C=S
			default:
				return 0x021B;									// C=C=X  X!=N,O,S
				}
			}
		
		int nrofHets = getHeteros(carbonylC,hetAtm,atm);
		
		switch (mMol.getAtomicNo( atm )) {
		case 7:
			switch(getAtomZ( carbonylC )) {
			case 4:
				if (nrofHets == 1) {
					if (mMol.getAtomicNo( hetAtm[0] ) == 7)
						return 0x021E;							// ?-N=C=N-?
					if (mMol.getAtomicNo( hetAtm[0] ) == 8
					 || mMol.getAtomicNo( hetAtm[0] ) == 16)
						return -1;  							// ?-N=C=O, ?-N=C=S
					return 0x021F;								// ?-N=C=X  X!=N,O,S
					}
		
				if (mMol.getAtomicNo( hetAtm[0] )
				  > mMol.getAtomicNo( hetAtm[1] )) {
					int temp = hetAtm[0];
					hetAtm[0] = hetAtm[1];
					hetAtm[1] = temp;
					}
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 7) {
					switch (mMol.getAtomicNo( hetAtm[1] )) {
					case  7: return 0x0230;						// N-C(=N-?)-N
					case  8: return 0x0231;						// N-C(=N-?)-O
					case 16: return 0x0232;						// N-C(=N-?)-S
					default: return 0x0233;						// N-C(=N-?)-X  X!=S,O,N
						}
					}
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
					switch (mMol.getAtomicNo( hetAtm[1] )) {
					case  8: return 0x0234;						// O-C(=N-?)-O
					case 16: return 0x0235;						// O-C(=N-?)-S
					default: return 0x0236;						// O-C(=N-?)-X  X!=S,O,N
						}
					}
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 16
				 && mMol.getAtomicNo( hetAtm[1] ) == 16)
					return 0x0238;								// S-C(=N-?)-S
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 16
				 || mMol.getAtomicNo( hetAtm[1] ) == 16)
					return 0x0239;								// S-C(=N-?)-X  X!=S,O,N
		
				return 0x023A;									// Y-C(=N-?)-X  X,Y!=S,O,N
			case 3:
				if (checkAnhydride(hetAtm[0], carbonylC)) {
					if (getAtomSigma( carbonylC ) != 0) {
						switch (mMol.getAtomicNo( hetAtm[0] )) {
						case  7: return 0x0220;					// C-C(=N-?)-N-?=X
						case  8: return 0x0221;					// C-C(=N-?)-O-?=X
						case 16: return 0x0222;					// C-C(=N-?)-S-?=X
						default: return 0x0223;					// C-C(=N-?)-Y-?=X  Y!=S,O,N
							}
						}
					else
						{
						switch (mMol.getAtomicNo( hetAtm[0] )) {
						case  7: return 0x0228;					// HC(=N-?)-N-?=X
						case  8: return 0x0229;					// HC(=N-?)-O-?=X
						case 16: return 0x022A;					// HC(=N-?)-S-?=X
						default: return 0x022B;					// ?-N=CH-X-?=Y  X!=S,O,N
							}
						}
					}
				else
					{
					if (getAtomSigma( carbonylC ) != 0) {
						switch (mMol.getAtomicNo( hetAtm[0] )) {
						case  7: return 0x0224;					// C-C(=N-?)-N  !anhydride Type
						case  8: return 0x0225;					// C-C(=N-?)-O  !anhydride Type
						case 16: return 0x0226;					// C-C(=N-?)-S  !anhydride Type
						default: return 0x0227;					// C-C(=N-?)-X  X!=S,O,N  !aT
							}
						}
					else
						{
						switch (mMol.getAtomicNo( hetAtm[0] )) {
						case  7: return 0x022C;					// HC(=N-?)-N  !anhydride Type
						case  8: return 0x022D;					// HC(=N-?)-O  !anhydride Type
						case 16: return 0x022E;					// HC(=N-?)-S  !anhydride Type
						default: return 0x022F;					// ?-N=CH-X  X!=S,O,N  !anhydr Type
							}
						}
					}
			default:
				if (getAtomSigma( carbonylC ) == 2) {
					if (mMol.isAllylicAtom( carbonylC ))
						return 0x023C;							// C=C-C(=N-?)-C
					else
						return 0x023D;							// C-C(=N-?)(!allyl)-C
					}
				if (getAtomSigma( carbonylC ) == 1) {
					if (mMol.isAllylicAtom( carbonylC ))
						return 0x023E;							// C=C-CH=N-?
					else
						return 0x023F;							// C-CH(!allyl)=N-?
					}
				return 0x01D4;									// H2C=N-?
				}
		case 8:
			switch(getAtomZ( carbonylC )) {
			case 4:
				if (nrofHets == 1) {
					mMol.setAtomMarker(hetAtm[0], false);
					return 0X021C;								// -N=C=O type
					}
		
				if (mMol.getAtomicNo( hetAtm[0] )
				  > mMol.getAtomicNo( hetAtm[1] )) {
					int temp = hetAtm[0];
					hetAtm[0] = hetAtm[1];
					hetAtm[1] = temp;
					}
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 7) {
					switch (mMol.getAtomicNo( hetAtm[1] )) {
					case  7: return 0x0280;						// N-C(=O)-N
					case  8: return 0x0281;						// N-C(=O)-O
					case 16: return 0x0282;						// N-C(=O)-S
					default: return 0x0283;						// N-C(=O)-X  X!=S,O,N
						}
					}
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
					switch (mMol.getAtomicNo( hetAtm[1] )) {
					case  8: return 0x0284;							// O-C(=O)-O
					case 16: return 0x0285;							// O-C(=O)-S
					default: return 0x0286;				// O-C(=O)-X  X!=S,O,N
						}
					}
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 16
				 && mMol.getAtomicNo( hetAtm[1] ) == 16)
					return 0x0288;									// S-C(=O)-S
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 16
				 || mMol.getAtomicNo( hetAtm[1] ) == 16)
					return 0x0289;							// S-C(=O)-X  X!=S,O,N
		
				return 0x028A;							// Y-C(=O)-X  X,Y!=S,O,N
			case 3:
				if (checkAnhydride(hetAtm[0],carbonylC)) {
					if (getAtomSigma( carbonylC ) != 0) {
						switch (mMol.getAtomicNo( hetAtm[0] )) {
						case  7: return 0x0290;					// C-C(=O)-N-?=X
						case  8: return 0x0291;					// C-C(=O)-O-?=X
						case 16: return 0x0292;					// C-C(=O)-S-?=X
						default: return 0x0293;		// C-C(=O)-Y-?=X  Y!=S,O,N
							}
						}
					else {
						switch (mMol.getAtomicNo( hetAtm[0] )) {
						case  7: return 0x0298;					// HC(=O)-N-?=X
						case  8: return 0x0299;					// HC(=O)-O-?=X
						case 16: return 0x029A;					// HC(=O)-S-?=X
						default: return 0x029B;			// O=CH-X-?=X  X!=S,O,N
							}
						}
					}
				else {
					if (getAtomSigma( carbonylC ) != 0) {
						switch (mMol.getAtomicNo( hetAtm[0] )) {
						case  7: return 0x0294;	// C-C(=O)-N  !anhydride Type
						case  8: return 0x0295;	// C-C(=O)-O  !anhydride Type
						case 16: return 0x0296;	// C-C(=O)-S  !anhydride Type
						default: return 0x0297;		// C-C(=O)-X  X!=S,O,N  !aT
							}
						}
					else {
						switch (mMol.getAtomicNo( hetAtm[0] )) {
						case  7: return 0x029C;		// HC(=O)-N  !anhydride Type
						case  8: return 0x029D;		// HC(=O)-O  !anhydride Type
						case 16: return 0x029E;		// HC(=O)-S  !anhydride Type
						default: return 0x029F;  // O=CH-X  X!=S,O,N  !anhydr Type
							}
						}
					}
			default:
				if (getAtomSigma( carbonylC ) == 2) {
					if (isQuinone(carbonylC)) return 0x00BC;		// quinone
					if (getCarbons(carbonylC,carbon) != 2) return -2;
					if (mMol.isAromaticAtom( carbon[0] )
					 && mMol.isAromaticAtom( carbon[1] ))
						return 0x02B0;							// Ar-C(=O)-Ar
					int nextCarbon = -1;
					if (mMol.isAromaticAtom( carbon[0] ))
						nextCarbon = carbon[1];
					if (mMol.isAromaticAtom( carbon[1] ))
						nextCarbon = carbon[0];
					if (nextCarbon != -1) {						// Ar-C(=O)-C
						if (getAtomPi(nextCarbon) != 0)
							return 0x02A0;						// Ar-C(=O)-C=C
						switch (getAtomZ( nextCarbon )) {
						case 3:
							return 0x02AC;						// X3C-C(=O)-Ar
						case 2:
							return 0x02AD;						// R-CX2-C(=O)-Ar
						case 1:
							if (getAtomSigma( nextCarbon ) == 3)
								return 0x02AE;				// (C-)2CX-C(=O)-Ar
							else
								return 0x02AF;				// R-CH(-X)-C(=O)-Ar
						default:
							if (getAtomSigma( nextCarbon ) == 4)
								return 0x02B2;				// (C-)3C-C(=O)-Ar
							else
								return 0x02B3;					// R2CH-C(=O)-Ar
							}
						}
					if (getAtomPi(carbon[0]) != 0
					 && getAtomPi(carbon[1]) != 0)
						return 0x02A1;							// C=C-C(=O)-C=C
					nextCarbon = -1;
					if (getAtomPi(carbon[0]) != 0) nextCarbon = carbon[1];
					if (getAtomPi(carbon[1]) != 0) nextCarbon = carbon[0];
					if (nextCarbon != -1) {						// C=C-C(=O)-C(pi=0)
						switch (getAtomZ( nextCarbon )) {
						case 3:
							return 0x02A4;						// X3C-C(=O)-C=C
						case 2:
							return 0x02A5;					// R-CX2-C(=O)-C=C
						case 1:
							if (getAtomSigma( nextCarbon ) == 3)
								return 0x02A6;				// (C-)2CX-C(=O)-C=C
							else
								return 0x02A7;				// R-CH(-X)-C(=O)-C=C
						default:
							if (getAtomSigma( nextCarbon ) == 4)
								return 0x02A2;				// (C-)3C-C(=O)-C=C
							else
								return 0x02A3;					// R2CH-C(=O)-C=C
							}
						}
					int allSigma = getAtomSigma( carbon[0] ) + getAtomSigma( carbon[1] );
					int allZ = getAtomZ( carbon[0] ) + getAtomZ( carbon[1] );
					if (allZ > 4)
						return 0x02A8;						// ?-CX2-C(=O)-CX3
					if (allZ > 2)
						return 0x02A9;					// R-CX2-C(=O)-CX2-R type
					if (allZ != 0) {
						if (allSigma + allZ > 6)
							return 0x02AA;			// C-CHX-C(=O)-CX(-C)2 type
						else
							return 0x02AB;			// R-CH(-X)-C(=O)-CH2-R type
						}
					if (allSigma > 6)
						return 0x02B4;				// (C-)3C-C(=O)-CH(-C)2 type
					else
						return 0x02B5;					// R-CH2-C(=O)-CH2-R type
					}
				if (getAtomSigma( carbonylC ) == 1) {
					int nextCarbon = getNextConn(carbonylC,atm);
					if (mMol.isAromaticAtom( nextCarbon ))
						return 0x01B0;								// Ar-CH=O
					if (getAtomPi(nextCarbon) == 2)
						return 0x01B3;								// C#C-CH=O
					if (getAtomPi(nextCarbon) != 0)
						return 0x01B2;								// C=C-CH=O
					switch (getAtomZ( nextCarbon ))
						{
					case 3:
						return 0x01B8;								// X3C-CH=O
					case 2:
						if (getAtomSigma( nextCarbon ) == 1)
							return 0x01BB;							// X2CH-CH=O
						else
							return 0x01BA;							// C-CX2-CH=O
					case 1:
						if (getAtomSigma( nextCarbon ) == 1)
							return 0x01BC;							// X-CH2-CH=O
						if (getAtomSigma( nextCarbon ) == 2)
							return 0x01BD;						// XCH(-C)-CH=O
						else
							return 0x01BE;						// XC(-C)2-CH=O
					default:
						if (getAtomSigma( nextCarbon ) < 3)
							return 0x01B7;							// RCH2-CH=O
						if (getAtomSigma( nextCarbon ) == 3)
							return 0x01B6;						// (C-)2CH-CH=O
						else
							return 0x01B4;						// (C-)3C-CH=O
						}
					}
				return 0x01D5;											// H2C=O
				}
		case 15:
			return 0x0066;												// C=PH
		case 16:
			switch(getAtomZ( carbonylC ))
				{
			case 4:
				if (nrofHets == 1) {
					mMol.setAtomMarker(hetAtm[0], false);
					return 0x021D;								// -N=C=S type
					}
		
				if (mMol.getAtomicNo( hetAtm[0] )
				  > mMol.getAtomicNo( hetAtm[1] )) {
					int temp = hetAtm[0];
					hetAtm[0] = hetAtm[1];
					hetAtm[1] = temp;
					}
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 7) {
					switch (mMol.getAtomicNo( hetAtm[1] )) {
					case  7: return 0x02C0;							// N-C(=S)-N
					case  8: return 0x02C1;							// N-C(=S)-O
					case 16: return 0x02C2;							// N-C(=S)-S
					default: return 0x02C3;				// N-C(=S)-X  X!=S,O,N
						}
					}
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
					switch (mMol.getAtomicNo( hetAtm[1] )) {
					case  8: return 0x02C4;							// O-C(=S)-O
					case 16: return 0x02C5;							// O-C(=S)-S
					default: return 0x02C6;				// O-C(=S)-X  X!=S,O,N
						}
					}
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 16
				 && mMol.getAtomicNo( hetAtm[1] ) == 16)
					return 0x02C8;									// S-C(=S)-S
		
				if (mMol.getAtomicNo( hetAtm[0] ) == 16
				 || mMol.getAtomicNo( hetAtm[1] ) == 16)
					return 0x02C9;							// S-C(=S)-X  X!=S,O,N
		
				return 0x02CA;							// Y-C(=S)-X  X,Y!=S,O,N
			case 3:
				if (getAtomSigma( carbonylC ) != 0) {
					switch (mMol.getAtomicNo( hetAtm[0] ))
						{
					case  7: return 0x02D0;							// C-C(=S)-N
					case  8: return 0x02D1;							// C-C(=S)-O
					case 16: return 0x02D2;							// C-C(=S)-S
					default: return 0x02D3;				// C-C(=S)-X  X!=S,O,N
						}
					}
				else {
					switch (mMol.getAtomicNo( hetAtm[0] )) {
					case  7: return 0x02D4;							// HC(=S)-N
					case  8: return 0x02D5;							// HC(=S)-O
					case 16: return 0x02D6;							// HC(=S)-S
					default: return 0x02D7;					// S=CH-X  X!=S,O,N
						}
					}
			default:
				if (getAtomSigma( carbonylC ) == 2) {
					if (mMol.isAllylicAtom( carbonylC ))
						return 0x02CC;							// C=C-C(=S)-C
					else
						return 0x02CD;					// C(pi=0)-C(=S)-C(pi=0)
					}
				if (getAtomSigma( carbonylC ) == 1) {
					if (mMol.isAllylicAtom( carbonylC ))
						return 0x02CE;								// C=C-CH=S
					else
						return 0x02CF;							// C(pi=0)-CH=S
					}
				return 0x01D6;											// H2C=S
				}
		case 34:
			switch(getAtomZ( carbonylC )) {
			case 4:
				return 0x019C;										// X-C(=Se)-X
			case 3:
				if (getAtomSigma( carbonylC ) != 0)
					return 0x019A;									// C-C(=Se)-X
				else
					return 0x019B;									// HC(=Se)-X
			default:
				if (getAtomSigma( carbonylC ) == 2)
					return 0x0198;									// C-C(=Se)-C
				else
					return 0x0199;									// R-CH(=Se)
				}
		default: return -2;
			}
		}

	private int classAs(int atm) {
		mMol.setAtomMarker(atm, false);
		
		if (getAtomSigma( atm ) == 0) {
			switch (getAtomZ( atm )) {
			case 3: return 0x0078;										// AsX3
			case 5: return 0x0079;										// AsX5
			default: return 0x007A;									// rare AsXn
				}
			}
		
		if (getAtomSigma( atm ) == 1) {
			switch (getAtomZ( atm )) {
			case 0: return 0x007C;								// C-AsH2 type
			case 2: return 0x007D;										// C-AsX2
			case 4: return 0x007E;										// C-AsX4
			default: return 0x007F;								// rare C-AsXn
				}
			}
		
		return 0x007B;													// C-As-C
		}

	private int classTe(int atm) {
		mMol.setAtomMarker(atm, false);
		
		if (getAtomSigma( atm ) == 0) {
			switch (getAtomZ( atm )) {
			case 2: return 0x00F8;										// TeX2
			case 4: return 0x00F9;										// TeX4
			case 6: return 0x00FA;										// TeX6
			default: return 0x00FB;									// rare TeXn
				}
			}
		
		if (getAtomSigma( atm ) == 1) {
			switch (getAtomZ( atm )) {
			case 1: return 0x00F4;										// C-Te-X
			case 3: return 0x00F5;										// C-TeX3
			case 5: return 0x00F6;										// C-TeX5
			default: return 0x00F7;								// rare C-TeXn
				}
			}
		
		if (getAtomSigma( atm ) == 2) {
			switch (getAtomZ( atm )) {
			case 0: return 0x00FC;										// C-Te-C
			case 2: return 0x00FD;									// C-TeX2-C
			case 4: return 0x00FE;									// C-TeX4-C
			default: return 0x00FF;								// rare C-TeXn-C
				}
			}
		
		return -1;
		}

	private int classSe(int atm) {
		int[] carbon = new int[ExtendedMolecule.cMaxConnAtoms];
		int[] hetAtm = new int[ExtendedMolecule.cMaxConnAtoms];

		mMol.setAtomMarker(atm, false);
		
		int nrofCarbs = getCarbons(atm,carbon);
		for (int i=0; i<nrofCarbs; i++)
			classCSe(carbon[i]);
		
		if (getAtomZ( atm ) == 0) {
			if (getAtomSigma( atm ) == 2)
				return 0x019E;											// C-Se-C
			else
				return 0x019F;											// C-SeH
			}

		int[] tripleBnds = new int[1];
		int[] doubleBnds = new int[1];
		int nrofHets = getSortedHeteros(atm,hetAtm,tripleBnds,doubleBnds,-1);
		if (tripleBnds[0] != 0) return -2;
		
		if (getAtomZ( atm ) == 1) {
			if (getAtomSigma( atm ) == 1) {
				switch (mMol.getAtomicNo( hetAtm[0] )) {
				case 8:
					if (mMol.getConnAtoms( hetAtm[0] ) == 1) {
						mMol.setAtomMarker(hetAtm[0], false);
						return 0x0180;									// C-SeOH
						}
					else
						return 0x0181;							// C-Se-O-?  ?!=O
				case 17:
					mMol.setAtomMarker(hetAtm[0], false);
					return 0x0182;									// C-Se-Cl
				case 34:
					if (mMol.getConnAtoms( hetAtm[0] ) == 1) {
						mMol.setAtomMarker(hetAtm[0], false);
						return 0x0184;								// C-SeSeH
						}
					else {
						if (getAtomZ( hetAtm[0] ) == 1)
							return 0x0185;							// C-SeSe-C
						else
							return 0x0186;						// C-SeSe(z>1)
						}
				default:
					return 0x0183;							// C-SeX  X!=O,Cl,Se
					}
				}
			else
				return 0x0187;											// X-SeH
			}
		
		if (getAtomZ( atm ) == 2) {
			if (nrofHets == 1
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && getAtomSigma( atm ) == 2) {
				mMol.setAtomMarker(hetAtm[0], false);
				return 0x0188;										// C-Se(=O)-C
				}
			else
				return 0x0189;									// rare Se(z=2)
			}
		
		if (getAtomZ( atm ) == 3) {
			if (nrofHets == 2
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && getAtomSigma( atm ) == 1) {
				mMol.setAtomMarker(hetAtm[0], false);
				if (mMol.getAtomicNo( hetAtm[1] ) == 8) {
					if (mMol.getConnAtoms( hetAtm[1] ) == 1) {
						mMol.setAtomMarker(hetAtm[1], false);
						return 0x0190;							// C-Se(=O)-OH
						}
					else
						return 0x0191;						// C-Se(=O)-O-?  ?!=H
					}
				else
					return 0x0192;							// C-Se(=O)-X  X!=O
				}
			else
				return 0x0193;										// rare SeX3
			}
		
		if (getAtomZ( atm ) == 4) {
			if (nrofHets == 2
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && mMol.getAtomicNo( hetAtm[1] ) == 8
			 && getAtomSigma( atm ) == 2) {
				mMol.setAtomMarker(hetAtm[0], false);
				mMol.setAtomMarker(hetAtm[1], false);
				return 0x018C;										// C-SeO2-C
				}
			if (nrofHets == 3
			&& mMol.getAtomicNo( hetAtm[0] ) == 8) {
				mMol.setAtomMarker(hetAtm[0], false);
				if (mMol.getAtomicNo( hetAtm[1] ) == 8
				 && mMol.getAtomicNo( hetAtm[2] ) == 8)
					return 0x018D;									// O-Se(=O)-O
				else
					return 0x018E;							// X-Se(=O)-X  X!=O
				}
			return 0x018F;											// rare SeX4
			}
		
		if (getAtomZ( atm ) == 5) {
			if (nrofHets == 3
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && mMol.getAtomicNo( hetAtm[1] ) == 8
			 && getAtomSigma( atm ) == 1) {
				mMol.setAtomMarker(hetAtm[0], false);
				mMol.setAtomMarker(hetAtm[1], false);
				if (mMol.getAtomicNo( hetAtm[2] ) == 8) {
					if (mMol.getConnAtoms( hetAtm[2] ) == 1) {
						mMol.setAtomMarker(hetAtm[1], false);
						return 0x0194;								// C-SeO2-OH
						}
					else
						return 0x0195;						// C-SeO2-O-?  ?!=H
					}
				else
					return 0x0196;								// C-SeO2-X  X!=O
				}
			else
				return 0x0197;										// rare SeX5
			}
		
		if (getAtomZ( atm ) == 6) {
			if (nrofHets == 4
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && mMol.getAtomicNo( hetAtm[1] ) == 8
			 && mMol.getAtomicNo( hetAtm[2] ) == 8
			 && mMol.getAtomicNo( hetAtm[3] ) == 8) {
				mMol.setAtomMarker(hetAtm[0], false);
				mMol.setAtomMarker(hetAtm[1], false);
				return 0x018A;										// O-SeO2-O
				}
			return 0x018B;											// rare SeX6
			}
		
		return -1;
		}

	private void classCSe(int atm) {
		if (mMol.isAromaticAtom( atm )) {
			storeEClass(0x01A0);										// Ar-Se
			return;
			}
		if (hasDBondToHetero( atm )) {
			storeEClass(0x01A8);										// X=C-Se
			return;
			}
		if (getAtomZ( atm ) == 1) {
			if (getAtomPi(atm) != 0) {
				if (getAtomSigma( atm ) == 1)
					storeEClass(0x01A2);						// C#C-Se / C=CH-Se
				else
					storeEClass(0x01A3);							// C=C(-C)-Se
				return;
				}
		
			if (mMol.isAllylicAtom( atm ))
				storeEClass(0x01A7);								// C=C-CR2-Se
		
			switch (getAtomSigma( atm )) {
			case 3:
				storeEClass(0x01A4);								// (C-)3C-Se
				return;
			case 2:
				storeEClass(0x01A5);								// (C-)2CH-Se
				return;
			default:
				storeEClass(0x01A6);									// R-CH2-Se
				return;
				}
			}
		if (getAtomZ( atm ) == 2)
			{
			if (getAtomPi(atm) != 0) {
				storeEClass(0x01A1);								// C=C(-X)-Se
				return;
				}
			switch (getAtomSigma( atm )) {
			case 2:
				storeEClass(0x01A9);							// (C-)2C(-X)-Se
				return;
			case 1:
				storeEClass(0x01AA);								// C-CH(-X)-Se
				return;
			default:
				storeEClass(0x01AB);									// X-CH2-Se
				return;
				}
			}
		if (getAtomZ( atm ) == 3) {
			if (getAtomSigma( atm ) != 0)
				storeEClass(0x01AC);									// C-CX2-Se
			else
				storeEClass(0x01AD);									// H-CX2-Se
			return;
			}
		storeEClass(0x01AE);											// CX3-Se
		}

	private int classP(int atm) {
		int[] hetAtm = new int[ExtendedMolecule.cMaxConnAtoms];
		
		mMol.setAtomMarker(atm, false);
		
		int phenyl = 0;
		int vinyl = 0;
		int alkyl = 0;

		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
			int connAtm = mMol.getConnAtom( atm, i );
			if (mMol.getAtomicNo( connAtm ) != 6) continue;
			if (mMol.isAromaticAtom( connAtm )) {
				phenyl++;
				continue;
				}
			if (getAtomPi(connAtm) != 0) {
				vinyl++;
				continue;
				}
			alkyl++;
			}
		
		if (getAtomZ( atm ) == 0) {
			if (mMol.getAtomCharge( atm ) == 1)
				return 0x0060;											// R4P(+)
			if (getAtomPi(atm) == 0) {
				if (getAtomSigma( atm ) == 3) {
					if (phenyl == 3) return 0x0068;						// Ar3P
					return 0x0069;							// (C-)3P except Ar3P
					}
				if (getAtomSigma( atm ) == 2)
					return 0x0062;									// (C-)2PH
				else
					return 0x0063;									// C-PH2 type
				}
			if (getAtomPi(atm) == 1) {
				if (getAtomSigma( atm ) == 4)
					return 0x0064;									// (C-)3P=C
				if (getAtomSigma( atm ) == 2)
					return 0x0065;										// C-P=C
				else
					return 0x0067;									// rare P=C
				}
			return -2;
			}
		
		if (mMol.getAtomCharge( atm ) == 1)
			return 0x0061;										// RnP(+)X(4-n)
		
		if (getAtomZ( atm ) == 1)
			return 0x006E;											// R2PX type
		
		int[] tripleBnds = new int[1];
		int[] doubleBnds = new int[1];
		int nrofHets = getSortedHeteros(atm,hetAtm,tripleBnds,doubleBnds,-1);
		if (tripleBnds[0] != 0) return -2;
		
		if (getAtomZ( atm ) == 2) {
			if (nrofHets == 1) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
					mMol.setAtomMarker(hetAtm[0], false);
					return 0x006C;									// R3P=O type
					}
				return 0x006D;										// R3P=X,RP=X
				}
			return 0x006F;											// R3PX2,RPX2
			}
		
		if (getAtomZ( atm ) == 3) {
			if (nrofHets == 2) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
					mMol.setAtomMarker(hetAtm[0], false);
					return 0x0070;									// R2P(=O)-X
					}
				return 0x0071;									// R2P(=X)-Y  X!=O
				}
			if (mMol.getAtomicNo( hetAtm[0] ) == 8
			 && mMol.getAtomicNo( hetAtm[1] ) == 8
			 && mMol.getAtomicNo( hetAtm[2] ) == 8)
				return 0x006A;										// P(-O-?)3
			else
				return 0x006B;							// PX3 except P(-O-?)3
			}
		
		if (getAtomZ( atm ) == 4) {
			if (nrofHets == 3) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
					mMol.setAtomMarker(hetAtm[0], false);
					if (mMol.getAtomicNo( hetAtm[1] ) == 8
					 && mMol.getAtomicNo( hetAtm[2] ) == 8)
						return 0x0072;							// RP(=O)(OR)2
					else
						return 0x0073;							// RP(=O)X2  X!=O
					}
				return 0x0074;									// RP(=X)Y2  X!=O
				}
			return -2;
			}
		
		if (getAtomZ( atm ) == 5) {
			if (nrofHets == 4) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
					mMol.setAtomMarker(hetAtm[0], false);
					if (mMol.getAtomicNo( hetAtm[1] ) == 8
					 && mMol.getAtomicNo( hetAtm[2] ) == 8
					 && mMol.getAtomicNo( hetAtm[3] ) == 8)
						return 0x0076;								// O=P(OR)3
					else
						return 0x0077;								// O=PX3  X!=O
					}
				return 0x0075;										// X=PY3  X!=O
				}
			return -2;
			}
		return -2;
		}

	private int classS(int atm) {
		int[] carbon = new int[ExtendedMolecule.cMaxConnAtoms];
		int[] hetAtm = new int[ExtendedMolecule.cMaxConnAtoms];
		
		mMol.setAtomMarker(atm, false);
		
		if (!mMol.isAromaticAtom( atm )) {
			int nrofCarbs = getCarbons(atm,carbon);
			for (int i=0; i<nrofCarbs; i++)
				classCS(carbon[i],atm);
			}
		
		if (mMol.isAromaticAtom( atm )) return 0x00BB;	// S(arom)
		
		if (mMol.getAtomCharge( atm ) == 1)
			return 0x02DF;												// S(+)
		
		if (getAtomZ( atm ) == 0) {
			if (getAtomSigma( atm ) == 2)
				return 0x02DC;											// C-S-C
			else
				return 0x02DD;											// C-SH
			}
		
		int[] tripleBnds = new int[1];
		int[] doubleBnds = new int[1];
		int nrofHets = getSortedHeteros(atm,hetAtm,tripleBnds,doubleBnds,-1);
		if (tripleBnds[0] != 0) return -2;
		
		if (getAtomZ( atm ) == 1) {
			if (getAtomSigma( atm ) == 1) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 16) {
					if (mMol.getConnAtoms( hetAtm[0] ) == 1) {
						mMol.setAtomMarker(hetAtm[0], false);
						return 0x02D8;									// C-S-SH
						}
					else {
						if (getAtomZ( hetAtm[0] ) == 1)
							return 0x02D9;							// C-S-S-C
						else
							return 0x02DA;							// C-S-S(z>1)
						}
					}
				else
					return 0x02DB;									// C-S-X  X!=S
				}
			else {
				if (mMol.getConnAtoms( atm ) > 1)
					return -1;					// forget X-S-Met
				return 0x02DE;											// X-SH
				}
			}
		
		if (getAtomZ( atm ) == 2) {
			if (nrofHets == 1
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && getAtomSigma( atm ) == 2) {
				mMol.setAtomMarker(hetAtm[0], false);
				return 0x02E0;										// C-S(=O)-C
				}
			else
				return 0x02E1;									// rare S(z=2)
			}
		
		if (getAtomZ( atm ) == 3) {
			if (nrofHets == 2
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && getAtomSigma( atm ) == 1) {
				mMol.setAtomMarker(hetAtm[0], false);
				if (mMol.getAtomicNo( hetAtm[1] ) == 8)
					return 0x02E8;									// C-S(=O)-O
				else
					return 0x02E9;								// C-S(=O)-X  X!=O
				}
			else
				return 0x02EA;									// rare S(z=3)
			}
		
		if (getAtomZ( atm ) == 4) {
			if (nrofHets == 2
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && mMol.getAtomicNo( hetAtm[1] ) == 8
			 && getAtomSigma( atm ) == 2) {
				mMol.setAtomMarker(hetAtm[0], false);
				mMol.setAtomMarker(hetAtm[1], false);
				return 0x02E4;										// C-SO2-C
				}
			if (nrofHets == 3
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && getAtomSigma( atm ) == 0) {
				mMol.setAtomMarker(hetAtm[0], false);
				if (mMol.getAtomicNo( hetAtm[1] ) == 8
				 && mMol.getAtomicNo( hetAtm[2] ) == 8)
					return 0x02E5;									// O-S(=O)-O
				else
					return 0x02E6;								// X-S(=O)-X  X!=O
				}
			return 0x02E7;										// rare S(z=4)
			}
		
		if (getAtomZ( atm ) == 5) {
			if (nrofHets == 3
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && mMol.getAtomicNo( hetAtm[1] ) == 8
			 && getAtomSigma( atm ) == 1) {
				mMol.setAtomMarker(hetAtm[0], false);
				mMol.setAtomMarker(hetAtm[1], false);
				if (mMol.getAtomicNo( hetAtm[2] ) == 8)
					return 0x02EC;									// C-SO2-O
				else
					return 0x02ED;								// rare C-SO2-X
				}
			else
				return 0x02EE;									// rare S(z=5)
			}
		
		if (getAtomZ( atm ) == 6) {
			if (nrofHets == 4
			 && mMol.getAtomicNo( hetAtm[0] ) == 8
			 && mMol.getAtomicNo( hetAtm[1] ) == 8) {
				mMol.setAtomMarker(hetAtm[0], false);
				mMol.setAtomMarker(hetAtm[1], false);
				return 0x02E2;										// X-SO2-X
				}
			return 0x02E3;										// rare S(z=6)
			}
		
		return -1;
		}

	private void classCS(int carbon, int sulfur) {
		int[] hetAtm = new int[ExtendedMolecule.cMaxConnAtoms];
		
		if (bondToMet(carbon)) return;
		
		if (mMol.isAromaticAtom( carbon )) {
			storeEClass(0x01C0);									// C(arom)-S
			return;
			}
		
		if (getAtomPi(carbon) == 2) {
			storeEClass(0x01C1);										// C#C-S
			return;
			}
		
		int[] tripleBnds = new int[1];
		int[] doubleBnds = new int[1];
		int nrofHets = getSortedHeteros(carbon,hetAtm, tripleBnds,doubleBnds,sulfur);
		
		if (tripleBnds[0] != 0) {
			storeEClass(0x01D0);										// X#C-S
			return;
			}
		
		if (hasDBondToHetero( carbon )) {
			storeEClass(0x01D1);										// X=C-S
			return;
			}
		
		if (getAtomPi(carbon) == 1) {
			if (getAtomZ( carbon ) == 2) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 16)
					storeEClass(0x01C8);								// C=C(S)-S
				else
					storeEClass(0x01C9);							// C=CX-S  X!=S
				}
			else {
				if (getAtomSigma( carbon ) == 2)
					storeEClass(0x01C2);							// C=C(-C)-S
				else
					storeEClass(0x01C3);								// C=CH-S
				}
			return;
			}
		
		if (getAtomZ( carbon ) == 1) {
			if (mMol.isAllylicAtom( carbon ))
				storeEClass(0x01C7);								// C=C-CR2-S
		
			switch (getAtomSigma( carbon )) {
			case 3:
				storeEClass(0x01C6);									// (C-)3C-S
				return;
			case 2:
				storeEClass(0x01C5);								// (C-)2CH-S
				return;
			default:
				storeEClass(0x01C4);									// R-CH2-S
				return;
				}
			}
		
		if (getAtomZ( carbon ) == 2) {
			if (getAtomSigma( carbon ) == 0) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 16)
					storeEClass(0x01CE);								// S-CH2-S
				else
					storeEClass(0x01CF);							// X-CH2-S  X!=S
				return;
				}
			if (getAtomSigma( carbon ) == 1) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 16)
					storeEClass(0x01CC);							// S-CH(-C)-S
				else
					storeEClass(0x01CD);						// X-CH(-C)-S  X!=S
				return;
				}
			if (getAtomSigma( carbon ) == 2) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 16)
					storeEClass(0x01CA);							// S-C(-C)2-S
				else
					storeEClass(0x01CB);						// X-C(-C)2-S  X!=S
				return;
				}
			}
		
		if (getAtomZ( carbon ) == 3)
			{
			storeEClass(0x01D2);										// R-CX2-S
			return;
			}
		
		if (getAtomZ( carbon ) == 4)
			{
			storeEClass(0x01D3);										// S-CX3
			}
		
		return;
		}

	private int classI(int atm) {
		final int[] iodineClass = {
		 0x0160,						// Ar-I
		 0x0161,						// C#C-I
		 0x0163,					// C=CH-I
		 0x0162,					// C=C(-C)-I
		 0x0164,					// C=C-CR2-I
		 0x0167,		// non allylic (C-)3C-I
		 0x0166,		// non allylic (C-)2CH-I
		 0x0165,		// non allylic R-CH2-I
		 0x0168,					// C=CX-I
		 0x0169,					// C=C-CX-I
		 0x016A,		// non allylic R2CX-I
		 0x0170,					// R-CO-I
		 0x0172,					// R-C(=N)-I
		 0x0171,			// R-C(=X)-I  X!=O,N
		 0x016C,					// I3C-X=Y
		 0x016D,				// I3C-R(!stab)
		 0x016E,		// stabilized I-CX2-R
		 0x016F,	// non stabilized I-CX2-R
		 0x0173,						// I-C#N
		 0x0174,					// I-C(=O)-X
		 0x0175,					// I-C(=N)-X
		 0x0176,			// I-C(=X)-Y  X!=N,O
		 0x0177,						// I-CX3
		 0x017C,				// I at amide N
		 0x017D,				// I at amine N
		 0x0179,				// rare I-O-X
		 0x0178,						// I-O-C
		 0x017E,					// I at B
		 0x017F,					// I at Si
		 0x017A,					// I at P
		 0x017B						// I at S
		  };
		int eClass;
		
		mMol.setAtomMarker(atm, false);
		
		if (mMol.getConnAtoms( atm ) > 1)		// more than one neighbors to iodine
			return 0x016B;
		
		eClass = classHal(atm,53);
		if (eClass < 0) return eClass;
		return iodineClass[eClass];
		}

	private int classBr(int atm) {
		final int[] bromineClass = {
		 0x0140,						// Ar-Br
		 0x0141,						// C#C-Br
		 0x0143,					// C=CH-Br
		 0x0142,					// C=C(-C)-Br
		 0x0144,					// C=C-CR2-Br
		 0x0147,		// non allylic (C-)3C-Br
		 0x0146,		// non allylic (C-)2CH-Br
		 0x0145,		// non allylic R-CH2-Br
		 0x0148,					// C=CX-Br
		 0x0149,					// C=C-CX-Br
		 0x014A,		// non allylic R2CX-Br
		 0x0150,					// R-CO-Br
		 0x0152,					// R-C(=N)-Br
		 0x0151,			// R-C(=X)-Br  X!=O,N
		 0x014C,					// Br3C-X=Y
		 0x014D,				// Br3C-R(!stab)
		 0x014E,		// stabilized Br-CX2-R
		 0x014F,	// non stabilized Br-CX2-R
		 0x0153,						// Br-C#N
		 0x0154,					// Br-C(=O)-X
		 0x0155,					// Br-C(=N)-X
		 0x0156,			// Br-C(=X)-Y  X!=N,O
		 0x0157,						// Br-CX3
		 0x015C,				// Br at amide N
		 0x015D,				// Br at amine N
		 0x0159,				// rare Br-O-X
		 0x0158,						// Br-O-C
		 0x015E,					// Br at B
		 0x015F,					// Br at Si
		 0x015A,					// Br at P
		 0x015B						// Br at S
		  };
		int eClass;
		
		mMol.setAtomMarker(atm, false);
		
		if (mMol.getConnAtoms( atm ) > 1)	// more than one neighbors to bromine
			return -2;
		
		eClass = classHal(atm,35);
		if (eClass < 0) return eClass;
		return bromineClass[eClass];
		}

	private int classCl(int atm) {
		final int[] chlorineClass = {
		 0x0120,						// Ar-Cl
		 0x0121,						// C#C-Cl
		 0x0123,					// C=CH-Cl
		 0x0122,					// C=C(-C)-Cl
		 0x0124,					// C=C-CR2-Cl
		 0x0127,		// non allylic (C-)3C-Cl
		 0x0126,		// non allylic (C-)2CH-Cl
		 0x0125,		// non allylic R-CH2-Cl
		 0x0128,					// C=CX-Cl
		 0x0129,					// C=C-CX-Cl
		 0x012A,		// non allylic R2CX-Cl
		 0x0130,					// R-CO-Cl
		 0x0132,					// R-C(=N)-Cl
		 0x0131,			// R-C(=X)-Cl  X!=O,N
		 0x012C,					// Cl3C-X=Y
		 0x012D,				// Cl3C-R(!stab)
		 0x012E,		// stabilized Cl-CX2-R
		 0x012F,	// non stabilized Cl-CX2-R
		 0x0133,						// Cl-C#N
		 0x0134,					// Cl-C(=O)-X
		 0x0135,					// Cl-C(=N)-X
		 0x0136,			// Cl-C(=X)-Y  X!=N,O
		 0x0137,						// Cl-CX3
		 0x013C,				// Cl at amide N
		 0x013D,				// Cl at amine N
		 0x0139,				// rare Cl-O-X
		 0x0138,						// Cl-O-C
		 0x013E,					// Cl at B
		 0x013F,					// Cl at Si
		 0x013A,					// Cl at P
		 0x013B						// Cl at S
		  };
		int eClass;
		
		mMol.setAtomMarker(atm, false);
		
		if (mMol.getConnAtoms( atm ) > 1)	// more than one neighbors to chlorine
			return -2;
		
		eClass = classHal(atm,17);
		if (eClass < 0) return eClass;
		return chlorineClass[eClass];
		}

	private int classF(int atm) {
		final int[] fluorineClass = {
		 0x0100,						// Ar-F
		 0x0101,						// C#C-F
		 0x0103,					// C=CH-F
		 0x0102,					// C=C(-C)-F
		 0x0104,					// C=C-CR2-F
		 0x0107,		// non allylic (C-)3C-F
		 0x0106,		// non allylic (C-)2CH-F
		 0x0105,		// non allylic R-CH2-F
		 0x0108,					// C=CX-F
		 0x0109,					// C=C-CX-F
		 0x010A,		// non allylic R2CX-F
		 0x0110,					// R-CO-F
		 0x0112,					// R-C(=N)-F
		 0x0111,			// R-C(=X)-F  X!=O,N
		 0x010C,					// F3C-X=Y
		 0x010D,				// F3C-R(!stab)
		 0x010E,		// stabilized F-CX2-R
		 0x010F,	// non stabilized F-CX2-R
		 0x0113,						// F-C#N
		 0x0114,					// F-C(=O)-X
		 0x0115,					// F-C(=N)-X
		 0x0116,			// F-C(=X)-Y  X!=N,O
		 0x0117,						// F-CX3
		 0x011C,				// F at amide N
		 0x011D,				// F at amine N
		 0x0119,				// rare F-O-X
		 0x0118,						// F-O-C
		 0x011E,					// F at B
		 0x011F,					// F at Si
		 0x011A,					// F at P
		 0x011B						// F at S
		  };
		int eClass;
		
		mMol.setAtomMarker(atm, false);
		
		if (mMol.getConnAtoms( atm ) > 1)	// more than one neighbors to fluorine
			return -2;
		
		eClass = classHal(atm,9);
		if (eClass < 0) return eClass;
		return fluorineClass[eClass];
		}

	private int classHal(int atm, int halType) {
		int[] hetAtm = new int[3];
		
		if (mMol.getConnAtoms( atm ) != 1)	// HHal not classified
			return -2;
		
		int connAtm = mMol.getConnAtom( atm, 0 );
		if (mMol.getAtomicNo( connAtm ) == 6) {		// halogene at carbon
			if (bondToMet(connAtm)) return -1;
		
			if (mMol.isAromaticAtom( connAtm ))
				return 0;						// aryl halide
		
			switch(getAtomZ( connAtm )) {
			case 1:											// one hetero atom at carbon
				if (getAtomPi(connAtm) == 2)
					return 1;					// alkinyl halide
				if (getAtomPi(connAtm) == 1) {
					if (getAtomSigma( connAtm ) == 1)
						return 2;				// C=CH-Hal
					else
						return 3;				// C=C(-C)-Hal
					}
				if (mMol.isAllylicAtom( connAtm )) {
					return 4;					// allylic halide
					}
				else {
					switch(getAtomSigma( connAtm )) {
					case 3: return 5;			// !allylic tert. alkyl halide
					case 2: return 6;			// !allylic sec.  alkyl halide
					default: return 7;			// !allylic prim. alkyl halide
						}
					}
			case 2:											// two hetero atoms at carbon
				if (getHeteros(connAtm,hetAtm,atm) != 1) return -2;
		
				if (getAtomPi(connAtm) != 0)
					return 8;					// 1-halo-1-hetero-1-alkene
				if (mMol.isAllylicAtom( connAtm ))
					return 9;					// allylic 1-halo-1-hetero-alkane
				else
					return 10;					// !allylic 1-halo-1-hetero-alkane
			case 3:											// three hetero atoms at carbon
				int nrofHets = getHeteros(connAtm,hetAtm,atm);
				if (nrofHets == 1) {						// double bond to hetero atom
					if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
						mMol.setAtomMarker(hetAtm[0], false);
						return 11;				// R-CO-Hal
						}
					if (mMol.getAtomicNo( hetAtm[0] ) == 7) {
						return 12;				// Hal-C(=N)-R
						}
					else return 13;				// Hal-C(z=3)=X  X!=O,N
					}
				else {										// three single bonded hetero atoms
					if (mMol.getAtomicNo( hetAtm[0] ) == halType
					 && mMol.getAtomicNo( hetAtm[1] ) == halType) {
						mMol.setAtomMarker(hetAtm[0], false);
						mMol.setAtomMarker(hetAtm[1], false);
						if (mMol.isStabilizedAtom( connAtm ))
							return 14;			// Hal3C-X=Y
						else
							return 15;			// Hal3C(!stab)
						}
					if (mMol.isStabilizedAtom( connAtm ))
						return 16;				// stabilized Hal-CX2
					else
						return 17;				// !stabilized Hal-CX2
					}
			case 4:											// four hetero atoms at carbon
				int[] tripleBnds = new int[1];
				int[] doubleBnds = new int[1];
				nrofHets = getSortedHeteros(connAtm,hetAtm, tripleBnds,doubleBnds,atm);
				if (nrofHets == 1) {						// triple bond to hetero atom
					if (mMol.getConnAtoms( hetAtm[0] ) == 2)
						mMol.setAtomMarker(hetAtm[0], false);
					return 18;					// Hal-C#N
					}
				if (nrofHets == 2) {						// double bond to hetero atom
					if (mMol.getAtomicNo( hetAtm[0] ) == 8)
						return 19;				// Hal-C(=O)-Het
					if (mMol.getAtomicNo( hetAtm[0] ) == 7)
						return 20;				// Hal-C(=N)-Het
					else
						return 21;				// Hal-C(=X)-Het  X!=N,O
					}
				if (nrofHets == 3)				// three single bonded hetero atoms
					return 22;					// Hal-CHet3
			default: return -2;
				}
			}
		
		if (mMol.getAtomicNo( connAtm ) == 7) {		// halogene at nitrogen
			if (mMol.isStabilizedAtom( connAtm ))
				return 23;						// Hal-N(stab)
			else
				return 24;						// Hal-N(!stab)
			}
		
		if (mMol.getAtomicNo( connAtm ) == 8) {		// halogene at oxygen
			int nextConn = getNextConn(connAtm,atm);

			if (nextConn == -1 || mMol.getAtomicNo( nextConn ) == 6)
				return 26;						// Hal-O-C, Hal-OH
			else
				return 25;						// Hal-O-X  X=Met,Het
			}
		
		if (mMol.getAtomicNo( connAtm ) == 5)
			return 27;							// halogene at boron
		
		if (mMol.getAtomicNo( connAtm ) == 14)
			return 28;							// halogene at silicon
		
		if (mMol.getAtomicNo( connAtm ) == 15)
			return 29;							// halogene at phosphorous
		
		if (mMol.getAtomicNo( connAtm ) == 16)
			return 30;							// halogene at sulfur
		
		return -2;
		}

	private int classN(int atm) {
		int[] conn = new int[ExtendedMolecule.cMaxConnAtoms];
		int[] hetAtm = new int[ExtendedMolecule.cMaxConnAtoms];
		
		mMol.setAtomMarker(atm, false);
		
		if (mMol.getConnAtoms( atm ) == 0)	// don't handle NH3
			return -2;
		
		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
			int connAtm = mMol.getConnAtom( atm, i );
			if (mMol.isElectropositive( connAtm )) {
				if (mMol.getAtomicNo( connAtm ) == 14) {
					if (mMol.getConnAtoms( atm ) == 3) {
						if (mMol.isStabilizedAtom( atm ))
							return 0x009D;					// Si-N(-C,X)-C=X
						else
							return 0x009F;			// Si-N(-C,X)2  !stblzd N
						}
					if (mMol.isStabilizedAtom( atm ))
						return 0x009C;						// Si-NH-C=X type
					else
						return 0x009E;					// Si-NH-?  !stblzd N
					}
				if (mMol.getAtomicNo( connAtm ) == 5)
					return 0x00BE;									// B at N
				return 0x00BF;											// N-Met
				}
			}
		
		if (mMol.isAromaticAtom( atm ))
			{
			if (getAtomPi(atm) != 0) return 0x00B8;				// pyridine type
			else					return 0x00B9;				// pyrrole type
			}
		
		if (isThreeRingAtom( atm )) {
			int count = 0;
			for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
				conn[count] = mMol.getConnAtom( atm, i );
				if (isThreeRingAtom( conn[count] )) count++;
				}
			if (count > 1) {
				if (mMol.getAtomicNo( conn[0] ) == 6						// aziridine
				 && mMol.getAtomicNo( conn[1] ) == 6) {
					if (getAtomPi(conn[0]) != 0
					 || getAtomPi(conn[1]) != 0)
						return 0x02BC;					// alkylidene aziridine
					if (getAtomZ( conn[0] ) == 1
					 && getAtomZ( conn[1] ) == 1)
						{
						if (mMol.isAllylicAtom( conn[0] )
						 || mMol.isAllylicAtom( conn[1] ))
							return 0x02BD;					// alkenyl aziridine
						else
							return 0x02BF;					// alkyl aziridine
						}
					return 0x02BE;							// hetero aziridine
					}
				return 0x00BD;								// dihetero 3-ring
				}
			}
		
		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {	// classify all C-N carbons
			if (mMol.getConnBondOrder( atm, i ) != 1) continue;
			int connAtm = mMol.getConnAtom( atm, i );
			if (mMol.getAtomicNo( connAtm ) == 6) classCN(connAtm);
			}
		
												// These values may be changed from classCN()	*/
		
		int[] tripleBnds = new int[1];
		int[] doubleBnds = new int[1];
		int nrofHets = getSortedHeteros(atm,hetAtm,tripleBnds,doubleBnds,-1);
		
		if (mMol.getAtomCharge( atm ) == -1) return -1;
		
		if (mMol.getAtomCharge( atm ) == 1) {
			switch (getAtomZ( atm )) {
			case 4:
				if (mMol.getConnAtoms( atm ) == 2) {
					if (mMol.getAtomicNo( mMol.getConnAtom( atm, 0 ) ) == 7
					 && mMol.getAtomicNo( mMol.getConnAtom( atm, 1 ) ) == 7) {
						int connAtm = getNextConn(mMol.getConnAtom( atm, 0 ),atm);
						if (connAtm == -1)
							connAtm = getNextConn(mMol.getConnAtom( atm, 1 ),atm);
						if (connAtm == -1) return -2;
						if (hasDBondToHetero( connAtm ))
							return 0x0200;						// X=?-N3 type
						if (mMol.getAtomicNo( connAtm ) != 6)
							return 0x0201;								// X-N3
						if (getAtomPi(connAtm) != 0)
							return 0x0202;						// C=C-N3 type
						else
							return 0x0203;					// aliph C-N3 type
						}
					}
				if (mMol.getConnAtoms( atm ) == 3) {
					if (mMol.getAtomicNo( mMol.getConnAtom( atm, 0 ) ) == 8
					 && mMol.getAtomicNo( mMol.getConnAtom( atm, 1 ) ) == 8
					 && mMol.getAtomicNo( mMol.getConnAtom( atm, 1 ) ) == 8) {
						for (int i=0; i<3; i++) {
							int connAtm = mMol.getConnAtom( atm, i );
							if (mMol.getConnAtoms( connAtm ) == 1)
								mMol.setAtomMarker(connAtm, false);
							}
						return 0x0205;									// O-NO2
						}
					}
				return -1;
			case 3:
				if (nrofHets == 3) return 0x0213;					// R-N(+)X3
				if (nrofHets == 2) {
					if (mMol.getAtomicNo( hetAtm[0] ) == 8
					 && mMol.getAtomicNo( hetAtm[1] ) == 8
					 && mMol.getConnAtoms( hetAtm[1] ) == 1) {
						mMol.setAtomMarker(hetAtm[0], false);
						mMol.setAtomMarker(hetAtm[1], false);
						return 0x0209;									// R-NO2
						}
					if (mMol.getAtomicNo( hetAtm[0] ) == 7
					 && mMol.getAtomicNo( hetAtm[1] ) == 8) {
						if (mMol.getConnAtoms( hetAtm[1] ) == 1)
							mMol.setAtomMarker(hetAtm[1], false);
						return 0x020A;								// N=NR(+)-O
						}
					return 0x020B;								// rare X=NR(+)-Y
					}
				return 0x0207;										// R-N(+)#N
			case 2:
				return 0x0212;										// R2N(+)X2
			case 1:
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
					if (mMol.getConnAtoms( hetAtm[0] ) == 1)
						mMol.setAtomMarker(hetAtm[0], false);
					switch (getAtomPi(atm)) {
					case 2: return 0x020F;							// -C#N(+)-O
					default:return 0x020E;							// R3N(+)-O
						}
					}
				return 0x0211;									// R3N(+)-X  X!=O
			default:
				switch (getAtomPi(atm)) {
				case 2: return 0x0206;							// C-N(+)#C type
				default:return 0x0210;							// R4N(+) type
					}
				}
			}
		
		switch (getAtomZ( atm )) {
		case 3:
			if (tripleBnds[0] != 0) return -1;
			if (doubleBnds[0] != 0) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 7) {
					if (mMol.getAtomicNo( hetAtm[1] ) == 7)
						return 0x0267;									// N-N=N-
					if (mMol.getAtomicNo( hetAtm[1] ) == 8)
						return 0x0266;								// -O-N=N-
					return 0x0265;								// X-N=N-  X!=O,N
					}
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
					mMol.setAtomMarker(hetAtm[0], false);
					if (mMol.getAtomicNo( hetAtm[1] ) == 7)
						return 0x0260;									// N-N=O
					if (mMol.getAtomicNo( hetAtm[1] ) == 8)
						return 0x0261;									// -O-N=O
					return 0x0262;								// X-N=O  X!=O,N
					}
				return 0x0263;									// Y-N=X  X!=O,N
				}
			return 0x0264;												// NX3
		case 2:
			if (doubleBnds[0] != 0) {
				if (mMol.getAtomicNo( hetAtm[0] ) == 7) {
					if (mMol.getAtomCharge( hetAtm[0] ) == 1)
						return -1;
					return 0x0268;										// R-N=N-
					}
				if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
					mMol.setAtomMarker(hetAtm[0], false);
					return 0x0269;										// R-N=O
					}
				return 0x026A;									// R-N=X  X!=N,O
				}
			return 0x026B;												// R-NX2
		case 1:
			if (mMol.getAtomicNo( hetAtm[0] ) == 7) {
				switch (getAtomSigma( atm )) {
				case 2: return 0x026C;								// (C-)2N-N
				case 1: return 0x0274;									// C-NH-N
				default:return 0x0270;									// NH2-N
					}
				}
			if (mMol.getAtomicNo( hetAtm[0] ) == 8) {
				switch (getAtomSigma( atm )) {
				case 2:
					if (mMol.getConnAtoms( hetAtm[0] ) == 1) {
						mMol.setAtomMarker(hetAtm[0], false);
						return 0x026E;								// (C-)2N-OH
						}
					else
						return 0x026F;						// (C-)2N-O-?  ?!=H
				case 1:
					if (mMol.getConnAtoms( hetAtm[0] ) == 1) {
						mMol.setAtomMarker(hetAtm[0], false);
						return 0x0276;								// C-NH-OH
						}
					else
						return 0x0277;							// C-NH-O-?  ?!=H
				default:return 0x0271;									// -O-NH2
					}
				}
			switch (getAtomSigma( atm )) {
			case 2: return 0x026D;							// (C-)2N-X  X!=N,O
			case 1: return 0x0275;								// C-NH-X  X!=N,O
			default:return 0x0272;								// NH2-X  X!=N,O
				}
		default:
			if (mMol.getConnBondOrder( atm, 0 ) == 3) {
				int connAtm = mMol.getConnAtom( atm, 0 );
				int nextConn = getNextConn(connAtm, atm);
				if (nextConn == -1)
					return 0x01DC;										// H-CN
				if (mMol.isElectronegative( nextConn ))
					return 0x01DD;										// X-CN
				if (mMol.isElectropositive( nextConn ))
					{
					if (mMol.getAtomicNo( nextConn ) == 5
					 || mMol.getAtomicNo( nextConn ) == 14)
						return 0x01DE;							// B-CN; Si-CN
					return 0x01DF;							// Met-CN, Met!=B,Si
					}
				if (getAtomPi(nextConn) != 0)
					return 0x01D8;										// C=C-CN
				switch(getAtomZ( nextConn ))
					{
				case 0:
					return 0x01D9;										// R3C-CN
				case 1:
					return 0x01DA;									// R2CX-CN
				default:
					return 0x01DB;							// R-CX2-CN; CX3-CN
					}
				}

			int acyl = 0;
			for (int i=0; i<mMol.getConnAtoms( atm ); i++)
				if (hasDBondToHetero( mMol.getConnAtom( atm, i ) ))
					acyl++;
		
			switch (getAtomSigma( atm )) {
			case 3:
				switch (acyl) {
				case 0: return 0x0278;						// amine type (C-)3N
				case 1: return 0x027D;						// amide type (C-)3N
				default:return 0x027F;						// imide type (C-)3N
					}
			case 2:
				switch (acyl) {
				case 0: return 0x0279;						// amine type (C-)2NH
				case 1: return 0x027C;						// amide type (C-)2NH
				default:return 0x027E;						// imide type (C-)2NH
					}
			default:
				if (acyl != 0) return 0x027B;						// amide type C-NH2
				else		return 0x027A;						// amine type C-NH2
				}
			}
		}

	private void classCN(int atm) {
		if (mMol.isAromaticAtom( atm )) {
			storeEClass(0x0240);											// Ar-N
			return;
			}
		if (hasDBondToHetero( atm )) {
			storeEClass(0x0248);										// X=C-N
			return;
			}
		if (getAtomZ( atm ) == 1) {
			if (getAtomPi(atm) != 0) {
				if (getAtomSigma( atm ) == 1)
					storeEClass(0x0242);						// C#C-N / C=CH-N
				else
					storeEClass(0x0243);							// C=C(-C)-N
				return;
				}
		
			if (mMol.isAllylicAtom( atm ))
				storeEClass(0x0247);								// C=C-CR2-N
		
			switch (getAtomSigma( atm )) {
			case 3:
				storeEClass(0x0244);									// (C-)3C-N
				return;
			case 2:
				storeEClass(0x0245);								// (C-)2CH-N
				return;
			default:
				storeEClass(0x0246);									// R-CH2-N
				return;
				}
			}
		if (getAtomZ( atm ) == 2) {
			if (getAtomPi(atm) != 0) {
				storeEClass(0x0241);								// C=C(-X)-N
				return;
				}
			switch (getAtomSigma( atm )) {
			case 2:
				storeEClass(0x0249);								// (C-)2C(-X)-N
				return;
			case 1:
				storeEClass(0x024A);								// C-CH(-X)-N
				return;
			default:
				storeEClass(0x024B);									// X-CH2-N
				return;
				}
			}
		if (getAtomZ( atm ) == 3) {
			if (getAtomSigma( atm ) != 0)
				storeEClass(0x024C);									// C-CX2-N
			else
				storeEClass(0x024D);									// H-CX2-N
			return;
			}
		storeEClass(0x024E);											// CX3-N
		}

	private int classO(int atm) {
		final int[] base = {
		  0, 21, 41, 60, 78, 95,111,126,140,153,165,
		176,186,195,203,210,216,221,225,228,230,231  };
		final int[] alcClass = {
		0x0350,0x035F,0x0352,0x0353,0x0354,0x0355,0x0356,0x0357,0x0358,0x0359,0x035A,
		0x035B,0x035C,0x035D,0x035E,0x038F,0x03DF,0x03CF,0x039F,0x03EF,0x03FF,0x03AF  };
		final int[] etherClass = {
		0x0344,0x0300,0x0334,0x0335,0x0340,0x0341,0x0342,0x0343,0x0314,0x01E4,0x01E6,
		0x01E7,0x02F0,0x02F1,0x03B0,0x0380,0x03D0,0x03C0,0x0390,0x03E0,0x03F0,0x03A0,
		0x0301,0x0325,0x0324,0x0304,0x0305,0x0306,0x0307,0x0318,0x0308,0x030A,0x030B,
		0x030C,0x030D,0x03B1,0x0381,0x03D1,0x03C1,0x0391,0x03E1,0x03F1,0x03A1,
		0x0330,0x0331,0x0338,0x0339,0x033A,0x033B,0x0316,0x0329,0x032D,0x032F,0x0321,
		0x0323,0x0327,0x0382,0x03D2,0x03C2,0x0392,0x03E2,0x03F2,0x03A2,
		0x0332,0x033C,0x033D,0x033E,0x033F,0x0317,0x0328,0x032C,0x032E,0x0320,0x0322,
		0x0326,0x0383,0x03D3,0x03C3,0x0393,0x03E3,0x03F3,0x03A3,
		0x0345,0x0346,0x0348,0x0349,0x0310,0x01E0,0x01E8,0x01EC,0x01F8,0x01FC,0x03B4,
		0x0384,0x03D4,0x03C4,0x0394,0x03E4,0x03F4,0x03A4,
		0x0347,0x034A,0x034C,0x0311,0x01E1,0x01E9,0x01ED,0x01F9,0x01FD,0x03B5,0x0385,
		0x03D5,0x03C5,0x0395,0x03E5,0x03F5,0x03A5,
		0x034D,0x034E,0x0312,0x01E2,0x01EA,0x01EE,0x01FA,0x01FE,0x03B6,0x0386,0x03D6,
		0x03C6,0x0396,0x03E6,0x03F6,0x03A6,
		0x034F,0x0313,0x01E3,0x01EB,0x01EF,0x01FB,0x01FF,0x03B7,0x0387,0x03D7,0x03C7,
		0x0397,0x03E7,0x03F7,0x03A7,
		0x0315,0x031E,0x031D,0x031C,0x031A,0x031B,0x03B8,0x0388,0x03D8,0x03C8,0x0398,
		0x03E8,0x03F8,0x03A8,
		0x01F0,0x01F2,0x01F4,0x02F8,0x02FC,0x03B9,0x0389,0x03D9,0x03C9,0x0399,0x03E9,
		0x03F9,0x03A9,
		0x01F3,0x01F6,0x02FA,0x02FE,0x03BA,0x038A,0x03DA,0x03CA,0x039A,0x03EA,0x03FA,
		0x03AA,
		0x01F7,0x02FB,0x02FF,0x03BB,0x038B,0x03DB,0x03CB,0x039B,0x03EB,0x03FB,0x03AB,
		0x02F4,0x02F6,0x03BC,0x038C,0x03DC,0x03CC,0x039C,0x03EC,0x03FC,0x03AC,
		0x02F7,0x03BD,0x038D,0x03DD,0x03CD,0x039D,0x03ED,0x03FD,0x03AD,
		0x03BE,0x038E,0x03DE,0x03CE,0x039E,0x03EE,0x03FE,0x03AE,
		0x036A,0x0370,0x0378,0x0373,0x0360,0x0368,0x037C,
		0x0375,0x0374,0x0371,0x0365,0x036D,0x0376,
		0x037A,0x0379,0x0364,0x036C,0x037B,
		0x0372,0x0361,0x0369,0x037D,
		0x0366,0x0367,0x0362,
		0x036E,0x036F,
		0x037E  };
		
		mMol.setAtomMarker(atm, false);
		
		if (mMol.getConnAtoms( atm ) == 0)				// don't handle H2O
			return -2;
		
		if (mMol.isAromaticAtom( atm )) return 0x00BA;		// O(arom)
		
		if (mMol.getConnAtoms( atm ) == 1) {
			if (mMol.getConnBondOrder( atm, 0 ) == 2) return -1;
										// carbonyl classified elsewhere
			int class1 = classAtmO(mMol.getConnAtom( atm, 0 ));
			return alcClass[class1];								// 22 alcohols
			}
		
		int connAtm1 = mMol.getConnAtom( atm, 0 );
		int connAtm2 = mMol.getConnAtom( atm, 1 );
		
		if (isThreeRingAtom( atm )) {
			if (mMol.getAtomicNo( connAtm1 ) == 6							// oxirane
			 && mMol.getAtomicNo( connAtm2 ) == 6) {
				if (getAtomPi(connAtm1) != 0
				 || getAtomPi(connAtm2) != 0)
					return 0x02B8;								// allene oxide
				if (getAtomZ( connAtm1 ) == 1
				 && getAtomZ( connAtm2 ) == 1)
					{
					if (mMol.isAllylicAtom( connAtm1 )
					 || mMol.isAllylicAtom( connAtm2 ))
						return 0x02B9;						// alkenyl oxirane
					else
						return 0x02BB;							// alkyl oxirane
					}
				return 0x02BA;							// hetero subst. oxirane
				}
			return 0x00BD;									// dihetero 3-ring
			}
		
		int class1 = classAtmO(connAtm1);
		int class2 = classAtmO(connAtm2);
		if (class1 > class2) {
			int temp = class1;
			class1 = class2;
			class2 = temp;
			}
		
		return etherClass[base[class1]+class2];
		}

	private int classAtmO(int atm) {
		if (mMol.getAtomicNo( atm ) == 6) {
			if (mMol.isAromaticAtom( atm )) return 0;			// Ar-O
			if (hasDBondToHetero( atm )) return 1;	// X=C-O
			if (getAtomZ( atm ) == 1) {
				if (getAtomPi(atm) != 0) {
					if (getAtomSigma( atm ) == 1)
						return 2;								// C#C-O / C=CH-O
					else
						return 3;									// C=C(-C)-O
					}
				switch (getAtomSigma( atm )) {
				case 3:  return 4;									// (C-)3C-O
				case 2:  return 5;									// (C-)2CH-O
				case 1:  return 6;									// C-CH2-O
				default: return 7;										// CH3-O
					}
				}
			if (getAtomZ( atm ) == 2) {
				if (getAtomPi(atm) != 0)
					return 8;										// C=C(-X)-O
				switch (getAtomSigma( atm )) {
				case 2:  return  9;								// (C-)2C(-X)-O
				case 1:  return 10;								// C-CH(-X)-O
				default: return 11;									// X-CH2-O
					}
				}
			if (getAtomZ( atm ) == 3) {
				if (getAtomSigma( atm ) != 0)
					return 12;										// C-CX2-O
				else
					return 13;										// H-CX2-O
				}
			return 14;													// CX3-O
			}
		else {
			switch (mMol.getAtomicNo( atm )) {
			case  5: return 15;											// B-O
			case  7: return 16;											// N-O
			case  8: return 17;											// O-O
			case 14: return 18;											// Si-O
			case 15: return 19;											// P-O
			case 16: return 20;											// S-O
			default: return 21;										// rare X-O
				}
			}
		}

	private int classC(int atm) {
		final int[][] alkyneClass = {	{ -1,0x00C0,0x00C2,0x00C1,0x00C3 },
										{ 0x00C0,0x00C8,0x00CA,0x00C9,0x00CB },
										{ 0x00C2,0x00CA,0x00C4,0x00C5,0x00C6 },
										{ 0x00C1,0x00C9,0x00C5,0x00CC,0x00CD },
										{ 0x00C3,0x00CB,0x00C6,0x00CD,0x00C7 } };
		
		mMol.setAtomMarker(atm, false);
		
		if (isThreeRingAtom( atm )) {
			int[] member = new int[3];
			getRingMembers( atm, 3, false, member );
			if (mMol.getAtomicNo( member[1] ) != 6
			 || mMol.getAtomicNo( member[2] ) != 6) return -2;
			boolean vinyl = false;
			boolean stab = false;
			boolean allyl = false;
			for (int i=0; i<3; i++) {
				mMol.setAtomMarker(member[i], false);
				if (getAtomPi(member[i]) != 0) vinyl = true;
				if (mMol.isStabilizedAtom( member[i] )) stab = true;
				if (mMol.isAllylicAtom( member[i] )) allyl = true;
				}
			if (vinyl) return 0x00E4;				// C3-ring with double bond
			if (stab)  return 0x00E6;						// stabilized C3-ring
			if (allyl) return 0x00E7;							// allylic C3-ring
			return 0x00E5;								// featureless C3-ring
			}
		
		if (getAtomPi(atm) == 0) {
			if (mMol.getAtomParity( atm ) == 0)
				return -2;
		
			if (getAtomSigma( atm ) + getAtomZ( atm ) == 4)
				return 0x00F0;									// XnC*(-C)4-n
			if (mMol.isStabilizedAtom( atm ))
				return 0x00F1;										// X=?-C*H
			if (mMol.isAllylicAtom( atm ))
				return 0x00F2;										// C=C-C*H
			return 0x00F3;									// featureless C*H
			}

		int[] member = new int[ExtendedMolecule.cMaxConnAtoms];
		int nrofPiConns = getPiConnCarbs(atm,member);
		if (getAtomPi(atm) == 2 && nrofPiConns == 1) {				// C#C
			mMol.setAtomMarker(member[0], false);
			int class1 = classAcetylen(atm,member[0]);
			int class2 = classAcetylen(member[0],atm);
			if (class1 == 5 || class2 == 5) return 0x00CE;			// ?-C#C-B
			if (class1 == 6 || class2 == 6) return 0x00CF;	/*?-C#C-Met  Met,?!=B*/
			return alkyneClass[class1][class2];
			}
		
		if (bondToMet(atm)) return -1;
		for (int i=0; i<nrofPiConns; i++)
			if (bondToMet(member[i])) return -1;
		
		if (getAtomPi(atm) == 2) {
			boolean stab = false;
			for (int i=0; i<2; i++) {
				if (mMol.isStabilizedAtom( member[i] )) stab = true;
				mMol.setAtomMarker(member[i], false);
				}
			if (getAtomZ( member[0] ) != 0
			 || getAtomZ( member[1] ) != 0) {
				if (stab) return 0x00E8;						// X-C=C=C-C=X type
				else		return 0x00E9;						// X-C=C=CR2 type
				}
			else {
				if (stab) return 0x00EA;						// C=C=C-C=X type
				else		return 0x00EB;							// dull C=C=C
				}
			}
		
		if (getAtomPi(member[0]) == 2) return -1;	// don't classify allene here
		mMol.setAtomMarker(member[0], false);
		
		if (mMol.isAromaticAtom( atm )) {
			int[] aromAtms = new int[16];
			int aroms = getRingMembers( atm, 0, true, aromAtms );
			if (aroms == 0)
				return -2;
		
			boolean stab = false;
			boolean hetero = false;
			for (int i=0; i<aroms; i++) {
				mMol.setAtomMarker(aromAtms[i], false);
				for (int j=0; j<mMol.getConnAtoms( aromAtms[i] ); j++) {
					int connAtm = mMol.getConnAtom( aromAtms[i], j );
					if (mMol.isAromaticAtom( connAtm )) continue;
					if (hasDBondToHetero( connAtm )) {
						stab = true;
						continue;
						}
					if (mMol.isElectronegative( connAtm )) {
						hetero = true;
						}
					}
				}
			if (stab)   storeEClass(0x00EC);					// X=?-Ar
			if (hetero) storeEClass(0x00ED);					// X-Ar
			if (!stab && !hetero) return 0x00EE;				// featureless Ar
			return -1;
			}
		
		member[1] = atm;
		int vinyl = 0;
		int stab = 0;
		int hetero = 0;
		int alkyl = 0;
		int phenyl = 0;

		for (int i=0; i<2; i++) {
			for (int j=0; j<mMol.getConnAtoms( member[i] ); j++) {
				int connAtm = mMol.getConnAtom( member[i], j );
				if (mMol.isElectropositive( connAtm )) continue;
				if (connAtm == member[1-i]) continue;
				if (mMol.isAromaticAtom( connAtm )) {
					phenyl++;
					continue;
					}
				if (hasDBondToHetero( connAtm )) {
					stab++;
					continue;
					}
				if (mMol.isElectronegative( connAtm )) {
					hetero++;
					continue;
					}
				if (getAtomPi(connAtm) != 0) {
					vinyl++;
					continue;
					}
				alkyl++;
				}
			}
		if (hetero > 2) {
			if (stab != 0) return 0x00DF;								// X=?-XC=CX2
			else		return 0x00E0;									// X2C=CX-?
			}
		if (hetero == 2) {
			if (stab != 0) return 0x00DE;							// X=?-XC=CX-? type
			else		return 0x00E1;							// ?-XC=CX-? type
			}
		if (stab > 2) return 0x00D8;						// (X=?-)2C=C-?=X type
		if (stab == 2) {
			if (hetero != 0) return 0x00DC;						// X=?-XC=C-?=X type
			else		return 0x00D9;						// X=?-RC=CR-?=X type
			}
		if (hetero == 1) {
			if (stab == 1) return 0x00DD;						// X=?-C=CX type
			if (vinyl + alkyl + phenyl >1) return 0x00E3;		// C-C=CX-C type
			return 0x00E2;										// R-HC=CH-X type
			}
		if (stab == 1) {
			if (vinyl + alkyl + phenyl >1) return 0x00DA;   // C-C=C(-C)-?=X type
			return 0x00DB;									// R-HC=CH-C=X type
			}
		switch (vinyl + alkyl + phenyl) {
		case 4:
			if (alkyl == 4) return 0x00D6;						// (C-)2C=C(-C)2
			return 0x00D7;									// C=C-C(-C)=C(-C)2
		case 3:
			if (alkyl == 3) return 0x00D4;						// C-CH=C(-C)2
			return 0x00D5;									// C=C-CH=C(-C)2 type
		case 2:
			if (alkyl == 2) return 0x00D2;						// C-CH=CH-C type
			return 0x00D3;									// C=C-CH=CH-C type
		default:
			if (alkyl == 1) return 0x00D0;							// C-CH=CH2
			return 0x00D1;											// C=C-CH=CH2
			}
		}

	private int classAcetylen(int atm, int other) {
		int nextConn;
		
		if (mMol.getConnAtoms( atm ) == 1) return 0;			// C#C-H
		nextConn = getNextConn(atm,other);
		
		if (mMol.isAromaticAtom( nextConn )) return 1;			// C#C-C=X
		if (mMol.isElectronegative( nextConn )) return 2;			// C#C-X
		
		if (mMol.isElectropositive( nextConn )) {
			if (mMol.getAtomicNo( nextConn ) == 5) return 5;	// C#C-B
			else								return 6;   	// ?-C#C-Met  Met,?!=B
			}
		
		if (getAtomPi(nextConn) != 0) return 3;						// C#C-C=C
			return 4;											// C#C-C(aliphatic)
		}

	private void storeEClass( int eClass ) {
		if (eClass == -1) return;					// FG was ealier classified by another modul
		
		if (eClass == -2) return;					// don't care about FG recognition
													// error here; just don't store FG

		Integer count = mFunctionalGroupCountMap.get(eClass);
		mFunctionalGroupCountMap.put(eClass, count == null ? 1 : count + 1);

		return;
		}

	private int getPiConnCarbs( int atm, int[] member) {
		int count = 0;
		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
			if (mMol.getConnBondOrder( atm, i ) > 1) {
				int connAtm = mMol.getConnAtom( atm, i );
				if (mMol.getAtomicNo( connAtm ) == 6)
					member[count++] = connAtm;
				}
			}
		return count;
		}

	/**
	 * Returns all ring member atoms
	 * @param atom
	 * @param ringSize if 0 then return ring of any size
	 * @param aromatic consider aromatic rings only
	 * @param member
	 * @return
	 */
	private int getRingMembers( int atom, int ringSize, boolean aromatic, int member[] ) {
		RingCollection ringSet = mMol.getRingSet();
		
		member[0] = atom;
		for (int i=0; i<ringSet.getSize(); i++) {
			if ((ringSet.getRingSize(i) == ringSize || ringSize == 0)
			 && (ringSet.isAromatic(i) || !aromatic)) {
				int[] ringAtom = ringSet.getRingAtoms(i);
				for (int j=0; j<ringAtom.length; j++) {
					if (ringAtom[j] == atom) {
						int no = 1;
						for (int k=0; k<ringAtom.length; k++)
							if (ringAtom[k] != atom)
								member[no++] = ringAtom[k];
						return ringAtom.length;
						}
					}
				}
			}
		
		return 0;
		}

	private boolean hasDBondToHetero(int atom) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (mMol.getConnBondOrder(atom, i) == 2
			 && mMol.isElectronegative(mMol.getConnAtom(atom, i)))
				return true;
		return false;
		}

	private int getCarbons(int atm, int[] carbon) {
		int carbons = 0;
		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
			int connAtm = mMol.getConnAtom( atm, i );
			if (mMol.getAtomicNo( connAtm ) == 6)
				carbon[carbons++] = connAtm;
			}
		return carbons;
		}

	private int getHeteros(int atm, int[] hetAtm, int notThisOne) {
		int heteros = 0;
		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
			int connAtm = mMol.getConnAtom( atm, i );
			if (connAtm == notThisOne) continue;
			if (mMol.isElectronegative( connAtm ))
				hetAtm[heteros++] = connAtm;
			}
		return heteros;
		}

	private int getSortedHeteros(int atm, int[] hetAtm, int[] tripleBnds, int[] doubleBnds, int notThisOne) {
		int heteros = 0;
		
		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
			if (mMol.getConnBondOrder( atm, i ) != 3) continue;
			int connAtm = mMol.getConnAtom( atm, i );
			if (connAtm == notThisOne) continue;
			if (mMol.isElectronegative( connAtm ))
				hetAtm[heteros++] = connAtm;
			}
		
		tripleBnds[0] = heteros;
		
		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
			if (mMol.getConnBondOrder( atm, i ) != 2) continue;
			int connAtm = mMol.getConnAtom( atm, i );
			if (connAtm == notThisOne) continue;
			if (mMol.isElectronegative( connAtm ))
				hetAtm[heteros++] = connAtm;
			}
		
		doubleBnds[0] = heteros - tripleBnds[0];
		
		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
			if (mMol.getConnBondOrder( atm, i ) != 1) continue;
			int connAtm = mMol.getConnAtom( atm, i );
			if (connAtm == notThisOne) continue;
			if (mMol.isElectronegative( connAtm ))
				hetAtm[heteros++] = connAtm;
			}
		
		return heteros;
		}

	private int getNextConn(int atm, int notThisOne) {
		for (int i=0; i<mMol.getConnAtoms( atm ); i++) {
			if (mMol.getConnAtom( atm, i ) == notThisOne) continue;
			return mMol.getConnAtom( atm, i );
			}
		return -1;
		}

	/**
	 * Determine whether rear side of hetero atoms are backed by a carbonyl- or similar group.
	 * @param hetAtm
	 * @param notThisOne
	 * @return
	 */
	private boolean checkAnhydride(int hetAtm, int notThisOne) {
		for (int i=0; i<mMol.getConnAtoms( hetAtm ); i++) {
			int connAtm = mMol.getConnAtom( hetAtm, i );
			if (connAtm == notThisOne) continue;
			if (hasDBondToHetero( connAtm )) return true;
			}
		return false;
		}

	private boolean bondToMet(int atm) {
		for (int i=0; i<mMol.getConnAtoms(atm); i++)
			if (mMol.isMetalAtom( mMol.getConnAtom(atm, i)))
				return true;
		
		return false;
		}
		
		
	private boolean isQuinone(int atm) {
		int[] members = new int[6];
		int[] hetAtm = new int[ExtendedMolecule.cMaxConnAtoms];
		
		if (getRingMembers( atm, 6, false, members ) != 6) return false;
		boolean secCOfound = false;
		for (int i=1; i<6; i++) {
			if (mMol.getAtomicNo( members[i] ) != 6) return false;
			if (getAtomPi(members[i]) == 0) {
				if ((i & 1) == 0) return false;
				if (secCOfound) return false;
				if (getAtomZ( members[i] ) != 2) return false;
				if (getHeteros(members[i],hetAtm,-1) != 1) return false;
				if (mMol.getAtomicNo( hetAtm[0] ) != 8) return false;
				secCOfound = true;
				}
			}
		if (!secCOfound) return false;
		for (int i=1; i<6; i++)
			mMol.setAtomMarker(members[i], false);
		return true;
		}

	private boolean isThreeRingAtom(int atom) {
		return mMol.getAtomRingSize(atom) == 3;
		}

	private int getAtomPi(int atom) {
		int pi = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			if (mMol.getAtomicNo(mMol.getConnAtom(atom, i)) == 6)
				pi += mMol.getConnBondOrder(atom, i) - 1;
			}
		return pi;
		}

	private int getAtomSigma(int atom) {
		int sigma = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			if (mMol.getAtomicNo(mMol.getConnAtom(atom, i)) == 6)
				sigma++;
			}
		return sigma;
		}

	private int getAtomZ(int atom) {
		int z = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			if (mMol.isElectronegative(mMol.getConnAtom(atom, i)))
				z += mMol.getConnBondOrder(atom, i);
			}
		return z;
		}
	}
