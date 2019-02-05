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
*/

package com.actelion.research.chem;

public class MolecularFormula {
	private static final double[] sRelativeMass = { 0.0,
	1.00794, 4.0026, 6.9410, 9.0122, 10.811, 12.011,    //  H , He, Li, Be, B , C , 
	 14.007, 15.999, 18.998, 20.180, 22.990, 24.305,    //  N , O , F , Ne, Na, Mg, 
	 26.982, 28.086, 30.974, 32.066, 35.453, 39.948,    //  Al, Si, P , S , Cl, Ar, 
	 39.098, 40.078, 44.956, 47.867, 50.942, 51.996,    //  K , Ca, Sc, Ti, V , Cr, 
	 54.938, 55.845, 58.933, 58.693, 63.546, 65.390,    //  Mn, Fe, Co, Ni, Cu, Zn, 
	 69.723, 72.610, 74.922, 78.960, 79.904, 83.800,    //  Ga, Ge, As, Se, Br, Kr, 
	 85.468, 87.620, 88.906, 91.224, 92.906, 95.940,    //  Rb, Sr, Y , Zr, Nb, Mo, 
	 98.906, 101.07, 102.91, 106.42, 107.87, 112.41,    //  Tc, Ru, Rh, Pd, Ag, Cd, 
	 114.82, 118.71, 121.76, 127.60, 126.90, 131.29,    //  In, Sn, Sb, Te, I , Xe, 
	 132.91, 137.33, 138.91, 140.12, 140.91, 144.24,    //  Cs, Ba, La, Ce, Pr, Nd, 
	 146.92, 150.36, 151.96, 157.25, 158.93, 162.50,    //  Pm, Sm, Eu, Gd, Tb, Dy, 
	 164.93, 167.26, 168.93, 173.04, 174.97, 178.49,    //  Ho, Er, Tm, Yb, Lu, Hf, 
	 180.95, 183.84, 186.21, 190.23, 192.22, 195.08,    //  Ta, W,  Re, Os, Ir, Pt, 
	 196.97, 200.59, 204.38, 207.20, 208.98, 209.98,    //  Au, Hg, Tl, Pb, Bi, Po, 
	 209.99, 222.02, 223.02, 226.03, 227.03, 232.04,    //  At, Rn, Fr, Ra, Ac, Th, 
	 231.04, 238.03, 237.05, 239.05, 241.06, 244.06,    //  Pa, U,  Np, Pu, Am, Cm, 
	 249.08, 252.08, 252.08, 257.10, 258.10, 259.10,    //  Bk, C,  Es, Fm, Md, No, 
	 262.11, 267.12, 268.13, 271.13, 270.13, 277.15,    //  Lr, Rf, Db, Sg, Bh, Hs,
	 276.15, 281.17, 281.17, 283.17, 285.18, 289.19,    //  Mt ,Ds ,Rg ,Cn ,Nh ,Fl,
	 289.19, 293.20, 294.21, 294.21,    0.0,    0.0,    //  Mc ,Lv ,Ts ,Og ,?? ,??,
	    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    //  ??, ??, ??, ??, ??, ??, 
	    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    //  ??, ??, R4, R5, R6, R7, 
	    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    //  R8, R9, R10,R11,R12,R13, 
	    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    //  R14,R15,R16,R1, R2, R3, 
	    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    //  A , A1, A2, A3, ??, ??, 
	 2.0141, 3.0160,    0.0,    0.0,    0.0,    0.0,    //  D , T , X , R , H2, H+, 
	    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    //  Nnn,HYD,Pol,??, ??, ??, 
	    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    //  ??, ??, ??, ??, ??, ??, 
	    0.0,    0.0, 71.0787, 156.18828, 114.10364, 115.0877,  // ??, ??, Ala,Arg,Asn,Asp,    relative masses of amino acids reduced by one H2O
	103.1447,  128.13052, 129.11458,  57.05182, 137.14158, 113.15934, //  Cys,Gln,Glu,Gly,His,Ile,
	113.15934, 128.17428, 131.19846, 147.17646,  97.11658,  87.0777,  //  Leu,Lys,Met,Phe,Pro,Ser,
	101.10458, 186.2134 , 163.17546,  99.13246 };                     //  Thr,Trp,Tyr,Val,

	private static final double[] sAbsoluteMass = { 0.0, 
	   1.007825,   4.00260 ,   7.016003,   9.012182,	//  H, He,Li,Be
	  11.009305,  12.000000,  14.003074,  15.994915,	//  B ,C ,N, O
	  18.998403,  19.992435,  22.989767,  23.985042,	//  F ,Ne,Na,Mg
	  26.98153 ,  27.976927,  30.973762,  31.972070,	//  Al,Si,P ,S
	  34.968852,  39.962384,  38.963707,  39.962591,	//  Cl,Ar,K ,Ca
	  44.955910,  47.947947,  50.943962,  51.940509,	//  Sc,Ti,V ,Cr
	  54.938047,  55.934939,  58.933198,  57.935346,	//  Mn,Fe,Co,Ni
	  62.939598,  63.929145,  68.925580,  73.921177,	//  Cu,Zn,Ga,Ge
	  74.921594,  79.916520,  78.918336,  83.911507,	//  As,Se,Br,Kr
	  84.911794,  87.905619,  88.905849,  89.904703,	//  Rb,Sr,Y ,Zr
	  92.906377,  97.905406,  89.923810, 101.904348,	//  Nb,Mo,Tc,Ru
	 102.905500, 105.903478, 106.905092, 113.903357,	//  Rh,Pd,Ag,Cd
	 114.903880, 119.902200, 120.903821, 129.906229,	//  In,Sn,Sb,Te
	 126.904473, 131.904144, 132.905429, 137.905232,	//  I ,Xe,Cs,Ba
	 138.906346, 139.905433, 140.907647, 141.907719,	//  La,Ce,Pr,Nd
	 135.923980, 151.919729, 152.921225, 157.924099,	//  Pm,Sm,Eu,Gd
	 158.925342, 163.929171, 164.930319, 165.930290,	//  Tb,Dy,Ho,Er
	 168.934212, 173.938859, 174.940770, 179.946545,	//  Tm,Yb,Lu,Hf
	 180.947992, 183.950928, 186.955744, 191.961467,	//  Ta,W, Re,Os
	 192.962917, 194.964766, 196.966543, 201.970617,	//  Ir,Pt,Au,Hg
	 204.974401, 207.976627, 208.980374, 193.988180,	//  Tl,Pb,Bi,Po
	 195.995730, 199.995700, 201.004110, 206.003800,	//  At,Rn,Fr,Ra
	 210.009230, 232.038054, 216.018960, 238.050784,	//  Ac,Th,Pa,U
	 229.036230, 232.041169, 237.050050, 238.053020,	//  Np,Pu,Am,Cm
	 242.061940, 240.062280, 243.069470, 243.074460,	//  Bk,C ,Es,Fm
	 248.082750, 251.088870, 253.095150, 257.102950,	//  Md,No,Lr,Rf
	 257.107770, 271.13,     270.13,     277.15,		//  Db,Sg,Bh,Hs
	 276.15,     281.17,     281.17,     283.17,		//  Mt,Ds,Rg,Cn
	 285.18,     289.19,     289.19,     291.20,		//  Nh,Fl,Mc,Lv
	 294.21,     294.21,       0.0f,       0.0f,		//  Ts,Og,??,??
	   0.0f,       0.0f,       0.0f,       0.0f,		//  ??,??,??,??
	   0.0f,       0.0f,       0.0f,       0.0f,		//  ??,??,??,??
	   0.0f,       0.0f,       0.0f,       0.0f,		//  R4,R5,R6,R7
	   0.0f,       0.0f,       0.0f,       0.0f,		//  R8,R9,R10,R11
	   0.0f,       0.0f,       0.0f,       0.0f,		//  R12,R13,R14,R15
	   0.0f,       0.0f,       0.0f,       0.0f,		//  R16,R1,R2,R3
	   0.0f,       0.0f,       0.0f,       0.0f,		//  A ,A1,A2,A3
	   0.0f,       0.0f,       2.0140,     3.01605,		//  ??,??,D ,T
	   0.0f,       0.0f,       0.0f,       0.0f,   	    //  X ,R ,H2,H+
	   // 12/18/2007 CxR Added the following
       0.0f,       0.0f,       0.0f,       0.0f,        //  Nnn,HYD,Pol,??
       0.0f,       0.0f,       0.0f,       0.0f,        //  ??,??,??,??
       0.0f,       0.0f,       0.0f,       0.0f,        //  ??,??,??,??
       0.0f,       0.0f,       0.0f,       0.0f,        //  ??, ??, Ala,Arg
       0.0f,       0.0f,       0.0f,       0.0f,        //  Asn,Asp,,Cys,Gln
       0.0f,       0.0f,       0.0f,       0.0f,        //  Glu,Gly,His,Ile,
       0.0f,       0.0f,       0.0f,       0.0f,        //  Leu,Lys,Met,Phe
       0.0f,       0.0f,       0.0f,       0.0f,        //  Pro,Ser,Thr,Trp,
       0.0f,       0.0f,                                //  Tyr,Val
	};

	private static final int[] sFirstFormulaAtomicNo = { 6, 1, 7, 8 };

	private int[] mAtomCount,mAtomicNo;
	private double mAbsoluteIsotopeWeightIncrement,mRelativeIsotopeWeightIncrement;

	public MolecularFormula(ExtendedMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		int[] atomCount = new int[Molecule.cMaxAtomicNo+1];
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			switch (mol.getAtomicNo(atom)) {
			case 171:	// Ala
				atomCount[1] += 5;
				atomCount[6] += 3;
				atomCount[7] += 1;
				atomCount[8] += 1;
				break;
			case 172:	// Arg
				atomCount[1] += 12;
				atomCount[6] += 6;
				atomCount[7] += 4;
				atomCount[8] += 1;
				break;
			case 173:	// Asn
				atomCount[1] += 6;
				atomCount[6] += 4;
				atomCount[7] += 2;
				atomCount[8] += 2;
				break;
			case 174:	// Asp
				atomCount[1] += 5;
				atomCount[6] += 4;
				atomCount[7] += 1;
				atomCount[8] += 3;
				break;
			case 175:	// Cys
				atomCount[1] += 5;
				atomCount[6] += 3;
				atomCount[7] += 1;
				atomCount[8] += 1;
				atomCount[16] += 1;
				break;
			case 176:	// Gln
				atomCount[1] += 8;
				atomCount[6] += 5;
				atomCount[7] += 2;
				atomCount[8] += 2;
				break;
			case 177:	// Glu
				atomCount[1] += 7;
				atomCount[6] += 5;
				atomCount[7] += 1;
				atomCount[8] += 3;
				break;
			case 178:	// Gly
				atomCount[1] += 3;
				atomCount[6] += 2;
				atomCount[7] += 1;
				atomCount[8] += 1;
				break;
			case 179:	// His
				atomCount[1] += 7;
				atomCount[6] += 6;
				atomCount[7] += 3;
				atomCount[8] += 1;
				break;
			case 180:	// Ile
				atomCount[1] += 11;
				atomCount[6] += 6;
				atomCount[7] += 1;
				atomCount[8] += 1;
				break;
			case 181:	// Leu
				atomCount[1] += 11;
				atomCount[6] += 6;
				atomCount[7] += 1;
				atomCount[8] += 1;
				break;
			case 182:	// Lys
				atomCount[1] += 12;
				atomCount[6] += 6;
				atomCount[7] += 2;
				atomCount[8] += 1;
				break;
			case 183:	// Met
				atomCount[1] += 9;
				atomCount[6] += 5;
				atomCount[7] += 1;
				atomCount[8] += 1;
				atomCount[16] += 1;
				break;
			case 184:	// Phe
				atomCount[1] += 9;
				atomCount[6] += 9;
				atomCount[7] += 1;
				atomCount[8] += 1;
				break;
			case 185:	// Pro
				atomCount[1] += 7;
				atomCount[6] += 5;
				atomCount[7] += 1;
				atomCount[8] += 1;
				break;
			case 186:	// Ser
				atomCount[1] += 5;
				atomCount[6] += 3;
				atomCount[7] += 1;
				atomCount[8] += 2;
				break;
			case 187:	// Thr
				atomCount[1] += 7;
				atomCount[6] += 4;
				atomCount[7] += 1;
				atomCount[8] += 2;
				break;
			case 188:	// Trp
				atomCount[1] += 10;
				atomCount[6] += 11;
				atomCount[7] += 2;
				atomCount[8] += 1;
				break;
			case 189:	// Tyr
				atomCount[1] += 9;
				atomCount[6] += 9;
				atomCount[7] += 1;
				atomCount[8] += 2;
				break;
			case 190:	// Val
				atomCount[1] += 9;
				atomCount[6] += 5;
				atomCount[7] += 1;
				atomCount[8] += 1;
				break;
			case 1:
				switch (mol.getAtomMass(atom)) {
				case 0:	// natural abundance
				case 1:	// exclusively H1
					atomCount[1]++;
					break;
				case 2:
					atomCount[151]++;
					break;
				case 3:
					atomCount[152]++;
					break;
					}
				break;
			default:
				atomCount[mol.getAtomicNo(atom)]++;
				break;
				}
			}

		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			if (mol.getAtomicNo(atom) >= 171 && mol.getAtomicNo(atom) <= 190) // amino acids
				atomCount[1] += 2 - mol.getOccupiedValence(atom);
			else {
				atomCount[1] += mol.getImplicitHydrogens(atom);
			}
			// count number of present atomicNos and initialize arrays
		int atomicNoCount = 0;
		for (int i=1; i<=Molecule.cMaxAtomicNo; i++)
			if (atomCount[i] != 0)
				atomicNoCount++;
		mAtomCount = new int[atomicNoCount];
		mAtomicNo = new int[atomicNoCount];

			// list elements in sFirstFormulaAtomicNo[] first
		atomicNoCount = 0;
		for (int i=0; i<sFirstFormulaAtomicNo.length; i++) {
			if (atomCount[sFirstFormulaAtomicNo[i]] != 0) {
				mAtomCount[atomicNoCount] = atomCount[sFirstFormulaAtomicNo[i]];
				mAtomicNo[atomicNoCount] = sFirstFormulaAtomicNo[i];
				atomicNoCount++;
				atomCount[sFirstFormulaAtomicNo[i]] = 0;
				}
			}

			// list remaining elements in alphabetical order
		while (true) {
			String lowestLabel = "zzz";
			int lowestAtomicNo = -1;

			for (int atomicNo=1; atomicNo<=Molecule.cMaxAtomicNo; atomicNo++)
				if (atomCount[atomicNo] > 0
				 && lowestLabel.compareTo(Molecule.cAtomLabel[atomicNo]) > 0) {
						lowestLabel = Molecule.cAtomLabel[atomicNo];
						lowestAtomicNo = atomicNo;
						}

			if (lowestAtomicNo == -1)
				break;

			mAtomCount[atomicNoCount] = atomCount[lowestAtomicNo];
			mAtomicNo[atomicNoCount] = lowestAtomicNo;
			atomicNoCount++;
			atomCount[lowestAtomicNo] = 0;
			}

		mAbsoluteIsotopeWeightIncrement = 0.0;
		mRelativeIsotopeWeightIncrement = 0.0;
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.getAtomicNo(atom) != 1
             && !mol.isNaturalAbundance(atom)) {
				int atomicNo = mol.getAtomicNo(atom);
				int atomMass = mol.getAtomMass(atom);
				mAbsoluteIsotopeWeightIncrement += IsotopeHelper.getAbsoluteMass(atomicNo, atomMass)
										         - sAbsoluteMass[atomicNo];
				mRelativeIsotopeWeightIncrement += IsotopeHelper.getAbsoluteMass(atomicNo, atomMass)
                                                 - sRelativeMass[atomicNo];
				}
			}
		}

	public double getRelativeWeight() {
		double weight = mRelativeIsotopeWeightIncrement;
		for (int i=0; i<mAtomCount.length; i++)
			weight += (double)mAtomCount[i] * sRelativeMass[mAtomicNo[i]];
		return weight;
		}

	public double getAbsoluteWeight() {
		double weight = mAbsoluteIsotopeWeightIncrement;
		for (int i=0; i<mAtomCount.length; i++)
			weight += (double)mAtomCount[i] * sAbsoluteMass[mAtomicNo[i]];
		return weight;
		}

	public String getFormula() {
		StringBuffer formula = new StringBuffer();
		for (int i=0; i<mAtomCount.length; i++) {
			formula.append(Molecule.cAtomLabel[mAtomicNo[i]]);
			if (mAtomCount[i] > 1)
				formula.append(mAtomCount[i]);
			}
		return formula.toString();
		}

	@Override
	public boolean equals(Object f) {
		if (f == this)
			return true;
		if (!(f instanceof MolecularFormula))
			return false;

		for (int i=0; i<mAtomCount.length; i++)
			if (mAtomCount[i] != ((MolecularFormula)f).mAtomCount[i])
				return false;

		return true;
		}
	}
