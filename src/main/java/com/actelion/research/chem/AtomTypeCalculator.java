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


public class AtomTypeCalculator {
	public static final int cPropertiesAll						= 0x00003FBE;
	public static final int cPropertiesForSolubility			= 0x00000860;
	public static final int cPropertiesForCLogPCharges			= 0x00001861;
	public static final int cPropertiesForCLogP					= 0x00000861;
	public static final int cPropertiesForMutator               = 0x000028BE;
	public static final int cPropertiesBasicType				= 0x0000003C;
	
	public static final int cPropertiesAtomSmallRing			= 0x00000001;
									// distinguish between ring and non-ring atoms
	public static final int cPropertiesAtomRingSize				= 0x00000002;
									// consider also the ring size
	public static final int cPropertiesAtomAromatic				= 0x00000004;
									// consider also atom aromaticity
	public static final int cPropertiesAtomAllylic				= 0x00000008;
									// consider also whether the atom is in allylic/benzylic position
	public static final int cPropertiesAtomStabilized			= 0x00000010;
									// consider also whether the atom is stabilized by a neighbour carbonyl or similar

	public static final int cPropertiesConnBondOrder			= 0x00000020;
									// bond order (arom=0) to neighbours
	public static final int cPropertiesConnAtomTypeSimple		= 0x00000040;
									// encoded atomicNo of neighbours (simple encoding)
	public static final int cPropertiesConnAtomType				= 0x00000080;
									// encoded atomicNo of neighbours
	public static final int cPropertiesConnAtomNeighbours		= 0x00000100;
									// has neighbour further substituents
	public static final int cPropertiesConnAtomNeighboursExact	= 0x00000200;
									// number of neighbours' further substituents
	public static final int cPropertiesConnAtomSmallRing		= 0x00000400;
									// is neighbour in small ring
	public static final int cPropertiesConnAtomAromatic			= 0x00000800;
									// is neighbour aromatic
	public static final int cPropertiesAtomCharged			    = 0x00001000;
									// is atom charged and atom is ampholytic
	public static final int cPropertiesAtomRingBondCount	    = 0x00002000;
	// is atom charged


	private static final short cAtomicNoCode[] = {-1,
      -1,     -1,      0,      0,      1,      2,   //  H  ,He ,Li ,Be ,B  ,C  ,
       3,      4,      5,     -1,      0,      0,   //  N , O  ,F  ,Ne ,Na ,Mg ,
       0,      6,      7,      8,      9,     -1,   //  Al ,Si ,P  ,S  ,Cl ,Ar ,
       0,      0,     10,     10,     10,     10,   //  K  ,Ca ,Sc ,Ti ,V  ,Cr ,
      10,     10,     10,     10,     10,     10,   //  Mn ,Fe ,Co ,Ni ,Cu ,Zn ,
       1,     11,     11,     12,     13,     -1,   //  Ga ,Ge ,As ,Se ,Br ,Kr ,
       0,      0,     10,     10,     10,     10,   //  Rb ,Sr ,Y  ,Zr ,Nb ,Mo ,
      10,     10,     10,     10,     10,     10,   //  Tc ,Ru ,Rh ,Pd ,Ag ,Cd ,
       0,      0,      0,     11,     14,     -1,   //  In ,Sn ,Sb ,Te ,I  ,Xe ,
       0,      0,     15,     15,     15,     15,   //  Cs ,Ba ,La ,Ce ,Pr ,Nd ,
      15,     15,     15,     15,     15,     15,   //  Pm ,Sm ,Eu ,Gd ,Tb ,Dy ,
      15,     15,     15,     15,     15,     15,   //  Ho ,Er ,Tm ,Yb ,Lu ,Hf ,
      10,     10,     10,     10,     10,     10,   //  Ta ,W , Re ,Os ,Ir ,Pt ,
      10,     10,      1,      1,      1,      1,   //  Au ,Hg ,Tl ,Pb ,Bi ,Po ,
      -1,     -1,     -1,     -1,     15,     15,   //  At ,Rn ,Fr ,Ra ,Ac ,Th ,
      15,     15,     15,     15,     15,     15,   //  Pa ,U , Np ,Pu ,Am ,Cm ,
      15,     15,     15,     15,     15,     15,   //  Bk ,Cf ,Es ,Fm ,Md ,No ,
      15,     -1,     -1,     -1,     -1,     -1,   //  Lr ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,R1 ,R2 ,R3 ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  A  ,A1 ,A2 ,A3 ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  D  ,T  ,X  ,R  ,H2 ,H+
      -1,     -1,     -1,     -1,     -1,     -1,   //  Nnn,HYD,Pol,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,Ala,Arg,Asn,Asp,
      -1,     -1,     -1,     -1,     -1,     -1,   //  Cys,Gln,Glu,Gly,His,Ile,
      -1,     -1,     -1,     -1,     -1,     -1,   //  Leu,Lys,Met,Phe,Pro,Ser,
      -1,     -1,     -1,     -1 };                 //  Thr,Trp,Tyr,Val,

	private static final short cSimpleAtomicNoCode[] = {-1,
      -1,     -1,      0,      0,      0,      2,   //  H  ,He ,Li ,Be ,B  ,C  ,
       5,      5,      5,     -1,      0,      0,   //  N , O  ,F  ,Ne ,Na ,Mg ,
       0,      0,      9,      9,      9,     -1,   //  Al ,Si ,P  ,S  ,Cl ,Ar ,
       0,      0,      0,      0,      0,      0,   //  K  ,Ca ,Sc ,Ti ,V  ,Cr ,
       0,      0,      0,      0,      0,      0,   //  Mn ,Fe ,Co ,Ni ,Cu ,Zn ,
       0,      0,      0,      9,      9,     -1,   //  Ga ,Ge ,As ,Se ,Br ,Kr ,
       0,      0,      0,      0,      0,      0,   //  Rb ,Sr ,Y  ,Zr ,Nb ,Mo ,
       0,      0,      0,      0,      0,      0,   //  Tc ,Ru ,Rh ,Pd ,Ag ,Cd ,
       0,      0,      0,      0,      9,     -1,   //  In ,Sn ,Sb ,Te ,I  ,Xe ,
       0,      0,      0,      0,      0,      0,   //  Cs ,Ba ,La ,Ce ,Pr ,Nd ,
       0,      0,      0,      0,      0,      0,   //  Pm ,Sm ,Eu ,Gd ,Tb ,Dy ,
       0,      0,      0,      0,      0,      0,   //  Ho ,Er ,Tm ,Yb ,Lu ,Hf ,
       0,      0,      0,      0,      0,      0,   //  Ta ,W , Re ,Os ,Ir ,Pt ,
       0,      0,      0,      0,      0,      0,   //  Au ,Hg ,Tl ,Pb ,Bi ,Po ,
      -1,     -1,     -1,     -1,      0,      0,   //  At ,Rn ,Fr ,Ra ,Ac ,Th ,
       0,      0,      0,      0,      0,      0,   //  Pa ,U , Np ,Pu ,Am ,Cm ,
       0,      0,      0,      0,      0,      0,   //  Bk ,Cf ,Es ,Fm ,Md ,No ,
       0,     -1,     -1,     -1,     -1,     -1,   //  Lr ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,R1 ,R2 ,R3 ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  A  ,A1 ,A2 ,A3 ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  D  ,T  ,X  ,R  ,H2 ,H+
      -1,     -1,     -1,     -1,     -1,     -1,   //  Nnn,HYD,Pol,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
      -1,     -1,     -1,     -1,     -1,     -1,   //  ?? ,?? ,Ala,Arg,Asn,Asp,
      -1,     -1,     -1,     -1,     -1,     -1,   //  Cys,Gln,Glu,Gly,His,Ile,
      -1,     -1,     -1,     -1,     -1,     -1,   //  Leu,Lys,Met,Phe,Pro,Ser,
      -1,     -1,     -1,     -1 };                 //  Thr,Trp,Tyr,Val,

	private static final String cAtomicNoCodeString[] = { 
		"MainGroupMetal",
		"Boron",
		"Carbon",
		"Nitrogen",
		"Oxygen",
		"Fluor",
		"Silicon",
		"Phosphorous",
		"Sulfur",
		"Chlorine",
		"Transition Metal",
		"MainGroupNonMetal",
		"Selene",
		"Bromine",
		"Iodine",
		"LanthanideOrActinide" };

	private static final String cSimpleAtomicNoCodeString[] = {
			"Metal",
			"-",
			"Carbon",
			"-",
			"-",
			"Small Hetero",
			"-",
			"-",
			"-",
			"Large Hetero",
			"-",
			"-",
			"-",
			"-",
			"-",
			"-" };


	public static String getAtomicNoCodeString(int atomicNoCode) {
		return cAtomicNoCodeString[atomicNoCode];
		}

	public static String getSimpleAtomicNoCodeString(int atomicNoCode) {
		return cSimpleAtomicNoCodeString[atomicNoCode];
		}

	public static void printAtomType(StereoMolecule mol, int atom) {
		try {
			printAtomType(getAtomType(mol, atom));
			}
		catch (Exception e) {
			System.out.println(e);
			}
		}

	public static String getHeaderString(int mode) {
		StringBuffer sb = new StringBuffer();
		sb.append("AtomicNo-Code");

		if ((mode & cPropertiesAtomRingSize) != 0)
			sb.append("\tRingsize");

		if ((mode & cPropertiesAtomSmallRing) != 0)
			sb.append("\tRingmember");

		if ((mode & cPropertiesAtomAromatic) != 0)
			sb.append("\tAromatic");

		if ((mode & cPropertiesAtomAllylic) != 0)
			sb.append("\tAllylic");

		if ((mode & cPropertiesAtomStabilized) != 0)
			sb.append("\tStabilized");

		if ((mode & cPropertiesAtomCharged) != 0)
			sb.append("\tCharged\tAmpholytic");

		if ((mode & cPropertiesAtomRingBondCount) != 0)
			sb.append("\tRingBonds");

		for (int i=1; i<=4; i++) {
			if ((mode & cPropertiesConnBondOrder) != 0)
				sb.append("\tConnBondType"+i);

			if ((mode & cPropertiesConnAtomType) != 0)
				sb.append("\tConnAtomType"+i);
			else if ((mode & cPropertiesConnAtomTypeSimple) != 0)
				sb.append("\tSimpleConnAtomType"+i);

			if ((mode & cPropertiesConnAtomNeighbours) != 0)
				sb.append("\tConnMoreNeighbours"+i);

			if ((mode & cPropertiesConnAtomSmallRing) != 0)
				sb.append("\tConnIsSmallRingMember"+i);

			if ((mode & cPropertiesConnAtomAromatic) != 0)
				sb.append("\tConnIsAromatic"+i);
			}

		return sb.toString();
		}

	public static String getTypeString(long type, int mode) {
		StringBuffer sb = new StringBuffer();
		sb.append(getAtomicNoCodeString((int)(type & 15)));

		if ((mode & cPropertiesAtomRingSize) != 0) {
			long ringSize = ((type & 112) >> 4);
			if (ringSize != 0) ringSize += 2;
			sb.append("\t"+ringSize);
			}
		if ((mode & cPropertiesAtomSmallRing) != 0)
			sb.append("\t"+(((type & 64) != 0) ? "yes" : "no"));

		if ((mode & cPropertiesAtomAromatic) != 0)
			sb.append("\t"+(((type & 128) != 0) ? "yes" : "no"));

		if ((mode & cPropertiesAtomAllylic) != 0)
			sb.append("\t"+(((type & 256) != 0) ? "yes" : "no"));

		if ((mode & cPropertiesAtomStabilized) != 0)
			sb.append("\t"+(((type & 512) != 0) ? "yes" : "no"));

		if ((mode & cPropertiesAtomCharged) != 0) {
			sb.append("\t" + (((type & 0x0004000000000000L) != 0) ? "yes" : "no"));
			sb.append("\t" + (((type & 0x0008000000000000L) != 0) ? "yes" : "no"));
			}

		if ((mode & cPropertiesAtomRingBondCount) != 0) {
			long ringBondCount = (type & 0x0030000000000000L) >> 52;
			if (ringBondCount != 0)
				ringBondCount++;
			sb.append("\t" + ringBondCount);
			}

		for (int i=1; i<=4; i++) {
			type >>= 10;
			long neighbourType = type & 1023;
			if (neighbourType != 0) {
				if ((mode & cPropertiesConnBondOrder) != 0)
					sb.append("\t"+(neighbourType >> 8));

				if ((mode & cPropertiesConnAtomType) != 0)
					sb.append("\t"+getAtomicNoCodeString((int)(neighbourType & 15)));
				else if ((mode & cPropertiesConnAtomTypeSimple) != 0)
					sb.append("\t"+getSimpleAtomicNoCodeString((int)(neighbourType & 15)));

				if ((mode & cPropertiesConnAtomNeighbours) != 0) {
					long otherNeighbours = 3 & (neighbourType >> 4);
					if ((mode & cPropertiesConnAtomNeighboursExact) == 0)
						sb.append("\t"+(otherNeighbours == 0 ? "no" : "yes"));
					else
						sb.append("\t"+otherNeighbours);
					}

				if ((mode & cPropertiesConnAtomSmallRing) != 0)
					sb.append("\t"+(((neighbourType & 64) != 0) ? "yes" : "no"));

				if ((mode & cPropertiesConnAtomAromatic) != 0)
					sb.append("\t"+(((neighbourType & 128) != 0) ? "yes" : "no"));
				}
			}

		return sb.toString();
		}


	public static String toString(long type) {
		return toString(type, cPropertiesAll);
		}

	public static String toString(long type, int mode) {
		StringBuffer sb = new StringBuffer();
		sb.append(getAtomicNoCodeString((int)(type & 15))+":");

		if ((mode & cPropertiesAtomCharged) != 0) {
			if((type & 0x0004000000000000L) !=0 )
				sb.append("Chg");
			if((type & 0x0008000000000000L) !=0 )
				sb.append("Amp");
			}

		if ((mode & cPropertiesAtomRingBondCount) != 0) {
			long ringBondCount = (type & 0x0030000000000000L) >> 52;
			if (ringBondCount != 0)
				ringBondCount++;
			sb.append("Rb"+ringBondCount);
			}

		if ((mode & cPropertiesAtomRingSize) != 0) {
			long ringSize = ((type & 112) >> 4);
			if (ringSize != 0) ringSize += 2;
			sb.append("R" + ringSize);
			}

		if ((mode & cPropertiesAtomAromatic) != 0 && (type & 128) != 0)
			sb.append("Ar");

		if ((mode & cPropertiesAtomAllylic) != 0 && (type & 256) != 0)
			sb.append("Al");

		if ((mode & cPropertiesAtomStabilized) != 0 && (type & 512) != 0)
			sb.append("St");

		for (int j=0; j<4; j++) {
			type >>= 10;
			if (type == 0)
				break;

			if ((mode & cPropertiesConnBondOrder) != 0) {
				int bondType = (int) (3 & (type >> 8));
				switch (bondType) {
					case 0:
						sb.append(" *");
						break;
					case 1:
						sb.append(" -");
						break;
					case 2:
						sb.append(" =");
						break;
					case 3:
						sb.append(" #");
						break;
					}
				}

			sb.append("{");
			if ((mode & cPropertiesConnAtomType) != 0)
				sb.append(getAtomicNoCodeString((int)(type & 15))+":");
			else if ((mode & cPropertiesConnAtomTypeSimple) != 0)
				sb.append(getSimpleAtomicNoCodeString((int)(type & 15))+":");

			long neighbours = ((type & 48) >> 4) + 1;
			if ((mode & cPropertiesConnAtomNeighbours) != 0)
				sb.append("N" + neighbours);

			if ((mode & cPropertiesConnAtomSmallRing) != 0 && (type & 64) != 0)
				sb.append("Ri");

			if ((mode & cPropertiesConnAtomAromatic) != 0 && (type & 128) != 0)
				sb.append("Ar");

			sb.append("}");
			}
		
		return sb.toString();
		}

	public static void printAtomType(long type) {
		System.out.println(toString(type));
		}


	public static long getAtomType(StereoMolecule mol, int atom) throws Exception {
		return getAtomType(mol, atom, cPropertiesAll);
		}


	public static long getAtomType(StereoMolecule mol, int atom, int mode) throws Exception {
		mol.ensureHelperArrays(Molecule.cHelperRings);

		long[] neighbourType = new long[mol.getConnAtoms(atom)];
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			long connAtomType = 0;
			long connBondType = 0;

			if ((mode & cPropertiesConnBondOrder) != 0) {
				long connBondOrder = mol.getConnBondOrder(atom, i);
				if (connBondOrder < 3 && mol.isAromaticBond(mol.getConnBond(atom, i)))
					connBondOrder = 0;
				connBondType += connBondOrder;
				}
				
			int connAtom = mol.getConnAtom(atom, i);

			if ((mode & cPropertiesConnAtomType) != 0) {
				if (cAtomicNoCode[mol.getAtomicNo(connAtom)] == -1)
					throw new Exception("unsupported atomicNo:"+mol.getAtomicNo(connAtom));
				connAtomType += cAtomicNoCode[mol.getAtomicNo(connAtom)];
				}
			else if ((mode & cPropertiesConnAtomTypeSimple) != 0) {
				if (cSimpleAtomicNoCode[mol.getAtomicNo(connAtom)] == -1)
					throw new Exception("unsupported atomicNo:"+mol.getAtomicNo(connAtom));
				connAtomType += cSimpleAtomicNoCode[mol.getAtomicNo(connAtom)];
				}

			if ((mode & cPropertiesConnAtomNeighbours) != 0) {
				int otherNeighbours = mol.getConnAtoms(connAtom) - 1;
				if (otherNeighbours > 3)
					otherNeighbours = 3;

				if ((mode & cPropertiesConnAtomNeighboursExact) == 0)
					if (otherNeighbours > 1)
						otherNeighbours = 1;

				connAtomType |= (otherNeighbours << 4);
				}

			if ((mode & cPropertiesConnAtomSmallRing) != 0)
				if (mol.isSmallRingAtom(connAtom))
					connAtomType |= 64;

			if ((mode & cPropertiesConnAtomAromatic) != 0)
				if (mol.isAromaticAtom(connAtom))
					connAtomType |= 128;

			long theType = connAtomType | (connBondType  << 8);

			int index=0;
			while (theType < neighbourType[index])
				index++;

			for (int j=i; j>index; j--)
				neighbourType[j] = neighbourType[j-1];

			neighbourType[index] = theType;
		}

		int neighbours = (mol.getConnAtoms(atom) < 4) ? mol.getConnAtoms(atom) : 4;
		long atomType = 0;
		for (int i=0; i<neighbours; i++) {
			atomType <<= 10;
			atomType |= neighbourType[i];
		}

		atomType <<= 10;
		if (cAtomicNoCode[mol.getAtomicNo(atom)] == -1)
			throw new Exception("unsupported atomicNo:"+mol.getAtomicNo(atom));
  		atomType |= cAtomicNoCode[mol.getAtomicNo(atom)];

		if ((mode & cPropertiesAtomRingSize) != 0) {
			int ringSize = mol.getAtomRingSize(atom);
			if (ringSize > 9)
				ringSize = 9;
			if (ringSize > 2)
				ringSize -= 2;
	  		atomType |= (ringSize << 4);
			}
		else if ((mode & cPropertiesAtomSmallRing) != 0)
			if (mol.isSmallRingAtom(atom))
				atomType |= 64;

		if ((mode & cPropertiesAtomAromatic) != 0)
			if (mol.isAromaticAtom(atom))
				atomType |= 128;

		if ((mode & cPropertiesAtomAllylic) != 0)
			if (mol.isAllylicAtom(atom))
				atomType |= 256;

		if ((mode & cPropertiesAtomStabilized) != 0)
			if (mol.isStabilizedAtom(atom))
				atomType |= 512;

		if ((mode & cPropertiesAtomCharged) != 0) {
			if(AtomFunctionAnalyzer.hasUnbalancedAtomCharge(mol, atom))
				atomType |= 0x0004000000000000L;
						
			if(AtomFunctionAnalyzer.isBasicNitrogen(mol, atom)){
				for (int i = 0; i < mol.getAtoms(); i++) {
					if(AtomFunctionAnalyzer.isAcidicOxygen(mol, i)){
						atomType |= 0x0008000000000000L;	// ampholytic
						break;
						}
					}
				}
			}

		if ((mode & cPropertiesAtomRingBondCount) != 0) {
			long ringBondCount = mol.getAtomRingBondCount(atom);
			if (ringBondCount != 0)
				ringBondCount--;
			atomType |= (ringBondCount << 52);	// ampholytic
			}

			return atomType;
		}
	}