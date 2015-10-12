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

package com.actelion.research.chem;


public class AtomTypeCalculator {
	public static final int cPropertiesAll						= 0x00001FBE;
	public static final int cPropertiesForSolubility			= 0x00000860;
	public static final int cPropertiesForCLogPCharges			= 0x00001861;
	public static final int cPropertiesForCLogP					= 0x00000861;
    public static final int cPropertiesForMutator               = 0x0000087E;
	public static final int cPropertiesBasicType				= 0x0000003C;
	
	public static final int cPropertiesAtomSmallRing			= 0x00000001;
									// distinguish between ring and non-ring atoms
	public static final int cPropertiesAtomRingSize				= 0x00000002;
									// consider also the ring size
	public static final int cPropertiesAtomAromatic				= 0x00000004;
									// consider also the ring size
	public static final int cPropertiesAtomAllylic				= 0x00000008;
									// consider also the ring size
	public static final int cPropertiesAtomStabilized			= 0x00000010;
									// consider also the ring size

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
      10,      0,      0,      0,      0,      0,   //  Tc ,Ru ,Rh ,Pd ,Ag ,Cd ,
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

	
	public static String getAtomicNoCodeString(int atomicNoCode) {
		return cAtomicNoCodeString[atomicNoCode];
		}


	public static void printAtomType(StereoMolecule mol, int atom) {
		try {
			printAtomType(getAtomType(mol, atom));
			}
		catch (Exception e) {
			System.out.println(e);
			}
		}

	public static String toString(long type) {
		StringBuffer sb = new StringBuffer();
		sb.append(getAtomicNoCodeString((int)(type & 15))+":");
		
		if((type & 0x0004000000000000L)!=0){
			sb.append("(charged)");
		}
		if((type & 0x0008000000000000L)!=0){
			sb.append("(ampholytic)");
		}
		
		long ringSize = ((type & 112) >> 4);
		if (ringSize != 0) ringSize += 2;
		sb.append("R"+ringSize);
		if ((type & 128) != 0) sb.append("Ar");
		if ((type & 256) != 0) sb.append("Al");
		if ((type & 512) != 0) sb.append("St");
		type >>= 10;
		for (int j=0; j<4; j++) {
			if (type == 0)
				break;

			int bondType = (int)(3 & (type >> 8));
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

			sb.append("{");
			sb.append(AtomTypeCalculator.getAtomicNoCodeString((int)(type & 15))+":");
			long neighbours = ((type & 48) >> 4) + 1;
			sb.append("N"+neighbours);
			if ((type & 64) != 0) sb.append("Ri");
			if ((type & 128) != 0) sb.append("Ar");

			sb.append("}");

			type >>= 10;
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

				connAtomType += (otherNeighbours << 4);
				}

			if ((mode & cPropertiesConnAtomSmallRing) != 0)
				if (mol.isSmallRingAtom(connAtom))
					connAtomType += 64;

			if ((mode & cPropertiesConnAtomAromatic) != 0)
				if (mol.isAromaticAtom(connAtom))
					connAtomType += 128;

			long theType = connAtomType + (connBondType  << 8);

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
			atomType += neighbourType[i];
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
				atomType += 128;

		if ((mode & cPropertiesAtomAllylic) != 0)
			if (mol.isAllylicAtom(atom))
				atomType += 256;

		if ((mode & cPropertiesAtomStabilized) != 0)
			if (mol.isStabilizedAtom(atom))
				atomType += 512;
		
		
		if((atomType & 0x0004000000000000L)!=0){
			RuntimeException ex = new RuntimeException("Bit already set!");
			ex.printStackTrace();
		}
		
		if((atomType & 0x0008000000000000L)!=0){
			RuntimeException ex = new RuntimeException("Bit already set!");
			ex.printStackTrace();
		}
		
		
		if ((mode & cPropertiesAtomCharged) != 0) {
			
			if(AtomFunctionAnalyzer.hasUnbalancedAtomCharge(mol, atom))
				atomType += 0x0004000000000000L;
						
			boolean ampholytic=false;
			
			if(AtomFunctionAnalyzer.isBasicNitrogen(mol, atom)){
				for (int i = 0; i < mol.getAtoms(); i++) {
					if(AtomFunctionAnalyzer.isAcidicOxygen(mol, i)){
						ampholytic=true;
						break;
					}
				}
			} 
//			else if(AtomFunctionAnalyzer.isAcidicOxygen(mol, atom)){
//				for (int i = 0; i < mol.getAtoms(); i++) {
//					if(AtomFunctionAnalyzer.isBasicNitrogen(mol, i)){
//						ampholytic=true;
//						break;
//					}
//				}
//			}
			
			if(ampholytic)
				atomType += 0x0008000000000000L;
			
		}

		return atomType;
	}
}//end_of_class