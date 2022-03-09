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

package com.actelion.research.chem;


public class AtomTypeCalculator {
	public static final int cPropertiesAll						= 0x00007FBE;
	public static final int cPropertiesForSolubility			= 0x00000860;
	public static final int cPropertiesForCLogPCharges			= 0x00001861;
	public static final int cPropertiesForCLogP					= 0x00000861;
	public static final int cPropertiesForMutator               = 0x00007DBE;
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
	public static final int cPropertiesAtomCharged			    = 0x00001000;
									// is atom charged and atom is ampholytic
	public static final int cPropertiesAtomRingCount		    = 0x00002000;
									// ring count of atoms

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
									// is small ring neighbour atom
	public static final int cPropertiesConnAtomAromatic			= 0x00000800;
									// is aromatic neighbour atom
	public static final int cPropertiesConnAtomStabilized		= 0x00004000;
									// is stabilized neighbour atom

	private static final long ATOM_FLAG_COUNT = 15;
	private static final long CONN_FLAG_COUNT = 11;

	private static final long ATOM_FLAGS_ATOMICNO  = 0x0000000F;
	private static final long ATOM_FLAGS_RINGSIZE  = 0x00000070;
	private static final long ATOM_SHIFT_RINGSIZE  = 4;
	private static final long ATOM_FLAGS_RINGCOUNT = 0x00000380;
	private static final long ATOM_SHIFT_RINGCOUNT = 7;
	private static final long ATOM_FLAG_SMALLRING  = 0x00000040;	// uses one of the ringsize flags
	private static final long ATOM_FLAG_AROMATIC   = 0x00000400;
	private static final long ATOM_FLAG_ALLYLIC    = 0x00000800;
	private static final long ATOM_FLAG_STABILIZED = 0x00001000;
	private static final long ATOM_FLAG_CHARGED    = 0x00002000;
	private static final long ATOM_FLAG_AMPHOLYTIC = 0x00004000;

	private static final long CONN_FLAGS_ALL        = 0x000007FF;
	private static final long CONN_FLAGS_ATOMICNO   = 0x0000000F;
	private static final long CONN_FLAGS_BONDORDER  = 0x00000030;
	private static final long CONN_SHIFT_BONDORDER  = 4;
	private static final long CONN_FLAGS_NEIGHBOURS = 0x000000C0;
	private static final long CONN_SHIFT_NEIGHBOURS = 6;
	private static final long CONN_FLAG_SMALLRING   = 0x00000100;
	private static final long CONN_FLAG_AROMATIC    = 0x00000200;
	private static final long CONN_FLAG_STABILIZED  = 0x00000400;

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

		if ((mode & cPropertiesAtomRingCount) != 0)
			sb.append("\tRingCount");

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

			if ((mode & cPropertiesConnAtomStabilized) != 0)
				sb.append("\tConnIsStabilized"+i);
			}

		return sb.toString();
		}

	public static String getTypeString(long type, int mode) {
		StringBuffer sb = new StringBuffer();
		sb.append(getAtomicNoCodeString((int)(type & ATOM_FLAGS_ATOMICNO)));

		if ((mode & cPropertiesAtomRingSize) != 0) {
			long ringSize = ((type & ATOM_FLAGS_RINGSIZE) >> ATOM_SHIFT_RINGSIZE);
			if (ringSize != 0) ringSize += 2;
			sb.append("\t"+ringSize);
			}
		if ((mode & cPropertiesAtomSmallRing) != 0)
			sb.append("\t"+(((type & ATOM_FLAG_SMALLRING) != 0) ? "yes" : "no"));

		if ((mode & cPropertiesAtomAromatic) != 0)
			sb.append("\t"+(((type & ATOM_FLAG_AROMATIC) != 0) ? "yes" : "no"));

		if ((mode & cPropertiesAtomAllylic) != 0)
			sb.append("\t"+(((type & ATOM_FLAG_ALLYLIC) != 0) ? "yes" : "no"));

		if ((mode & cPropertiesAtomStabilized) != 0)
			sb.append("\t"+(((type & ATOM_FLAG_STABILIZED) != 0) ? "yes" : "no"));

		if ((mode & cPropertiesAtomCharged) != 0) {
			sb.append("\t" + (((type & ATOM_FLAG_CHARGED) != 0) ? "yes" : "no"));
			sb.append("\t" + (((type & ATOM_FLAG_AMPHOLYTIC) != 0) ? "yes" : "no"));
			}

		if ((mode & cPropertiesAtomRingCount) != 0) {
			long ringCount = (type & ATOM_FLAGS_RINGCOUNT) >> ATOM_SHIFT_RINGCOUNT;
			sb.append("\t" + ringCount);
			}

		type >>= ATOM_FLAG_COUNT;
		for (int i=1; i<=4; i++) {
			long neighbourType = type & CONN_FLAGS_ALL;
			if (neighbourType != 0) {
				if ((mode & cPropertiesConnBondOrder) != 0)
					sb.append("\t"+((neighbourType & CONN_FLAGS_BONDORDER) >> CONN_SHIFT_BONDORDER));

				if ((mode & cPropertiesConnAtomType) != 0)
					sb.append("\t"+getAtomicNoCodeString((int)(neighbourType & CONN_FLAGS_ATOMICNO)));
				else if ((mode & cPropertiesConnAtomTypeSimple) != 0)
					sb.append("\t"+getSimpleAtomicNoCodeString((int)(neighbourType & CONN_FLAGS_ATOMICNO)));

				if ((mode & cPropertiesConnAtomNeighbours) != 0) {
					long otherNeighbours = ((neighbourType & CONN_FLAGS_NEIGHBOURS) >> CONN_SHIFT_NEIGHBOURS);
					if ((mode & cPropertiesConnAtomNeighboursExact) == 0)
						sb.append("\t"+(otherNeighbours == 0 ? "no" : "yes"));
					else
						sb.append("\t"+otherNeighbours);
					}

				if ((mode & cPropertiesConnAtomSmallRing) != 0)
					sb.append("\t"+(((neighbourType & CONN_FLAG_SMALLRING) != 0) ? "yes" : "no"));

				if ((mode & cPropertiesConnAtomAromatic) != 0)
					sb.append("\t"+(((neighbourType & CONN_FLAG_AROMATIC) != 0) ? "yes" : "no"));

				if ((mode & cPropertiesConnAtomStabilized) != 0)
					sb.append("\t"+(((neighbourType & CONN_FLAG_STABILIZED) != 0) ? "yes" : "no"));
				}
			type >>= CONN_FLAG_COUNT;
			}

		return sb.toString();
		}


	public static String toString(long type) {
		return toString(type, cPropertiesAll);
		}

	public static String toString(long type, int mode) {
		StringBuffer sb = new StringBuffer();
		sb.append(getAtomicNoCodeString((int)(type & ATOM_FLAGS_ATOMICNO))+":");

		if ((mode & cPropertiesAtomCharged) != 0) {
			if((type & ATOM_FLAG_CHARGED) !=0 )
				sb.append("Chg");
			if((type & ATOM_FLAG_AMPHOLYTIC) !=0 )
				sb.append("Amp");
			}

		if ((mode & cPropertiesAtomRingCount) != 0) {
			long ringCount = (type & ATOM_FLAGS_RINGCOUNT) >> ATOM_SHIFT_RINGCOUNT;
			sb.append("Rc"+ringCount);
			}

		if ((mode & cPropertiesAtomRingSize) != 0) {
			long ringSize = ((type & ATOM_FLAGS_RINGSIZE) >> ATOM_SHIFT_RINGSIZE);
			if (ringSize != 0) ringSize += 2;
			sb.append("Rs" + ringSize);
			}

		if ((mode & cPropertiesAtomAromatic) != 0 && (type & ATOM_FLAG_AROMATIC) != 0)
			sb.append("Ar");

		if ((mode & cPropertiesAtomAllylic) != 0 && (type & ATOM_FLAG_ALLYLIC) != 0)
			sb.append("Al");

		if ((mode & cPropertiesAtomStabilized) != 0 && (type & ATOM_FLAG_STABILIZED) != 0)
			sb.append("St");

		type >>= ATOM_FLAG_COUNT;
		for (int j=0; j<4; j++) {
			if (type == 0)
				break;

			if ((mode & cPropertiesConnBondOrder) != 0) {
				int bondType = (int) ((type & CONN_FLAGS_BONDORDER) >> CONN_SHIFT_BONDORDER);
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
				sb.append(getAtomicNoCodeString((int)(type & CONN_FLAGS_ATOMICNO))+":");
			else if ((mode & cPropertiesConnAtomTypeSimple) != 0)
				sb.append(getSimpleAtomicNoCodeString((int)(type & CONN_FLAGS_ATOMICNO))+":");

			long neighbours = ((type & CONN_FLAGS_NEIGHBOURS) >> CONN_SHIFT_NEIGHBOURS) + 1;
			if ((mode & cPropertiesConnAtomNeighbours) != 0)
				sb.append("N" + neighbours);

			if ((mode & cPropertiesConnAtomSmallRing) != 0 && (type & CONN_FLAG_SMALLRING) != 0)
				sb.append("Ri");

			if ((mode & cPropertiesConnAtomAromatic) != 0 && (type & CONN_FLAG_AROMATIC) != 0)
				sb.append("Ar");

			if ((mode & cPropertiesConnAtomStabilized) != 0 && (type & CONN_FLAG_STABILIZED) != 0)
				sb.append("St");

			sb.append("}");

			type >>= CONN_FLAG_COUNT;
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
		int neighbourCount = 0;

		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int connAtom = mol.getConnAtom(atom, i);

			if (mol.getAtomicNo(connAtom) == 1)	// skip explicit hydrogen
				continue;

			long connType = 0;

			if ((mode & cPropertiesConnBondOrder) != 0) {
				long connBondOrder = mol.getConnBondOrder(atom, i);
				if (mode == cPropertiesForMutator) {	// for the Mutator this distinction is more appropriate
					if (connBondOrder < 3 && mol.isDelocalizedBond(mol.getConnBond(atom, i)) && mol.getAtomPi(atom) == 1)
						connBondOrder = 0;
					}
				else {
					if (connBondOrder < 3 && mol.isAromaticBond(mol.getConnBond(atom, i)))
						connBondOrder = 0;
				}

				connType |= (connBondOrder << CONN_SHIFT_BONDORDER);
				}
				
			if ((mode & cPropertiesConnAtomType) != 0) {
				if (cAtomicNoCode[mol.getAtomicNo(connAtom)] == -1)
					throw new Exception("unsupported atomicNo:"+mol.getAtomicNo(connAtom));
				connType += cAtomicNoCode[mol.getAtomicNo(connAtom)];
				}
			else if ((mode & cPropertiesConnAtomTypeSimple) != 0) {
				if (cSimpleAtomicNoCode[mol.getAtomicNo(connAtom)] == -1)
					throw new Exception("unsupported atomicNo:"+mol.getAtomicNo(connAtom));
				connType += cSimpleAtomicNoCode[mol.getAtomicNo(connAtom)];
				}

			if ((mode & cPropertiesConnAtomNeighbours) != 0) {
				int otherNeighbours = mol.getConnAtoms(connAtom) - 1;
				if (otherNeighbours > 3)
					otherNeighbours = 3;

				if ((mode & cPropertiesConnAtomNeighboursExact) == 0)
					if (otherNeighbours > 1)
						otherNeighbours = 1;

				connType |= (otherNeighbours << CONN_SHIFT_NEIGHBOURS);
				}

			if ((mode & cPropertiesConnAtomSmallRing) != 0)
				if (mol.isSmallRingAtom(connAtom))
					connType |= CONN_FLAG_SMALLRING;

			if ((mode & cPropertiesConnAtomAromatic) != 0)
				if (mol.isAromaticAtom(connAtom))
					connType |= CONN_FLAG_AROMATIC;

			if ((mode & cPropertiesConnAtomStabilized) != 0)
				if (mol.isStabilizedAtom(connAtom))
					connType |= CONN_FLAG_STABILIZED;

			int index=0;
			while (connType < neighbourType[index])
				index++;

			for (int j=i; j>index; j--)
				neighbourType[j] = neighbourType[j-1];

			neighbourType[index] = connType;

			neighbourCount++;
			}

		if (neighbourCount > 4)
			neighbourCount = 4;

		long atomType = 0;
		for (int i=0; i<neighbourCount; i++) {
			atomType <<= CONN_FLAG_COUNT;
			atomType |= neighbourType[i];
			}

		atomType <<= ATOM_FLAG_COUNT;
		if (cAtomicNoCode[mol.getAtomicNo(atom)] == -1)
			throw new Exception("unsupported atomicNo:"+mol.getAtomicNo(atom));
  		atomType |= cAtomicNoCode[mol.getAtomicNo(atom)];

		if ((mode & cPropertiesAtomRingSize) != 0) {
			int ringSize = mol.getAtomRingSize(atom);
			if (ringSize > 9)
				ringSize = 9;
			if (ringSize > 2)
				ringSize -= 2;
	  		atomType |= (ringSize << ATOM_SHIFT_RINGSIZE);
			}
		else if ((mode & cPropertiesAtomSmallRing) != 0)
			if (mol.isSmallRingAtom(atom))
				atomType |= ATOM_FLAG_SMALLRING;

		if ((mode & cPropertiesAtomAromatic) != 0)
			if (mol.isAromaticAtom(atom))
				atomType |= ATOM_FLAG_AROMATIC;

		if ((mode & cPropertiesAtomAllylic) != 0)
			if (mol.isAllylicAtom(atom))
				atomType |= ATOM_FLAG_ALLYLIC;

		if ((mode & cPropertiesAtomStabilized) != 0)
			if (mol.isStabilizedAtom(atom))
				atomType |= ATOM_FLAG_STABILIZED;

		if ((mode & cPropertiesAtomCharged) != 0) {
			if(AtomFunctionAnalyzer.hasUnbalancedAtomCharge(mol, atom))
				atomType |= ATOM_FLAG_CHARGED;
						
			if(AtomFunctionAnalyzer.isBasicNitrogen(mol, atom)){
				for (int i = 0; i < mol.getAtoms(); i++) {
					if(AtomFunctionAnalyzer.isAcidicOxygen(mol, i)){
						atomType |= ATOM_FLAG_AMPHOLYTIC;	// ampholytic
						break;
						}
					}
				}
			}

		if ((mode & cPropertiesAtomRingCount) != 0) {
			long ringCount = mol.getAtomRingCount(atom, 10);
			atomType |= (ringCount << ATOM_SHIFT_RINGCOUNT);	// ampholytic
			}

		return atomType;
		}
	}