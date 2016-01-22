/*
* Copyright (c) 1997 - 2015
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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Arrays;


public class Molecule implements Serializable {
	// We need the serialVersionUID in order to make sure, that after obfuscation the Molecule keeps
	// the same serialVersionUID since the inter-JVM drag & drop checks the Serializable version ID
	// which is either generated at via hash or can explicitly set via serialVersionUID. Due to the
	// fact, that the obfuscation process alters member/type names the generated hash is going to
	// be changed as well.
	// In addition to the above requirement, the class name should not be obfuscated at all!

	static final long serialVersionUID = 0x20100310;	// after splitting bond flags and query features

	public static final int cMaxAtomicNo = 190;

		// parity based on atom positions in atom table (as MDL parity)
	private static final int cAtomFlagsParity		= 0x000003;
	public static final int cAtomParityNone			= 0x000000;
	public static final int cAtomParity1			= 0x000001;
	public static final int cAtomParity2			= 0x000002;
	public static final int cAtomParityUnknown	  = 0x000003;

	public static final int cAtomParityIsPseudo	 = 0x000004;
	protected static final int cAtomFlagSmallRing	= 0x000008;

	public static final int cAtomRadicalState		= 0x000030;
	public static final int cAtomRadicalStateShift	= 4;
	public static final int cAtomRadicalStateNone	= 0x000000;
	public static final int cAtomRadicalStateS		= 0x000010;
	public static final int cAtomRadicalStateD		= 0x000020;
	public static final int cAtomRadicalStateT		= 0x000030;

	private static final int cAtomFlagsColor		= 0x0001C0;
	public static final int cAtomColorNone = 0x000000;
	public static final int cAtomColorBlue			= 0x000040;
	public static final int cAtomColorRed			= 0x000080;
	public static final int cAtomColorGreen			= 0x0000C0;
	public static final int cAtomColorMagenta		= 0x000100;
	public static final int cAtomColorOrange		= 0x000140;
	public static final int cAtomColorDarkGreen		= 0x000180;
	public static final int cAtomColorDarkRed		= 0x0001C0;
	private static final int cAtomFlagSelected 		= 0x000200;

	protected static final int cAtomFlagsHelper		= 0x0003FC0F;
	protected static final int cAtomFlagsHelper2	= 0x00007C08;
	protected static final int cAtomFlagsHelper3	= 0x08038007;

	protected static final int cAtomFlagsRingBonds	= 0x000C00;
	protected static final int cAtomFlags2RingBonds = 0x000400;
	protected static final int cAtomFlags3RingBonds = 0x000800;
	protected static final int cAtomFlags4RingBonds = 0x000C00;
	protected static final int cAtomFlagAromatic	= 0x001000;
	protected static final int cAtomFlagAllylic		= 0x002000;
	protected static final int cAtomFlagStabilized	= 0x004000;

	private static final int cAtomFlagsCIPParity	= 0x018000;
	private static final int cAtomFlagsCIPParityShift = 15;
	public static final int cAtomCIPParityNone	  = 0x000000;
	public static final int cAtomCIPParityRorM	  = 0x000001;
	public static final int cAtomCIPParitySorP		= 0x000002;
	public static final int cAtomCIPParityProblem   = 0x000003;

	protected static final int cAtomFlagStereoProblem = 0x020000;
	protected static final int cAtomFlagMarked		= 0x040000;

		// MDL's enhanced stereochemical representation (ESR group and type may be assigned
		// to TH and allene stereo centers as well as to BINAP kind of stereo bonds)
	public static final int cESRTypeAbs				= 0;
	public static final int cESRTypeAnd				= 1;
	public static final int cESRTypeOr				= 2;
	public static final int cESRMaxGroups			= 32;
	public static final int cESRGroupBits			= 5;

	protected static final int cAtomFlagsESR		= 0x03F80000;
	private static final int cAtomFlagsESRType		= 0x00180000;
	private static final int cAtomFlagsESRTypeShift = 19;
	private static final int cAtomFlagsESRGroup		= 0x03E00000;
	private static final int cAtomFlagsESRGroupShift = 21;

	protected static final int cAtomFlagConfigurationUnknown = 0x04000000;
	private static final int cAtomFlagIsStereoCenter = 0x08000000;

	protected static final int cAtomFlagsValence	= 0xF0000000;
	private static final int cAtomFlagsValenceShift = 28;

	public static final int cAtomQFNoOfBits			= 29;
	public static final int cAtomQFAromStateBits	= 2;
	public static final int cAtomQFAromStateShift	= 1;
	public static final int cAtomQFRingStateBits	= 4;
	public static final int cAtomQFRingStateShift	= 3;
	public static final int cAtomQFHydrogenBits		= 4;
	public static final int cAtomQFHydrogenShift	= 7;
	public static final int cAtomQFPiElectronBits   = 3;
	public static final int cAtomQFPiElectronShift  = 14;
	public static final int cAtomQFNeighbourBits	= 5;
	public static final int cAtomQFNeighbourShift   = 17;
	public static final int cAtomQFRingSizeBits	 = 3;
	public static final int cAtomQFRingSizeShift	= 22;
	public static final int cAtomQFChargeBits		= 3;
	public static final int cAtomQFChargeShift		= 25;
	public static final int cAtomQFSimpleFeatures	= 0x0E3FC7FE;
	public static final int cAtomQFNarrowing		= 0x0E3FC7FE;
	public static final int cAtomQFAny				= 0x00000001;
	public static final int cAtomQFAromState		= 0x00000006;
	public static final int cAtomQFAromatic			= 0x00000002;
	public static final int cAtomQFNotAromatic		= 0x00000004;
	public static final int cAtomQFRingState		= 0x00000078;
	public static final int cAtomQFNotChain			= 0x00000008;
	public static final int cAtomQFNot2RingBonds	= 0x00000010;
	public static final int cAtomQFNot3RingBonds	= 0x00000020;
	public static final int cAtomQFNot4RingBonds	= 0x00000040;
	public static final int cAtomQFHydrogen			= 0x00000780;
	public static final int cAtomQFNot0Hydrogen		= 0x00000080;
	public static final int cAtomQFNot1Hydrogen		= 0x00000100;
	public static final int cAtomQFNot2Hydrogen		= 0x00000200;
	public static final int cAtomQFNot3Hydrogen	 = 0x00000400;
	public static final int cAtomQFNoMoreNeighbours	= 0x00000800;
	public static final int cAtomQFMoreNeighbours	= 0x00001000;
	public static final int cAtomQFMatchStereo		= 0x00002000;
	public static final int cAtomQFPiElectrons	  = 0x0001C000;
	public static final int cAtomQFNot0PiElectrons  = 0x00004000;
	public static final int cAtomQFNot1PiElectron   = 0x00008000;
	public static final int cAtomQFNot2PiElectrons  = 0x00010000;
	public static final int cAtomQFNeighbours	   = 0x003E0000;  // these QF refer to non-H neighbours
	public static final int cAtomQFNot0Neighbours   = 0x00020000;
	public static final int cAtomQFNot1Neighbour	= 0x00040000;
	public static final int cAtomQFNot2Neighbours   = 0x00080000;
	public static final int cAtomQFNot3Neighbours   = 0x00100000;
	public static final int cAtomQFNot4Neighbours   = 0x00200000;  // this is not 4 or more neighbours
	public static final int cAtomQFRingSize		 = 0x01C00000;
	public static final int cAtomQFCharge			= 0x0E000000;
	public static final int cAtomQFNotChargeNeg		= 0x02000000;
	public static final int cAtomQFNotCharge0		= 0x04000000;
	public static final int cAtomQFNotChargePos		= 0x08000000;
	public static final int cAtomQFFlatNitrogen		= 0x10000000;  // currently only used in TorsionDetail

	public static final int cBondTypeSingle			= 0x00000001;
	public static final int cBondTypeDouble			= 0x00000002;
	public static final int cBondTypeTriple			= 0x00000004;
	public static final int cBondTypeDown			= 0x00000009;
	public static final int cBondTypeUp				= 0x00000011;
	public static final int cBondTypeCross			= 0x0000001A;
	public static final int cBondTypeMetalLigand	= 0x00000020;
	public static final int cBondTypeDelocalized	= 0x00000040;
	public static final int cBondTypeDeleted		= 0x00000080;
	public static final int cBondTypeIncreaseOrder  = 0x0000007F;

	private static final int cBondTypeMaskSimple	= 0x00000067;	// masks
	protected static final int cBondTypeMaskStereo	= 0x00000018;

	protected static final int cBondFlagsHelper2	= 0x000003C0;
	protected static final int cBondFlagsHelper3	= 0x0000003F;

		// double bond E/Z parities based on atom positions in atom table
	private static final int cBondFlagsParity		= 0x00000003;
	public static final int cBondParityNone			= 0x00000000;
	public static final int cBondParityEor1			= 0x00000001;
	public static final int cBondParityZor2			= 0x00000002;
	public static final int cBondParityUnknown		= 0x00000003;
	private static final int cBondParityIsPseudo	= 0x00000004;

	private static final int cBondFlagsCIPParity	= 0x00000030;
	protected static final int cBondFlagsCIPParityShift = 4;
	public static final int cBondCIPParityNone		= 0x00000000;
	public static final int cBondCIPParityEorP		= 0x00000001;
	public static final int cBondCIPParityZorM	  = 0x00000002;
	public static final int cBondCIPParityProblem   = 0x00000003;

	protected static final int cBondFlagRing		= 0x00000040;
	protected static final int cBondFlagSmallRing	= 0x00000080;
	protected static final int cBondFlagAromatic	= 0x00000100;
	protected static final int cBondFlagDelocalized	= 0x00000200;
	protected static final int cBondFlagsESR		= 0x0001FC00;
	private static final int cBondFlagsESRType	  = 0x00000C00;
	private static final int cBondFlagsESRTypeShift = 10;
	private static final int cBondFlagsESRGroup	 = 0x0001F000;
	private static final int cBondFlagsESRGroupShift = 12;
	private static final int cBondFlagBGHilited		= 0x00020000;
	private static final int cBondFlagFGHilited		= 0x00040000;

	private static final int cBondParityUnknownOrNone	= 0x1000000;
	// This hint/flag is set by CoordinateInventor for double bonds without
	// given EZ-parity because coordinates may imply a not intended EZ-parity.
	// The setBondParity() method clears this flag. The Canonizer considers
	// this flag when calculating EZ-parities.

	public static final int cBondQFNoOfBits			= 20;
	public static final int cBondQFBondTypesBits	= 4;
	public static final int cBondQFBondTypesShift	= 0;
	public static final int cBondQFRingStateBits	= 2;
	public static final int cBondQFRingStateShift	= 4;
	public static final int cBondQFBridgeBits	   = 8;
	public static final int cBondQFBridgeShift	  = 6;
	public static final int cBondQFBridgeMinBits	= 4;
	public static final int cBondQFBridgeMinShift   = 6;
	public static final int cBondQFBridgeSpanBits   = 4;
	public static final int cBondQFBridgeSpanShift  = 10;
	public static final int cBondQFRingSizeBits	 = 3;
	public static final int cBondQFRingSizeShift	= 14;
	public static final int cBondQFAromStateBits	= 2;
	public static final int cBondQFAromStateShift	= 18;
	public static final int cBondQFAllFeatures		= 0x000FFFFF;
	public static final int cBondQFSimpleFeatures	= 0x000C003F;
	public static final int cBondQFNarrowing		= 0x000C0030;
	public static final int cBondQFBondTypes		= 0x0000000F;
	public static final int cBondQFSingle			= 0x00000001;
	public static final int cBondQFDouble			= 0x00000002;
	public static final int cBondQFTriple			= 0x00000004;
	public static final int cBondQFDelocalized		= 0x00000008;
	public static final int cBondQFRingState		= 0x00000030;
	public static final int cBondQFNotRing			= 0x00000010;
	public static final int cBondQFRing				= 0x00000020;
	public static final int cBondQFBridge		   = 0x00003FC0;
	public static final int cBondQFBridgeMin		= 0x000003C0;
	public static final int cBondQFBridgeSpan	   = 0x00003C00;
	public static final int cBondQFRingSize			= 0x0001C000;
	public static final int cBondQFMatchStereo		= 0x00020000;
	public static final int cBondQFAromState		= 0x000C0000;
	public static final int cBondQFAromatic			= 0x00040000;
	public static final int cBondQFNotAromatic		= 0x00080000;


	public static final int cHelperNone = 0x0000;
	public static final int cHelperBitNeighbours = 0x0001;
	public static final int cHelperBitRings = 0x0002;
	public static final int cHelperBitParities = 0x0004;
	public static final int cHelperBitCIP = 0x0008;

	public static final int cHelperBitSymmetrySimple = 0x0010;
	public static final int cHelperBitSymmetryDiastereotopic = 0x0020;
	public static final int cHelperBitSymmetryEnantiotopic = 0x0040;
	public static final int cHelperBitIncludeNitrogenParities = 0x0080;

	public static final int cHelperBitsStereo = 0x00FC;

	public static final int cHelperNeighbours = cHelperBitNeighbours;
	public static final int cHelperRings = cHelperNeighbours | cHelperBitRings;
	public static final int cHelperParities = cHelperRings | cHelperBitParities;
	public static final int cHelperCIP = cHelperParities | cHelperBitCIP;

	public static final int cHelperSymmetrySimple = cHelperCIP | cHelperBitSymmetrySimple;
	public static final int cHelperSymmetryDiastereotopic = cHelperCIP | cHelperBitSymmetryDiastereotopic;
	public static final int cHelperSymmetryEnantiotopic = cHelperCIP | cHelperBitSymmetryEnantiotopic;

	public static final int cChiralityIsomerCountMask   = 0x00FFFF;
	public static final int cChiralityUnknown		  	= 0x000000;
	public static final int cChiralityNotChiral			= 0x010000;
	public static final int cChiralityMeso			  = 0x020000; // this has added the number of meso isomers
	public static final int cChiralityRacemic			= 0x030000;
	public static final int cChiralityKnownEnantiomer   = 0x040000;
	public static final int cChiralityUnknownEnantiomer = 0x050000;
	public static final int cChiralityEpimers		   = 0x060000;
	public static final int cChiralityDiastereomers	 = 0x070000; // this has added the number of diastereomers

	public static final String cAtomLabel[] = { "?",
		"H"  ,"He" ,"Li" ,"Be" ,"B"  ,"C"  ,"N"  ,"O"  ,
		"F"  ,"Ne" ,"Na" ,"Mg" ,"Al" ,"Si" ,"P"  ,"S"  ,
		"Cl" ,"Ar" ,"K"  ,"Ca" ,"Sc" ,"Ti" ,"V"  ,"Cr" ,
		"Mn" ,"Fe" ,"Co" ,"Ni" ,"Cu" ,"Zn" ,"Ga" ,"Ge" ,
		"As" ,"Se" ,"Br" ,"Kr" ,"Rb" ,"Sr" ,"Y"  ,"Zr" ,
		"Nb" ,"Mo" ,"Tc" ,"Ru" ,"Rh" ,"Pd" ,"Ag" ,"Cd" ,
		"In" ,"Sn" ,"Sb" ,"Te" ,"I"  ,"Xe" ,"Cs" ,"Ba" ,
		"La" ,"Ce" ,"Pr" ,"Nd" ,"Pm" ,"Sm" ,"Eu" ,"Gd" ,
		"Tb" ,"Dy" ,"Ho" ,"Er" ,"Tm" ,"Yb" ,"Lu" ,"Hf" ,
		"Ta" ,"W"  ,"Re" ,"Os" ,"Ir" ,"Pt" ,"Au" ,"Hg" ,
		"Tl" ,"Pb" ,"Bi" ,"Po" ,"At" ,"Rn" ,"Fr" ,"Ra" ,
		"Ac" ,"Th" ,"Pa" ,"U"  ,"Np" ,"Pu" ,"Am" ,"Cm" ,
		"Bk" ,"Cf" ,"Es" ,"Fm" ,"Md" ,"No" ,"Lr" ,"??" ,
		"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,
		"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,
		"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,
		"R4" ,"R5" ,"R6" ,"R7" ,"R8" ,"R9" ,"R10","R11",	// R4 to R16 do not belong to the MDL set
		"R12","R13","R14","R15","R16","R1" ,"R2" ,"R3" ,
		"A"  ,"A1" ,"A2" ,"A3" ,"??" ,"??" ,"D"  ,"T"  ,
		"X"  ,"R"  ,"H2" ,"H+" ,"Nnn","HYD","Pol","??" ,
		"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,"??" ,
		"??" ,"??" ,"Ala","Arg","Asn","Asp","Cys","Gln",
		"Glu","Gly","His","Ile","Leu","Lys","Met","Phe",
		"Pro","Ser","Thr","Trp","Tyr","Val" };

	public static final short cRoundedMass[] = { 0,
	   1,	  4,	  7,	  9,	 11,	 12,   //  H  ,He ,Li ,Be ,B  ,C  ,
	  14,	 16,	 19,	 20,	 23,	 24,   //  N , O  ,F  ,Ne ,Na ,Mg ,
	  27,	 28,	 31,	 32,	 35,	 40,   //  Al ,Si ,P  ,S  ,Cl ,Ar ,
	  39,	 40,	 45,	 48,	 51,	 52,   //  K  ,Ca ,Sc ,Ti ,V  ,Cr ,
	  55,	 56,	 59,	 58,	 63,	 64,   //  Mn ,Fe ,Co ,Ni ,Cu ,Zn ,
	  69,	 74,	 75,	 80,	 79,	 84,   //  Ga ,Ge ,As ,Se ,Br ,Kr ,
	  85,	 88,	 89,	 90,	 93,	 98,   //  Rb ,Sr ,Y  ,Zr ,Nb ,Mo ,
	   0,	102,	103,	106,	107,	114,   //  Tc ,Ru ,Rh ,Pd ,Ag ,Cd ,
	 115,	120,	121,	130,	127,	132,   //  In ,Sn ,Sb ,Te ,I  ,Xe ,
	 133,	138,	139,	140,	141,	142,   //  Cs ,Ba ,La ,Ce ,Pr ,Nd ,
	   0,	152,	153,	158,	159,	164,   //  Pm ,Sm ,Eu ,Gd ,Tb ,Dy ,
	 165,	166,	169,	174,	175,	180,   //  Ho ,Er ,Tm ,Yb ,Lu ,Hf ,
	 181,	184,	187,	192,	193,	195,   //  Ta ,W , Re ,Os ,Ir ,Pt ,
	 197,	202,	205,	208,	209,	  0,   //  Au ,Hg ,Tl ,Pb ,Bi ,Po ,
	   0,	  0,	  0,	  0,	  0,	232,   //  At ,Rn ,Fr ,Ra ,Ac ,Th ,
	   0,	238,	  0,	  0,	  0,	  0,   //  Pa ,U , Np ,Pu ,Am ,Cm ,
	   0,	  0,	  0,	  0,	  0,	  0,   //  Bk ,Cf ,Es ,Fm ,Md ,No ,
	   0,	  0,	  0,	  0,	  0,	  0,   //  Lr ,?? ,?? ,?? ,?? ,?? ,
	   0,	  0,	  0,	  0,	  0,	  0,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
	   0,	  0,	  0,	  0,	  0,	  0,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
	   0,	  0,	  0,	  0,	  0,	  0,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
	   0,	  0,	  0,	  0,	  0,	  0,   //  ?? ,?? ,R4 ,R5 ,R6 ,R7 ,
	   0,	  0,	  0,	  0,	  0,	  0,   //  R8 ,R9 ,R10,R11,R12,R13,
	   0,	  0,	  0,	  0,	  0,	  0,   //  R14,R15,R16,R1 ,R2 ,R3 ,
	   0,	  0,	  0,	  0,	  0,	  0,   //  A  ,A1 ,A2 ,A3 ,?? ,?? ,
	   2,	  3,	  0,	  0,	  0,	  0,   //  D  ,T  ,X  ,R  ,H2 ,H+
	   0,	  0,	  0,	  0,	  0,	  0,   //  Nnn,HYD,Pol,?? ,?? ,?? ,
	   0,	  0,	  0,	  0,	  0,	  0,   //  ?? ,?? ,?? ,?? ,?? ,?? ,
	   0,	  0,	 71,	156,	114,	115,   //  ?? ,?? ,Ala,Arg,Asn,Asp,
	 103,	128,	129,	 57,	137,	113,   //  Cys,Gln,Glu,Gly,His,Ile,
	 113,	128,	131,	147,	 97,	 87,   //  Leu,Lys,Met,Phe,Pro,Ser,
	 101,	186,	163,	 99 };					//  Thr,Trp,Tyr,Val,

	public static final int cDefaultAtomValence = 6;
	public static final byte cAtomValence[][] = { null,
		{ 1 }, { 0 }, { 1 }, { 2 }, { 3 }, { 4 }, { 3 }, { 2 }, { 1 }, { 0 },   // H to Ne
		{ 1 }, { 2 }, { 3 }, { 4 }, { 3,5 }, { 2,4,6 }, { 1,3,5,7 }, { 0 },	 // Na to Ar
		{ 1 }, { 2 }, null, null, null, null, null, null, null, null,		   // K to Ni
		null, null, { 2,3 }, { 2,4 }, { 3,5 }, { 2,4,6 }, { 1,3,5,7 }, { 0,2 }, // Cu to Kr
		{ 1,2,3,4 }, { 2 }, null, null, null, null, null, null, null, null,	 // Rb to Pd
		null, null, { 1,2,3 }, { 2,4 }, { 3,5 }, { 2,4,6 }, { 1,3,5,7 },		// Ag to I
		{ 0,2,4,6 }, { 1 }, { 2 } };														 // Xe to Ba

	transient protected int mMaxAtoms;
	transient protected int mMaxBonds;

	transient protected int mValidHelperArrays;
	transient protected int mAllAtoms;
	transient protected int mAllBonds;
	transient protected int mAtomicNo[];
	transient protected int mAtomCharge[];
	transient protected int mAtomMapNo[];
	transient protected int mAtomMass[];
	transient protected int mAtomFlags[];
	transient protected int mAtomQueryFeatures[];
	transient protected int mBondAtom[][];
	transient protected int mBondType[];
	transient protected int mBondFlags[];
	transient protected int mBondQueryFeatures[];
	transient protected float mAtomX[];
	transient protected float mAtomY[];
	transient protected float mAtomZ[];
	transient protected boolean mIsFragment;
	transient protected boolean mIsRacemate;	 	// to indicate a molfileV2's chiral flat to be 0
	transient protected boolean mProtectHydrogen;  	// protects hydrogens atoms from being converted to query features
	transient protected int mChirality;			// property set by Canonizer
	transient protected int[][] mAtomList;
	transient protected byte[][] mAtomCustomLabel;

	transient private float mZoomRotationX,mZoomRotationY;
	transient private float mOriginalAngle[];
	transient private float mOriginalDistance[];
	transient private String mName;

    private float defaultBondLength = 24.0f;

    public static int getAtomicNoFromLabel(String atomLabel) {
		for (int i=1; i<cAtomLabel.length; i++)
			if (atomLabel.equalsIgnoreCase(cAtomLabel[i]))
				return i;
		return 0;
		}


	public static float getAngle(float x1, float y1, float x2, float y2) {
		float angle;
		float xdiff = x2 - x1;
		float ydiff = y2 - y1;

		if (ydiff != 0) {
			angle = (float)Math.atan(xdiff/ydiff);
			if (ydiff < 0) {
				if (xdiff < 0)
					angle -= Math.PI;
				else
					angle += Math.PI;
				}
			}
		else
			angle = (xdiff > 0f) ? (float)Math.PI/2 : -(float)Math.PI/2;

		return angle;
		}


	public static float getAngleDif(float angle1, float angle2) {
		float angleDif = angle1 - angle2;
		while (angleDif < -Math.PI)
			angleDif += 2 * Math.PI;
		while (angleDif > Math.PI)
			angleDif -= 2 * Math.PI;
		return angleDif;
		}


	public Molecule() {
		mMaxAtoms = mMaxBonds = 256;
		init();
		}


	public Molecule(int maxAtoms, int maxBonds) {
		mMaxAtoms = Math.max(1, maxAtoms);
		mMaxBonds = Math.max(1, maxBonds);
		init();
		}


	private void init() {
		mValidHelperArrays = cHelperNone;
		mAtomicNo = new int[mMaxAtoms];
		mAtomCharge = new int[mMaxAtoms];
		mAtomMapNo = new int[mMaxAtoms];
		mAtomX = new float[mMaxAtoms];
		mAtomY = new float[mMaxAtoms];
		mAtomZ = new float[mMaxAtoms];
		mAtomMass = new int[mMaxAtoms];
		mAtomFlags = new int[mMaxAtoms];
		mAtomQueryFeatures = new int[mMaxAtoms];
		mAtomList = null;
		mAtomCustomLabel = null;
		mBondAtom = new int[2][mMaxBonds];
		mBondType = new int[mMaxBonds];
		mBondFlags = new int[mMaxBonds];
		mBondQueryFeatures = new int[mMaxBonds];
		}


	/**
	 * High level function for constructing a molecule.
	 * @param x
	 * @param y
	 * @return
	 */
	public int addAtom(float x, float y) {
		return addAtom(x, y, 0);
		}


	/**
	 * High level function for constructing a molecule.
	 * @param x
	 * @param y
	 * @param z
	 * @return
	 */
	public int addAtom(float x, float y, float z) {
		int atom = addAtom(6);
		mAtomX[atom] = x;
		mAtomY[atom] = y;
		mAtomZ[atom] = z;
		return atom;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param atomLabel
	 * @return
	 */
	public int addAtom(String atomLabel) {
		int atomicNo = getAtomicNoFromLabel(atomLabel);
		return (atomicNo == 0) ? -1 : addAtom(atomicNo);
		}


	/**
	 * High level function for constructing a molecule.
	 * @param atomicNo
	 * @return
	 */
	public int addAtom(int atomicNo) {
		if (mAllAtoms >= mMaxAtoms)
			setMaxAtoms(mMaxAtoms*2);

		mAtomicNo[mAllAtoms] = 0;			// default
		setAtomicNo(mAllAtoms, atomicNo);	// sets atomicNo and mass

		mAtomCharge[mAllAtoms] = 0;
		mAtomFlags[mAllAtoms] = 0;
		mAtomQueryFeatures[mAllAtoms] = 0;
		mAtomMapNo[mAllAtoms] = 0;
		mAtomX[mAllAtoms] = 0;
		mAtomY[mAllAtoms] = 0;
		mAtomZ[mAllAtoms] = 0;

		if (mAtomList != null)
			mAtomList[mAllAtoms] = null;
		if (mAtomCustomLabel != null)
			mAtomCustomLabel[mAllAtoms] = null;

		mValidHelperArrays = cHelperNone;
		return mAllAtoms++;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param atm1
	 * @param atm2
	 * @param type
	 * @return
	 */
	public int addBond(int atm1,int atm2,int type) {
		if (atm1 == atm2)
			return -1;

		for (int bnd=0; bnd<mAllBonds; bnd++) {
			if (mBondAtom[0][bnd] == atm1 && mBondAtom[1][bnd] == atm2
			 || mBondAtom[0][bnd] == atm2 && mBondAtom[1][bnd] == atm1) {
				if (mBondType[bnd] < type)
					mBondType[bnd] = type;
				return bnd;
				}
			}

		if (mAllBonds >= mMaxBonds)
			setMaxBonds(mMaxBonds*2);

		mBondAtom[0][mAllBonds] = atm1;
		mBondAtom[1][mAllBonds] = atm2;
		mBondType[mAllBonds] = type;
		mBondFlags[mAllBonds] = 0;
		mBondQueryFeatures[mAllBonds] = 0;
		mValidHelperArrays = cHelperNone;
//		checkAtomParity(atm1);
//		checkAtomParity(atm2);
		return mAllBonds++;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param x
	 * @param y
	 * @param atomicNo
	 * @param mass
	 * @param abnormalValence
	 * @param radical
	 * @return
	 */
	public boolean addOrChangeAtom(float x, float y, int atomicNo, int mass, int abnormalValence, int radical) {
		int atom = findAtom(x,y);
		if (atom == -1) {
			if (mAllAtoms >= mMaxAtoms)
				setMaxAtoms(mMaxAtoms*2);

			atom = addAtom(atomicNo);
			mAtomX[atom] = x;
			mAtomY[atom] = y;
			mAtomMass[atom] = mass;
			setAtomAbnormalValence(atom, abnormalValence);
			setAtomRadical(atom, radical);
			return true;
			}

		return changeAtom(atom, atomicNo, mass, abnormalValence, radical);
		}


	/**
	 * High level function for constructing a molecule.
	 * @param atm1
	 * @param atm2
	 * @param type
	 * @return
	 */
	public int addOrChangeBond(int atm1,int atm2,int type) {
		for (int bnd=0; bnd<mAllBonds; bnd++) {
			if (mBondAtom[0][bnd] == atm1 && mBondAtom[1][bnd] == atm2
			 || mBondAtom[0][bnd] == atm2 && mBondAtom[1][bnd] == atm1) {
				changeBond(bnd,type);
				mValidHelperArrays = cHelperNone;
				return bnd;
				}
			}

		if(mAllBonds >= mMaxBonds)
			setMaxBonds(mMaxBonds*2);

		mBondAtom[0][mAllBonds] = atm1;
		mBondAtom[1][mAllBonds] = atm2;
		mBondType[mAllBonds] = type;
		mBondFlags[mAllBonds] = 0;
		mBondQueryFeatures[mAllBonds] = 0;
		mValidHelperArrays = cHelperNone;
		return mAllBonds++;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param x
	 * @param y
	 * @param ringSize
	 * @param aromatic
	 * @return
	 */
	public boolean addRing(float x, float y, int ringSize, boolean aromatic) {
		while(mAllAtoms + ringSize > mMaxAtoms)
			setMaxAtoms(mMaxAtoms*2);
		while(mAllBonds + ringSize > mMaxBonds)
			setMaxBonds(mMaxBonds*2);

		int atom = findAtom(x,y);
		if (atom != -1)
			return addRingToAtom(atom, ringSize, aromatic);

		int bond = findBond(x,y);
		if (bond != -1)
			return addRingToBond(bond, ringSize, aromatic);

		// new ring in empty space
		atom = addAtom(x,y);
		float cornerAngle = ((float)Math.PI * (ringSize-2))/ringSize;
		polygon(atom, ringSize, atom,aromatic, 0, (float)Math.PI - cornerAngle);
		mValidHelperArrays = cHelperNone;
		return true;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param atom
	 * @param ringSize
	 * @param aromatic
	 * @return
	 */
	public boolean addRingToAtom(int atom, int ringSize, boolean aromatic) {
		if ((aromatic && getOccupiedValence(atom) > 1)
		 || (!aromatic && getOccupiedValence(atom) > 2))
			return false;

		int angles = 0;
		float angle[] = new float[4];
		for (int i=0; i<mAllBonds; i++) {
			for (int j=0; j<2; j++) {
				if (mBondAtom[j][i] == atom) {
					if (angles == 2) {
						angles = 3;
						break;
						}
					angle[angles++] = getBondAngle(atom,mBondAtom[1-j][i]);
					}
				}
			if (angles == 3)
				break;
			}
		if (angles == 3)
			return false;

		float newAngle = (angles == 1) ? angle[0] + (float)Math.PI
				: (Math.abs(angle[0] - angle[1]) > (float)Math.PI) ? (angle[0] + angle[1])/2
				: (angle[0] + angle[1])/2 + (float)Math.PI;

				float cornerAngle = ((float)Math.PI * (ringSize-2))/ringSize;
		polygon(atom, ringSize, atom, aromatic, newAngle-cornerAngle/2, (float)Math.PI - cornerAngle);
		mValidHelperArrays = cHelperNone;
//				checkAtomParity(atom);
		return true;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param bond
	 * @param ringSize
	 * @param aromatic
	 * @return
	 */
	public boolean addRingToBond(int bond, int ringSize, boolean aromatic) {
		int bondAtom[] = new int[2];
		float bondAngle[] = new float[2];

		bondAtom[0] = mBondAtom[0][bond];
		bondAtom[1] = mBondAtom[1][bond];
		if (getOccupiedValence(bondAtom[0]) > 3) return false;
		if (getOccupiedValence(bondAtom[1]) > 3) return false;
		int angles = 0;
		float angle[] = new float[4];
		for (int i=0; i<mAllBonds; i++) {
			if (i == bond) continue;
			for (int j=0; j<2; j++) {
				for (int k=0; k<2; k++) {
					if (mBondAtom[j][i] == bondAtom[k]) {
						if (angles == 4) {
							angles = 5;
							break;
							}
						angle[angles++] = getBondAngle(bondAtom[k],
												   mBondAtom[1-j][i]);
						}
					}
				if (angles == 5) break;
				}
			if (angles == 5) break;
			}
		if (angles == 5) return false;

		bondAngle[0] = getBondAngle(bondAtom[0], bondAtom[1]);
		int atomNo;
		if (bondAngle[0] < 0) {
			bondAngle[1] = bondAngle[0] + (float)Math.PI;
			atomNo = 0;
			}
		else {
			bondAngle[1] = bondAngle[0];
			bondAngle[0] = bondAngle[1] - (float)Math.PI;
			atomNo = 1;
			}
		int side = 0;
		for (int i=0; i<angles; i++) {
			if ((angle[i] > bondAngle[0]) && (angle[i] < bondAngle[1]))
				side--;
			else
				side++;
			}

		atomNo = (side > 0) ? 1-atomNo : atomNo;
		float cornerAngle = ((float)Math.PI * (ringSize-2))/ringSize;
		polygon(bondAtom[atomNo], ringSize-1,
				bondAtom[1-atomNo], aromatic,
				bondAngle[(side > 0) ? 0 : 1] + (float)Math.PI - cornerAngle,
				(float)Math.PI - cornerAngle);

		mValidHelperArrays = cHelperNone;
//		checkAtomParity(bondAtom[0]);
//		checkAtomParity(bondAtom[1]);
		return true;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param atom
	 * @param atomicNo
	 * @param mass
	 * @param abnormalValence
	 * @param radical
	 * @return
	 */
	public boolean changeAtom(int atom, int atomicNo, int mass, int abnormalValence, int radical) {
		if ((atomicNo == 1 || atomicNo == 151 || atomicNo == 152)
		 && getOccupiedValence(atom) > 1)
			return false;

		mAtomQueryFeatures[atom] &= ~cAtomQFAny;
		if (mAtomList != null)
			mAtomList[atom] = null;
		if (mAtomCustomLabel != null)
			mAtomCustomLabel[atom] = null;

		if (atomicNo == mAtomicNo[atom]
		 && mass == mAtomMass[atom]
		 && abnormalValence == getAtomAbnormalValence(atom)
		 && radical == getAtomRadical(atom))
			return false;

		if (atomicNo == 151 || atomicNo == 152) {	// 'D' or 'T'
			mass = atomicNo - 149;
			atomicNo = 1;
			}

		mAtomFlags[atom] &= (cAtomFlagsColor | cAtomFlagSelected);
		mAtomicNo[atom] = atomicNo;
		mAtomMass[atom] = mass;
		mAtomCharge[atom] = 0;
		mAtomQueryFeatures[atom] = 0;
		setAtomAbnormalValence(atom, abnormalValence);
		setAtomRadical(atom, radical);
		removeMappingNo(mAtomMapNo[atom]);

		mValidHelperArrays = cHelperNone;
//		checkAtomParity(atom);
		return true;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param x
	 * @param y
	 * @param positive
	 * @return
	 */
	public boolean changeAtomCharge(float x, float y, boolean positive) {
		int atom = findAtom(x,y);
		return (atom == -1) ? false : changeAtomCharge(atom, positive);
		}


	/**
	 * High level function for constructing a molecule.
	 * @param atom
	 * @param positive
	 * @return
	 */
	public boolean changeAtomCharge(int atom, boolean positive) {
		if (positive) {
			if (mAtomCharge[atom] > 3) return false;
			mAtomCharge[atom]++;
			}
		else {
			if (mAtomCharge[atom] < -3) return false;
			mAtomCharge[atom]--;
			}
		mValidHelperArrays = cHelperNone;
//		checkAtomParity(atom);
		return true;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param bnd
	 * @param type
	 * @return
	 */
	public boolean changeBond(int bnd, int type) {
		boolean bondWasChanged = false;
		int oldType = mBondType[bnd];

		if (type == cBondTypeIncreaseOrder) {
			bondWasChanged = incrementBondOrder(bnd);
			}
		else if (validateBondType(bnd, type)) {
			if (type == cBondTypeUp || type == cBondTypeDown) {
				boolean bondAtAtom1Qualifies = qualifiesAsStereoBond(bnd, mBondAtom[0][bnd]);
				boolean bondAtAtom2Qualifies = qualifiesAsStereoBond(bnd, mBondAtom[1][bnd]);
				if (type == oldType) {
					// If both atoms are stereocenters (or none is recognized yet as stereo center)
					// or if only one atom is a stereo center and the bond points to the other one
					// then we can invert the bond direction.
					if (bondAtAtom1Qualifies == bondAtAtom2Qualifies
					 || bondAtAtom2Qualifies) {
						int temp = mBondAtom[0][bnd];
						mBondAtom[0][bnd] = mBondAtom[1][bnd];
						mBondAtom[1][bnd] = temp;
						bondWasChanged = true;
						}
					}
				else {
						// if stereo center information available assure proper bond direction
					if (!bondAtAtom1Qualifies && bondAtAtom2Qualifies) {
						int temp = mBondAtom[0][bnd];
						mBondAtom[0][bnd] = mBondAtom[1][bnd];
						mBondAtom[1][bnd] = temp;
						}
					mBondType[bnd] = type;
					bondWasChanged = true;
					}
				}
			else {
				mBondType[bnd] = type;
				bondWasChanged = true;
				}
			}

		if (bondWasChanged) {
			mValidHelperArrays = (oldType & cBondTypeMaskSimple)
								 == (type & cBondTypeMaskSimple) ?
						mValidHelperArrays & cHelperRings
					  : cHelperNone;
			mBondQueryFeatures[bnd] = 0;
			}

		return bondWasChanged;
		}


	/**
	 * Checks whether this bond may rightfully be an up/down stereo bond with its pointed end
	 * connected to atom. 
	 * i.e. whether it is a stereo center, an allene end atom, or an atom of a BINAP bond.
	 * @return true if an attached stereo bond aids in stereo recognition
	 */
	private boolean qualifiesAsStereoBond(int bond, int atom) {
		if (getBondOrder(bond) != 1)
			return false;

		if ((mAtomFlags[atom] & cAtomFlagsParity) != cAtomParityNone)
			return true;

		// stereo bond at allene 
		for (int i=0; i<mAllBonds; i++)
			if (i != bond
			 && mBondType[i] == cBondTypeDouble
			 &&	((mBondAtom[0][i] == atom && (mAtomFlags[mBondAtom[1][i]] & cAtomFlagsParity) != cAtomParityNone)
			  || (mBondAtom[1][i] == atom && (mAtomFlags[mBondAtom[0][i]] & cAtomFlagsParity) != cAtomParityNone)))
				return true;

		// stereo bond indicating BINAP axial chirality
		for (int i=0; i<mAllBonds; i++)
			if (i != bond
			 && mBondType[i] == cBondTypeSingle
			 && (mBondAtom[0][i] == atom || mBondAtom[1][i] == atom)
			 && (mBondFlags[i] & cBondFlagsParity) != cBondParityNone)
				return true;

		return false;
		}


	/**
	 * Copies all atoms and bonds of mol to the end of this Molecule's atom and bond
	 * tables. If mol is a fragment then this Molecule's fragment flag is set to true
	 * and all query features of mol are also copied.
	 * High level function for constructing a molecule.
	 * @param mol
	 * @return atom mapping from original mol to this molecule after incorporation of mol
	 */
	public int[] addMolecule(Molecule mol) {
		mIsFragment |= mol.mIsFragment;

		int[] atomMap = new int[mol.mAllAtoms];
		int esrGroupCountAND = renumberESRGroups(cESRTypeAnd);
		int esrGroupCountOR = renumberESRGroups(cESRTypeOr);
		for (int atom=0; atom<mol.mAllAtoms; atom++) {
			atomMap[atom] = mol.copyAtom(this, atom, esrGroupCountAND, esrGroupCountOR);
			}
		for (int bond=0; bond<mol.mAllBonds; bond++) {
			mol.copyBond(this, bond, esrGroupCountAND, esrGroupCountOR, atomMap, false);
			}

		mIsRacemate = (mIsRacemate && mol.mIsRacemate);
		mChirality = cChiralityUnknown;
		mValidHelperArrays = cHelperNone;
		return atomMap;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param substituent
	 * @param connectionAtom
	 * @return atom mapping from substituent to this molecule after addition of substituent
	 */
	public int[] addSubstituent(Molecule substituent, int connectionAtom) {
		int[] atomMap = new int[substituent.mAllAtoms];
		int esrGroupCountAND = renumberESRGroups(cESRTypeAnd);
		int esrGroupCountOR = renumberESRGroups(cESRTypeOr);
		for (int atom=0; atom<substituent.mAllAtoms; atom++) {
			if (substituent.getAtomicNo(atom) != 0)
				atomMap[atom] = substituent.copyAtom(this, atom, esrGroupCountAND, esrGroupCountOR);
			else
				atomMap[atom] = connectionAtom;
			}
		for (int bond=0; bond<substituent.mAllBonds; bond++) {
			substituent.copyBond(this, bond, esrGroupCountAND, esrGroupCountOR, atomMap, false);
			}

		mIsRacemate = (mIsRacemate && substituent.mIsRacemate);
		mChirality = cChiralityUnknown;
		mValidHelperArrays = cHelperNone;
		return atomMap;
		}


	/**
	 * Copies this molecule including parity settings, if valid.
	 * The original content of destMol is replaced.
	 * Helper arrays are not copied and need to be recalculated if needed.
	 * @param destMol
	 */
	public void copyMolecule(Molecule destMol) {
		destMol.mAtomList = null;
		destMol.mAtomCustomLabel = null;

		destMol.mIsFragment = mIsFragment;	// to allow copying of atom/bond query features

		destMol.mAllAtoms = 0;
		for (int atom=0; atom<mAllAtoms;atom++)
			copyAtom(destMol, atom, 0, 0);

		destMol.mAllBonds = 0;
		for (int bnd=0; bnd<mAllBonds;bnd++)
			copyBond(destMol, bnd, 0, 0, null, false);

		copyMoleculeProperties(destMol);
		}


	/**
	 * Creates a new atom in destMol and copies all source atom properties
	 * including atom list, custom label, flags, and mapNo to it.
	 * @param destMol
	 * @param sourceAtom
	 * @param esrGroupOffsetAND -1 to create new ESR group or destMol ESR group count from esrGroupCountAND()
	 * @param esrGroupOffsetOR -1 to create new ESR group or destMol ESR group count from esrGroupCountOR()
	 * @return index of new atom in destMol
	 */
	public int copyAtom(Molecule destMol, int sourceAtom, int esrGroupOffsetAND, int esrGroupOffsetOR) {
		int destAtom = destMol.mAllAtoms;
		if (destAtom >= destMol.mMaxAtoms)
			destMol.setMaxAtoms(destMol.mMaxAtoms*2);

		int esrType = getAtomESRType(sourceAtom);
		int esrGroup = -1;
		if (esrType == cESRTypeAnd) {
			if (esrGroupOffsetAND == -1)   // create a new ESR group for this atom
				esrGroup = destMol.renumberESRGroups(esrType);
			else	// take existing group and add offset that should be the
					// ESR group member count of destMol before starting to add atoms
				esrGroup = Math.min(cESRMaxGroups, esrGroupOffsetAND + getAtomESRGroup(sourceAtom));
			}
		else if (esrType == cESRTypeOr) {
			if (esrGroupOffsetOR == -1)   // create a new ESR group for this atom
				esrGroup = destMol.renumberESRGroups(esrType);
			else	// take existing group and add offset that should be the
					// ESR group member count of destMol before starting to add atoms
				esrGroup = Math.min(cESRMaxGroups, esrGroupOffsetOR + getAtomESRGroup(sourceAtom));
			}

		destMol.mAtomicNo[destAtom] = mAtomicNo[sourceAtom];
		destMol.mAtomCharge[destAtom] = mAtomCharge[sourceAtom];
		destMol.mAtomMass[destAtom] = mAtomMass[sourceAtom];
		destMol.mAtomFlags[destAtom] = mAtomFlags[sourceAtom];
		destMol.mAtomQueryFeatures[destAtom] = destMol.mIsFragment ? mAtomQueryFeatures[sourceAtom] : 0;
		destMol.mAtomX[destAtom] = mAtomX[sourceAtom];
		destMol.mAtomY[destAtom] = mAtomY[sourceAtom];
		destMol.mAtomZ[destAtom] = mAtomZ[sourceAtom];
		destMol.mAtomMapNo[destAtom] = mAtomMapNo[sourceAtom];

		if (destMol.mAtomList != null)
			destMol.mAtomList[destAtom] = null;
		if (mAtomList != null && mAtomList[sourceAtom] != null && destMol.mIsFragment) {
			if (destMol.mAtomList == null)
				destMol.mAtomList = new int[destMol.mAtomicNo.length][];

			destMol.mAtomList[destAtom] = new int[mAtomList[sourceAtom].length];
			for (int i=0; i<mAtomList[sourceAtom].length; i++)
				destMol.mAtomList[destAtom][i] = mAtomList[sourceAtom][i];
			}

		if (destMol.mAtomCustomLabel != null)
			destMol.mAtomCustomLabel[destAtom] = null;
		if (mAtomCustomLabel != null && mAtomCustomLabel[sourceAtom] != null) {
			if (destMol.mAtomCustomLabel == null)
				destMol.mAtomCustomLabel = new byte[destMol.mAtomicNo.length][];

			destMol.mAtomCustomLabel[destAtom] = new byte[mAtomCustomLabel[sourceAtom].length];
			for (int i=0; i<mAtomCustomLabel[sourceAtom].length; i++)
				destMol.mAtomCustomLabel[destAtom][i] = mAtomCustomLabel[sourceAtom][i];
			}

		if (esrGroup != -1) {
			destMol.mAtomFlags[destAtom] &= ~cAtomFlagsESRGroup;
			destMol.mAtomFlags[destAtom] |= (esrGroup << cAtomFlagsESRGroupShift);
			}

		destMol.mAllAtoms++;
		destMol.mValidHelperArrays = cHelperNone;

		return destAtom;
		}


	/**
	 * @param destMol
	 * @param sourceBond
	 * @param esrGroupOffsetAND -1 to create new ESR group or destMol ESR group count from esrGroupCountAND()
	 * @param esrGroupOffsetOR -1 to create new ESR group or destMol ESR group count from esrGroupCountOR()
	 * @param atomMap
	 * @param useBondTypeDelocalized
	 * @return
	 */
	public int copyBond(Molecule destMol, int sourceBond, int esrGroupOffsetAND, int esrGroupOffsetOR,
						int[] atomMap, boolean useBondTypeDelocalized) {
		int destBond = destMol.mAllBonds;
		if (destBond >= destMol.mMaxBonds)
			destMol.setMaxBonds(destMol.mMaxBonds * 2);

		int esrType = getBondESRType(sourceBond);
		int esrGroup = -1;
		if (esrType == cESRTypeAnd) {
			if (esrGroupOffsetAND == -1)   // create a new ESR group for this atom
				esrGroup = destMol.renumberESRGroups(esrType);
			else	// take existing group and add offset that should be the
					// ESR group member count of destMol before starting to add atoms
				esrGroup = Math.min(cESRMaxGroups, esrGroupOffsetAND + getBondESRGroup(sourceBond));
			}
		if (esrType == cESRTypeOr) {
			if (esrGroupOffsetOR == -1)   // create a new ESR group for this atom
				esrGroup = destMol.renumberESRGroups(esrType);
			else	// take existing group and add offset that should be the
					// ESR group member count of destMol before starting to add atoms
				esrGroup = Math.min(cESRMaxGroups, esrGroupOffsetOR + getBondESRGroup(sourceBond));
			}

		for (int i=0; i<2; i++)
			destMol.mBondAtom[i][destBond] = (atomMap == null) ?
					mBondAtom[i][sourceBond] : atomMap[mBondAtom[i][sourceBond]];

		int bondType = (useBondTypeDelocalized
					 && (mBondFlags[sourceBond] & cBondFlagDelocalized) != 0) ?
									cBondTypeDelocalized : mBondType[sourceBond];
		destMol.mBondType[destBond] = bondType;
		destMol.mBondFlags[destBond] = mBondFlags[sourceBond];
		destMol.mBondQueryFeatures[destBond] = destMol.mIsFragment ? mBondQueryFeatures[sourceBond] : 0;

		if (esrGroup != -1) {
			destMol.mBondFlags[destBond] &= ~cBondFlagsESRGroup;
			destMol.mBondFlags[destBond] |= (esrGroup << cBondFlagsESRGroupShift);
			}

		destMol.mAllBonds++;
		destMol.mValidHelperArrays = cHelperNone;

		return destBond;
		}


	/**
	 * Copies name,isFragment,chirality and validity of parity & CIP flags.
	 * When copying molecules parts only or when changing the atom order during copy,
	 * then atom parities or CIP parities may not be valid anymore and
	 * invalidateHelperArrays([affected bits]) should be called in these cases.
	 * @param destMol
	 */
	public void copyMoleculeProperties(Molecule destMol) {
		destMol.mIsFragment = mIsFragment;
		destMol.mIsRacemate = mIsRacemate;
		destMol.mProtectHydrogen = mProtectHydrogen;
		destMol.mChirality = mChirality;
		destMol.mName = (mName == null) ? null : new String(mName);
		destMol.mValidHelperArrays = (mValidHelperArrays & (cHelperBitParities | cHelperBitCIP));
		}

	/**
	 * Clears helperBits from mValidHelperArrays.
	 * @param helperBits
	 */
	public void invalidateHelperArrays(int helperBits) {
		mValidHelperArrays &= ~helperBits;
		}


	/**
	 * For the given ESR type (AND or OR) renumbers all group indexes starting from 0.
	 * Use this, if stereo center deletion or other operations caused an inconsisten ESR
	 * number state. Molecule and derived methods do this automatically.
	 * @param type cESRTypeAnd or cESRTypeOr
	 * @return number of ESR groups
	 */
	public int renumberESRGroups(int type) {
		if (type == cESRTypeAbs)
			return 0;

		// reassign group numbers from bottom up
		boolean[] groupUsed = null;
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (getAtomESRType(atom) == type) {
				if (groupUsed == null)
					groupUsed = new boolean[cESRMaxGroups];
				groupUsed[getAtomESRGroup(atom)] = true;
				}
			}
		for (int bond=0; bond<mAllBonds; bond++) {
			if (getBondESRType(bond) == type) {
				if (groupUsed == null)
					groupUsed = new boolean[cESRMaxGroups];
				groupUsed[getBondESRGroup(bond)] = true;
				}
			}
	
		int newIndex = 0;
		if (groupUsed != null) {
			int[] newGroup = new int[cESRMaxGroups];
			for (int i=0; i<cESRMaxGroups; i++)
				if (groupUsed[i])
					newGroup[i] = newIndex++;
	
			for (int atom=0; atom<mAllAtoms; atom++) {
				if (getAtomESRType(atom) == type) {
					int group = newGroup[getAtomESRGroup(atom)];
					mAtomFlags[atom] &= ~cAtomFlagsESRGroup;
					mAtomFlags[atom] |= (group << cAtomFlagsESRGroupShift);
					}
				}
			for (int bond=0; bond<mAllBonds; bond++) {
				if (getBondESRType(bond) == type) {
					int group = newGroup[getBondESRGroup(bond)];
					mBondFlags[bond] &= ~cBondFlagsESRGroup;
					mBondFlags[bond] |= (group << cBondFlagsESRGroupShift);
					}
				}
			}
	
		return newIndex;
		}

	/**
	 * Swaps two atoms' indexes/locations in the atom table. This is used to move hydrogen atoms
	 * to the end of the table and for some testsing purposes.
	 * @param atom1
	 * @param atom2
	 */
	protected void swapAtoms(int atom1, int atom2) {
		int tempInt = mAtomicNo[atom1];
		mAtomicNo[atom1] = mAtomicNo[atom2];
		mAtomicNo[atom2] = tempInt;
		tempInt = mAtomCharge[atom1];
		mAtomCharge[atom1] = mAtomCharge[atom2];
		mAtomCharge[atom2] = tempInt;
		tempInt = mAtomMass[atom1];
		mAtomMass[atom1] = mAtomMass[atom2];
		mAtomMass[atom2] = tempInt;
		tempInt = mAtomFlags[atom1];
		mAtomFlags[atom1] = mAtomFlags[atom2];
		mAtomFlags[atom2] = tempInt;
		tempInt = mAtomQueryFeatures[atom1];
		mAtomQueryFeatures[atom1] = mAtomQueryFeatures[atom2];
		mAtomQueryFeatures[atom2] = tempInt;
		tempInt = mAtomMapNo[atom1];
		mAtomMapNo[atom1] = mAtomMapNo[atom2];
		mAtomMapNo[atom2] = tempInt;
		float tempFloat = mAtomX[atom1];
		mAtomX[atom1] = mAtomX[atom2];
		mAtomX[atom2] = tempFloat;
		tempFloat = mAtomY[atom1];
		mAtomY[atom1] = mAtomY[atom2];
		mAtomY[atom2] = tempFloat;
		tempFloat = mAtomZ[atom1];
		mAtomZ[atom1] = mAtomZ[atom2];
		mAtomZ[atom2] = tempFloat;
		if (mAtomList != null) {
			int[] tempList = mAtomList[atom1];
			mAtomList[atom1] = mAtomList[atom2];
			mAtomList[atom2] = tempList;
			}
		if (mAtomCustomLabel != null) {
			byte[] tempList = mAtomCustomLabel[atom1];
			mAtomCustomLabel[atom1] = mAtomCustomLabel[atom2];
			mAtomCustomLabel[atom2] = tempList;
			}

		for (int bond=0; bond<mAllBonds; bond++) {
			for (int i=0; i<2; i++) {
				if (mBondAtom[i][bond] == atom1)
					mBondAtom[i][bond] = atom2;
				else if (mBondAtom[i][bond] == atom2)
					mBondAtom[i][bond] = atom1;
				}
			}

		mValidHelperArrays = cHelperNone;
		}


	/**
	 * Swaps two bonds' indexes/locations in the atom table. This is used to move hydrogen atoms
	 * to the end of the table and for some testsing purposes.
	 * @param bond1
	 * @param bond2
	 */
	protected void swapBonds(int bond1, int bond2) {
		int temp = mBondAtom[0][bond1];
		mBondAtom[0][bond1] = mBondAtom[0][bond2];
		mBondAtom[0][bond2] = temp;

		temp = mBondAtom[1][bond1];
		mBondAtom[1][bond1] = mBondAtom[1][bond2];
		mBondAtom[1][bond2] = temp;

		temp = mBondType[bond1];
		mBondType[bond1] = mBondType[bond2];
		mBondType[bond2] = temp;

		temp = mBondFlags[bond1];
		mBondFlags[bond1] = mBondFlags[bond2];
		mBondFlags[bond2] = temp;

		temp = mBondQueryFeatures[bond1];
		mBondQueryFeatures[bond1] = mBondQueryFeatures[bond2];
		mBondQueryFeatures[bond2] = temp;

		mValidHelperArrays = cHelperNone;
		}


	/**
	 * High level function for constructing a molecule.
	 * After the deletion the original order of atom and bond indexes is retained.
	 * @param atom
	 */
	public void deleteAtom(int atom) {
		for (int bnd=0; bnd<mAllBonds; bnd++) {
			for (int i=0; i<2; i++) {
				if (mBondAtom[i][bnd] == atom) {
					mBondType[bnd] = cBondTypeDeleted;	// mark for delete
					int bonds = 0;
					for (int j=0; j<mAllBonds; j++) {
						if (j == bnd) continue;
						if ((mBondAtom[0][j] == mBondAtom[1-i][bnd])
						 || (mBondAtom[1][j] == mBondAtom[1-i][bnd]))
							bonds++;
						}
					if (bonds == 0) {
						removeMappingNo(mAtomMapNo[mBondAtom[1-i][bnd]]);
						mAtomicNo[mBondAtom[1-i][bnd]] = -1;
						}			// mark for delete
					}
				}
			}
		removeMappingNo(mAtomMapNo[atom]);
		mAtomicNo[atom] = -1;		// mark for delete
		if (mAtomList != null)
			mAtomList[atom] = null;
		if (mAtomCustomLabel != null)
			mAtomCustomLabel[atom] = null;
		compressMolTable();
		mValidHelperArrays = cHelperNone;
		}


	/**
	 * High level function for constructing a molecule.
	 * @param x
	 * @param y
	 * @return
	 */
	public boolean deleteAtomOrBond(float x, float y) {
		int atom = findAtom(x,y);
		if (atom != -1) {
			if ((mAtomFlags[atom] & cAtomFlagSelected) != 0)
				deleteSelectedAtoms();
			else
				deleteAtom(atom);
			mValidHelperArrays = cHelperNone;
			return true;
			}
		int bnd = findBond(x,y);
		if (bnd != -1) {
			if (((mAtomFlags[mBondAtom[0][bnd]]
				& mAtomFlags[mBondAtom[1][bnd]])
			   & cAtomFlagSelected) != 0)
				deleteSelectedAtoms();
			else
				deleteBondAndSurrounding(bnd);
			mValidHelperArrays = cHelperNone;
			return true;
			}
		return false;
		}


	/**
	 * High level function for constructing a molecule.
	 * After the deletion the original order of atom and bond indexes is retained.
	 * @param bond
	 */
	public void deleteBond(int bond) {
		mBondType[bond] = cBondTypeDeleted;
		compressMolTable();
		mValidHelperArrays = cHelperNone;
		}


	/**
	 * High level function for constructing a molecule.
	 * After the deletion the original order of atom and bond indexes is retained.
	 * @param bond
	 */
	public void deleteBondAndSurrounding(int bond) {
		for (int i=0; i<2; i++) {
			int bonds = 0;
			for (int j=0; j<mAllBonds; j++) {
				if (j == bond) continue;
				if ((mBondAtom[0][j] == mBondAtom[i][bond])
				 || (mBondAtom[1][j] == mBondAtom[i][bond]))
					bonds++;
				}
			if (bonds == 0) {
				removeMappingNo(mAtomMapNo[mBondAtom[i][bond]]);
				mAtomicNo[mBondAtom[i][bond]] = -1;
				}							// mark for delete
			}
		mBondType[bond] = cBondTypeDeleted;
		compressMolTable();
		mValidHelperArrays = cHelperNone;
		}


	/**
	 * High level function for constructing a molecule.
	 * This method deletes atoms flagged in deleteAtom and all bonds leading to them.
	 * After the deletion the original order of atom and bond indexes is retained.
	 * @param deleteAtom
	 * @return
	 */
	public int[] deleteAtoms(boolean[] deleteAtom) {
		boolean found = false;
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (deleteAtom[atom]) {
				markAtomForDeletion(atom);
				found = true;
				}
			}
		if (!found)
			return null;

		return deleteMarkedAtomsAndBonds();
		}


	/**
	 * High level function for constructing a molecule.
	 * After the deletion the original order of atom and bond indexes is retained.
	 * @param atomList
	 * @return
	 */
	public int[] deleteAtoms(int[] atomList) {
		if (atomList.length == 0)
			return null;

		for (int i=0; i<atomList.length; i++)
			markAtomForDeletion(atomList[i]);

		return deleteMarkedAtomsAndBonds();
		}


	/**
	 * High level function for constructing a molecule.
	 * Delete all selected atoms and all bonds attached to them.
	 * After the deletion the original order of atom and bond indexes is retained.
	 * @return
	 */
	public boolean deleteSelectedAtoms() {
		boolean found = false;
		for (int atom=0; atom<mAllAtoms; atom++) {
			if ((mAtomFlags[atom] & cAtomFlagSelected) != 0) {
				markAtomForDeletion(atom);
				found = true;
				}
			}

		return found ? deleteMarkedAtomsAndBonds() != null : false;
		}


	/**
	 * Marks this atom to be deleted in a later call to deleteMarkedAtomsAndBonds().
	 * @param atom
	 */
	public void markAtomForDeletion(int atom) {
		mAtomicNo[atom] = -1;		// mark for delete
		}


	/**
	 * Marks this bond to be deleted in a later call to deleteMarkedAtomsAndBonds().
	 * @param bond
	 */
	public void markBondForDeletion(int bond) {
		mBondType[bond] = cBondTypeDeleted;		// mark for delete
		}


	/**
	 * Checks whether this atom was marked to be deleted and not deleted yet.
	 * @param atom
	 * @return
	 */
	public boolean isAtomMarkedForDeletion(int atom) {
		return (mAtomicNo[atom] == -1);
		}


	/**
	 * Checks whether this bond was marked to be deleted and not deleted yet.
	 * @param bond
	 * @return
	 */
	public boolean isBondMarkedForDeletion(int bond) {
		return (mBondType[bond] == cBondTypeDeleted);
		}


	/**
	 * High level function for constructing a molecule.
	 * Deletes all atoms and bonds from the molecule, which were marked before for deletion
	 * by calling markAtomForDeletion() or markBondForDeletion(). Bonds connecting atoms
	 * of which at least one is marked for deletion, are deleted automatically and don't
	 * require to be explicitly marked.<br>
	 * When multiple atoms and/or bonds need to be deleted, marking them and calling
	 * this method is more efficient than deleting them individually with deleteAtom() and
	 * deleteBond().
	 * Bonds, whose atoms carry opposite charges are treated in the following manner: If only one of
	 * the two bond atoms is kept, then its absolute charge will be reduced by 1.
	 * After the deletion the original order of atom and bond indexes is retained.
	 * @return mapping from old to new atom indices; null if no atoms nor bonds were deleted
	 */
	public int[] deleteMarkedAtomsAndBonds() {
		boolean found = false;
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (mAtomicNo[atom] == -1) {
				found = true;
				removeMappingNo(mAtomMapNo[atom]);
				}
			}
		for (int bond=0; bond<mAllBonds; bond++) {
			if (mBondType[bond] == cBondTypeDeleted) {
				found = true;
				}
			else if (mAtomicNo[mBondAtom[0][bond]] == -1
				  || mAtomicNo[mBondAtom[1][bond]] == -1) {
				mBondType[bond] = cBondTypeDeleted;  // mark for delete
				found = true;
				}
			}

		if (found) {
			mValidHelperArrays = cHelperNone;
			return compressMolTable();
			}

		return null;
		}


	/**
	 * Empties the molecule to serve as container for constructing a new molecule,
	 * e.g. by multiply calling addAtom(...), addBond(...) and other high level methods.
	 */
	public void deleteMolecule() {
		mAllAtoms = 0;
		mAllBonds = 0;
		mIsFragment = false;
		mIsRacemate = false;
		mChirality = cChiralityUnknown;
		mAtomList = null;
		mAtomCustomLabel = null;
		mName = null;
		mValidHelperArrays = cHelperNone;
		}


	public void removeAtomSelection() {
		for (int i=0; i<mAllAtoms; i++)
			mAtomFlags[i] &= ~cAtomFlagSelected;
		}


	public void removeAtomColors() {
		for (int i=0; i<mAllAtoms; i++)
			mAtomFlags[i] &= ~cAtomFlagsColor;
		}


	/**
	 * This removes all custom labels from the atoms.
	 */
	public void removeAtomCustomLabels() {
		mAtomCustomLabel = null;
		}


	public void removeAtomMarkers() {
		for (int i=0; i<mAllAtoms; i++)
			mAtomFlags[i] &= ~cAtomFlagMarked;
		}


	public void removeBondHiliting() {
		for (int i=0; i<mAllBonds; i++)
			mBondFlags[i] &= ~(cBondFlagBGHilited | cBondFlagFGHilited);
		}


	/**
	 * @param pickx
	 * @param picky
	 * @return index of closest of nearby atoms or -1, if no atom is close
	 */
	public int findAtom(float pickx, float picky) {
		int foundAtom = -1;
		float avbl = this.getAverageBondLength();
		float foundDistanceSquare = Float.MAX_VALUE;
		float maxDistanceSquare = avbl * avbl / 12.0f;
		for (int atom=0; atom<mAllAtoms; atom++) {
			float distanceSquare = (pickx-mAtomX[atom]) * (pickx-mAtomX[atom])
								  + (picky-mAtomY[atom]) * (picky-mAtomY[atom]);
			if (distanceSquare < maxDistanceSquare
			 && distanceSquare < foundDistanceSquare) {
				foundDistanceSquare = distanceSquare;
				foundAtom = atom;
				}
			}
		return foundAtom;
		}


	/**
	 * @param pickx
	 * @param picky
	 * @return index of closest of nearby bonds or -1, if no bond is close
	 */
	public int findBond(float pickx, float picky) {
		int foundBond = -1;
		float maxDistance = this.getAverageBondLength();
		float foundDistance = Float.MAX_VALUE;
		for (int bond=0; bond<mAllBonds; bond++) {
			float x1 = mAtomX[mBondAtom[0][bond]];
			float y1 = mAtomY[mBondAtom[0][bond]];
			float x2 = mAtomX[mBondAtom[1][bond]];
			float y2 = mAtomY[mBondAtom[1][bond]];

			float dx = x2 - x1;// if pick position not within bond circle continue
			float dy = y2 - y1;
			float bondLength = (float)Math.sqrt((dx*dx + dy*dy));
			float centralX = (x1+x2) / 2.0f;
			float centralY = (y1+y2) / 2.0f;
			dx = pickx - centralX;
			dy = picky - centralY;
			if (Math.sqrt(dx*dx + dy*dy) > bondLength / 2.0) continue;

			float distance;
			if (x2 == x1)// if pick posn closer to bond than cPickRange return bnd
				distance = Math.abs(x1 - pickx);
			else {
				float constA = (y2-y1)/(x1-x2);
				float constC = - constA * x1 - y1;
				distance = Math.abs((constA * pickx + picky + constC)
						 / (float)Math.sqrt(constA * constA + 1f));
				}
			if (distance < maxDistance
			 && distance < foundDistance) {
				foundDistance = distance;
				foundBond = bond;
				}
			}
		return foundBond;
		}


	/**
	 * @return the number of all atoms, which includes hydrogen atoms.
	 */
	public int getAllAtoms() {
		return mAllAtoms;
		}


	/**
	 * @return the number of all bonds, which includes those connecting hydrogen atoms.
	 */
	public int getAllBonds() {
		return mAllBonds;
		}


	/**
	 * Get an atom's defined maximum valance if different from the default one.
	 * @param atom
	 * @return valence 0-14: new maximum valence; -1: use default
	 */
	public int getAtomAbnormalValence(int atom) {
		return ((mAtomFlags[atom] & cAtomFlagsValence) >>> cAtomFlagsValenceShift) -1;
		}


	/**
	 * @param atom
	 * @return the formal atom charge
	 */
	public int getAtomCharge(int atom) {
		return mAtomCharge[atom];
		}


	/**
	 * The atom Cahn-Ingold-Prelog parity is a calculated property available above/equal helper level cHelperCIP.
	 * It encodes the stereo configuration of an atom with its neighbors using up/down-bonds
	 * or 3D-atom-coordinates, whatever is available. It depends on the atom indices of the neighbor
	 * atoms and their orientation is space. This method is called by the Canonizer and usually should not
	 * be called otherwise.
	 * @param atom
	 * @return one of cAtomCIPParityNone,cAtomCIPParityRorM,cAtomCIPParitySorP,cAtomCIPParityProblem
	 */
	public int getAtomCIPParity(int atom) {
		return (mAtomFlags[atom] & cAtomFlagsCIPParity) >> cAtomFlagsCIPParityShift;
		}


	public int getAtomColor(int atom) {
		return mAtomFlags[atom] & cAtomFlagsColor;
		}


	/**
	 * This is MDL's enhanced stereo representation (ESR).
	 * Stereo atoms and bonds with the same ESR type (AND or OR) and the same ESR group number
	 * are in the same group, i.e. within this group they have the defined (relative) stereo configuration.
	 * @param atom
	 * @return group index starting with 0
	 */
	public int getAtomESRGroup(int atom) {
		if (getAtomESRType(atom) != cESRTypeAnd
		 && getAtomESRType(atom) != cESRTypeOr)
			return -1;
		else
			return (mAtomFlags[atom] & cAtomFlagsESRGroup) >> cAtomFlagsESRGroupShift;
		}


	/**
	 * This is MDL's enhanced stereo representation (ESR).
	 * Stereo atoms and bonds with the same ESR type (AND or OR) and the same ESR group number
	 * are in the same group, i.e. within this group they have the defined (relative) stereo configuration.
	 * @param atom
	 * @return one of cESRTypeAbs,cESRTypeAnd,cESRTypeOr
	 */
	public int getAtomESRType(int atom) {
		return (mAtomFlags[atom] & cAtomFlagsESRType) >> cAtomFlagsESRTypeShift;
		}


	/**
	 * In addition to the natural atomic numbers, we support additional pseudo atomic numbers.
	 * Most of these are for compatibility with the MDL atom table, e.g. for amino acids and R-groups.
	 * D and T are accepted for setting, but are on-the-fly translated to H with the proper atom mass.
	 * @param atom
	 * @return
	 */
	public int getAtomicNo(int atom) {
		return mAtomicNo[atom];
		}


	/**
	 * If a custom atom label is set, a molecule depiction displays
	 * the custom label instead of the original one.
	 * @param atom
	 * @return null or previously defined atom custom label
	 */
	public String getAtomCustomLabel(int atom) {
		return (mAtomCustomLabel == null) ? null
			 : (mAtomCustomLabel[atom] == null) ? null
			 : new String(mAtomCustomLabel[atom]);
		}


	/**
	 * This method is more efficient than getAtomCustomLabel(),
	 * because internally atom custom labels are stored as a byte[].
	 * Use this if you can work with bytes and don't need a String.
	 * @param atom
	 * @return null or previously defined atom custom label as byte[]
	 */
	public byte[] getAtomCustomLabelBytes(int atom) {
		return (mAtomCustomLabel == null) ? null : mAtomCustomLabel[atom];
		}


	/**
	 * @param atom
	 * @return standard atom label of the atom: C,Li,Sc,...
	 */
	public String getAtomLabel(int atom) {
		return cAtomLabel[mAtomicNo[atom]];
		}


	/**
	 * The list of atoms that are allowed at this position during sub-structure search.
	 * (or refused atoms, if atom query feature cAtomQFAny is set).
	 * @param atom
	 * @return null or sorted list of unique atomic numbers, if defined
	 */
	public int[] getAtomList(int atom) {
		return (mAtomList == null) ? null : mAtomList[atom];
		}


	public String getAtomListString(int atom) {
		if (mAtomList == null || mAtomList[atom] == null)
			return ((mAtomQueryFeatures[atom] & cAtomQFAny) != 0) ?
						"" : cAtomLabel[mAtomicNo[atom]];

		String listString = "";
		for (int i=0; i<mAtomList[atom].length; i++) {
			if (i > 0)
				listString = listString.concat(",");
			int atomicNo = mAtomList[atom][i];
			listString = listString.concat(Molecule.cAtomLabel[atomicNo]);
			}

		return listString;
		}


	/**
	 * Returns an atom mapping number within the context of a reaction.
	 * Atoms that that share the same mapping number on the reactant and product side
	 * are considered to be the same atom.
	 * @param atom
	 * @return
	 */
	public int getAtomMapNo(int atom) {
		return Math.abs(mAtomMapNo[atom]);
		}


	/**
	 * @param atom
	 * @return atom mass, if is specific isotop, otherwise 0 for natural abundance
	 */
	public int getAtomMass(int atom) {
		return mAtomMass[atom];
		}


	/**
	 * The atom parity is a calculated property available above/equal helper level cHelperParities.
	 * It describes the stereo configuration of a chiral atom and is calculated either from
	 * 2D-atom-coordinates and up/down-bonds or from 3D-atom-coordinates, whatever is available.
	 * It depends on the atom indexes of the neighbor atoms and their orientation in space.<br>
	 * The parity is defined as follows: Look at the chiral atom such that its neighbor atom with the
	 * highest atom index (or the hydrogen atom if it is implicit) is oriented to the back.
	 * If the remaining three neighbors are in clockwise order (considering ascending atom indexes)
	 * than the parity is 1. If they are in anti-clockwise order, then the parity is 2.<br>
	 * For linear chirality (allenes): Look along the straight line of double bonds such that the
	 * rear neighbor with the lower atom index points to the top. If the front neighbor with the
	 * lower atom index points to the right than the parity is 1.<br>
	 * @param atom
	 * @return one of cAtomParity1,cAtomParity2,cAtomParityNone,cAtomParityUnknown
	 */
	public int getAtomParity(int atom) {
		return mAtomFlags[atom] & cAtomFlagsParity;
		}


	/**
	 * Returns all set query features for this atom. In order to get all features related to a certain subject
	 * use something like this: <i>getAtomQueryFeatures() & cAtomQFHydrogen</i>
	 * @param atom
	 * @return
	 */
	public int getAtomQueryFeatures(int atom) {
		return mAtomQueryFeatures[atom];
		}


	/**
	 * Gets an atom's radical state as singulet,dublet,triplet or none
	 * @param atom
	 * @return one of cAtomRadicalStateNone,cAtomRadicalStateS,cAtomRadicalStateD,cAtomRadicalStateT
	 */
	public int getAtomRadical(int atom) {
		return mAtomFlags[atom] & cAtomRadicalState;
		}


	public float getAtomX(int atom) {
		return mAtomX[atom];
		}


	public float getAtomY(int atom) {
		return mAtomY[atom];
		}


	public float getAtomZ(int atom) {
		return mAtomZ[atom];
		}


	public float getAverageBondLength() {
		return getAverageBondLength(mAllAtoms, mAllBonds);
		}



	public float getAverageBondLength(int atoms, int bonds) {
		for (int bond=0; bond<bonds; bond++)
			if ((mBondQueryFeatures[bond] & cBondQFBridge) != 0)
				bonds--;

		if (bonds == 0) {
				// since this function is used to get an idea about the scale
				// of the molecule return as approximation a mean atom distance
			if (mAllAtoms < 2)
				return defaultBondLength;

			float sum = 0.0f;
			int count = 0;
			for (int atom1=1; atom1<mAllAtoms; atom1++) {
				for (int atom2=0; atom2<atom1; atom2++) {
					float dx = mAtomX[atom1] - mAtomX[atom2];
					float dy = mAtomY[atom1] - mAtomY[atom2];
					float dz = mAtomZ[atom1] - mAtomZ[atom2];
					sum += Math.sqrt(dx*dx+dy*dy+dz*dz);
					count++;
					}
				}
			return Math.min(defaultBondLength, (float)Math.sqrt(mAllAtoms) * sum / (2f * count));
			}

		float avblSum = 0.0f;
		for (int bond=0; bond<bonds; bond++) {
			if ((mBondQueryFeatures[bond] & cBondQFBridge) == 0) {
				float dx = mAtomX[mBondAtom[1][bond]] - mAtomX[mBondAtom[0][bond]];
				float dy = mAtomY[mBondAtom[1][bond]] - mAtomY[mBondAtom[0][bond]];
				float dz = mAtomZ[mBondAtom[1][bond]] - mAtomZ[mBondAtom[0][bond]];
				avblSum += Math.sqrt(dx*dx+dy*dy+dz*dz);
				}
			}
		return avblSum / bonds;
		}


	public float getBondAngle(int atom1, int atom2) {
		return getAngle(mAtomX[atom1], mAtomY[atom1], mAtomX[atom2], mAtomY[atom2]);
		}


	/**
	 * @param no 0 or 1
	 * @param bond
	 * @return atom index
	 */
	public int getBondAtom(int no,int bond) {
		return mBondAtom[no][bond];
		}


	/**
	 * The bond Cahn-Ingold-Prelog parity is a calculated property available above/equal helper level cHelperCIP.
	 * It encodes the stereo configuration of a bond with its neighbors using 2D-coordinates and up/down-bonds
	 * or 3D-atom-coordinates, whatever is available. It depends on the atom indices of the neighbor
	 * atoms and their orientation is space. This method is called by the Canonizer and usually should not
	 * be called otherwise. Considered are E/Z-double bonds and M/P-BINAP type single bonds.
	 * @param bond
	 * @return one of cBondCIPParityNone,cBondCIPParityEorP,cBondCIPParityZorM,cBondCIPParityProblem
	 */
	public int getBondCIPParity(int bond) {
		return (mBondFlags[bond] & cBondFlagsCIPParity) >> cBondFlagsCIPParityShift;
		}


	/**
	 * This is MDL's enhanced stereo representation (ESR).
	 * Stereo atoms and bonds with the same ESR type (AND or OR) and the same ESR group number
	 * are in the same group, i.e. within this group they have the defined (relative) stereo configuration.
	 * @param bond
	 * @return group index starting with 0
	 */
	public int getBondESRGroup(int bond) {
		if (getBondESRType(bond) != cESRTypeAnd
		 && getBondESRType(bond) != cESRTypeOr)
			return -1;
		else
			return (mBondFlags[bond] & cBondFlagsESRGroup) >> cBondFlagsESRGroupShift;
		}


	/**
	 * This is MDL's enhanced stereo representation (ESR).
	 * Stereo atoms and bonds with the same ESR type (AND or OR) and the same ESR group number
	 * are in the same group, i.e. within this group they have the defined (relative) stereo configuration.
	 * @param bond
	 * @return one of cESRTypeAbs,cESRTypeAnd,cESRTypeOr
	 */
	public int getBondESRType(int bond) {
		return (mBondFlags[bond] & cBondFlagsESRType) >> cBondFlagsESRTypeShift;
		}


	/**
	 * @param bond
	 * @return bond length calculated from atom 2D-coordinates.
	 */
	public float getBondLength(int bond) {
		int atom1 = mBondAtom[0][bond];
		int atom2 = mBondAtom[1][bond];
		float xdif = mAtomX[atom2] - mAtomX[atom1];
		float ydif = mAtomY[atom2] - mAtomY[atom1];
		return (float)Math.sqrt(xdif * xdif + ydif * ydif);
		}


	/**
	 * Delocalized bonds, i.e. bonds in an aromatic 6-membered ring, are returned as 1.
	 * Ligand field bonds are returned as 0.
	 * @param bond
	 * @return for organic molecules 1,2, or 3
	 */
	public int getBondOrder(int bond) {
		switch (mBondType[bond] & cBondTypeMaskSimple) {
		case cBondTypeSingle:
		case cBondTypeDelocalized: return 1;
		case cBondTypeDouble: return 2;
		case cBondTypeTriple: return 3;
		default: return 0;	// metal bond
			}
		}


	/**
	 * Returns the pre-calculated bond parity, e.g. cBondParityEor1.
	 * To distinguish double bond parities (E/Z) from parities of axial
	 * chirality, e.g. BINAP type (1/2) simply check with getBondOrder(bond):
	 * If the order is 2, then the parity describes E/Z, otherwise an axial parity. 
	 * @param bnd
	 * @return one of cBondParity???
	 */
	public int getBondParity(int bnd) {
		return mBondFlags[bnd] & cBondFlagsParity;
		}


	public int getBondQueryFeatures(int bnd) {
		return mBondQueryFeatures[bnd];
		}


	public boolean isBondBridge(int bond) {
		return (mBondQueryFeatures[bond] & cBondQFBridge) != 0;
		}


	public int getBondBridgeMinSize(int bond) {
		return (mBondQueryFeatures[bond] & cBondQFBridgeMin) >> cBondQFBridgeMinShift;
		}


	public int getBondBridgeMaxSize(int bond) {
		return ((mBondQueryFeatures[bond] & cBondQFBridgeMin) >> cBondQFBridgeMinShift)
			 + ((mBondQueryFeatures[bond] & cBondQFBridgeSpan) >> cBondQFBridgeSpanShift);
		}


	/**
	 * Returns bond type combining bond order and stereo orientation.
	 * @param bond
	 * @return one of cBondTypeSingle,cBondTypeDouble,cBondTypeUp,cBondTypeCross,...
	 */
	public int getBondType(int bond) {
		return mBondType[bond];
		}


	/**
	 * This is the bond type without stereo information.
	 * @param bond
	 * @return cBondTypeSingle,cBondTypeDouble,cBondTypeTriple,cBondTypeDelocalized
	 */
	public int getBondTypeSimple(int bond) {
		return mBondType[bond] & cBondTypeMaskSimple;
		}

	
	/**
	 * Gets the overall chirality of the molecule, which is a calculated information considering:
	 * Recognition of stereo centers and stereo bonds, defined ESR features, meso detection.
	 * The chirality combines the knowledge about how many stereo isomers are represented,
	 * whether all of these are meso, whether we have one defined stereo isomer, a mixture
	 * of racemates, epimers, or other diastereomers.
	 * The information is used during depiction.
	 */
	public int getChirality() {
		return mChirality;
		}


	/**
	 * The currently defined maximum of atoms, which increases automatically when using high level
	 * construction methods and new atoms exceed the current maximum.
	 * @return
	 */
	public int getMaxAtoms() {
		return mMaxAtoms;
		}


	/**
	 * Usually called automatically and hardly needed to be called.
	 * @param v
	 */
	public void setMaxAtoms(int v) {
		mAtomicNo = Arrays.copyOf(mAtomicNo, v);
		mAtomCharge = Arrays.copyOf(mAtomCharge, v);
		mAtomMapNo = Arrays.copyOf(mAtomMapNo, v);
		mAtomX = Arrays.copyOf(mAtomX, v);
		mAtomY = Arrays.copyOf(mAtomY, v);
		mAtomZ = Arrays.copyOf(mAtomZ, v);
		mAtomMass = Arrays.copyOf(mAtomMass, v);
		mAtomFlags = Arrays.copyOf(mAtomFlags, v);
		mAtomQueryFeatures = Arrays.copyOf(mAtomQueryFeatures, v);
		if (mAtomList != null)
			mAtomList = Arrays.copyOf(mAtomList, v);
		if (mAtomCustomLabel != null)
			mAtomCustomLabel = Arrays.copyOf(mAtomCustomLabel, v);
		mMaxAtoms = v;
		}


	/**
	 * The currently defined maximum of bonds, which increases automatically when using high level
	 * construction methods and new bonds exceed the current maximum.
	 * @return
	 */
	public int getMaxBonds() {
		return mMaxBonds;
		}

	/**
	 * Usually called automatically and hardly needed to be called.
	 * @param v
	 */
	public void setMaxBonds(int v) {
		mBondAtom[0] = Arrays.copyOf(mBondAtom[0], v);
		mBondAtom[1] = Arrays.copyOf(mBondAtom[1], v);
		mBondType = Arrays.copyOf(mBondType, v);
		mBondFlags = Arrays.copyOf(mBondFlags, v);
		mBondQueryFeatures = Arrays.copyOf(mBondQueryFeatures, v);
		mMaxBonds = v;
		}


	/**
	 * Allows to set a molecule name or identifier, that is, for instance, written to or read from molfiles.
	 * @return
	 */
	public String getName() {
		return mName;
		}


	/**
	 * The stereo problem flag is set by the stereo recognition (available equal/above helper level cHelperParities)
	 * if an atom has over- or under-specified stereo bonds attached, i.e. a stereo center with less or more than one
	 * up/down-bond, an non-stereo-center atom carrying (a) stereo bond(s), or a stereo center with neighbors coordinates
	 * such that the stereo configuration cannot be deduced. This flag is used by the depiction and causes affected atoms
	 * to be drawn in margenta.
	 * @param atom
	 * @return
	 */
	public boolean getStereoProblem(int atom) {
		return ((mAtomFlags[atom] & cAtomFlagStereoProblem) != 0);
		}


	/**
	 * @param atom
	 * @return whether the atom's stereo configuration was explicitly deplared unknown
	 */
	public boolean isAtomConfigurationUnknown(int atom) {
		return ((mAtomFlags[atom] & cAtomFlagConfigurationUnknown) != 0);
		}


	/**
	 * Pseudo paries are parities that indicate a relative configuration.
	 * It always needs at least 2 pseudo parities (atom or bond) within
	 * a part of a molecule to be meaningful.
	 * This information is calculated by ensureHelperArrays(Molecule.cHelperCIP).
	 * Molecules extracted from IDCode don't know about pseudo parities.
	 * @param atom
	 * @return wether this atom's parity is a relative configuration
	 */
	public boolean isAtomParityPseudo(int atom) {
		return ((mAtomFlags[atom] & cAtomParityIsPseudo) != 0);
		}


	/**
	 * Atoms with pseudo parities are not considered stereo centers.
	 * While parities are canonized and always refer to the full set
	 * of molecules (in case ESR groups are defined), this method
	 * returns true if this atom is a stereo center in any(!) of the
	 * individual molecules described by the ESR settings.
	 * @param atom
	 * @return true if atom is stereo center in at least one molecule after ESR resolution
	 */
	public boolean isAtomStereoCenter(int atom) {
		return ((mAtomFlags[atom] & cAtomFlagIsStereoCenter) != 0);
		}

	
	public boolean isBondParityPseudo(int bond) {
		return ((mBondFlags[bond] & cBondParityIsPseudo) != 0);
		}

	
	/**
	 * This hint/flag is set by CoordinateInventor for double bonds without given EZ-parity,
	 * because the new coordinates may imply a not intended EZ-parity. If parities are calculated
	 * later by the Canonizer is can correctly assign cBondParityUnknown if the bond is a stereo bond. 
	 * The setBondParity() method clears this flag.
	 * This method usually should not be called for other purposes.
	 * @return whether the bond parity was unknown when 2D- atom coordinates were created
	 */
	public boolean isBondParityUnknownOrNone(int bond) {
		return ((mBondFlags[bond] & cBondParityUnknownOrNone) != 0);
		}


	/**
	 * Molecule objects may represent complete molecules or sub-structure fragments,
	 * depending on, whether they are flagges as being a fragment or not. Both representations
	 * have much in common, but in certain aspects behave differently. Thus, complete molecules
	 * are considered to carry implicit hydrogens to fill unoccupied atom valences.
	 * Sub-structure fragments on the other hand may carry atom or bond query features.
	 * Depiction, sub-structure search, and other algorithms treat fragments and complete molecules differerently.
	 * @return
	 */
	public boolean isFragment() {
		return mIsFragment;
		}


	/**
	 * @param atom
	 * @return whether the atom has the natural isotop distribution
	 */
	public boolean isNaturalAbundance(int atom) {
		return (mAtomMass[atom] == 0);
		}


	/**
	 * @return true if atom is one of H,B,C,N,O,F,Si,P,S,Cl,As,Se,Br,Te,I
	 */
	public boolean isPurelyOrganic() {
		for (int atom=0; atom<mAllAtoms; atom++) {
			switch (mAtomicNo[atom]) {
			case  1:	// H
			case  5:	// B
			case  6:	// C
			case  7:	// N
			case  8:	// O
			case  9:	// F
			case 14:	// Si
			case 15:	// P
			case 16:	// S
			case 17:	// Cl
			case 33:	// As
			case 34:	// Se
			case 35:	// Br
			case 52:	// Te
			case 53:	// I
				continue;
			default:
				return false;
				}
			}
		return true;
		}


	public boolean isSelectedAtom(int atom) {
		return ((mAtomFlags[atom] & cAtomFlagSelected) != 0) ? true : false;
		}


	/**
	 * Atom marking may be used for any external purpose
	 */
	public boolean isMarkedAtom(int atom) {
		return ((mAtomFlags[atom] & cAtomFlagMarked) != 0) ? true : false;
		}


	/**
	 * Used for depiction only.
	 * @param bond
	 */
	public boolean isBondBackgroundHilited(int bond) {
		return ((mBondFlags[bond] & cBondFlagBGHilited) != 0) ? true : false;
		}


	/**
	 * Used for depiction only.
	 * @param bond
	 */
	public boolean isBondForegroundHilited(int bond) {
		return ((mBondFlags[bond] & cBondFlagFGHilited) != 0) ? true : false;
		}


	public boolean isSelectedBond(int bond) {
		return ((mAtomFlags[mBondAtom[0][bond]]
			   & mAtomFlags[mBondAtom[1][bond]]
			   & cAtomFlagSelected) != 0) ? true : false;
		}


	public boolean isAutoMappedAtom(int atom) {
		return (mAtomMapNo[atom] < 0);
		}


	/**
	 * Checks whether bond is drawn as up/down single bond
	 * @param bond
	 * @return true if bond is a stereo bond
	 */
	public boolean isStereoBond(int bond) {
		return mBondType[bond] == cBondTypeUp || mBondType[bond] == cBondTypeDown;
		}

	
	/**
	 * Checks whether bond is drawn as up/down single bond and is connected to atom with its pointed tip
	 * @param bond
	 * @param atom
	 * @return true if bond is a stereo bond referring to atom
	 */
	public boolean isStereoBond(int bond, int atom) {
		return (mBondType[bond] == cBondTypeUp || mBondType[bond] == cBondTypeDown) && mBondAtom[0][bond] == atom;
		}


	/**
	 * Low level method for constructing/modifying a molecule from scratch.
	 * Use setAtomicNo(), possibly setAtomX(), setAtomY() and other setAtomXXX() methods for new atoms.
	 * @param no
	 */
	public void setAllAtoms(int no) {
		mAllAtoms = no;
		mValidHelperArrays = cHelperNone;
		}


	/**
	 * Low level method for constructing/modifying a molecule from scratch.
	 * Use setBondType() and setBondAtom() if you increase the number of bonds with this method.
	 * @param no
	 */
	public void setAllBonds(int no) {
		mAllBonds = no;
		mValidHelperArrays = cHelperNone;
		}

	/**
	 * Set an atom's maximum valance to be different from the default one.
	 * If a carbon atom's valence is set to -1,0 or 4 its radical state is removed.
	 * If a carbon atom's valence is set to 2, a singulet carbene state is assumed.
	 * @param atom
	 * @param valence 0-14: new maximum valence; -1: use default
	 */
	public void setAtomAbnormalValence(int atom, int valence) {
		if (valence >= -1 && valence <= 14) {
			mAtomFlags[atom] &= ~cAtomFlagsValence;
			if (valence != getDefaultMaxValenceUncharged(atom))
				mAtomFlags[atom] |= ((1+valence) << cAtomFlagsValenceShift);

			if (mAtomicNo[atom] == 6) {
				if (valence == -1 || valence == 0 || valence == 2 || valence == 4) {
					mAtomFlags[atom] &= ~cAtomRadicalState;
					if (valence == 2)
						mAtomFlags[atom] |= cAtomRadicalStateS;
					}
				}
			}
		}


	public void setAtomCharge(int atom, int charge) {
		mAtomCharge[atom] = charge;
		mValidHelperArrays = cHelperNone;
		}


	public void setAtomColor(int atom,int color) {
		mAtomFlags[atom] &= ~cAtomFlagsColor;
		mAtomFlags[atom] |= color;
		}


	/**
	 * This is a user applied information, rather than a calculated value.
	 * The stereo center configuration is declared to be unknown.
	 * If the atom is recognized a stereo center, then its parity will be cAtomParityUnknown.
	 * @param atom
	 * @param u
	 */
	public void setAtomConfigurationUnknown(int atom, boolean u) {
		if (u)
			mAtomFlags[atom] |= cAtomFlagConfigurationUnknown;
		else
			mAtomFlags[atom] &= ~cAtomFlagConfigurationUnknown;

		mValidHelperArrays &= cHelperRings;
		}


	public void setAtomSelection(int atom,boolean s) {
		if (s)
			mAtomFlags[atom] |= cAtomFlagSelected;
		else
			mAtomFlags[atom] &= ~cAtomFlagSelected;
		}


	/**
	 * Atom marking may be used for any external purpose
	 */
	public void setAtomMarker(int atom,boolean s) {
		if (s)
			mAtomFlags[atom] |= cAtomFlagMarked;
		else
			mAtomFlags[atom] &= ~cAtomFlagMarked;
		}


	/**
	 * Set an atom's atomic number and defines the isotop to be natural abundance.
	 * @param atom
	 * @param no
	 */
	public void setAtomicNo(int atom,int no) {
		if ((no >= 0) && (no <= cMaxAtomicNo)) {
			if (no == 151 || no == 152) {	// 'D' or 'T'
				mAtomicNo[atom] = 1;
				mAtomMass[atom] = no - 149;
				}
			else {
				mAtomicNo[atom] = no;
				mAtomMass[atom] = 0;
				}
			mAtomFlags[atom] &= ~cAtomFlagsValence;

			mValidHelperArrays = cHelperNone;
			}
		}


	/**
	 * Defines a list of allowed/excluded atomic numbers for sub-structure matching.
	 * If this atom's query feature cAtomQFAny (any atom) is set, then the list is considered to be a NOT-list.
	 * Depending on cAtomQFAny the list must contain at least 1 or 2 members.
	 * @param atom
	 * @param list null or int[] of valid unique, but not sorted, atomic numbers
	 */
	public void setAtomList(int atom, int[] list) {
		if (mAtomList == null)
			mAtomList = new int[mMaxAtoms][];

		if (list != null)
			Arrays.sort(list);

		mAtomList[atom] = list;

		mValidHelperArrays = cHelperNone;
		mIsFragment = true;
		}


	/**
	 * Defines an atom list as query feature for substructure search
	 * @param atom
	 * @param list is null or a sorted int[] of valid atomic numbers
	 * @param isExcludeList true if atom is a wild card and list contains atoms to be excluded 
	 */
	public void setAtomList(int atom, int[] list, boolean isExcludeList) {
		if (list == null) {
			if (mAtomList != null)
				mAtomList[atom] = null;
			return;
			}

		if (list.length == 1 && !isExcludeList) {
			int atomicNo = list[0];
			if (mAtomicNo[atom] != atomicNo)
				changeAtom(atom, atomicNo, 0, -1, 0);
			if (mAtomList != null)
				mAtomList[atom] = null;
			return;
			}

		if (mAtomList == null)
			mAtomList = new int[mMaxAtoms][];

		mAtomList[atom] = list;

		if (isExcludeList)
			mAtomQueryFeatures[atom] |= cAtomQFAny;

		mValidHelperArrays = cHelperNone;
		mIsFragment = true;
		}


	/**
	 * Defines an atom mapping number within the context of a reaction.
	 * Atoms that that share the same mapping number on the reactant and product side
	 * are considered to be the same atom.
	 * @param atom
	 * @param mapNo
	 * @param autoMapped
	 */
	public void setAtomMapNo(int atom, int mapNo, boolean autoMapped) {
		mAtomMapNo[atom] = (autoMapped) ? -mapNo : mapNo;
		}


	/**
	 * Set atom to specific isotop or to have a natural isotop distribution
	 * @param atom
	 * @param mass rounded atom mass or 0 (default) for natural abundance
	 */
	public void setAtomMass(int atom, int mass) {
		mAtomMass[atom] = mass;
		mValidHelperArrays &= cHelperRings;
		}


	/**
	 * The atom parity is a calculated property available above/equal helper level cHelperParities.
	 * It describes the stereo configuration of a chiral atom and is calculated either from
	 * 2D-atom-coordinates and up/down-bonds or from 3D-atom-coordinates, whatever is available.
	 * It depends on the atom indices of the neighbor atoms and their orientation is space.<br>
	 * The parity is defined as follows: Look at the chiral atom such that its neighbor atom with the
	 * highest atom index (or the hydrogen atom if it is implicit) is oriented to the back.
	 * If the remaining three neighbors are in clockwise order (considering ascending atom indexes)
	 * than the parity is 1. If they are in anti-clockwise order, then the parity is 2.<br>
	 * For linear chirality (allenes): Look along the straight line of double bonds such that the
	 * rear neighbor with the lower atom index points to the top. If the front neighbor with the
	 * lower atom index points to the right than the parity is 1.<br>
	 * This method is called by the Canonizer and usually should not be called otherwise.
	 * @param atom
	 * @param parity one of cAtomParity1,cAtomParity2,cAtomParityNone,cAtomParityUnknown
	 * @param isPseudo true if the configuration is only meaningful relative to another one
	 */
	public void setAtomParity(int atom, int parity, boolean isPseudo) {
		mAtomFlags[atom] &= ~(cAtomFlagsParity | cAtomParityIsPseudo);
		mAtomFlags[atom] |= parity;
		if (isPseudo)
			mAtomFlags[atom] |= cAtomParityIsPseudo;
		}


	/**
	 * An atom is considered a stereo center, if it is a stereo center in at least in one of the
	 * molecule configurations represented by the ESR definitions. Pseudo stereo centers are not(!)
	 * considered to be a stereo center.
	 * This method is called by the Canonizer and usually should not be called otherwise.
	 * @param atom
	 * @param isStereoCenter
	 */
	protected void setAtomStereoCenter(int atom, boolean isStereoCenter) {
		mAtomFlags[atom] &= ~cAtomFlagIsStereoCenter;
		if (isStereoCenter)
			mAtomFlags[atom] |= cAtomFlagIsStereoCenter;
		}


	/**
	 * Introduce or remove an atom query feature and make sure, the molecule is flagged
	 * to be a sub-structure fragment (see setFragment()).
	 * A query feature is usually a flag, which if set, poses an additional atom/bond matching constraint
	 * for the sub-structure search and, thus, reduces the number of matching atoms and therefore also
	 * the number of molecules found. Often multiple query feature flags are related and grouped, e.g.
	 * to define the number of hydrogens atoms. These are the flags related to hydrogen neighbors:<br><br>
	 * public static final int cAtomQFHydrogen		= 0x00000780;<br>
	 * public static final int cAtomQFNot0Hydrogen	= 0x00000080;<br>
	 * public static final int cAtomQFNot1Hydrogen	= 0x00000100;<br>
	 * public static final int cAtomQFNot2Hydrogen	= 0x00000200;<br>
	 * public static final int cAtomQFNot3Hydrogen	= 0x00000400;<br>
	 * <p>An inverse logic needs to be applied to translate a user request to the bits needed. For example,
	 * to only accept atoms that have 1 or 2 hydrogen neighbors, we need to filter out all others. Thus, we
	 * would call<br>setAtomQueryFeature(atom, cAtomQFNot0Hydrogen | cAtomQFNot3Hydrogen, true);</p>
	 * <p>To match only atoms without hydrogen neighbors, call<br>setAtomQueryFeature(atom, cAtomQFHydrogen & ~cAtomQFNot3Hydrogen, true);<br>
	 * This mechanism allows a very efficient atom matching and therefore very fast sub-structure search.</p>
	 * @param atom
	 * @param feature one of cAtomQF...
	 * @param value if true, the feature is set, otherwise it is removed
	 */
	public void setAtomQueryFeature(int atom, int feature, boolean value) {
		if (value)
			mAtomQueryFeatures[atom] |= feature;
		else
			mAtomQueryFeatures[atom] &= ~feature;
		mValidHelperArrays = 0;	// there is an influence on occipied valence, bond order, etc.
		mIsFragment = true;
		}


	/**
	 * Sets an atom's radical state as singulet,dublet,triplet or none
	 * @param atom
	 * @param radical one of cAtomRadicalStateNone,cAtomRadicalStateS,cAtomRadicalStateD,cAtomRadicalStateT
	 */
	public void setAtomRadical(int atom, int radical) {
		mAtomFlags[atom] &= ~cAtomRadicalState;
		mAtomFlags[atom] |= radical;
		mValidHelperArrays &= cHelperRings;
		}


	/**
	 * The atom Cahn-Ingold-Prelog parity is a calculated property available above/equal helper level cHelperCIP.
	 * It encodes the stereo configuration of an atom with its neighbors using up/down-bonds
	 * or 3D-atom-coordinates, whatever is available. It depends on the atom indices of the neighbor
	 * atoms and their orientation is space. This method is called by the Canonizer and usually should not
	 * be called otherwise.
	 * @param atom
	 * @param parity one of cAtomCIPParityRorM,cAtomCIPParitySorP,cAtomCIPParityProblem
	 */
	public void setAtomCIPParity(int atom, int parity) {
		mAtomFlags[atom] &= ~cAtomFlagsCIPParity;
		mAtomFlags[atom] |= (parity << cAtomFlagsCIPParityShift);
		}


	public void setAtomX(int atom, float x) {
		mAtomX[atom] = x;
		mValidHelperArrays &= cHelperRings;
		}


	public void setAtomY(int atom, float y) {
		mAtomY[atom] = y;
		mValidHelperArrays &= cHelperRings;
		}


	public void setAtomZ(int atom, float z) {
		mAtomZ[atom] = z;
		mValidHelperArrays &= cHelperRings;
		}


	public void setAtomX(int atom, double x) {
		setAtomX(atom, (float) x);
		}


	public void setAtomY(int atom, double y) {
		setAtomY(atom, (float) y);
		}


	public void setAtomZ(int atom, double z) {
		setAtomZ(atom, (float) z);
		}


	public void setBondAtom(int no, int bond, int atom) {
		mBondAtom[no][bond] = atom;
		mValidHelperArrays = cHelperNone;
		}


	/**
	 * The bond Cahn-Ingold-Prelog parity is a calculated property available above/equal helper level cHelperCIP.
	 * It encodes the stereo configuration of a bond with its neighbors using 2D-coordinates and up/down-bonds
	 * or 3D-atom-coordinates, whatever is available. It depends on the atom indices of the neighbor
	 * atoms and their orientation is space. This method is called by the Canonizer and usually should not
	 * be called otherwise. Considered are E/Z-double bonds and M/P-BINAP type single bonds.
	 * @param bond
	 * @param parity one of cBondCIPParityEorP,cBondCIPParityZorM,cBondCIPParityProblem
	 */
	public void setBondCIPParity(int bond, int parity) {
		mBondFlags[bond] &= ~cBondFlagsCIPParity;
		mBondFlags[bond] |= (parity << cBondFlagsCIPParityShift);
		}


	/**
	 * Used for depiction only.
	 * @param bond
	 * @param s
	 */
	public void setBondBackgroundHiliting(int bond, boolean s) {
		if (s)
			mBondFlags[bond] |= cBondFlagBGHilited;
		else
			mBondFlags[bond] &= ~cBondFlagBGHilited;
		}


	/**
	 * Used for depiction only.
	 * @param bond
	 * @param s
	 */
	public void setBondForegroundHiliting(int bond, boolean s) {
		if (s)
			mBondFlags[bond] |= cBondFlagFGHilited;
		else
			mBondFlags[bond] &= ~cBondFlagFGHilited;
		}


	/**
	 * The bond parity is a calculated property available above/equal helper level cHelperParities.
	 * It encodes the stereo configuration of a double bond or BINAP type single bond from up/down-bonds
	 * and 2D-coordinates or 3D-atom-coordinates, whatever is available. It depends on the atom indices
	 * of the neighbor atoms and their orientation is space. This method is called by the Canonizer and
	 * usually should not be called otherwise.
	 * @param bond
	 * @param parity one of cBondParityEor1,cBondParityZor2,cBondParityNone,cBondParityUnknown
	 * @param isPseudo true if the configuration is only meaningful relative to another one
	 */
	public void setBondParity(int bond, int parity, boolean isPseudo) {
		mBondFlags[bond] &= ~(cBondFlagsParity | cBondParityIsPseudo | cBondParityUnknownOrNone);
		mBondFlags[bond] |= parity;
		if (isPseudo)
			mBondFlags[bond] |= cBondParityIsPseudo;
		}


	/**
	 * This hint/flag is set by CoordinateInventor for double bonds without given EZ-parity,
	 * because the new coordinates may imply a not intended EZ-parity. If parities are calculated
	 * later by the Canonizer is can correctly assign cBondParityUnknown if the bond is a stereo bond. 
	 * The setBondParity() method clears this flag.
	 * This method usually should not be called for other purposes.
	 * @param bond
	 */
	public void setBondParityUnknownOrNone(int bond) {
		mBondFlags[bond] |= cBondParityUnknownOrNone;
		}


	public void setBondQueryFeature(int bond, int feature, boolean value) {
		if (value)
			mBondQueryFeatures[bond] |= feature;
		else
			mBondQueryFeatures[bond] &= ~feature;
		mValidHelperArrays = 0;	// there is an influence on occipied valence, bond order, etc.
		mIsFragment = true;
		}


	/**
	 * Defines a bond type combining bod order and stereo orientation.
	 * @param bond
	 * @param type one of cBondTypeSingle,cBondTypeDouble,cBondTypeUp,cBondTypeCross,...
	 */
	public void setBondType(int bond,int type) {
		mBondType[bond] = type;
		mValidHelperArrays = cHelperNone;
		}


	/**
	 * Sets the overall chirality of the molecule taking into account:
	 * Recognition of stereo centers and stereo bonds, defined ESR features, meso detection.
	 * The chirality combines the knowledge about how many stereo isomers are represented,
	 * whether all of these are meso, whether we have one defined stereo isomer, a mixture
	 * of racemates, epimers, or other diastereomers.
	 * The information is used during depiction.
	 * This method is called by the Canonizer and usually should not be called otherwise.
	 * @param c
	 */
	public void setChirality(int c) {
		mChirality = c;
		}


	/**
	 * Fragment's query features are checked for consistency and normalized
	 * during helper array creation. As part of this, simple hydrogen atoms
	 * are converted into hydrogen-count query features. If hydrogen protection
	 * is enabled, explicit hydrogens are not touched.
	 * @param protectHydrogen
	 */
	public void setHydrogenProtection(boolean protectHydrogen) {
		mProtectHydrogen = protectHydrogen;
		}


	/**
	 * This is for compatibility with old MDL stereo representation
	 * that contained a 'chiral' flag to indicate that the molecule
	 * is not a racemate. If a molecule is constructed from a source
	 * format (e.g. a molfile version 2) that contains a 'chiral' flag
	 * then setToRacemate() needs to be called if the chiral flag is
	 * not(!) set. This causes after stereo center recognition to
	 * turn all absolute stereo centers into racemic ones.
	 */
	public void setToRacemate() {
		mIsRacemate = true;
		}

	/**
	 * If a custom atom label is set, a molecule depiction displays
	 * the custom label instead of the original one. Custom labels
	 * are not interpreted otherwise. However, they may optionally
	 * be encoded into idcodes; see Canonizer.encodeAtomCustomLabels().
	 * @param atom
	 * @param label null to remove custom label
	 */
	public void setAtomCustomLabel(int atom, byte[] label) {
		if (label != null && label.length == 0)
			label = null;

		if (label == null) {
			if (mAtomCustomLabel != null)
				mAtomCustomLabel[atom] = null;
			}
		else {
			if (mAtomCustomLabel == null)
				mAtomCustomLabel = new byte[mMaxAtoms][];
			mAtomCustomLabel[atom] = label;
			}
		}

	/**
	 * If a custom atom label is set, a molecule depiction displays
	 * the custom label instead of the original one. Custom labels
	 * are not interpreted otherwise. However, they may optionally
	 * be encoded into idcodes; see Canonizer.encodeAtomCustomLabels().
	 * This label equals the normal atom label, then the custom label
	 * is removed. This method is less efficient than the byte[] version:
	 * setAtomCustomLabel(int, byte[])
	 * @param atom
	 * @param label null to remove custom label
	 */
	public void setAtomCustomLabel(int atom, String label) {
		if (label != null) {
			if (label.length() == 0)
				label = null;
			else {
				int atomicNo = getAtomicNoFromLabel(label);
				if ((atomicNo != 0 && label.equals(cAtomLabel[atomicNo]))
				 || label.equals("?")) {
					setAtomicNo(atom, atomicNo);
					label = null;
					}
				}
			}

		if (label == null) {
			if (mAtomCustomLabel != null)
				mAtomCustomLabel[atom] = null;
			}
		else {
			if (mAtomCustomLabel == null)
				mAtomCustomLabel = new byte[mMaxAtoms][];
			mAtomCustomLabel[atom] = label.getBytes();
			}
		}

	/**
	 * This is MDL's enhanced stereo representation (ESR).
	 * Stereo atoms and bonds with the same ESR type (AND or OR) and the same ESR group number
	 * are in the same group, i.e. within this group they have the defined (relative) stereo configuration.
	 * @param atom
	 * @param type one of cESRTypeAbs,cESRTypeAnd,cESRTypeOr
	 * @param group index starting with 0
	 */
	public void setAtomESR(int atom, int type, int group) {
		if (type == cESRTypeAbs) {
			mAtomFlags[atom] &= ~cAtomFlagsESR;
			mAtomFlags[atom] |= (type << cAtomFlagsESRTypeShift);
			}
		else {
			if (group >= cESRMaxGroups)
				return;

			if (group == -1) {  // find unused group
				int maxGroup = -1;
				for (int i=0; i<mAllAtoms; i++)
					if (i != atom
					 && type == getAtomESRType(i)
					 && maxGroup < getAtomESRGroup(i))
						maxGroup = getAtomESRGroup(i);
				for (int i=0; i<mAllBonds; i++)
					if (type == getBondESRType(i)
					 && maxGroup < getBondESRGroup(i))
						maxGroup = getBondESRGroup(i);
				group = maxGroup + 1;
				if (group >= cESRMaxGroups)
					return;
				}
			mAtomFlags[atom] &= ~cAtomFlagsESR;
			mAtomFlags[atom] |= ((type << cAtomFlagsESRTypeShift)
							   | (group << cAtomFlagsESRGroupShift));
			}

		mValidHelperArrays &= cHelperRings;
		}


	/**
	 * MDL's enhanced stereo representation for BINAP type of stereo bonds.
	 * Stereo atoms and bonds with the same ESR type (AND or OR) and the same ESR group number
	 * are in the same group, i.e. within this group they have the defined (relative) stereo configuration.
	 * @param bond
	 * @param type one of cESRTypeAbs,cESRTypeAnd,cESRTypeOr
	 * @param group index starting with 0
	 */
	public void setBondESR(int bond, int type, int group) {
		if (type == cESRTypeAbs) {
			mBondFlags[bond] &= ~cBondFlagsESR;
			mBondFlags[bond] |= (type << cBondFlagsESRTypeShift);
			}
		else {
			if (group >= cESRMaxGroups)
				return;

			if (group == -1) {  // find unused group
				int maxGroup = -1;
				for (int i=0; i<mAllAtoms; i++)
					if (type == getAtomESRType(i)
					 && maxGroup < getAtomESRGroup(i))
						maxGroup = getAtomESRGroup(i);
				for (int i=0; i<mAllBonds; i++)
					if (i != bond
					 && type == getBondESRType(i)
					 && maxGroup < getBondESRGroup(i))
						maxGroup = getBondESRGroup(i);
				group = maxGroup + 1;
				if (group >= cESRMaxGroups)
					return;
				}
			mBondFlags[bond] &= ~cBondFlagsESR;
			mBondFlags[bond] |= ((type << cBondFlagsESRTypeShift)
							   | (group << cBondFlagsESRGroupShift));
			}

		mValidHelperArrays &= cHelperRings;
		}


	/**
	 * Molecule objects may represent complete molecules or sub-structure fragments,
	 * depending on, whether they are flagges as being a fragment or not. Both representations
	 * have much in common, but in certain aspects behave differently. Thus, complete molecules
	 * are considered to carry implicit hydrogens to fill unoccupied atom valences.
	 * Sub-structure fragments on the other hand may carry atom or bond query features.
	 * Depiction, sub-structure search, and other algorithms treat fragments and complete molecules differerently.
	 * @param isFragment if false, then all query features are removed
	 */
	public void setFragment(boolean isFragment) {
		mIsFragment = isFragment;

		if (!isFragment) {
			mAtomList = null;
			for (int atom=0; atom<mAllAtoms; atom++)
				mAtomQueryFeatures[atom] = 0;
			for (int bond=0; bond<mAllBonds; bond++) {
				mBondQueryFeatures[bond] = 0;
				if (mBondType[bond] == cBondTypeDelocalized)
					mBondType[bond] = cBondTypeSingle;
				}
			}

		mValidHelperArrays = 0;
		}


	public void setName(String name) {
		mName = name;
		}


	/**
	 * The stereo problem flag is set by the stereo recognition (available equal/above helper level cHelperParities)
	 * if an atom has over- or under-specified stereo bonds attached, i.e. a stereo center with less or more than one
	 * up/down-bond, an non-stereo-center atom carrying (a) stereo bond(s), or a stereo center with neighbors coordinates
	 * such that the stereo configuration cannot be deduced. This flag is used by the depiction and causes affected atoms
	 * to be drawn in margenta.
	 * This method usually is not used for other purposees.
	 * @param atom
	 */
	protected void setStereoProblem(int atom) {
		mAtomFlags[atom] |= cAtomFlagStereoProblem;
		}


	/**
	 * Removes all isotop information, i.e. sets all atoms to the natural isotop abundance.
	 * @return true if something was changed
	 */
	public boolean stripIsotopInfo() {
		boolean found = false;
		boolean hydrogenIsotopFound = false;
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (mAtomMass[atom] != 0) {
				mAtomMass[atom] = 0;
				found = true;
				if (mAtomicNo[atom] == 1)
					hydrogenIsotopFound = true;
				}
			}

		if (hydrogenIsotopFound)	// D and T are treated like non hydrogen atoms
			mValidHelperArrays = cHelperNone;

		return found;
		}


	/**
	 * Removes all user supplied stereo information: up/down-bonds,
	 * atom configuration unknown flags, atom and bond ESR assignments.
	 * Removes parity and CIP flags from helper state, such that stereo calculation
	 * is redone, if needed.
	 */
	public void stripStereoInformation() {
		mIsRacemate = false;
		for (int atom=0; atom<mAllAtoms; atom++)
			mAtomFlags[atom] &= ~(cAtomFlagConfigurationUnknown | cAtomFlagsESR);
		for (int bond=0; bond<mAllBonds; bond++)
			mBondType[bond] &= ~cBondTypeMaskStereo;

		mValidHelperArrays &= ~cHelperBitsStereo;
		}


	public void translateCoords(float dx, float dy) {
		for (int i=0; i<mAllAtoms; i++) {
			mAtomX[i] += dx;
			mAtomY[i] += dy;
			}
		mZoomRotationX += dx;
		mZoomRotationY += dy;
		}


	public void scaleCoords(float f) {
		for (int i=0; i<mAllAtoms; i++) {
			mAtomX[i] *= f;
			mAtomY[i] *= f;
			}
		}


	public void zoomAndRotateInit(float x, float y) {
		mZoomRotationX = x;
		mZoomRotationY = y;
		mOriginalAngle = new float[mAllAtoms];
		mOriginalDistance = new float[mAllAtoms];
		for (int atom=0; atom<mAllAtoms; atom++) {
			float dx = x - mAtomX[atom];
			float dy = y - mAtomY[atom];
			mOriginalDistance[atom] = (float)Math.sqrt(dx*dx+dy*dy);	// distance to center of gravity
			mOriginalAngle[atom] = getAngle(x, y, mAtomX[atom], mAtomY[atom]);
			}
		}


	public void zoomAndRotate(float zoom, float angle, boolean selected) {
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (!selected || isSelectedAtom(atom)) {
				float newDistance = mOriginalDistance[atom] * zoom;
				float newAngle = mOriginalAngle[atom] - angle;
				mAtomX[atom] = mZoomRotationX + newDistance*(float)Math.sin(newAngle);
				mAtomY[atom] = mZoomRotationY + newDistance*(float)Math.cos(newAngle);
				}
			}

		if (selected)	// if only parts of the molecule were rotated
						// the stereo calculations may not be valid anymore
			mValidHelperArrays &= cHelperRings;
		}


	private int getMaximumBondOrder(int bond) {
		int maxBondOrder = 3;
		for (int i=0; i<2; i++) {
			int atom = mBondAtom[i][bond];
			int max = getBondOrder(bond) + getMaxValence(atom) - getOccupiedValence(atom);
			if (maxBondOrder > max)
				maxBondOrder = max;
			}
		return maxBondOrder;
		}


	private boolean incrementBondOrder(int bond) {
		int maxBondOrder = getMaximumBondOrder(bond);

		if (mBondType[bond] == cBondTypeTriple) {
			mBondType[bond] = cBondTypeSingle;
			mValidHelperArrays = cHelperNone;
			return true;
			}

		if (mBondType[bond] == cBondTypeDouble) {
			mBondType[bond] = cBondTypeCross;
			mValidHelperArrays &= cHelperRings;
			if ((mBondFlags[bond] & cBondFlagSmallRing) == 0)
				return true;
			}

		if (mBondType[bond] == cBondTypeCross) {
			if (maxBondOrder == 3)
				mBondType[bond] = cBondTypeTriple;
			else
				mBondType[bond] = cBondTypeSingle;
			mValidHelperArrays = cHelperNone;
			return true;
			}

		if ((cBondTypeMaskStereo & mBondType[bond]) != 0) {
			mBondType[bond] = cBondTypeSingle;
			mValidHelperArrays &= cHelperRings;
			return true;
			}

		if (maxBondOrder < 2)
			return false;

		if (mBondType[bond] == cBondTypeSingle) {
			mBondType[bond] = cBondTypeDouble;
			mValidHelperArrays = cHelperNone;
			return true;
			}

		return false;
		}


	protected boolean validateBondType(int bond, int type) {
		int simpleType = type & Molecule.cBondTypeMaskSimple;

		int maxBondOrder = getMaximumBondOrder(bond);
		switch (simpleType) {
		case cBondTypeSingle:
		case cBondTypeDelocalized:
			return maxBondOrder >= 1;
		case cBondTypeDouble:
			return maxBondOrder >= 2;
		case cBondTypeTriple:
			return maxBondOrder >= 3;
		case cBondTypeMetalLigand:
			return isMetalAtom(mBondAtom[0][bond]) ^ isMetalAtom(mBondAtom[1][bond]);
		default:
			return false;
			}
		}


	/**
	 * The sum of bond orders of explicitly connected neighbour atoms.
	 * @param atom
	 * @return explicitly used valence
	 */
	protected int getOccupiedValence(int atom) {
			// ExtendedMolecule overwrites this function by updateing and using neighbour arrays
		return simpleGetValence(atom);
		}


	private int simpleGetValence(int atom) {
		int val=0;		// this version doesn't change helper arrays
		for (int bnd=0; bnd<mAllBonds; bnd++)
			if (mBondAtom[0][bnd] == atom
			 || mBondAtom[1][bnd] == atom)
				val += getBondOrder(bnd);
		return val;
		}


	/**
	 * This is the defined maximum valence (or set abnormal valence)
	 * neglecting atom charge or radical influences, e.g. N or N(+) -> 3.
	 * @param atom
	 * @return
	 */
	public int getMaxValenceUncharged(int atom) {
		int valence = getAtomAbnormalValence(atom);

		if (valence == -1)
			valence = getDefaultMaxValenceUncharged(atom);

		return valence;
		}

	
	/**
	 * This is the default maximum valence of the atom
	 * neglecting atom charge or radical influences, e.g. N or N(+) -> 3.
	 * If the atomic no has multiple valid max valences, it is the highest one.
	 * @param atom
	 * @return
	 */
	public int getDefaultMaxValenceUncharged(int atom) {
		byte[] valenceList = (mAtomicNo[atom] < cAtomValence.length) ?
							 cAtomValence[mAtomicNo[atom]] : null;
		return (valenceList == null) ? cDefaultAtomValence : valenceList[valenceList.length-1];
		}


	/**
	 * This is the defined maximum valence (or set abnormal valence)
	 * corrected by atom charge or radical influences, e.g. N(+) -> 4.
	 * @param atom
	 * @return
	 */
	public int getMaxValence(int atom) {
		return getMaxValenceUncharged(atom) + getElectronValenceCorrection(atom);
		}


	/**
	 * This is the maximum valence correction caused by atom charge
	 * or radical status, e.g. N+ -> 1; N- -> -1; Al+ -> -1; C+,C- -> -1
	 * @param atom
	 * @return
	 */
	public int getElectronValenceCorrection(int atom) {
		if (mAtomicNo[atom] >= 171 && mAtomicNo[atom] <= 190)
			return 0;

		int correction = 0;

		if ((mAtomFlags[atom] & cAtomRadicalState) == cAtomRadicalStateD)
			correction -= 1;
		if ((mAtomFlags[atom] & cAtomRadicalState) == cAtomRadicalStateS
		 || (mAtomFlags[atom] & cAtomRadicalState) == cAtomRadicalStateT)
			correction -= 2;

		int charge = mAtomCharge[atom];
		if (charge == 0 && mIsFragment) {
			if ((mAtomQueryFeatures[atom] & cAtomQFCharge) == cAtomQFNotCharge0+cAtomQFNotChargePos)
				charge = -1;
			if ((mAtomQueryFeatures[atom] & cAtomQFCharge) == cAtomQFNotCharge0+cAtomQFNotChargeNeg)
				charge = 1;
			}
		if (mAtomicNo[atom] == 6)
			correction -= Math.abs(charge);
		else if (isElectronegative(atom))
			correction += charge;
		else
			correction -= charge;

		return correction;
		}


	public static boolean isAtomicNoElectronegative(int atomicNo) {
		switch (atomicNo) {
		case  7: 	// N
		case  8: 	// O
		case  9: 	// F
		case 15: 	// P
		case 16: 	// S
		case 17: 	// Cl
		case 33: 	// As
		case 34: 	// Se
		case 35: 	// Br
		case 53: 	// I
			return true;
			}
		return false;
		}

	/**
	 * @param atom
	 * @return whether atom is an electronegative one
	 */
	public boolean isElectronegative(int atom) {
		return isAtomicNoElectronegative(mAtomicNo[atom]);
		}


	public static boolean isAtomicNoElectropositive(int atomicNo) {
		if (atomicNo == 1
		 || atomicNo == 6)
			return false;

		if (isAtomicNoElectronegative(atomicNo))
			return false;

		if (atomicNo == 2	// He
		 || atomicNo == 10	// Ne
		 || atomicNo == 18	// Ar
		 || atomicNo == 36	// Kr
		 || atomicNo == 54)	// Xe
			return false;

		if (atomicNo > 103)	// amino acids etc.
			return false;

		return true;
		}


	/**
	 * @param atom
	 * @return whether atom is an electropositive one
	 */
	public boolean isElectropositive(int atom) {
		return isAtomicNoElectropositive(mAtomicNo[atom]);
		}


	/**
	 * @param atom
	 * @return whether atom is any metal atom
	 */
	public boolean isMetalAtom(int atom) {
		int atomicNo = mAtomicNo[atom];
		return (atomicNo >=  3 && atomicNo <=  4)
			|| (atomicNo >= 11 && atomicNo <= 13)
			|| (atomicNo >= 19 && atomicNo <= 31)
			|| (atomicNo >= 37 && atomicNo <= 51)
			|| (atomicNo >= 55 && atomicNo <= 84)
			|| (atomicNo >= 87 && atomicNo <= 103);
		}


	/**
	 * @param atom
	 * @return true if this atom is not a metal and not a nobel gas
	 */
	public boolean isOrganicAtom(int atom) {
		int atomicNo = mAtomicNo[atom];
		return atomicNo == 1
			|| (atomicNo >=  5 && atomicNo <=  9)	// B,C,N,O,F
			|| (atomicNo >= 14 && atomicNo <= 17)	// Si,P,S,Cl
			|| (atomicNo >= 33 && atomicNo <= 35)	// As,Se,Br
			|| (atomicNo >= 52 && atomicNo <= 53);	// Te,I
		}

	
	protected void removeMappingNo(int mapNo) {
		for (int atom=0; atom<mAllAtoms; atom++)
			if (Math.abs(mAtomMapNo[atom]) == Math.abs(mapNo))
				mAtomMapNo[atom] = 0;
		}


	protected int[] compressMolTable() {
		// neutralize charges, if after deletion of a bi-polar bond one atom remains
		for (int bnd=0; bnd<mAllBonds; bnd++) {
			if (mBondType[bnd] == cBondTypeDeleted) {
				int atom1 = mBondAtom[0][bnd];
				int atom2 = mBondAtom[1][bnd];
				if (mAtomicNo[atom1] == -1
				  ^ mAtomicNo[atom2] == -1) {
					if (mAtomCharge[atom1] != 0
					 && mAtomCharge[atom2] != 0) {
						if (mAtomCharge[atom1] < 0
						  ^ mAtomCharge[atom2] < 0) {
							if (mAtomCharge[atom1] < 0) {
								mAtomCharge[atom1]++;
								mAtomCharge[atom2]--;
								}
							else {
								mAtomCharge[atom1]--;
								mAtomCharge[atom2]++;
								}
							}
						}
					}
				}
			}
		
		int newAtmNo[] = new int[mAllAtoms];
		int atomDest = 0;
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (mAtomicNo[atom] == -1) {
				newAtmNo[atom] = -1;
				continue;
				}
			if (atomDest < atom) {
				mAtomicNo[atomDest] = mAtomicNo[atom];
				mAtomCharge[atomDest] = mAtomCharge[atom];
				mAtomMass[atomDest] = mAtomMass[atom];
				mAtomFlags[atomDest] = mAtomFlags[atom];
				mAtomQueryFeatures[atomDest] = mAtomQueryFeatures[atom];
				mAtomMapNo[atomDest] = mAtomMapNo[atom];
				mAtomX[atomDest] = mAtomX[atom];
				mAtomY[atomDest] = mAtomY[atom];
				mAtomZ[atomDest] = mAtomZ[atom];
				if (mAtomList != null)
					mAtomList[atomDest] = mAtomList[atom];
				if (mAtomCustomLabel != null)
					mAtomCustomLabel[atomDest] = mAtomCustomLabel[atom];
				}
			newAtmNo[atom] = atomDest;
			atomDest++;
			}
		mAllAtoms = atomDest;

		int bondDest = 0;
		for (int bnd=0; bnd<mAllBonds; bnd++) {
			if (mBondType[bnd] == cBondTypeDeleted) continue;
			mBondType[bondDest] = mBondType[bnd];
			mBondFlags[bondDest] = mBondFlags[bnd];
			mBondQueryFeatures[bondDest] = mBondQueryFeatures[bnd];
			mBondAtom[0][bondDest] = newAtmNo[mBondAtom[0][bnd]];
			mBondAtom[1][bondDest] = newAtmNo[mBondAtom[1][bnd]];
			bondDest++;
			}
		mAllBonds = bondDest;

		return newAtmNo;
		}


	private void polygon(int atom, int bonds, int endAtm, boolean aromatic,
						 float actlAngle, float angleChange) {
		boolean dblBnd;
		int actlAtm,remoteAtm,bnd;
		float bondLength,xdiff,ydiff,newx,newy;

		bondLength = 0.0f;
		if (atom == endAtm) {
			bondLength = getAverageBondLength();
			}
		else {
			xdiff = mAtomX[atom] - mAtomX[endAtm];
			ydiff = mAtomY[atom] - mAtomY[endAtm];
			bondLength = (float)Math.sqrt(xdiff * xdiff + ydiff * ydiff);
			}

		actlAtm = atom;
		dblBnd = !(simpleGetValence(atom) == 3);

		for (int step=1; step<bonds; step++) {
			newx = mAtomX[actlAtm] + bondLength * (float)Math.sin(actlAngle);
			newy = mAtomY[actlAtm] + bondLength * (float)Math.cos(actlAngle);
			remoteAtm = -1;
			for (int i=0; i<mAllAtoms; i++) {
				if ((Math.abs(newx - mAtomX[i]) < 4)
				 && (Math.abs(newy - mAtomY[i]) < 4)) {
					remoteAtm = i;
					break;
					}
				}
			if (remoteAtm == -1) {
				remoteAtm = addAtom(newx,newy);
				mAtomX[remoteAtm] = newx;
				mAtomY[remoteAtm] = newy;
				mAtomZ[remoteAtm] = 0;
				}
			bnd = getBondNo(actlAtm,remoteAtm);
			if (bnd == -1) {
				bnd = addBond(actlAtm,remoteAtm,1);
				if (aromatic) {
					if (dblBnd) {
						if ((simpleGetValence(mBondAtom[0][bnd]) < 4)
						 && (simpleGetValence(mBondAtom[1][bnd]) < 3))
							mBondType[bnd] = cBondTypeDouble;
						}
					dblBnd = !dblBnd;
					}
				}
			actlAtm = remoteAtm;
			actlAngle += angleChange;
			}
		bnd = getBondNo(actlAtm,endAtm);
		if (bnd == -1)
			bnd = addBond(actlAtm,endAtm,1);

		if (aromatic)
			if (dblBnd)
				if ((simpleGetValence(mBondAtom[0][bnd]) < 4)
				 && (simpleGetValence(mBondAtom[1][bnd]) < 4))
					mBondType[bnd] = cBondTypeDouble;
		}


	private int getBondNo(int atm1,int atm2) {
		for (int bnd=0; bnd<mAllBonds; bnd++)
			if (((mBondAtom[0][bnd] == atm1)
			  && (mBondAtom[1][bnd] == atm2))
			 || ((mBondAtom[0][bnd] == atm2)
			  && (mBondAtom[1][bnd] == atm1)))
				if (mBondType[bnd] != cBondTypeDeleted)
					return bnd;
		return -1;
		}


	private void writeObject(ObjectOutputStream stream) throws IOException {		
		stream.writeInt(mAllAtoms);
		stream.writeInt(mAllBonds);
		stream.writeBoolean(mIsFragment);
		for (int atom=0; atom<mAllAtoms; atom++) {
			stream.writeInt(mAtomicNo[atom]);
			stream.writeInt(mAtomCharge[atom]);
			stream.writeInt(mAtomMass[atom]);
			stream.writeInt(mAtomFlags[atom] & ~cAtomFlagsHelper);
			stream.writeInt(mAtomQueryFeatures[atom]);
			stream.writeDouble(mAtomX[atom]);	// for compatibility with earlier double based coords
			stream.writeDouble(mAtomY[atom]);
			stream.writeDouble(mAtomZ[atom]);
			stream.writeInt(mAtomMapNo[atom]);

			if (mAtomList != null && mAtomList[atom] != null) {
				stream.writeInt(mAtomList[atom].length);
				for (int i=0; i<mAtomList[atom].length; i++)
					stream.writeInt(mAtomList[atom][i]);
				}
			else
				stream.writeInt(0);

			if (mAtomCustomLabel != null && mAtomCustomLabel[atom] != null) {
				stream.writeInt(mAtomCustomLabel[atom].length);
				for (int i=0; i<mAtomCustomLabel[atom].length; i++)
					stream.writeByte(mAtomCustomLabel[atom][i]);
				}
			else
				stream.writeInt(0);
			}
		for (int bond=0; bond<mAllBonds;bond++) {
			stream.writeInt(mBondAtom[0][bond]);
			stream.writeInt(mBondAtom[1][bond]);
			stream.writeInt(mBondType[bond]);
			stream.writeInt(mBondFlags[bond]);
			stream.writeInt(mBondQueryFeatures[bond]);
			}
		
		stream.writeObject(mName);
		}


	private void readObject(ObjectInputStream stream) throws IOException {
		mAllAtoms = stream.readInt();
		mAllBonds = stream.readInt();

		mMaxAtoms = mAllAtoms;
		mMaxBonds = mAllBonds;
		init();

		mIsFragment = stream.readBoolean();

		for (int atom=0; atom<mAllAtoms; atom++) {
			mAtomicNo[atom] = stream.readInt();
			mAtomCharge[atom] = stream.readInt();
			mAtomMass[atom] = stream.readInt();
			mAtomFlags[atom] = stream.readInt();
			mAtomQueryFeatures[atom] = stream.readInt();
			mAtomX[atom] = (float)stream.readDouble();	// for compatibility with earlier double based coords
			mAtomY[atom] = (float)stream.readDouble();
			mAtomZ[atom] = (float)stream.readDouble();
			mAtomMapNo[atom] = stream.readInt();

			int count = stream.readInt();
			if (count != 0) {
				if (mAtomList == null)
					mAtomList = new int[mMaxAtoms][];
				mAtomList[atom] = new int[count];
				for (int i=0; i<count; i++)
					mAtomList[atom][i] = stream.readInt();
				}

			count = stream.readInt();
			if (count != 0) {
				if (mAtomCustomLabel == null)
					mAtomCustomLabel = new byte[mMaxAtoms][];
				mAtomCustomLabel[atom] = new byte[count];
				for (int i=0; i<count; i++)
					mAtomCustomLabel[atom][i] = stream.readByte();
				}
			}

		for (int bond=0; bond<mAllBonds; bond++) {
			mBondAtom[0][bond] = stream.readInt();
			mBondAtom[1][bond] = stream.readInt();
			mBondType[bond] = stream.readInt();
			mBondFlags[bond] = stream.readInt();
			mBondQueryFeatures[bond] = stream.readInt();
			}
		
		try {
			mName = (String) stream.readObject();
			} catch(Exception e) {}

		mValidHelperArrays = cHelperNone;
		}

}
