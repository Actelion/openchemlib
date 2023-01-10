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

import com.actelion.research.util.DoubleFormat;

/**
 * Typically you should use IDCodeParser instead of this class. You may instantiate this class
 * if you need to avoid a dependency to the CoordinateInventor and if you pass encoded coordinates
 * together with any idcode for parsing.
 * We needed to introduce this class to avoid a cyclic dependency between the IDCodeParser and
 * the CoordinateInventor: If encoded atom coords are not given, then the IDcodeParser needs
 * to invent then in order to assign proper up-/down-bonds. The CoordinateInventor needs the
 * IDCodeParser to unpack its default template list.
 */
public class IDCodeParserWithoutCoordinateInvention {
	private StereoMolecule	mMol;
	private byte[]			mDecodingBytes;
	private	int				mIDCodeBitsAvail,mIDCodeTempData,mIDCodeBufferIndex;
	private boolean         mNeglectSpaceDelimitedCoordinates;

	protected boolean ensure2DCoordinates() {
		return false;
		}

	/**
	 * IDCodeParsers allow passing idcode and coordinates as one String with a space
	 * as separator in between. If an idcode is followed by a space and more, and if
	 * the following shall not be interpreted as encoded coordinates, then call this
	 * method after instantiation.
	 */
	public void neglectSpaceDelimitedCoordinates() {
		mNeglectSpaceDelimitedCoordinates = true;
		}

	/**
	 * Creates and returns a molecule from the idcode with its atom and bond arrays being
	 * just as large as needed to hold the molecule. Use this to conserve memory if no
	 * atoms or bonds are added to the molecule afterwards. This version of the method
	 * allows to pass idcode and atom coordinates in one String object.
	 * @param idcode null or idcode, which may contain coordinates separated by a space character
	 * @return
	 */
	public StereoMolecule getCompactMolecule(String idcode) {
		return (idcode == null || idcode.length() == 0) ? null : getCompactMolecule(idcode.getBytes(), null);
		}

	/**
	 * Creates and returns a molecule from the idcode with its atom and bond arrays being
	 * just as large as needed to hold the molecule. Use this to conserve memory if no
	 * atoms or bonds are added to the molecule afterwards.
	 * @param idcode may be null
	 * @return
	 */
	public StereoMolecule getCompactMolecule(byte[] idcode) {
		if (idcode == null || idcode.length == 0)
			return null;

		for (int i=2; i<idcode.length-2; i++)
			if (idcode[i] == ' ')
				return getCompactMolecule(idcode, idcode, 0, i+1);

		return getCompactMolecule(idcode, null);
		}

	/**
	 * Creates and returns a molecule from the idcode with its atom and bond arrays being
	 * just as large as needed to hold the molecule. Use this to conserve memory if no
	 * atoms or bonds are added to the molecule afterwards.
	 * @param idcode may be null
	 * @param coordinates may be null
	 * @return
	 */
	public StereoMolecule getCompactMolecule(String idcode, String coordinates) {
		return (idcode == null) ? null : getCompactMolecule(idcode.getBytes(),
							(coordinates == null) ? null : coordinates.getBytes());
		}

	/**
	 * Creates and returns a molecule from the idcode with its atom and bond arrays being
	 * just as large as needed to hold the molecule. Use this to conserve memory if no
	 * atoms or bonds are added to the molecule afterwards.
	 * @param idcode may be null
	 * @param coordinates may be null
	 * @return
	 */
	public StereoMolecule getCompactMolecule(byte[] idcode, byte[] coordinates) {
		return getCompactMolecule(idcode, coordinates, 0, 0);
		}

	public StereoMolecule getCompactMolecule(byte[] idcode, int idcodeStart) {
		return getCompactMolecule(idcode, null, idcodeStart, -1);
		}

	public StereoMolecule getCompactMolecule(byte[] idcode, byte[] coordinates, int idcodeStart, int coordsStart) {
		if (idcode == null)
			return null;

		decodeBitsStart(idcode, idcodeStart);
		int abits = decodeBits(4);
		int bbits = decodeBits(4);

		if (abits > 8)	// abits is the version number
			abits = bbits;

		int allAtoms = decodeBits(abits);
		int allBonds = decodeBits(bbits);

		StereoMolecule mol = new StereoMolecule(allAtoms, allBonds);
		parse(mol, idcode, coordinates, idcodeStart, coordsStart);
		return mol;
		}

	/**
	 * Parses the idcode and populates the given molecule to represent the passed idcode.
	 * This version of the method allows to pass idcode and atom coordinates in one String object.
	 * @param mol molecule object to be filled with the idcode content
	 * @param idcode null or idcode, which may contain coordinates separated by a space character
	 */
	public void parse(StereoMolecule mol, String idcode) {
		if (idcode == null || idcode.length() == 0) {
			parse(mol, (byte[])null, (byte[])null);
			return;
			}

		int index = idcode.indexOf(' ');
		if (index > 0 && index < idcode.length()-1)
			parse(mol, idcode.substring(0, index).getBytes(), idcode.substring(index+1).getBytes());
		else
			parse(mol, idcode.getBytes(), null);
		}

	/**
	 * Parses the idcode and populates the given molecule to represent the passed idcode.
	 * @param mol molecule object to be filled with the idcode content
	 * @param idcode null or valid idcode optionally concatenates with SPACE and encoded coordinates
	 */
	public void parse(StereoMolecule mol, byte[] idcode) {
		parse(mol, idcode, null);
		}

	/**
	 * Parses the idcode and populates the given molecule to represent the passed idcode.
	 * @param mol molecule object to be filled with the idcode content
	 * @param idcode may be null
	 * @param coordinates may be null
	 */
	public void parse(StereoMolecule mol, String idcode, String coordinates) {
		byte[] idcodeBytes = (idcode == null) ? null : idcode.getBytes();
		byte[] coordinateBytes = (coordinates == null) ? null : coordinates.getBytes();
		parse(mol, idcodeBytes, coordinateBytes);
		}

	/**
	 * Parses the idcode and populates the given molecule to represent the passed idcode.
	 * @param mol molecule object to be filled with the idcode content
	 * @param idcode may be null
	 * @param coordinates may be null
	 */
	public void parse(StereoMolecule mol, byte[] idcode, byte[] coordinates) {
		if (idcode == null || idcode.length == 0) {
			mol.clear();
			return;
			}

		parse(mol, idcode, coordinates, 0, 0);
		}

	/**
	 * Parses the idcode and populates the given molecule to represent the passed idcode.
	 * @param mol molecule object to be filled with the idcode content
	 * @param idcode may be null
	 * @param idcodeStart first byte index of idcode
	 */
	public void parse(StereoMolecule mol, byte[] idcode, int idcodeStart) {
		parse(mol, idcode, null, idcodeStart, -1);
		}

	/**
	 * Parses the idcode and populates the given molecule to represent the passed idcode.
	 * @param mol molecule object to be filled with the idcode content
	 * @param idcode may be null
	 * @param coordinates may be null
	 * @param idcodeStart first byte index of idcode
	 * @param coordsStart first byte indexif coordinates
	 */
	public void parse(StereoMolecule mol, byte[] idcode, byte[] coordinates, int idcodeStart, int coordsStart) {
		mol.clear();

		if (idcode==null || idcodeStart < 0 || idcodeStart >= idcode.length)
			return;

		mMol = mol;
		int version = Canonizer.cIDCodeVersion2;

		if (coordinates != null && (coordsStart < 0 || coordsStart >= coordinates.length))
			coordinates = null;

		decodeBitsStart(idcode, idcodeStart);
		int abits = decodeBits(4);
		int bbits = decodeBits(4);

		if (abits > 8) {	// abits is the version number
			version = abits;
			abits = bbits;
			}

		if (abits == 0) {
			mMol.setFragment(decodeBits(1) == 1);
			return;
			}

		int allAtoms = decodeBits(abits);
		int allBonds = decodeBits(bbits);
		int nitrogens = decodeBits(abits);
		int oxygens = decodeBits(abits);
		int otherAtoms = decodeBits(abits);
		int chargedAtoms = decodeBits(abits);
		for (int atom=0; atom<allAtoms; atom++)
			mMol.addAtom(6);
		for (int i=0; i<nitrogens; i++)
			mMol.setAtomicNo(decodeBits(abits), 7);
		for (int i=0; i<oxygens; i++)
			mMol.setAtomicNo(decodeBits(abits), 8);
		for (int i=0; i<otherAtoms; i++)
			mMol.setAtomicNo(decodeBits(abits),
							 decodeBits(8));
		for (int i=0; i<chargedAtoms; i++)
			mMol.setAtomCharge(decodeBits(abits),
							   decodeBits(4) - 8);

		int closureBonds = 1 + allBonds - allAtoms;
		int dbits = decodeBits(4);
		int base = 0;

		mMol.setAtomX(0, 0.0);
		mMol.setAtomY(0, 0.0);
		mMol.setAtomZ(0, 0.0);

		boolean decodeOldCoordinates = (coordinates != null && coordinates[coordsStart] >= '\'');
		double targetAVBL = 0.0;
		double xOffset = 0.0;
		double yOffset = 0.0;
		double zOffset = 0.0;
		boolean coordsAre3D = false;
		boolean coordsAreAbsolute = false;

		if (decodeOldCoordinates) {	// old coordinate encoding
			if ((coordinates.length > 2*allAtoms-2 && coordinates[2*allAtoms-2] == '\'')
			 || (coordinates.length > 3*allAtoms-3 && coordinates[3*allAtoms-3] == '\'')) {	// old faulty encoding
				coordsAreAbsolute = true;
				coordsAre3D = (coordinates.length == 3*allAtoms-3+9);
				int index = coordsAre3D ? 3*allAtoms-3 : 2*allAtoms-2;
				int avblInt = 86*((int)coordinates[index+1]-40)+(int)coordinates[index+2]-40;
				targetAVBL = Math.pow(10.0, avblInt/2000.0-1.0);
				index += 2;
				int xInt = 86*((int)coordinates[index+1]-40)+(int)coordinates[index+2]-40;
				xOffset = Math.pow(10.0, xInt/1500.0-1.0);
				index += 2;
				int yInt = 86*((int)coordinates[index+1]-40)+(int)coordinates[index+2]-40;
				yOffset = Math.pow(10.0, yInt/1500.0-1.0);
				if (coordsAre3D) {
					index += 2;
					int zInt = 86*((int)coordinates[index+1]-40)+(int)coordinates[index+2]-40;
					zOffset = Math.pow(10.0, zInt/1500.0-1.0);
					}
				}
			else {
				coordsAre3D = (coordinates.length == 3*allAtoms-3);
				}
			}

		// don't use 3D coordinates, if we need 2D
		if (ensure2DCoordinates() && coordsAre3D) {
			coordinates = null;
			decodeOldCoordinates = false;
			}

		for (int i=1; i<allAtoms; i++) {
			int dif = decodeBits(dbits);
			if (dif == 0) {
				if (decodeOldCoordinates) {
					mMol.setAtomX(i, mMol.getAtomX(0) + 8 * (coordinates[i*2-2]-83));
					mMol.setAtomY(i, mMol.getAtomY(0) + 8 * (coordinates[i*2-1]-83));
					if (coordsAre3D)
						mMol.setAtomZ(i, mMol.getAtomZ(0) + 8 * (coordinates[2*allAtoms-3+i]-83));
					}

				closureBonds++;
				continue;
				}

			base += dif - 1;

			if (decodeOldCoordinates) {
				mMol.setAtomX(i, mMol.getAtomX(base) + coordinates[i*2-2] - 83);
				mMol.setAtomY(i, mMol.getAtomY(base) + coordinates[i*2-1] - 83);
				if (coordsAre3D)
					mMol.setAtomZ(i, mMol.getAtomZ(base) + (coordinates[2*allAtoms-3+i]-83));
				}
			mMol.addBond(base, i, Molecule.cBondTypeSingle);
			}

		for (int i=0; i<closureBonds; i++)
			mMol.addBond(decodeBits(abits),
						 decodeBits(abits), Molecule.cBondTypeSingle);

		boolean[] isDelocalizedBond = new boolean[allBonds];

		for (int bond=0; bond<allBonds; bond++) {
			int bondOrder = decodeBits(2);
			switch (bondOrder) {
			case 0:
				isDelocalizedBond[bond] = true;
				break;
			case 2:
				mMol.setBondType(bond, Molecule.cBondTypeDouble);
				break;
			case 3:
				mMol.setBondType(bond, Molecule.cBondTypeTriple);
				break;
				}
			}

		int THCount = decodeBits(abits);
		for (int i=0; i<THCount; i++) {
			int atom = decodeBits(abits);
			if (version == Canonizer.cIDCodeVersion2) {
				int parity = decodeBits(2);
				if (parity == 3) {
					// this was the old discontinued Molecule.cAtomParityMix
					// version2 idcodes had never more than one center with parityMix
					mMol.setAtomESR(atom, Molecule.cESRTypeAnd, 0);
					mMol.setAtomParity(atom, Molecule.cAtomParity1, false);
					}
				else {
					mMol.setAtomParity(atom, parity, false);
					}
				}
			else {
				int parity = decodeBits(3);
				switch (parity) {
				case Canonizer.cParity1And:
					mMol.setAtomParity(atom, Molecule.cAtomParity1, false);
					mMol.setAtomESR(atom, Molecule.cESRTypeAnd, decodeBits(3));
					break;
				case Canonizer.cParity2And:
					mMol.setAtomParity(atom, Molecule.cAtomParity2, false);
					mMol.setAtomESR(atom, Molecule.cESRTypeAnd, decodeBits(3));
					break;
				case Canonizer.cParity1Or:
					mMol.setAtomParity(atom, Molecule.cAtomParity1, false);
					mMol.setAtomESR(atom, Molecule.cESRTypeOr, decodeBits(3));
					break;
				case Canonizer.cParity2Or:
					mMol.setAtomParity(atom, Molecule.cAtomParity2, false);
					mMol.setAtomESR(atom, Molecule.cESRTypeOr, decodeBits(3));
					break;
				default:
					mMol.setAtomParity(atom, parity, false);
					}
				}
			}

		if (version == Canonizer.cIDCodeVersion2)
			if ((decodeBits(1) == 0))   // translate chiral flag
				mMol.setToRacemate();

		int EZCount = decodeBits(bbits);
		for (int i=0; i<EZCount; i++) {
			int bond = decodeBits(bbits);
			if (mMol.getBondType(bond) == Molecule.cBondTypeSingle) {	// BINAP type of axial chirality
				int parity = decodeBits(3);
				switch (parity) {
				case Canonizer.cParity1And:
					mMol.setBondParity(bond, Molecule.cBondParityEor1, false);
					mMol.setBondESR(bond, Molecule.cESRTypeAnd, decodeBits(3));
					break;
				case Canonizer.cParity2And:
					mMol.setBondParity(bond, Molecule.cBondParityZor2, false);
					mMol.setBondESR(bond, Molecule.cESRTypeAnd, decodeBits(3));
					break;
				case Canonizer.cParity1Or:
					mMol.setBondParity(bond, Molecule.cBondParityEor1, false);
					mMol.setBondESR(bond, Molecule.cESRTypeOr, decodeBits(3));
					break;
				case Canonizer.cParity2Or:
					mMol.setBondParity(bond, Molecule.cBondParityZor2, false);
					mMol.setBondESR(bond, Molecule.cESRTypeOr, decodeBits(3));
					break;
				default:
					mMol.setBondParity(bond, parity, false);
					}
				}
			else {
				mMol.setBondParity(bond, decodeBits(2), false);	// double bond
				}
			}

		mMol.setFragment(decodeBits(1) == 1);

		int[] aromaticSPBond = null;

		int offset = 0;
		while (decodeBits(1) == 1) {
			int dataType = offset + decodeBits(4);
			switch (dataType) {
			case 0:	//	datatype 'AtomQFNoMoreNeighbours'
				int no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNoMoreNeighbours, true);
					}
				break;
			case 1:	//	datatype 'isotop'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					int mass = decodeBits(8);
					mMol.setAtomMass(atom, mass);
					}
				break;
			case 2:	//	datatype 'bond defined to be delocalized'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++) {
					int bond = decodeBits(bbits);
// This used to be Molecule.cBondTypeDelocalized, which is redundant to the bond order encoding
// Then it was wrongly fixed to cBondQFDelocalized, which is part of cBondQFBondTypes and encoded as type 10
// We can take it out entirely without sacrifycing idcode compatibility
//					mMol.setBondQueryFeature(bond, Molecule.cBondQFDelocalized, true);
//System.out.println("wrong outdated 'delocalized bond'; idcode:"+new String(mDecodingBytes));
					}
				break;
			case 3:	//	datatype 'AtomQFMoreNeighbours'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFMoreNeighbours, true);
					}
				break;
			case 4:	//	datatype 'AtomQFRingState'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long ringState = (long)decodeBits(Molecule.cAtomQFRingStateBits) << Molecule.cAtomQFRingStateShift;
					mMol.setAtomQueryFeature(atom, ringState, true);
					}
				break;
			case 5:	//	datatype 'AtomQFAromState'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long aromState = (long)decodeBits(Molecule.cAtomQFAromStateBits) << Molecule.cAtomQFAromStateShift;
					mMol.setAtomQueryFeature(atom, aromState, true);
					}
				break;
			case 6:	//	datatype 'AtomQFAny'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFAny, true);
					}
				break;
			case 7:	//	datatype 'AtomQFHydrogen'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long hydrogen = (long)decodeBits(Molecule.cAtomQFHydrogenBits) << Molecule.cAtomQFHydrogenShift;
					mMol.setAtomQueryFeature(atom, hydrogen, true);
					}
				break;
			case 8:	//	datatype 'AtomList'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					int atoms = decodeBits(4);
					int[] atomList = new int[atoms];
					for (int j=0; j<atoms; j++) {
						int atomicNo = decodeBits(8);
						atomList[j] = atomicNo;
						}
					mMol.setAtomList(atom, atomList);
					}
				break;
			case 9:	//	datatype 'BondQFRingState'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++) {
					int bond = decodeBits(bbits);
					int ringState = decodeBits(Molecule.cBondQFRingStateBits) << Molecule.cBondQFRingStateShift;
					mMol.setBondQueryFeature(bond, ringState, true);
					}
				break;
			case 10://	datatype 'BondQFBondTypes'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++) {
					int bond = decodeBits(bbits);
					int bondTypes = decodeBits(Molecule.cBondQFBondTypesBits) << Molecule.cBondQFBondTypesShift;
					mMol.setBondQueryFeature(bond, bondTypes, true);
					}
				break;
			case 11:	//	datatype 'AtomQFMatchStereo'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFMatchStereo, true);
					}
				break;
			case 12:	//  datatype 'bond defined to be a bridge from n1 to n2 atoms'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++) {
					int bond = decodeBits(bbits);
					int bridgeData = decodeBits(Molecule.cBondQFBridgeBits) << Molecule.cBondQFBridgeShift;
					mMol.setBondQueryFeature(bond, bridgeData, true);
					}
				break;
			case 13: //  datatype 'AtomQFPiElectrons'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long piElectrons = (long)decodeBits(Molecule.cAtomQFPiElectronBits) << Molecule.cAtomQFPiElectronShift;
					mMol.setAtomQueryFeature(atom, piElectrons, true);
					}
				break;
			case 14: //  datatype 'AtomQFNeighbours'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long neighbours = (long)decodeBits(Molecule.cAtomQFNeighbourBits) << Molecule.cAtomQFNeighbourShift;
					mMol.setAtomQueryFeature(atom, neighbours, true);
					}
				break;
			case 15: //  datatype 'start next feature set'
			case 31:
				offset += 16;
				break;
			case 16: //  datatype 'AtomQFSmallRingSize'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long ringSize = (long)decodeBits(Molecule.cAtomQFSmallRingSizeBits) << Molecule.cAtomQFSmallRingSizeShift;
					mMol.setAtomQueryFeature(atom, ringSize, true);
					}
				break;
			case 17: //  datatype 'AtomAbnormalValence'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					mMol.setAtomAbnormalValence(atom, decodeBits(4));
					}
				break;
			case 18: //  datatype 'AtomCustomLabel'
				no = decodeBits(abits);
				int lbits = decodeBits(4);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					int count = decodeBits(lbits);
					byte[] label = new byte[count];
					for (int j=0; j<count; j++)
						label[j] = (byte)decodeBits(7);
					mMol.setAtomCustomLabel(atom, new String(label));
					}
				break;
			case 19: //  datatype 'AtomQFCharge'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long charge = (long)decodeBits(Molecule.cAtomQFChargeBits) << Molecule.cAtomQFChargeShift;
					mMol.setAtomQueryFeature(atom, charge, true);
					}
				break;
			case 20: //  datatype 'BondQFRingSize'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++) {
					int bond = decodeBits(bbits);
					int ringSize = decodeBits(Molecule.cBondQFRingSizeBits) << Molecule.cBondQFRingSizeShift;
					mMol.setBondQueryFeature(bond, ringSize, true);
					}
				break;
			case 21: //  datatype 'AtomRadicalState'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					mMol.setAtomRadical(atom, decodeBits(2) << Molecule.cAtomRadicalStateShift);
					}
				break;
			case 22:	//	datatype 'flat nitrogen'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFFlatNitrogen, true);
					}
				break;
			case 23:	//	datatype 'BondQFMatchStereo'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++) {
					int bond = decodeBits(bbits);
					mMol.setBondQueryFeature(bond, Molecule.cBondQFMatchStereo, true);
					}
				break;
			case 24:	//	datatype 'cBondQFAromState'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++) {
					int bond = decodeBits(bbits);
					int aromState = decodeBits(Molecule.cBondQFAromStateBits) << Molecule.cBondQFAromStateShift;
					mMol.setBondQueryFeature(bond, aromState, true);
					}
				break;
			case 25:	//	datatype 'atom selection'
				for (int i=0; i<allAtoms; i++)
					if (decodeBits(1) == 1)
						mMol.setAtomSelection(i, true);
				break;
			case 26:	//	datatype 'delocalized high order bond'
				no = decodeBits(bbits);
				aromaticSPBond = new int[no];
				for (int i=0; i<no; i++)
					aromaticSPBond[i] = decodeBits(bbits);
				break;
			case 27:	//	datatype 'part of an exclude group'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFExcludeGroup, true);
					}
				break;
			case 28: //  datatype 'coordinate bond'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++)
					mMol.setBondType(decodeBits(bbits), Molecule.cBondTypeMetalLigand);
				break;
			case 29:	//	datatype 'reaction parity hint'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long hint = (long)decodeBits(Molecule.cAtomQFRxnParityBits) << Molecule.cAtomQFRxnParityShift;
					mMol.setAtomQueryFeature(atom, hint, true);
					}
				break;
			case 30: //  datatype 'AtomQFNewRingSize'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long ringSize = (long)decodeBits(Molecule.cAtomQFNewRingSizeBits) << Molecule.cAtomQFNewRingSizeShift;
					mMol.setAtomQueryFeature(atom, ringSize, true);
					}
				break;
			case 32: //  datatype 'AtomQFNewRingSize'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long stereoState = (long)decodeBits(Molecule.cAtomQFStereoStateBits) << Molecule.cAtomQFStereoStateShift;
					mMol.setAtomQueryFeature(atom, stereoState, true);
					}
				break;
			case 33: //  datatype 'AtomQFENeighbours'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					long eNeighbours = (long)decodeBits(Molecule.cAtomQFENeighbourBits) << Molecule.cAtomQFENeighbourShift;
					mMol.setAtomQueryFeature(atom, eNeighbours, true);
					}
				break;
			case 34:	//	datatype 'AtomQFHetereoAromatic'
				no = decodeBits(abits);
				for (int i=0; i<no; i++) {
					int atom = decodeBits(abits);
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFHeteroAromatic, true);
					}
				break;
			case 35:	//	datatype 'BondQFMatchFormalOrder'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++) {
					int bond = decodeBits(bbits);
					mMol.setBondQueryFeature(bond, Molecule.cBondQFMatchFormalOrder, true);
					}
				break;
			case 36:	//	datatype 'cBondQFRareBondType'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++) {
					int bond = decodeBits(bbits);
					int bondType = decodeBits(Molecule.cBondQFRareBondTypesBits) << Molecule.cBondQFRareBondTypesShift;
					mMol.setBondQueryFeature(bond, bondType, true);
					}
				break;
			case 37:	//	datatype 'rare order bond'
				no = decodeBits(bbits);
				for (int i=0; i<no; i++) {
					int bond = decodeBits(bbits);
					int bondType = decodeBits(1) == 0 ? Molecule.cBondTypeQuadruple : Molecule.cBondTypeQuintuple;
					mMol.setBondType(bond, bondType);
					}
				break;
				}
			}

		new AromaticityResolver(mMol).locateDelocalizedDoubleBonds(isDelocalizedBond);

		if (aromaticSPBond != null)
			for (int bond:aromaticSPBond)
				mMol.setBondType(bond, mMol.getBondType(bond) == Molecule.cBondTypeDouble ?
						Molecule.cBondTypeTriple : Molecule.cBondTypeDouble);

		if (coordinates == null
		 && !mNeglectSpaceDelimitedCoordinates
		 && idcode.length > mIDCodeBufferIndex+1
		 && (idcode[mIDCodeBufferIndex+1] == ' ' || idcode[mIDCodeBufferIndex+1] == '\t')) {
			coordinates = idcode;
			coordsStart = mIDCodeBufferIndex+2;
			}

		if (coordinates != null) {
			try {
				if (coordinates[coordsStart] == '!' || coordinates[coordsStart] == '#') {    // new coordinate format
					decodeBitsStart(coordinates, coordsStart + 1);
					coordsAre3D = (decodeBits(1) == 1);
					coordsAreAbsolute = (decodeBits(1) == 1);
					int resolutionBits = 2 * decodeBits(4);
					int binCount = (1 << resolutionBits);

					double factor;
					int from = 0;
					int bond = 0;
					for (int atom = 1; atom < allAtoms; atom++) {
						if (bond < allBonds && mMol.getBondAtom(1, bond) == atom) {
							from = mMol.getBondAtom(0, bond++);
							factor = 1.0;
							}
						else {
							from = 0;
							factor = 8.0;
							}
						mMol.setAtomX(atom, mMol.getAtomX(from) + factor * (decodeBits(resolutionBits) - binCount / 2));
						mMol.setAtomY(atom, mMol.getAtomY(from) + factor * (decodeBits(resolutionBits) - binCount / 2));
						if (coordsAre3D)
							mMol.setAtomZ(atom, mMol.getAtomZ(from) + factor * (decodeBits(resolutionBits) - binCount / 2));
						}

					if (coordinates[coordsStart] == '#') {    // we have 3D-coordinates that include implicit hydrogen coordinates
						int hydrogenCount = 0;

						// we need to cache hCount, because otherwise getImplicitHydrogens() would create helper arrays with every call
						int[] hCount = new int[allAtoms];
						for (int atom = 0; atom < allAtoms; atom++)
							hydrogenCount += (hCount[atom] = mMol.getImplicitHydrogens(atom));

						for (int atom = 0; atom < allAtoms; atom++) {
							for (int i = 0; i < hCount[atom]; i++) {
								int hydrogen = mMol.addAtom(1);
								mMol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);

								mMol.setAtomX(hydrogen, mMol.getAtomX(atom) + (decodeBits(resolutionBits) - binCount / 2));
								mMol.setAtomY(hydrogen, mMol.getAtomY(atom) + (decodeBits(resolutionBits) - binCount / 2));
								if (coordsAre3D)
									mMol.setAtomZ(hydrogen, mMol.getAtomZ(atom) + (decodeBits(resolutionBits) - binCount / 2));
								}
							}

						allAtoms += hydrogenCount;
						allBonds += hydrogenCount;
						}

					double avblDefault = coordsAre3D ? 1.5 : Molecule.getDefaultAverageBondLength();
					double avbl = mMol.getAverageBondLength(allAtoms, allBonds, avblDefault);

					if (coordsAreAbsolute) {
						targetAVBL = decodeAVBL(decodeBits(resolutionBits), binCount);
						xOffset = targetAVBL * decodeShift(decodeBits(resolutionBits), binCount);
						yOffset = targetAVBL * decodeShift(decodeBits(resolutionBits), binCount);
						if (coordsAre3D)
							zOffset = targetAVBL * decodeShift(decodeBits(resolutionBits), binCount);

						factor = targetAVBL / avbl;
						for (int atom = 0; atom < allAtoms; atom++) {
							mMol.setAtomX(atom, xOffset + factor * mMol.getAtomX(atom));
							mMol.setAtomY(atom, yOffset + factor * mMol.getAtomY(atom));
							if (coordsAre3D)
								mMol.setAtomZ(atom, zOffset + factor * mMol.getAtomZ(atom));
							}
						}
					else {    // with new format 2D and 3D coordinates are scaled to average bond lengths of 1.5 Angstrom
						targetAVBL = 1.5;
						factor = targetAVBL / avbl;
						for (int atom = 0; atom < allAtoms; atom++) {
							mMol.setAtomX(atom, factor * mMol.getAtomX(atom));
							mMol.setAtomY(atom, factor * mMol.getAtomY(atom));
							if (coordsAre3D)
								mMol.setAtomZ(atom, factor * mMol.getAtomZ(atom));
							}
						}
					}
				else {    // old coordinate format
					if (coordsAre3D && !coordsAreAbsolute && targetAVBL == 0.0) // if no scaling factor is given, then scale to mean bond length = 1.5
						targetAVBL = 1.5;

					if (targetAVBL != 0.0 && mMol.getAllBonds() != 0) {
						double avbl = 0.0;
						for (int bond = 0; bond < mMol.getAllBonds(); bond++) {
							double dx = mMol.getAtomX(mMol.getBondAtom(0, bond)) - mMol.getAtomX(mMol.getBondAtom(1, bond));
							double dy = mMol.getAtomY(mMol.getBondAtom(0, bond)) - mMol.getAtomY(mMol.getBondAtom(1, bond));
							double dz = coordsAre3D ? mMol.getAtomZ(mMol.getBondAtom(0, bond)) - mMol.getAtomZ(mMol.getBondAtom(1, bond)) : 0.0f;
							avbl += Math.sqrt(dx * dx + dy * dy + dz * dz);
							}
						avbl /= mMol.getAllBonds();
						double f = targetAVBL / avbl;
						for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
							mMol.setAtomX(atom, mMol.getAtomX(atom) * f + xOffset);
							mMol.setAtomY(atom, mMol.getAtomY(atom) * f + yOffset);
							if (coordsAre3D)
								mMol.setAtomZ(atom, mMol.getAtomZ(atom) * f + zOffset);
							}
						}
					}
				}
			catch (Exception e) {
				e.printStackTrace();
				System.err.println("Faulty id-coordinates:"+e.toString()+" "+new String(idcode)+" "+new String(coordinates));
				coordinates = null;
				coordsAre3D = false;
				}
			}

		boolean coords2DAvailable = (coordinates != null && !coordsAre3D);

		// If we have or create 2D-coordinates, then we need to set all double bonds to a cross bond, which
		// - have distinguishable substituents on both ends, i.e. is a stereo double bond
		// - are not in a small ring
		// Here we don't know, whether a double bond without E/Z parity is a stereo bond with unknown
		// configuration or not a stereo bond. Therefore, we need to set a flag, that causes the Canonizer
		// during the next stereo recognition with atom coordinates to assign an unknown configuration rather
		// than E or Z based on created or given coordinates.
		// In a next step these double bonds are converted into cross bonds by
		if (coords2DAvailable || ensure2DCoordinates()) {
			mMol.ensureHelperArrays(Molecule.cHelperRings);
			for (int bond=0; bond<mMol.getBonds(); bond++)
				if (mMol.getBondOrder(bond) == 2
				 && !mMol.isSmallRingBond(bond)
				 && mMol.getBondParity(bond) == Molecule.cBondParityNone)
					mMol.setBondParityUnknownOrNone(bond);
			}

		mMol.setParitiesValid(0);

		if (!coords2DAvailable && ensure2DCoordinates()) {
			try {
				inventCoordinates(mMol);
				coords2DAvailable = true;
				}
			catch (Exception e) {
				e.printStackTrace();
				System.err.println("2D-coordinate creation failed:"+e.toString()+" "+new String(idcode));
				}
			}

		if (coords2DAvailable) {
			mMol.setStereoBondsFromParity();
			mMol.setUnknownParitiesToExplicitlyUnknown();
			}
		else if (!coordsAre3D) {
			mMol.setParitiesValid(0);
			}
		}

	protected void inventCoordinates(StereoMolecule mol) throws Exception {
		throw new Exception("Unexpected request to invent coordinates. Check source code logic!");
		}

	/**
	 * This method parses an id-coordinate string (new format only) and writes the coordinates into a Coordinates array.
	 * If the id-coordinates contain implicit hydrogen coordinates, then this method does not(!!!) add these hydrogen atoms
	 * to the Molecule. Thus, for 3D-coordinates with implicit hydrogen coordinates, you need to make sure that all
	 * of the Molecule's hydrogen atoms are explicit and that the Coordinates array's size covers all hydrogens atoms.
	 * For instance, if parsing idcodes and coordinates of a conformer set, you may parse the first conformer with one
	 * of the getCompactMolecule() or parse() methods.
	 * This adds all implicit hydrogens as explicit ones to the Molecule and conformer object. All subsequent conformers
	 * may be generated by instantiating a new Conformer from the molecule and using this method to populate its coordinates.
	 * @param encodedCoords id-coordinates string in new format
	 * @param coordsStart offset
	 * @param mol Molecule with all hydrogens being explicit
	 * @param coords Coordinates array to be populated by coordinates taken from encoded string
	 * @throws Exception
	 */
	public void parseCoordinates(byte[] encodedCoords, int coordsStart, StereoMolecule mol, Coordinates[] coords) throws Exception {
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		int atomCount = mol.getAtoms();
		int bondCount = mol.getBonds();

		decodeBitsStart(encodedCoords, coordsStart + 1);
		boolean coordsAre3D = (decodeBits(1) == 1);
		boolean coordsAreAbsolute = (decodeBits(1) == 1);
		int resolutionBits = 2 * decodeBits(4);
		int binCount = (1 << resolutionBits);

		double factor;
		int from = 0;
		int bond = 0;
		for (int atom = 1; atom < atomCount; atom++) {
			if (bond < bondCount && mol.getBondAtom(1, bond) == atom) {
				from = mol.getBondAtom(0, bond++);
				factor = 1.0;
				}
			else {
				from = 0;
				factor = 8.0;
				}
			coords[atom].x = coords[from].x + factor * (decodeBits(resolutionBits) - binCount / 2);
			coords[atom].y = coords[from].y + factor * (decodeBits(resolutionBits) - binCount / 2);
			if (coordsAre3D)
				coords[atom].z = coords[from].z + factor * (decodeBits(resolutionBits) - binCount / 2);
			}

		double avbl = coordsAre3D ? 1.5 : Molecule.getDefaultAverageBondLength();
		if (bondCount != 0)
			for (int b=0; b<bondCount; b++)
				avbl += coords[mol.getBondAtom(0, b)].distance(coords[mol.getBondAtom(1, b)]);
		avbl /= bondCount;

		if (encodedCoords[coordsStart] == '#') {    // we have 3D-coordinates that include implicit hydrogen coordinates
			int hydrogenCount = 0;

			int hydrogen = atomCount;
			for (int atom = 0; atom < atomCount; atom++) {
				int hCount = mol.getAllConnAtoms(atom) - mol.getConnAtoms(atom);
				for (int i = 0; i < hCount; i++) {
					coords[hydrogen].x = coords[atom].x + (decodeBits(resolutionBits) - binCount / 2);
					coords[hydrogen].y = coords[atom].y + (decodeBits(resolutionBits) - binCount / 2);
					if (coordsAre3D)
						coords[hydrogen].z = coords[atom].z + (decodeBits(resolutionBits) - binCount / 2);

					hydrogen++;
					}
				hydrogenCount += hCount;
				}

			atomCount += hydrogenCount;
			bondCount += hydrogenCount;
			}

		if (coordsAreAbsolute) {
			double targetAVBL = decodeAVBL(decodeBits(resolutionBits), binCount);
			double xOffset = targetAVBL * decodeShift(decodeBits(resolutionBits), binCount);
			double yOffset = targetAVBL * decodeShift(decodeBits(resolutionBits), binCount);
			double zOffset = 0;
			if (coordsAre3D)
				zOffset = targetAVBL * decodeShift(decodeBits(resolutionBits), binCount);

			factor = targetAVBL / avbl;
			for (int atom = 0; atom < atomCount; atom++) {
				coords[atom].x = xOffset + factor * coords[atom].x;
				coords[atom].y = yOffset + factor * coords[atom].y;
				if (coordsAre3D)
					coords[atom].z = zOffset + factor * coords[atom].z;
				}
			}
		else {    // with new format 2D and 3D coordinates are scaled to average bond lengths of 1.5 Angstrom
			double targetAVBL = 1.5;
			factor = targetAVBL / avbl;
			for (int atom = 0; atom < atomCount; atom++) {
				coords[atom].x = factor * coords[atom].x;
				coords[atom].y = factor * coords[atom].y;
				if (coordsAre3D)
					coords[atom].z = factor * coords[atom].z;
				}
			}
		}

	public void parseMapping(byte[] mapping) {
		parseMapping(mapping, 0);
		}

	public void parseMapping(byte[] mapping, int mappingStart) {
		if (mapping == null || mapping.length <= mappingStart || mapping[mappingStart] < 0x40)
			return;

		decodeBitsStart(mapping, mappingStart);
		int nbits = decodeBits(4);
		boolean autoMappingFound = (decodeBits(1) == 1);
		boolean manualMappingFound = (decodeBits(1) == 1);
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			int mapNo = decodeBits(nbits);
			boolean autoMapped = autoMappingFound;
			if (autoMappingFound && manualMappingFound)
				autoMapped = (decodeBits(1) == 1);
			mMol.setAtomMapNo(atom, mapNo, autoMapped);
			}
		}

	public boolean coordinatesAre3D(String idcode, String coordinates) {
		return coordinates != null && coordinatesAre3D(idcode.getBytes(), coordinates.getBytes());
		}

	public boolean coordinatesAre3D(byte[] idcode, byte[] coordinates) {
		return coordinatesAre3D(idcode, coordinates, 0, 0);
		}

	public boolean coordinatesAre3D(byte[] idcode, byte[] coordinates, int idcodeStart, int coordsStart) {
		if (coordinates == null || coordinates.length <= coordsStart)
			return false;

		if (coordinates[coordsStart] == '!' || coordinates[coordsStart] == '#') {
			// current version starts with '!' (ASC 33) or '#' (ASC 35) (includs implicit hydrogen coordinates)
			// further versions may start with ASC 36 to 38
			decodeBitsStart(coordinates, coordsStart+1);
			return (decodeBits(1) == 1);
			}
		else {	// old format uses ACSII 39 and higher
			int allAtoms = getAtomCount(idcode, idcodeStart);
			return (allAtoms != 0
				 && coordinates.length >= coordsStart+3*allAtoms-3
				 && coordinates[coordsStart+2*allAtoms-2] != '\'');
			}
		}

	public boolean coordinatesAreAbsolute(String coordinates) {
		return coordinates != null && coordinatesAreAbsolute(coordinates.getBytes());
		}

	public boolean coordinatesAreAbsolute(byte[] coordinates) {
		return coordinatesAreAbsolute(coordinates, 0);
		}

	public boolean coordinatesAreAbsolute(byte[] coordinates, int coordStart) {
		if (coordinates == null || coordinates.length <= coordStart)
			return false;

		if (coordinates[coordStart] >= '\'') {	// old format uses ACSII 39 and higher
			for (int i=coordStart; i<coordinates.length; i++)
				if (coordinates[i] == '\'' || coordinates[i] == '&')
					return true;
			}
		else if (coordinates[coordStart] == '!' || coordinates[coordStart] == '#') {
			// current version starts with '!' (ASC 33) or '#' (ASC 35) (includs implicit hydrogen coordinates)
			// further versions may start with ASC 36 to 38
			decodeBitsStart(coordinates, coordStart+1);
			decodeBits(1);	// skip 3D information
			return (decodeBits(1) == 1);
			}

		return false;
		}

	public int getIDCodeVersion(String idcode) {
		if (idcode == null || idcode.length() == 0)
			return -1;

		return getIDCodeVersion(idcode.getBytes());
		}

	public int getIDCodeVersion(byte[] idcode) {
		int version = Canonizer.cIDCodeVersion2;

		decodeBitsStart(idcode, 0);
		int abits = decodeBits(4);
		if (abits > 8)	// abits is the version number
			version = abits;

		return version;
		}

	public int getAtomCount(String idcode) {
		if (idcode == null || idcode.length() == 0)
			return 0;

		return getAtomCount(idcode.getBytes(), 0);
		}

	public int getAtomCount(byte[] idcode, int offset) {
		if (idcode == null || idcode.length <= offset)
			return 0;

		decodeBitsStart(idcode, offset);
		int abits = decodeBits(4);
		int bbits = decodeBits(4);

		if (abits > 8)	// abits is the version number
			abits = bbits;

		if (abits == 0)
			return 0;

		return decodeBits(abits);
		}

	/**
	 * Determines atom and bond counts of the given idcode
	 * @param idcode
	 * @param count null or int[2], which is filled and returned
	 * @return int[] with atom and bond count as first and second values
	 */
	public int[] getAtomAndBondCounts(String idcode, int[] count) {
		if (idcode == null || idcode.length() == 0)
			return null;

		return getAtomAndBondCounts(idcode.getBytes(), 0, count);
		}

	/**
	 * Determines atom and bond counts of the given idcode
	 * @param idcode
	 * @param offset
	 * @param count null or int[2], which is filled and returned
     * @return int[] with atom and bond count as first and second values
     */
	public int[] getAtomAndBondCounts(byte[] idcode, int offset, int[] count) {
		if (idcode == null || idcode.length == 0)
			return null;

		decodeBitsStart(idcode, 0);
		int abits = decodeBits(4);
		int bbits = decodeBits(4);

		if (abits > 8)   // abits is the version number
			abits = bbits;

		if (count == null)
			count = new int[2];

		if (abits == 0) {
			count[0] = 0;
			count[1] = 0;
			}
		else {
			count[0] = decodeBits(abits);
			count[1] = decodeBits(bbits);
			}

		return count;
		}

	private void decodeBitsStart(byte[] bytes, int offset) {
		mIDCodeBitsAvail = 6;
		mIDCodeBufferIndex = offset;
		mDecodingBytes = bytes;
		mIDCodeTempData = (bytes[mIDCodeBufferIndex] & 0x3F) << 11;
		}

	private int decodeBits(int bits) {
		int allBits = bits;

		int data = 0;
		while (bits != 0) {
			if (mIDCodeBitsAvail == 0) {
				mIDCodeTempData = (mDecodingBytes[++mIDCodeBufferIndex] & 0x3F) << 11;
				mIDCodeBitsAvail = 6;
				}
			data |= ((0x00010000 & mIDCodeTempData) >> (16 - allBits + bits));
			mIDCodeTempData <<= 1;
			bits--;
			mIDCodeBitsAvail--;
			}
		return data;
		}

	private double decodeAVBL(int value, int binCount) {
		return Math.pow(10, Math.log10(200/0.1) * value / (binCount - 1) - 1);
		}

	private double decodeShift(int value, int binCount) {
		int halfBinCount = binCount / 2;
		boolean isNegative = (value >= halfBinCount);
		if (isNegative)
			value -= halfBinCount;
		double steepness = binCount/32;
		double doubleValue = steepness * value / (halfBinCount - value);
		return isNegative ? -doubleValue : doubleValue;
		}

	public void printContent(byte[] idcode, byte[] coordinates) {
		try {
			int version = Canonizer.cIDCodeVersion2;

			if (idcode == null || idcode.length == 0)
				return;

			if (coordinates != null && coordinates.length == 0)
				coordinates = null;

			System.out.println("idcode: " + new String(idcode));
			if (coordinates != null)
				System.out.println("coords: " + new String(coordinates));

			decodeBitsStart(idcode, 0);
			int abits = decodeBits(4);
			int bbits = decodeBits(4);

			if (abits > 8) {    // abits is the version number
				version = abits;
				abits = bbits;
			}

			System.out.println("version:" + version);

			int allAtoms = decodeBits(abits);
			if (allAtoms == 0)
				return;

			int allBonds = decodeBits(bbits);
			int nitrogens = decodeBits(abits);
			int oxygens = decodeBits(abits);
			int otherAtoms = decodeBits(abits);
			int chargedAtoms = decodeBits(abits);

			System.out.println("allAtoms:" + allAtoms + " allBonds:" + allBonds);
			if (nitrogens != 0) {
				System.out.print("nitrogens:");
				for (int i = 0; i < nitrogens; i++)
					System.out.print(" " + decodeBits(abits));
				System.out.println();
			}
			if (oxygens != 0) {
				System.out.print("oxygens:");
				for (int i = 0; i < oxygens; i++)
					System.out.print(" " + decodeBits(abits));
				System.out.println();
			}
			if (otherAtoms != 0) {
				System.out.print("otherAtoms:");
				for (int i = 0; i < otherAtoms; i++)
					System.out.print(" " + decodeBits(abits) + ":" + decodeBits(8));
				System.out.println();
			}
			if (chargedAtoms != 0) {
				System.out.print("chargedAtoms:");
				for (int i = 0; i < chargedAtoms; i++)
					System.out.print(" " + decodeBits(abits) + ":" + (decodeBits(4) - 8));
				System.out.println();
			}

			int closureBonds = 1 + allBonds - allAtoms;
			int dbits = decodeBits(4);
			int base = 0;

			int[][] bondAtom = new int[2][allBonds];
			int bondCount = 0;
			for (int i = 1; i < allAtoms; i++) {
				int dif = decodeBits(dbits);
				if (dif == 0) {
					closureBonds++;
					continue;
				}
				base += dif - 1;
				bondAtom[0][bondCount] = base;
				bondAtom[1][bondCount++] = i;
			}

			for (int i = 0; i < closureBonds; i++) {
				bondAtom[0][bondCount] = decodeBits(abits);
				bondAtom[1][bondCount++] = decodeBits(abits);
			}

			int[] bondOrder = new int[allBonds];
			System.out.print("bonds:");
			for (int bond = 0; bond < allBonds; bond++) {
				System.out.print(" " + bondAtom[0][bond]);
				bondOrder[bond] = decodeBits(2);
				System.out.print(bondOrder[bond] == 0 ? "." : bondOrder[bond] == 1 ? "-" : bondOrder[bond] == 2 ? "=" : "#");
				System.out.print("" + bondAtom[1][bond]);
			}
			System.out.println();

			int THCount = decodeBits(abits);
			if (THCount != 0) {
				System.out.print("parities:");
				for (int i = 0; i < THCount; i++) {
					int atom = decodeBits(abits);
					if (version == Canonizer.cIDCodeVersion2) {
						int parity = decodeBits(2);
						if (parity == 3) {
							// this was the old discontinued Molecule.cAtomParityMix
							// version2 idcodes had never more than one center with parityMix
							System.out.print(" " + atom + ":1&0");
						} else {
							System.out.print(" " + atom + ":" + parity);
						}
					} else {
						int parity = decodeBits(3);
						switch (parity) {
							case Canonizer.cParity1And:
								System.out.print(" " + atom + ":1&" + decodeBits(3));
								break;
							case Canonizer.cParity2And:
								System.out.print(" " + atom + ":2&" + decodeBits(3));
								break;
							case Canonizer.cParity1Or:
								System.out.print(" " + atom + ":1|" + decodeBits(3));
								break;
							case Canonizer.cParity2Or:
								System.out.print(" " + atom + ":2|" + decodeBits(3));
								break;
							default:
								System.out.print(" " + atom + ":" + parity);
						}
					}
				}
				System.out.println();
			}

			if (version == Canonizer.cIDCodeVersion2)
				if ((decodeBits(1) == 0))   // translate chiral flag
					System.out.println("isRacemate");

			int EZCount = decodeBits(bbits);
			if (EZCount != 0) {
				System.out.print("EZ:");
				for (int i = 0; i < EZCount; i++) {
					int bond = decodeBits(bbits);
					if (bondOrder[bond] == 1) {    // BINAP type of axial chirality
						int parity = decodeBits(3);
						switch (parity) {
							case Canonizer.cParity1And:
								System.out.print(" " + bond + ":1&" + decodeBits(3));
								break;
							case Canonizer.cParity2And:
								System.out.print(" " + bond + ":2&" + decodeBits(3));
								break;
							case Canonizer.cParity1Or:
								System.out.print(" " + bond + ":1|" + decodeBits(3));
								break;
							case Canonizer.cParity2Or:
								System.out.print(" " + bond + ":2|" + decodeBits(3));
								break;
							default:
								System.out.print(" " + bond + ":" + parity);
						}
					} else
						System.out.print(" " + bond + ":" + decodeBits(2));
				}
				System.out.println();
			}

			if (decodeBits(1) == 1)
				System.out.println("isFragment = true");

			int offset = 0;
			while (decodeBits(1) == 1) {
				int dataType = offset + decodeBits(4);
				switch (dataType) {
					case 0: //  datatype 'AtomQFNoMoreNeighbours'
						int no = decodeBits(abits);
						System.out.print("noMoreNeighbours:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits));
						System.out.println();
						break;
					case 1: //  datatype 'isotop'
						no = decodeBits(abits);
						System.out.print("mass:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(8));
						System.out.println();
						break;
					case 2: //  datatype 'bond defined to be delocalized'
						no = decodeBits(bbits);
						System.out.print("delocalizedBonds (outdated, redundant and wrong):");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(bbits));
						System.out.println();
						break;
					case 3: //  datatype 'AtomQFMoreNeighbours'
						no = decodeBits(abits);
						System.out.print("moreNeighbours:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits));
						System.out.println();
						break;
					case 4: //  datatype 'AtomQFRingState'
						no = decodeBits(abits);
						System.out.print("atomRingState:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFRingStateBits));
						System.out.println();
						break;
					case 5: //  datatype 'AtomQFAromState'
						no = decodeBits(abits);
						System.out.print("atomAromState:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFAromStateBits));
						System.out.println();
						break;
					case 6: //  datatype 'AtomQFAny'
						no = decodeBits(abits);
						System.out.print("atomAny:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits));
						System.out.println();
						break;
					case 7: //  datatype 'AtomQFHydrogen'
						no = decodeBits(abits);
						System.out.print("atomHydrogen:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFHydrogenBits));
						System.out.println();
						break;
					case 8: //  datatype 'AtomList'
						no = decodeBits(abits);
						System.out.print("atomList:");
						for (int i = 0; i < no; i++) {
							int atom = decodeBits(abits);
							int atoms = decodeBits(4);
							System.out.print(" " + atom);
							for (int j = 0; j < atoms; j++) {
								System.out.print(j == 0 ? ":" : ",");
								System.out.print("" + decodeBits(8));
							}
						}
						System.out.println();
						break;
					case 9: //  datatype 'BondQFRingState'
						no = decodeBits(bbits);
						System.out.print("bondRingState:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(bbits) + ":" + decodeBits(Molecule.cBondQFRingStateBits));
						System.out.println();
						break;
					case 10://  datatype 'BondQFBondTypes'
						no = decodeBits(bbits);
						System.out.print("bondTypes:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(bbits) + ":" + decodeBits(Molecule.cBondQFBondTypesBits));
						System.out.println();
						break;
					case 11:    //  datatype 'AtomQFMatchStereo'
						no = decodeBits(abits);
						System.out.print("atomMatchStereo:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits));
						System.out.println();
						break;
					case 12:    //  datatype 'bond defined to be a bridge from n1 to n2 atoms'
						no = decodeBits(bbits);
						for (int i = 0; i < no; i++) {
							System.out.print("bridgeBond:" + decodeBits(bbits));
							int min = decodeBits(Molecule.cBondQFBridgeMinBits);
							int max = min + decodeBits(Molecule.cBondQFBridgeSpanBits);
							System.out.println("(" + min + "-" + max + ")");
						}
						break;
					case 13: //  datatype 'AtomQFPiElectrons'
						no = decodeBits(abits);
						System.out.print("atomPiElectrons:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFPiElectronBits));
						System.out.println();
						break;
					case 14: //  datatype 'AtomQFNeighbours'
						no = decodeBits(abits);
						System.out.print("AtomQFNeighbours:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFNeighbourBits));
						System.out.println();
						break;
					case 15: //  datatype 'start second feature set'
					case 31:
						offset += 16;
						System.out.println("<start next feature set>");
						break;
					case 16: //  datatype 'AtomQFSmallRingSize'
						no = decodeBits(abits);
						System.out.print("AtomQFSmallRingSize:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFSmallRingSizeBits));
						System.out.println();
						break;
					case 17: //  datatype 'AtomAbnormalValence'
						no = decodeBits(abits);
						System.out.print("AtomAbnormalValence:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(4));
						System.out.println();
						break;
					case 18: //  datatype 'AtomCustomLabel'
						no = decodeBits(abits);
						System.out.print("AtomCustomLabel:");
						int lbits = decodeBits(4);
						for (int i = 0; i < no; i++) {
							int atom = decodeBits(abits);
							int count = decodeBits(lbits);
							byte[] label = new byte[count];
							for (int j = 0; j < count; j++)
								label[j] = (byte) decodeBits(7);
							System.out.print(" " + atom + ":" + new String(label));
						}
						System.out.println();
						break;
					case 19: //  datatype 'AtomQFCharge'
						no = decodeBits(abits);
						System.out.print("AtomQFCharge:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFChargeBits));
						System.out.println();
						break;
					case 20: //  datatype 'BondQFRingSize'
						no = decodeBits(bbits);
						System.out.print("BondQFRingSize:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(bbits) + ":" + decodeBits(Molecule.cBondQFRingSizeBits));
						System.out.println();
						break;
					case 21: //  datatype 'AtomRadicalState'
						no = decodeBits(abits);
						System.out.print("AtomRadicalState:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(2));
						System.out.println();
						break;
					case 22:    //	datatype 'flat nitrogen'
						no = decodeBits(abits);
						System.out.print("AtomQFFlatNitrogen:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits));
						System.out.println();
						break;
					case 23:    //	datatype 'cBondQFMatchStereo'
						no = decodeBits(bbits);
						System.out.print("cBondQFMatchStereo:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits));
						System.out.println();
						break;
					case 24:    //	datatype 'cBondQFAromatic'
						no = decodeBits(bbits);
						System.out.print("BondQFAromState:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(bbits) + ":" + decodeBits(Molecule.cBondQFAromStateBits));
						System.out.println();
						break;
					case 25:    //	datatype 'atom selection'
						System.out.print("AtomSelection:");
						for (int i = 0; i < allAtoms; i++)
							if (decodeBits(1) == 1)
								System.out.print(" " + i);
						System.out.println();
						break;
					case 26:    //	datatype 'delocalized high order bond'
						System.out.print("DelocalizedHigherOrderBonds:");
						no = decodeBits(bbits);
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(bbits));
						break;
					case 27:    //	datatype 'part of an exclude group'
						no = decodeBits(abits);
						System.out.print("AtomQFExcludeGroup:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits));
						System.out.println();
						break;
					case 28:    //	datatype 'coordinate bond'
						no = decodeBits(bbits);
						System.out.print("Coordinate Bonds:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(bbits));
						System.out.println();
						break;
					case 29:
						no = decodeBits(abits);
						System.out.print("ReactionParityHint:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFRxnParityBits));
						System.out.println();
						break;
					case 30: //  datatype 'AtomQFNewRingSize'
						no = decodeBits(abits);
						System.out.print("AtomQFNewRingSize:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFNewRingSizeBits));
						System.out.println();
						break;
					case 32: //  datatype 'AtomQFStereoStateBits'
						no = decodeBits(abits);
						System.out.print("AtomQFStereoState:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFStereoStateBits));
						System.out.println();
						break;
					case 33: //  datatype 'AtomQFENeighbours'
						no = decodeBits(abits);
						System.out.print("AtomQFENeighbours:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits) + ":" + decodeBits(Molecule.cAtomQFENeighbourBits));
						System.out.println();
						break;
					case 34:    //	datatype 'in hetero aromatic ring'
						no = decodeBits(abits);
						System.out.print("AtomQFHeteroAromatic:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits));
						System.out.println();
						break;
					case 35:    //	datatype 'cBondQFMatchFormalOrder'
						no = decodeBits(bbits);
						System.out.print("BondQFMatchFormalOrder:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(abits));
						System.out.println();
						break;
					case 36:    //	datatype 'cBondQFRareBondType'
						no = decodeBits(bbits);
						System.out.print("BondQFRareBondType:");
						for (int i = 0; i < no; i++)
							System.out.print(" " + decodeBits(bbits) + ":" + decodeBits(Molecule.cBondQFRareBondTypesBits));
						System.out.println();
						break;
					case 37:    // datatype 'rare bond type'
						no = decodeBits(bbits);
						System.out.print("Rare Bond Type:");
						for (int i=0; i<no; i++)
							System.out.print(" " + decodeBits(bbits) + ":" + (decodeBits(1) == 0 ? "quadruple" : "quintuple"));
						break;
					}
				}

			if (coordinates != null) {
				if (coordinates[0] == '!' || coordinates[0] == '#') {    // new coordinate format
					decodeBitsStart(coordinates, 1);
					boolean coordsAre3D = (decodeBits(1) == 1);
					boolean coordsAreAbsolute = (decodeBits(1) == 1);
					int resolutionBits = 2 * decodeBits(4);
					int binCount = (1 << resolutionBits);

					double factor;

					int hydrogenCount = 0;
					int[] hCount = null;
					if (coordinates[0] == '#') {    // we have 3D-coordinates that include implicit hydrogen coordinates
						StereoMolecule mol = new IDCodeParserWithoutCoordinateInvention().getCompactMolecule(idcode);

						// we need to cache hCount, because otherwise getImplicitHydrogens() would create helper arrays with every call
						hCount = new int[allAtoms];
						for (int atom = 0; atom < allAtoms; atom++)
							hydrogenCount += (hCount[atom] = mol.getImplicitHydrogens(atom));
					}

					double[][] coords = new double[coordsAre3D ? 3 : 2][allAtoms + hydrogenCount];
					int from = 0;
					int bond = 0;
					System.out.print("Raw coords:");
					for (int atom = 1; atom < allAtoms; atom++) {
						if (bond < allBonds && bondAtom[1][bond] == atom) {
							from = bondAtom[0][bond++];
							factor = 1.0;
						} else {
							from = 0;
							factor = 8.0;
						}
						System.out.print(atom + " (");
						coords[0][atom] = coords[0][from] + factor * (decodeBits(resolutionBits) - binCount / 2);
						System.out.print((int) coords[0][atom] + ",");
						coords[1][atom] = coords[1][from] + factor * (decodeBits(resolutionBits) - binCount / 2);
						System.out.print((int) coords[1][atom]);
						if (coordsAre3D) {
							coords[2][atom] = coords[2][from] + factor * (decodeBits(resolutionBits) - binCount / 2);
							System.out.print("," + (int) coords[0][atom]);
						}
						System.out.print("), ");
						if ((atom & 3) == 3 || atom == allAtoms - 1)
							System.out.println();
					}

					// with new format 2D and 3D coordinates are scaled to average bond lengths of 1.5 Angstrom
					double avbl = 0;
					if (allBonds != 0) {
						for (bond = 0; bond < allBonds; bond++)
							avbl += getDistance(coords, bondAtom[0][bond], bondAtom[1][bond], coordsAre3D);
						avbl /= allBonds;    // avbl without hydrogen atoms
					} else {
						double defaultAVBL = coordsAre3D ? 1.5 : Molecule.getDefaultAverageBondLength();
						if (allAtoms < 2) {
							avbl = defaultAVBL;
						} else {
							double lowDistance = Double.MAX_VALUE;
							for (int atom1 = 1; atom1 < allAtoms; atom1++) {
								for (int atom2 = 0; atom2 < atom1; atom2++) {
									double distance = getDistance(coords, atom1, atom2, coordsAre3D);
									if (distance > 0 && distance < lowDistance)
										lowDistance = distance;
								}
							}
							avbl = (lowDistance == Double.MAX_VALUE) ? defaultAVBL : lowDistance;
						}
					}

					if (coordinates[0] == '#') {    // we have 3D-coordinates that include implicit hydrogen coordinates
						System.out.print("hydrogen coords (" + hydrogenCount + " expected): ");
						int hydrogen = allAtoms;
						for (int atom = 0; atom < allAtoms; atom++) {
							if (hCount[atom] != 0)
								System.out.print(atom);
							for (int i = 0; i < hCount[atom]; i++) {
								System.out.print(" (");
								coords[0][hydrogen] = coords[0][atom] + (decodeBits(resolutionBits) - binCount / 2);
								System.out.print((int) coords[0][hydrogen] + ",");
								coords[1][hydrogen] = coords[1][atom] + (decodeBits(resolutionBits) - binCount / 2);
								System.out.print((int) coords[1][hydrogen]);
								if (coordsAre3D) {
									coords[2][hydrogen] = coords[2][atom] + (decodeBits(resolutionBits) - binCount / 2);
									System.out.print("," + (int) coords[2][hydrogen]);
								}
								System.out.print("), ");
								hydrogen++;
							}
						}
						System.out.println();
					}

					System.out.print(coordsAreAbsolute ? "absolute coords:" : "relative coords:");
					if (hydrogenCount != 0)
						System.out.println("Coordinates contain " + hydrogenCount + " hydrogen atoms!");

					if (coordsAreAbsolute) {
						double targetAVBL = decodeAVBL(decodeBits(resolutionBits), binCount);
						double xOffset = targetAVBL * decodeShift(decodeBits(resolutionBits), binCount);
						double yOffset = targetAVBL * decodeShift(decodeBits(resolutionBits), binCount);
						double zOffset = 0;
						if (coordsAre3D)
							zOffset = targetAVBL * decodeShift(decodeBits(resolutionBits), binCount);
						System.out.println("Abs-coord transformation: targetAVBL:" + targetAVBL + " xOffset:" + xOffset + " yOffset:" + yOffset + " zOffset:" + zOffset);

						factor = targetAVBL / avbl;
						for (int atom = 0; atom < allAtoms; atom++) {
							coords[0][atom] = xOffset + factor * coords[0][atom];
							coords[1][atom] = xOffset + factor * coords[1][atom];
							if (coordsAre3D)
								coords[2][atom] = xOffset + factor * coords[2][atom];
						}
					} else {
						double targetAVBL = 1.5;
						factor = targetAVBL / avbl;
						for (int atom = 0; atom < allAtoms; atom++) {
							System.out.print(atom + " (");
							coords[0][atom] = coords[0][atom] * factor;
							System.out.print(DoubleFormat.toString(coords[0][atom]) + ",");
							coords[1][atom] = coords[1][atom] * factor;
							System.out.print(DoubleFormat.toString(coords[1][atom]));
							if (coordsAre3D) {
								coords[2][atom] = coords[2][atom] * factor;
								System.out.print("," + DoubleFormat.toString(coords[2][atom]));
							}
							System.out.print("), ");
							if ((atom & 3) == 3 || atom == allAtoms - 1)
								System.out.println();
						}
					}
				}
			}
			System.out.println();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}

	private double getDistance(double[][] coords, int atom1, int atom2, boolean coordsAre3D) {
		double dx = coords[0][atom1] - coords[0][atom2];
		double dy = coords[1][atom1] - coords[1][atom2];
		double dz = coordsAre3D ? coords[2][atom1] - coords[2][atom2] : 0;
		return Math.sqrt(dx*dx + dy*dy + dz*dz);
		}
	}