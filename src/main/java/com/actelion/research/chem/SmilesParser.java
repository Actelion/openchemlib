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

import java.util.TreeMap;


public class SmilesParser {
	private static final int MAX_BRACKET_LEVELS = 64;
	private static final int MAX_RE_CONNECTIONS = 64;
	private static final int MAX_AROMATIC_RING_SIZE = 15;
	private StereoMolecule mMol;
	private boolean[] mIsAromaticBond;

	/**
	 * Parses the given smiles into the molecule, creates proper atom coordinates
	 * to reflect correct double bond geometries and translates tetrahedral and allene
	 * parities into up/down-bonds.
	 * @param mol
	 * @param smiles
	 * @throws Exception
	 */
	public void parse(StereoMolecule mol, String smiles) throws Exception {
		parse(mol, smiles.getBytes(), true, true);
		}

	public void parse(StereoMolecule mol, byte[] smiles) throws Exception {
		parse(mol, smiles, true, true);
		}

	public void parse(StereoMolecule mol, byte[] smiles, boolean createCoordinates, boolean readStereoFeatures) throws Exception {
		mMol = mol;
		mMol.deleteMolecule();

		TreeMap<Integer,THParity> parityMap = null;

		int[] baseAtom = new int[MAX_BRACKET_LEVELS];
		baseAtom[0] = -1;

		int[] ringClosureAtom = new int[MAX_RE_CONNECTIONS];
		int[] ringClosurePosition = new int[MAX_RE_CONNECTIONS];
		boolean[] ringClosureInUse = new boolean[MAX_RE_CONNECTIONS];
		for (int i=0; i<MAX_RE_CONNECTIONS; i++)
			ringClosureAtom[i] = -1;

		int position = 0;
		int atomMass = 0;
		int fromAtom = -1;
		boolean squareBracketOpen = false;
		boolean percentFound = false;
		boolean smartsFeatureFound = false;
		int bracketLevel = 0;
		int smilesLength = smiles.length;
		int bondType = Molecule.cBondTypeSingle;

		while (smiles[position] <= 32)
			position++;

		while (position < smilesLength) {
			char theChar = (char)smiles[position++];

			if (Character.isLetter(theChar) || theChar == '*') {
				int atomicNo = 0;
				int explicitHydrogens = -1;
				boolean isWildCard = false;
				boolean parityFound = false;
				boolean isClockwise = false;
				if (squareBracketOpen) {
					if (theChar == 'R' && Character.isDigit(smiles[position])) {
						int noOfDigits = Character.isDigit(smiles[position+1]) ? 2 : 1;
						atomicNo = Molecule.getAtomicNoFromLabel(new String(smiles, position-1, 1+noOfDigits));
						position += noOfDigits;
						}
					else {
						int labelLength = Character.isLowerCase(smiles[position]) ? 2 : 1;
						atomicNo = Molecule.getAtomicNoFromLabel(new String(smiles, position-1, labelLength));
						position += labelLength-1;
						explicitHydrogens = 0;
						}

					if (smiles[position] == '@') {
						position++;
						if (smiles[position] == '@') {
							isClockwise = true;
							position++;
							}
						parityFound = true;
						}

					if (smiles[position] == 'H') {
						position++;
						explicitHydrogens = 1;
						if (Character.isDigit(smiles[position])) {
							explicitHydrogens = smiles[position] - '0';
							position++;
							}
						}
					}
				else if (theChar == '*') {
					atomicNo = 6;
					isWildCard = true;
					}
				else {
					switch (Character.toUpperCase(theChar)) {
					case 'B':
						if (position < smilesLength && smiles[position] == 'r') {
							atomicNo = 35;
							position++;
							}
						else
							atomicNo = 5;
						break;
					case 'C':
						if (position < smilesLength && smiles[position] == 'l') {
							atomicNo = 17;
							position++;
							}
						else
							atomicNo = 6;
						break;
					case 'F':
						atomicNo = 9;
						break;
					case 'I':
						atomicNo = 53;
						break;
					case 'N':
						atomicNo = 7;
						break;
					case 'O':
						atomicNo = 8;
						break;
					case 'P':
						atomicNo = 15;
						break;
					case 'S':
						atomicNo = 16;
						break;
						}
					}

				if (atomicNo == 0)
					throw new Exception("SmilesParser: unknown element label found");

				int atom = mMol.addAtom(atomicNo);	// this may be a hydrogen, if defined as [H]
				if (isWildCard) {
					smartsFeatureFound = true;
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFAny, true);
					}
				else {	// mark aromatic atoms
					mMol.setAtomMarker(atom, Character.isLowerCase(theChar));
					}

				// put explicitHydrogen into atomCustomLabel to keep atom-relation when hydrogens move to end of atom list in handleHydrogen()
				if (explicitHydrogens != -1 && atomicNo != 1) {	// no custom labels for hydrogen to get useful results in getHandleHydrogenMap()
					byte[] bytes = new byte[1];
					bytes[0] = (byte)explicitHydrogens;
					mMol.setAtomCustomLabel(atom, bytes);
					}

				fromAtom = baseAtom[bracketLevel];
				if (baseAtom[bracketLevel] != -1 && bondType != Molecule.cBondTypeDeleted) {
					mMol.addBond(atom, baseAtom[bracketLevel], bondType);
					}
				bondType = Molecule.cBondTypeSingle;
				baseAtom[bracketLevel] = atom;
				if (atomMass != 0) {
					mMol.setAtomMass(atom, atomMass);
					atomMass = 0;
					}

				if (readStereoFeatures) {
					THParity parity = (parityMap == null) ? null : parityMap.get(fromAtom);
					if (parity != null)	// if previous atom is a stereo center
						parity.addNeighbor(atom, position, atomicNo==1 && atomMass==0);

					if (parityFound) {	// if this atom is a stereo center
						if (parityMap == null)
							parityMap = new TreeMap<Integer,THParity>();
	
						// using position as hydrogenPosition is close enough
						parityMap.put(atom, new THParity(atom, fromAtom, explicitHydrogens, position, isClockwise));
						}
					}

				continue;
				}

			if (theChar == '.') {
				bondType = Molecule.cBondTypeDeleted;
				continue;
				}

			if (theChar == '=') {
				bondType = Molecule.cBondTypeDouble;
				continue;
				}

			if (theChar == '#') {
				bondType = Molecule.cBondTypeTriple;
				continue;
				}

			if (Character.isDigit(theChar)) {
				int number = theChar - '0';
				if (squareBracketOpen) {
					while (position < smilesLength
					 && Character.isDigit(smiles[position])) {
						number = 10 * number + smiles[position] - '0';
						position++;
						}
					atomMass = number;
					}
				else {
					if (percentFound
					 && position < smilesLength
					 && Character.isDigit(smiles[position])) {
						number = 10 * number + smiles[position] - '0';
						position++;
						}
					percentFound = false;
					if (number >= MAX_RE_CONNECTIONS)
						throw new Exception("SmilesParser: ringClosureAtom number out of range");
					if (ringClosureAtom[number] == -1) {
						ringClosureAtom[number] = baseAtom[bracketLevel];
						ringClosurePosition[number] = position-1;
						}
					else {
						if (ringClosureAtom[number] == baseAtom[bracketLevel])
							throw new Exception("SmilesParser: ring closure to same atom");

						if (readStereoFeatures && parityMap != null) {
							THParity parity = parityMap.get(ringClosureAtom[number]);
							if (parity != null)
								parity.addNeighbor(baseAtom[bracketLevel], ringClosurePosition[number], false);
							parity = parityMap.get(baseAtom[bracketLevel]);
							if (parity != null)
								parity.addNeighbor(ringClosureAtom[number], position-1, false);
							}

						mMol.addBond(baseAtom[bracketLevel], ringClosureAtom[number], bondType);
						ringClosureAtom[number] = -1;	// for number re-usage
						}
					bondType = Molecule.cBondTypeSingle;
					}
				continue;
				}

			if (theChar == '+') {
				if (!squareBracketOpen)
					throw new Exception("SmilesParser: '+' found outside brackets");
				int charge = 1;
				while (smiles[position] == '+') {
					charge++;
					position++;
					}
				if (charge == 1 && Character.isDigit(smiles[position])) {
					charge = smiles[position] - '0';
					position++;
					}
				mMol.setAtomCharge(baseAtom[bracketLevel], charge);
				continue;
				}

			if (theChar == '-') {
				if (!squareBracketOpen)
					continue;	// single bond

				int charge = -1;
				while (smiles[position] == '-') {
					charge--;
					position++;
					}
				if (charge == -1 && Character.isDigit(smiles[position])) {
					charge = '0' - smiles[position];
					position++;
					}
				mMol.setAtomCharge(baseAtom[bracketLevel], charge);
				continue;
				}

			if (theChar == '(') {
				if (baseAtom[bracketLevel] == -1)
					throw new Exception("Smiles with leading parenthesis are not supported");
				baseAtom[bracketLevel+1] = baseAtom[bracketLevel];
				bracketLevel++;
				continue;
				}

			if (theChar == ')') {
				bracketLevel--;
				continue;
				}

			if (theChar == '[') {
				if (squareBracketOpen)
					throw new Exception("SmilesParser: nested square brackets found");
				squareBracketOpen = true;
				continue;
				}

			if (theChar == ']') {
				if (!squareBracketOpen)
					throw new Exception("SmilesParser: closing bracket without opening one");
				squareBracketOpen = false;
				continue;
				}

			if (theChar == '%') {
				percentFound = true;
				continue;
				}

/*			if (theChar == '.') {
				if (bracketLevel != 0)
					throw new Exception("SmilesParser: '.' found within brackets");
				baseAtom[0] = -1;
//				for (int i=0; i<ringClosureAtom.length; i++)	we allow ringClosures between fragments separated by '.'
//					ringClosureAtom[i] = -1;
				continue;
				}*/

			if (theChar == ':') {
				if (!squareBracketOpen) {
					bondType = Molecule.cBondTypeDelocalized;
					continue;
					}

				int mapNo = 0;
				while (Character.isDigit(smiles[position])) {
					mapNo = 10 * mapNo + smiles[position] - '0';
					position++;
					}
				mMol.setAtomMapNo(baseAtom[bracketLevel], mapNo, false);
				continue;
				}

			if (theChar == '/') {
				if (readStereoFeatures)
					bondType = Molecule.cBondTypeUp;
				continue;	// encode slash temporarily in bondType
				}
			if (theChar == '\\') {
				if (readStereoFeatures)
					bondType = Molecule.cBondTypeDown;
				continue;	// encode backslash temporarily in bondType
				}

			throw new Exception("SmilesParser: unexpected character found: '"+theChar+"'");
			}

		// Check for unsatisfied open bonds
		if (bondType != Molecule.cBondTypeSingle)
			throw new Exception("SmilesParser: dangling open bond");
		for (int i=0; i<MAX_RE_CONNECTIONS; i++)
			if (ringClosureAtom[i] != -1)
				throw new Exception("SmilesParser: dangling ring closure");

		int[] handleHydrogenAtomMap = mMol.getHandleHydrogenMap();

		// If the number of explicitly defined hydrogens conflicts with the occupied and default valence, then set an abnormal valence.
		mMol.setHydrogenProtection(true);	// We may have a fragment. Therefore, prevent conversion of explicit H into a query feature.
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			if (mol.getAtomCustomLabel(atom) != null) {	// if we have the exact number of hydrogens
				if (!mMol.isMarkedAtom(atom)) {	// don't correct aromatic atoms
					int explicitHydrogen = mMol.getAtomCustomLabelBytes(atom)[0];
					if (mMol.getAtomicNo(atom) < Molecule.cAtomValence.length
					 && Molecule.cAtomValence[mMol.getAtomicNo(atom)] != null) {
						boolean compatibleValenceFound = false;
						int usedValence = mMol.getOccupiedValence(atom) - mMol.getElectronValenceCorrection(atom);
						for (byte valence:Molecule.cAtomValence[mMol.getAtomicNo(atom)]) {
							if (usedValence <= valence) {
								compatibleValenceFound = true;
								if (valence != usedValence+explicitHydrogen)
									mMol.setAtomAbnormalValence(atom, usedValence+explicitHydrogen);
								break;
								}
							}
						if (!compatibleValenceFound)
							mMol.setAtomAbnormalValence(atom, usedValence+explicitHydrogen);
						}
					}
				}
			}

		correctValenceExceededNitrogen();	// convert pyridine oxides and nitro into polar structures with valid nitrogen valences

		locateAromaticDoubleBonds();

		mMol.removeAtomCustomLabels();
		mMol.setHydrogenProtection(false);

		if (readStereoFeatures) {
			if (resolveStereoBonds())
				mMol.setParitiesValid(0);
			}

		if (createCoordinates || readStereoFeatures) {
			new CoordinateInventor().invent(mMol);

			if (readStereoFeatures) {
				if (parityMap != null) {
					for (THParity parity:parityMap.values())
						mMol.setAtomParity(parity.mCentralAtom, parity.calculateParity(handleHydrogenAtomMap), false);

					mMol.setParitiesValid(0);
					}

				mMol.setStereoBondsFromParity();
				mMol.setUnknownParitiesToExplicitlyUnknown();
				}
			}

		if (smartsFeatureFound)
			mMol.setFragment(true);
		}


	private void locateAromaticDoubleBonds() throws Exception {
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
		mIsAromaticBond = new boolean[mMol.getBonds()];

		// all explicitly defined aromatic bonds are taken
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondType(bond) == Molecule.cBondTypeDelocalized) {
				mMol.setBondType(bond, Molecule.cBondTypeSingle);
				mIsAromaticBond[bond] = true;
				}
			}

			// assume all bonds of small rings to be aromatic if the ring consists of aromatic atoms only
		RingCollection ringSet = new RingCollection(mMol, RingCollection.MODE_SMALL_AND_LARGE_RINGS);
		boolean[] isAromaticRing = new boolean[ringSet.getSize()];
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			int[] ringAtom = ringSet.getRingAtoms(ring);
			isAromaticRing[ring] = true;
			for (int i=0; i<ringAtom.length; i++) {
				if (!mMol.isMarkedAtom(ringAtom[i])) {
					isAromaticRing[ring] = false;
					break;
					}
				}
			if (isAromaticRing[ring]) {
				int[] ringBond = ringSet.getRingBonds(ring);
				for (int i=0; i<ringBond.length; i++)
					mIsAromaticBond[ringBond[i]] = true;
				}
			}

			// if ring bonds with two aromaticity markers are left, check whether
			// these are a member of a large ring that has all atoms marked as aromatic.
			// If yes then assume all of its bonds aromatic.
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (!mIsAromaticBond[bond]
			 && ringSet.getBondRingSize(bond) != 0
			 && mMol.isMarkedAtom(mMol.getBondAtom(0, bond))
			 && mMol.isMarkedAtom(mMol.getBondAtom(1, bond))) {
				addLargeAromaticRing(bond);
				}
			}

		mMol.ensureHelperArrays(Molecule.cHelperRings);	// to accomodate for the structure changes

			// Some Smiles contain 'aromatic' rings with atoms not being compatible
			// with a PI-bond. These include: tertiary non-charged nitrogen, [nH],
			// sulfur, non-charged oxygen, charged carbon, etc...
			// All these atoms and attached bonds are marked as handled to avoid
			// attached bonds to be promoted (changed to double bond) later.
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			if (isAromaticRing[ring]) {
				int[] ringAtom = ringSet.getRingAtoms(ring);
				for (int i=0; i<ringAtom.length; i++) {
					if (!qualifiesForPi(ringAtom[i])) {
						mMol.setAtomMarker(ringAtom[i], false);// mark: atom aromaticity handled
						for (int j=0; j<mMol.getConnAtoms(ringAtom[i]); j++)
							mIsAromaticBond[mMol.getConnBond(ringAtom[i], j)] = false;
						}
					}
				}
			}
		promoteObviousBonds();

		// promote fully delocalized 6-membered rings
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			if (isAromaticRing[ring] && ringSet.getRingSize(ring) == 6) {
				int[] ringBond = ringSet.getRingBonds(ring);
				boolean isFullyDelocalized = true;
				for (int bond:ringBond) {
					if (!mIsAromaticBond[bond]) {
						isFullyDelocalized = false;
						break;
						}
					}
				if (isFullyDelocalized) {
					promoteBond(ringBond[0]);
					promoteBond(ringBond[2]);
					promoteBond(ringBond[4]);
					promoteObviousBonds();
					}
				}
			}

			// handle remaining annelated rings (naphtalines, azulenes, etc.) starting from bridge heads (qualifyingNo=5)
			// and then handle and simple rings (qualifyingNo=4)
		boolean qualifyingBondFound;
		for (int qualifyingNo=5; qualifyingNo>=4; qualifyingNo--) {
			do {
				qualifyingBondFound = false;
				for (int bond=0; bond<mMol.getBonds(); bond++) {
					if (mIsAromaticBond[bond]) {
						int aromaticConnBonds = 0;
						for (int i=0; i<2; i++) {
							int bondAtom = mMol.getBondAtom(i, bond);
							for (int j=0; j<mMol.getConnAtoms(bondAtom); j++)
								if (mIsAromaticBond[mMol.getConnBond(bondAtom, j)])
									aromaticConnBonds++;
							}

						if (aromaticConnBonds == qualifyingNo) {
							promoteBond(bond);
							promoteObviousBonds();
							qualifyingBondFound = true;
							break;
							}
						}
					}
				} while (qualifyingBondFound);
			}

		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mIsAromaticBond[bond])
				throw new Exception("Assignment of aromatic double bonds failed");
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.isMarkedAtom(atom))
				throw new Exception("Assignment of aromatic double bonds failed");
		}


	private void addLargeAromaticRing(int bond) {
		int[] graphLevel = new int[mMol.getAtoms()];
		int graphAtom[] = new int[mMol.getAtoms()];
		int graphBond[] = new int[mMol.getAtoms()];
		int graphParent[] = new int[mMol.getAtoms()];

		int atom1 = mMol.getBondAtom(0, bond);
		int atom2 = mMol.getBondAtom(1, bond);
		graphAtom[0] = atom1;
		graphAtom[1] = atom2;
		graphBond[0] = -1;
		graphBond[1] = bond;
		graphLevel[atom1] = 1;
		graphLevel[atom2] = 2;
		graphParent[atom1] = -1;
		graphParent[atom2] = atom1;

		int current = 1;
		int highest = 1;
		while (current <= highest && graphLevel[graphAtom[current]] < MAX_AROMATIC_RING_SIZE) {
			int parent = graphAtom[current];
			for (int i=0; i<mMol.getConnAtoms(parent); i++) {
				int candidate = mMol.getConnAtom(parent, i);
				if (candidate != graphParent[parent]) {
					int candidateBond = mMol.getConnBond(parent, i);
					if (candidate == atom1) {	// ring closure
						graphBond[0] = candidateBond;
						for (int j=0; j<=highest; j++)
							mIsAromaticBond[graphBond[i]] = true;
						return;
						}
	
					if (mMol.isMarkedAtom(candidate)
					 && graphLevel[candidate] == 0) {
						highest++;
						graphAtom[highest] = candidate;
						graphBond[highest] = candidateBond;
						graphLevel[candidate] = graphLevel[parent]+1;
						graphParent[candidate] = parent;
						}
					}
				}
			current++;
			}
		return;
		}


	private boolean qualifiesForPi(int atom) {
		if ((mMol.getAtomicNo(atom) == 16 && mMol.getAtomCharge(atom) <= 0)
		 || (mMol.getAtomicNo(atom) == 6 && mMol.getAtomCharge(atom) != 0)
		 || !mMol.isMarkedAtom(atom))	// already marked as hetero-atom of another ring
			return false;

		int explicitHydrogens = (mMol.getAtomCustomLabel(atom) == null) ? 0 : mMol.getAtomCustomLabelBytes(atom)[0];
		if (mMol.getFreeValence(atom) - explicitHydrogens < 1)
			return false;

		if (mMol.getAtomicNo(atom) != 5
		 && mMol.getAtomicNo(atom) != 6
		 && mMol.getAtomicNo(atom) != 7
		 && mMol.getAtomicNo(atom) != 8
		 && mMol.getAtomicNo(atom) != 15	// P
		 && mMol.getAtomicNo(atom) != 16	// S
		 && mMol.getAtomicNo(atom) != 33	// As
		 && mMol.getAtomicNo(atom) != 34)	// Se
			return false;

		return true;
		}


	private void promoteBond(int bond) {
		if (mMol.getBondType(bond) == Molecule.cBondTypeSingle)
			mMol.setBondType(bond, Molecule.cBondTypeDouble);

		for (int i=0; i<2; i++) {
			int bondAtom = mMol.getBondAtom(i, bond);
			mMol.setAtomMarker(bondAtom, false);
			for (int j=0; j<mMol.getConnAtoms(bondAtom); j++)
				mIsAromaticBond[mMol.getConnBond(bondAtom, j)] = false;
			}
		}


	private void promoteObviousBonds() {
			// handle bond orders of aromatic bonds along the chains attached to 5- or 7-membered ring
		boolean terminalAromaticBondFound;
		do {
			terminalAromaticBondFound = false;
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (mIsAromaticBond[bond]) {
					boolean isTerminalAromaticBond = false;
					for (int i=0; i<2; i++) {
						boolean aromaticNeighbourFound = false;
						int bondAtom = mMol.getBondAtom(i, bond);
						for (int j=0; j<mMol.getConnAtoms(bondAtom); j++) {
							if (bond != mMol.getConnBond(bondAtom, j)
							 && mIsAromaticBond[mMol.getConnBond(bondAtom, j)]) {
								aromaticNeighbourFound = true;
								break;
								}
							}
						if (!aromaticNeighbourFound) {
							isTerminalAromaticBond = true;
							break;
							}
						}

					if (isTerminalAromaticBond) {
						terminalAromaticBondFound = true;
						promoteBond(bond);
						}
					}
				}
			} while (terminalAromaticBondFound);
		}

	/**
	 * This corrects N=O double bonds where the nitrogen has an exceeded valence
	 * by converting to a single bond and introducing separated charges.
	 * (e.g. pyridinoxides and nitro groups)
	 */
	private void correctValenceExceededNitrogen() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getAtomicNo(atom) == 7
			 && mMol.getAtomCharge(atom) == 0
			 && mMol.getOccupiedValence(atom) > 3
			 && mMol.getAtomPi(atom) > 0) {
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					int connAtom = mMol.getConnAtom(atom, i);
					int connBond = mMol.getConnBond(atom, i);
					if ((mMol.getBondOrder(connBond) > 1)
					 && mMol.isElectronegative(connAtom)) {
						if (mMol.getBondType(connBond) == Molecule.cBondTypeTriple)
							mMol.setBondType(connBond, Molecule.cBondTypeDouble);
						else
							mMol.setBondType(connBond, Molecule.cBondTypeSingle);
	
						mMol.setAtomCharge(atom, mMol.getAtomCharge(atom) + 1);
						mMol.setAtomCharge(connAtom, mMol.getAtomCharge(connAtom) - 1);
						break;
						}
					}
				}
			}
		}

	private boolean resolveStereoBonds() {
		mMol.ensureHelperArrays(Molecule.cHelperRings);

		boolean paritiesFound = false;
		int[] refAtom = new int[2];
		int[] refBond = new int[2];
		int[] otherAtom = new int[2];
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (!mMol.isSmallRingBond(bond)
			 && mMol.getBondType(bond) == Molecule.cBondTypeDouble) {
				for (int i=0; i<2; i++) {
					refAtom[i] = -1;
					otherAtom[i] = -1;
					int atom = mMol.getBondAtom(i, bond);
					for (int j=0; j<mMol.getConnAtoms(atom); j++) {
						int connBond = mMol.getConnBond(atom, j);
						if (connBond != bond) {
							if (mMol.getBondType(connBond) == Molecule.cBondTypeUp
							 || mMol.getBondType(connBond) == Molecule.cBondTypeDown) {
								refAtom[i] = mMol.getConnAtom(atom, j);
								refBond[i] = connBond;
								}
							else {
								otherAtom[i] = mMol.getConnAtom(atom, j);
								}
							}
						}
					if (refAtom[i] == -1)
						break;
					}
				if (refAtom[0] != -1 && refAtom[1] != -1) {
					boolean isZ = mMol.getBondType(refBond[0]) != mMol.getBondType(refBond[1]);
					boolean inversion = false;
					for (int i=0; i<2; i++) {
						if (otherAtom[i] != -1
						 && otherAtom[i] < refAtom[i])
							inversion = !inversion;
						}

					mMol.setBondParity(bond, isZ ^ inversion ? Molecule.cBondParityZor2
															 : Molecule.cBondParityEor1, false);
					paritiesFound = true;
					}
				}
			}

		// convert temporary stereo bonds back to plain single bonds
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mMol.getBondType(bond) == Molecule.cBondTypeUp
			 || mMol.getBondType(bond) == Molecule.cBondTypeDown)
				mMol.setBondType(bond, Molecule.cBondTypeSingle);

		return paritiesFound;
		}

	private class THParity {
		int mCentralAtom,mImplicitHydrogen,mFromAtom,mNeighborCount;
		int[] mNeighborAtom,mNeighborPosition;
		boolean[] mNeighborIsHydrogen;
		boolean mIsClockwise,mError;

		/**
		 * Instantiates a new parity object during smiles traversal.
		 * @param centralAtom index of atoms processed
		 * @param fromAtom index of parent atom of centralAtom (-1 if centralAtom is first atom in smiles)
		 * @param implicitHydrogen Daylight syntax: hydrogen atoms defined within square bracket of other atom
		 * @param isClockwise true if central atom is marked with @@ rather than @
		 */
		public THParity(int centralAtom, int fromAtom, int implicitHydrogen, int hydrogenPosition, boolean isClockwise) {
			if (implicitHydrogen != 0 && implicitHydrogen != 1) {
				mError = true;
				}
			else {
				mCentralAtom = centralAtom;
				mFromAtom = fromAtom;
				mImplicitHydrogen = implicitHydrogen;
				mIsClockwise = isClockwise;
				mNeighborCount = 0;
				mNeighborIsHydrogen = new boolean[4];
				mNeighborAtom = new int[4];
				mNeighborPosition = new int[4];

				// If we have a fromAtom and we have an implicit hydrogen,
				// then make the implicit hydrogen a normal neighbor.
				if (fromAtom != -1 && implicitHydrogen == 1) {
					// We put it at the end of the atom list with MAX_VALUE
					addNeighbor(Integer.MAX_VALUE, hydrogenPosition, true);
					mImplicitHydrogen = 0;
					}
				}
			}

		/**
		 * Adds a currently traversed neighbor or ring closure to parity object,
		 * which belongs to the neighbor's parent atom.
		 * In case of a ring closure the bond closure digit's position in the smiles
		 * rather than the neighbor's position is the relevant position used for parity
		 * determination.
		 * We need to track the atom, because neighbors are not necessarily added in atom
		 * sequence (ring closure with connection back to stereo center).
		 * @param position
		 * @param isHydrogen
		 */
		public void addNeighbor(int atom, int position, boolean isHydrogen) {
			if (mError)
				return;

			if (mNeighborCount == 4
			 || (mNeighborCount == 3 && mFromAtom != -1)) {
				mError = true;
				return;
				}

			mNeighborIsHydrogen[mNeighborCount] = isHydrogen;
			mNeighborAtom[mNeighborCount] = atom;
			mNeighborPosition[mNeighborCount] = position;
			mNeighborCount++;
			}

		public int calculateParity(int[] handleHydrogenAtomMap) {
			if (mError)
				return Molecule.cAtomParityUnknown;

			// We need to translate smiles-parse-time atom indexes to those that the molecule
			// uses after calling handleHydrogens, which is called from ensureHelperArrays().
			if (mFromAtom != -1)
				mFromAtom = handleHydrogenAtomMap[mFromAtom];
			for (int i=0; i<mNeighborCount; i++)
				if (mNeighborAtom[i] != Integer.MAX_VALUE)
					mNeighborAtom[i] = handleHydrogenAtomMap[mNeighborAtom[i]];

			if (mFromAtom == -1 && mImplicitHydrogen == 0) {
				// If we have no implicit hydrogen and the central atom is the first atom in the smiles,
				// then we assume that we have to take the first neighbor as from-atom (not described in Daylight theory manual).
				// Assumption: take the first neighbor as front atom, i.e. skip it when comparing positions
				int minPosition = Integer.MAX_VALUE;
				int minIndex = -1;
				for (int i=0; i<mNeighborCount; i++) {
					if (minPosition > mNeighborPosition[i]) {
						minPosition = mNeighborPosition[i];
						minIndex = i;
						}
					}
				mFromAtom = mNeighborAtom[minIndex];
				for (int i=minIndex+1; i<mNeighborCount; i++) {
					mNeighborAtom[i-1] = mNeighborAtom[i];
					mNeighborPosition[i-1] = mNeighborPosition[i];
					mNeighborIsHydrogen[i-1] = mNeighborIsHydrogen[i];
					}
				mNeighborCount--;
				}

			int totalNeighborCount = (mFromAtom == -1? 0 : 1) + mImplicitHydrogen + mNeighborCount;
			if (totalNeighborCount > 4 || totalNeighborCount < 3)
				return Molecule.cAtomParityUnknown;

			// We look from the hydrogen towards the central carbon if the fromAtom is a hydrogen or
			// if there is no fromAtom but the central atom has an implicit hydrogen.
			boolean fromAtomIsHydrogen = (mFromAtom == -1 && mImplicitHydrogen == 1)
									  || (mFromAtom != -1 && mMol.isSimpleHydrogen(mFromAtom));

			int hydrogenNeighborIndex = -1;
			for (int i=0; i<mNeighborCount; i++) {
				if (mNeighborIsHydrogen[i]) {
					if (hydrogenNeighborIndex != -1 || fromAtomIsHydrogen)
						return Molecule.cAtomParityUnknown;
					hydrogenNeighborIndex = i;
					}
				}

			// hydrogens are moved to the end of the atom list. If the hydrogen passes an odd number of
			// neighbor atoms on its way to the list end, we are effectively inverting the atom order.
			boolean isHydrogenTraversalInversion = false;
			if (hydrogenNeighborIndex != -1)
				for (int i=0; i<mNeighborCount; i++)
					if (!mNeighborIsHydrogen[i]
					 && mNeighborAtom[hydrogenNeighborIndex] < mNeighborAtom[i])
						isHydrogenTraversalInversion = !isHydrogenTraversalInversion;

			// If fromAtom is not a hydrogen, we consider it moved to highest atom index,
			// because
			boolean fromAtomTraversalInversion = false;
			if (mFromAtom != -1 && !fromAtomIsHydrogen)
				for (int i=0; i<mNeighborCount; i++)
					if (mFromAtom < mNeighborAtom[i])
						fromAtomTraversalInversion = !fromAtomTraversalInversion;

			int parity = (mIsClockwise
						^ isInverseOrder(mNeighborAtom, mNeighborPosition, mNeighborCount)
						^ fromAtomTraversalInversion
						^ isHydrogenTraversalInversion) ?
								Molecule.cAtomParity2 : Molecule.cAtomParity1;
/*
System.out.println();
System.out.println("central:"+mCentralAtom+(mIsClockwise?" @@":" @")+" from:"
				+((mFromAtom == -1)?"none":Integer.toString(mFromAtom))+" with "+mImplicitHydrogen+" hydrogens");
System.out.print("neighbors: "+mNeighborAtom[0]+"("+mNeighborPosition[0]+(mNeighborIsHydrogen[0]?",H":",non-H")+")");
for (int i=1; i<mNeighborCount; i++)
	System.out.print(", "+mNeighborAtom[i]+"("+mNeighborPosition[i]+(mNeighborIsHydrogen[i]?",H":",non-H")+")");
System.out.println();
System.out.println("parity:"+parity);
*/
			return parity;
			}

		private boolean isInverseOrder(int[] atom, int[] position, int count) {
			boolean inversion = false;
			for (int i=1; i<count; i++) {
				for (int j=0; j<i; j++) {
					if (atom[j] > atom[i])
						inversion = !inversion;
					if (position[j] > position[i])
						inversion = !inversion;
					}
				}
			return inversion;
			}
		}

	public static void main(String[] args) {
		System.out.println("ID-code equivalence test:");
		final String[][] data = { {	"N[C@@]([H])(C)C(=O)O",	"S-alanine",		"gGX`BDdwMUM@@" },
								  { "N[C@@H](C)C(=O)O",		"S-alanine",		"gGX`BDdwMUM@@" },
								  { "N[C@H](C(=O)O)C",		"S-alanine",		"gGX`BDdwMUM@@" },
								  { "[H][C@](N)(C)C(=O)O",	"S-alanine",		"gGX`BDdwMUM@@" },
								  { "[C@H](N)(C)C(=O)O",	"S-alanine",		"gGX`BDdwMUM@@" },
								  { "N[C@]([H])(C)C(=O)O",	"R-alanine",		"gGX`BDdwMUL`@" },
								  { "N[C@H](C)C(=O)O",		"R-alanine",		"gGX`BDdwMUL`@" },
								  { "N[C@@H](C(=O)O)C",		"R-alanine",		"gGX`BDdwMUL`@" },
								  { "[H][C@@](N)(C)C(=O)O",	"R-alanine",		"gGX`BDdwMUL`@" },
								  { "[C@@H](N)(C)C(=O)O",	"R-alanine",		"gGX`BDdwMUL`@" },
								  { "C[C@H]1CCCCO1",		"S-Methyl-pyran",	"gOq@@eLm]UUH`@" },
								  { "O1CCCC[C@@H]1C",		"S-Methyl-pyran",	"gOq@@eLm]UUH`@" },
								  { "[C@H](F)(B)O",			"S-Methyl-oxetan",	"gCaDDICTBSURH@" },
								  { "C1CO[C@H]1C",			"S-Methyl-oxetan",	"gKQ@@eLmUTb@" },
								  { "C1CO[C@@H](C)1",		"S-Methyl-oxetan",	"gKQ@@eLmUTb@" },
								  { "[C@H]1(C)CCO1",		"S-Methyl-oxetan",	"gKQ@@eLmUTb@" },
								  { "[H][C@]1(C)CCO1",		"S-Methyl-oxetan",	"gKQ@@eLmUTb@" },
								  { "[H][C@@]1(CCO1)C",		"S-Methyl-oxetan",	"gKQ@@eLmUTb@" },
								  { "[C@@]1([H])(C)CCO1",	"S-Methyl-oxetan",	"gKQ@@eLmUTb@" },
								  { "[C@]1(C)([H])CCO1",	"S-Methyl-oxetan",	"gKQ@@eLmUTb@" },
								  { "C1[C@@H]2COC2=N1",		"oxetan-azetin",	"gGy@LDimDvfja`@" },
								  { "CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2", "alpha-thujone", "dmLH@@RYe~IfyjjjkDaIh@" },
								  { "CN1CCC[C@H]1c2cccnc2",	"Nicotine",			"dcm@@@{IDeCEDUSh@UUECP@" },
								  { "CC[C@H](O1)CC[C@@]12CCCO2", "2S,5R-Chalcogran", "dmLD@@qJZY|fFZjjjdbH`@" },
								  { "CCCC",					"butane",			"gC`@Dij@@" },
								  { "C1C.CC1",				"butane",			"gC`@Dij@@" },
								  { "[CH3][CH2][CH2][CH3]",	"butane",			"gC`@Dij@@" },
								  { "C-C-C-C",				"butane",			"gC`@Dij@@" },
								  { "C12.C1.CC2",			"butane",			"gC`@Dij@@" },
								  { "[Na+].[Cl-]",			"NaCl",				"eDARHm@zd@@" },
								  { "[Na+]-[Cl-]",			"NaCl",				"error" },
								  { "[Na+]1.[Cl-]1",		"NaCl",				"error" },
								  { "c1ccccc1",				"benzene",			"gFp@DiTt@@@" },
								  { "C1=C-C=C-C=C1",		"benzene",			"gFp@DiTt@@@" },
								  { "C1:C:C:C:C:C:1",		"benzene",			"gFp@DiTt@@@" },
								  { "c1ccncc1",				"pyridine",			"gFx@@eJf`@@@" },
								  { "[nH]1cccc1",			"pyrrole",			"gKX@@eKcRp@" },
								  { "N1C=C-C=C1",			"pyrrole",			"gKX@@eKcRp@" },
								  { "[H]n1cccc1",			"pyrrole",			"gKX@@eKcRp@" },
								  { "[H]n1cccc1",			"pyrrole",			"gKX@@eKcRp@" },
								  { "c1cncc1",				"pyrrole no [nH]",	"error" },
								  { "[13CH4]",				"C13-methane",		"fH@FJp@" },
								  { "[35ClH]",				"35-chlorane",		"fHdP@qX`" },
								  { "[35Cl-]",				"35-chloride",		"fHtPxAbq@" },
								  { "[Na+].[O-]c1ccccc1",	"Na-phenolate",		"daxHaHCPBXyAYUn`@@@" },
								  { "c1cc([O-].[Na+])ccc1",	"Na-phenolate",		"daxHaHCPBXyAYUn`@@@" },
								  { "C[C@@](C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO",
									"Cephalostatin-1",
									"gdKe@h@@K`H@XjKHuYlnoP\\bbdRbbVTLbTrJbRaQRRRbTJTRTrfrfTTOBPHtFODPhLNSMdIERYJmShLfs]aqy|uUMUUUUUUE@UUUUMUUUUUUTQUUTPR`nDdQQKB|RIFbiQeARuQt`rSSMNtGS\\ct@@" },
									};

		StereoMolecule mol = new StereoMolecule();
		for (String[] test:data) {
			try {
				new SmilesParser().parse(mol, test[0]);
				String idcode = new Canonizer(mol).getIDCode();
				if (test[2].equals("error"))
					System.out.println("Should create error! "+test[1]+" smiles:"+test[0]+" idcode:"+idcode);
				else if (!test[2].equals(idcode))
					System.out.println("ERROR! "+test[1]+" smiles:"+test[0]+" is:"+idcode+" must:"+test[2]);
				}
			catch (Exception e) {
				if (!test[2].equals("error"))
					System.out.println("ERROR! "+test[1]+" smiles:"+test[0]+" exception:"+e.getMessage());
				}
			}
		}
	}