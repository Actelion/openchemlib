package com.actelion.research.chem;

import com.actelion.research.util.SortedList;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

public class SmilesAtomParser {
	private static final int HYDROGEN_ANY = -1;

	boolean parityFound,isClockwise,smartsFeatureFound;
	;
	int atomicNo,abnormalValence,charge,mapNo,explicitHydrogens;
	private final boolean mAllowCactvs;
	private final int mMode;
	private long atomQueryFeatures;      // translated from obvious SMARTS features
	private final SmilesParser mParentParser;
	private SortedList<Integer> atomList;
	private ArrayList<StereoMolecule> recursiveSmartsList,excludeGroupList;

	public SmilesAtomParser(SmilesParser parser, int mode) {
		mParentParser = parser;
		mMode = mode;
		mAllowCactvs = (mode & SmilesParser.MODE_NO_CACTUS_SYNTAX) == 0;
		atomicNo = -1;
		charge = 0;
		mapNo = 0;
		abnormalValence = -1;
		explicitHydrogens = HYDROGEN_ANY;
		atomQueryFeatures = 0L;
	}

	private void addAtomToList(int atomicNo) {
		if (atomList == null)
			atomList = new SortedList<>();

		atomList.add(atomicNo);
	}

	protected int parseAtomOutsideBrackets(byte[] smiles, int position, int endIndex, boolean allowSmarts) {
		if (smiles[position-1] == '*') {
			atomicNo = 6;
			atomQueryFeatures |= Molecule.cAtomQFAny;
		}
		else if (smiles[position-1] == '?') {
			atomicNo = 0;
		}
		else if ((smiles[position-1] == 'A' || smiles[position-1] == 'a') && allowSmarts) {
			atomicNo = 6;
			atomQueryFeatures |= Molecule.cAtomQFAny;
			atomQueryFeatures |= smiles[position-1] == 'A' ? Molecule.cAtomQFNotAromatic : Molecule.cAtomQFAromatic;
			smartsFeatureFound = true;
		}
		else {
			switch (Character.toUpperCase(smiles[position-1])) {
				case 'B':
					if (position < endIndex && smiles[position] == 'r') {
						atomicNo = 35;
						position++;
					}
					else
						atomicNo = 5;
					break;
				case 'C':
					if (position < endIndex && smiles[position] == 'l') {
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
		return position;
	}

	private int advanceJustAfterClosingBracket(byte[] smiles, int position) throws Exception {
		int level = 0;
		while (position < smiles.length && (smiles[position] != ']' || level != 0)) {
			if (smiles[position] == '[')
				level++;
			else if (smiles[position] == ']')
				level--;

			position++;
		}

		if (position == smiles.length)
			throw new Exception("SmilesParser: No matching closing bracket found.");

		return position+1;
	}

	/**
	 * @param smiles
	 * @param position points to second character of atom description, e.g. of the atom label
	 * @param endIndex
	 * @param allowSmarts
	 * @return position of first character after closing ']' (or delimiting ',' if enumerating SMARTS)
	 * @throws Exception
	 */
	protected int parseAtomInsideBrackets(byte[] smiles, int position, int endIndex, boolean allowSmarts, boolean allowAtomOptions) throws Exception {
		if (smiles[position-1] == '$') {  // recursive SMARTS
			recursiveSmartsList = new ArrayList<>();
			position += parseRecursiveGroup(smiles, position-1, recursiveSmartsList) - 1;

			if (smiles[position++] != ']') {
				if (!allowAtomOptions)
					throw new Exception("SmilesParser: A positive recursive SMARTS followed by another one or by atom query features is not supported. Position:" + (position - 1));

				if ((mMode & SmilesParser.MODE_ENUMERATE_SMARTS) == 0)
					position = advanceJustAfterClosingBracket(smiles, position);
			}

			return position;
		}

		if (smiles[position-1] == '*') {
			atomicNo = 6;
			atomQueryFeatures |= Molecule.cAtomQFAny;
		}
		else if (smiles[position-1] == '?') {
			atomicNo = 0;
		}
		else {
			boolean isNotList = (smiles[position-1] == '!');
			if (isNotList) {
				smartsFeatureFound = true;
				atomQueryFeatures |= Molecule.cAtomQFAny;
				position++;
			}

			// Handle this before checking for atom symbols, because R<n> (ring count) takes precedence to R1 - R16 (substituent pseudo label)
			if (smiles[position-1] == 'R' && allowSmarts && (Character.isDigit(smiles[position]) || (mAllowCactvs && smiles[position] == '{'))) {
				atomicNo = 6;
				atomQueryFeatures |= Molecule.cAtomQFAny;
				position--;
				if (isNotList)
					position--;
			}
			else {
				AtomLabelInfo labelInfo = new AtomLabelInfo();
				if (!parseAtomLabelInBrackets(smiles, position-1, endIndex, labelInfo))
					throw new Exception("SmilesParser: Unexpected character in atom definition:'"+((char)smiles[position-1])+"' position:"+(position-1));

				atomicNo = labelInfo.atomicNo;
				position += labelInfo.labelLength - 1;
				if ((mMode & SmilesParser.SMARTS_MODE_MASK) != SmilesParser.SMARTS_MODE_IS_SMARTS)
					explicitHydrogens = SmilesParser.HYDROGEN_IMPLICIT_ZERO;  // in case we have SMILES; neglected, if we process a SMARTS, which we may learn later when hitting a query feature

				// If we have a comma after the first atom label, then we need to parse a (positive) atom list.
				// In this case we also have to set aromaticity query features from upper and lower case symbols.
				if (allowSmarts && (smiles[position] == ',' || isNotList)) {
					boolean mayBeAromatic = labelInfo.mayBeAromatic;
					boolean mayBeAliphatic = labelInfo.mayBeAliphatic;
					int start = position - labelInfo.labelLength;
					while (start < endIndex) {
						if (!parseAtomLabelInBrackets(smiles, start, endIndex, labelInfo)) {
							if (!isNotList)
								throw new Exception("SmilesParser: Unexpected character in atom list:'"+((char)smiles[start])+"'. Position:"+start);
							// a not-list may be followed by ';' and another atom condition, while a positive list must not end with ','
							break;
						}

						if (labelInfo.atomicNo == 1) {
							if (!isNotList) // in not-lists we are allowed to remove hydrogens!
								throw new Exception("SmilesParser: Hydrogen is not supported in positive atom lists:'"+new String(Arrays.copyOfRange(smiles, start, endIndex))+"'. Position:"+start);
						}
						else {
							addAtomToList(labelInfo.atomicNo);
							mayBeAromatic |= labelInfo.mayBeAromatic;
							mayBeAliphatic |= labelInfo.mayBeAliphatic;
						}
						start += labelInfo.labelLength;
						if (isNotList && smiles[start] != ';' && smiles[start] != '&')  // negative lists with ';' or '&', e.g. "!#7;!#8"
							break;
						if (!isNotList && smiles[start] != ',')   // positive list: ',' e.g. "N,O"
							break;
						if (isNotList && smiles[start+1] != '!')
							break;
						start++;
						if (smiles[start] == '!')
							start++;
					}

					if (atomList != null && atomList.size() > 1) {
						explicitHydrogens = HYDROGEN_ANY;   // don't use implicit zero with atom lists
						if (!mayBeAliphatic)
							atomQueryFeatures |= Molecule.cAtomQFAromatic;
						else if (!mayBeAromatic)
							atomQueryFeatures |= Molecule.cAtomQFNotAromatic;
					}

					position = start;
				}
			}
		}

		SmilesRange range = new SmilesRange(smiles);
		int[] skipCount = new int[1];
		boolean squareBracketOpen = true;

		while (squareBracketOpen) {
			if (smiles[position] == '@') {
				position++;
				if (smiles[position] == '@') {
					isClockwise = true;
					position++;
				}
				parityFound = true;
				continue;
			}

			if (smiles[position] == ':') {
				position++;
				while (Character.isDigit(smiles[position])) {
					mapNo = 10 * mapNo + smiles[position] - '0';
					position++;
				}
				continue;
			}

			if (smiles[position] == '[')
				throw new Exception("SmilesParser: nested square brackets found. Position:"+position);

			if (smiles[position] == ']') {
				position++;
				squareBracketOpen = false;
				continue;
			}

			charge = parseCharge(smiles, position, skipCount);
			if (skipCount[0] != 0) {
				position += skipCount[0];

				// explicit charge=0 is usually meant as query feature
				if (charge == 0)
					atomQueryFeatures |= Molecule.cAtomQFNotChargeNeg | Molecule.cAtomQFNotChargePos;
				continue;
			}

			boolean isNot = (smiles[position] == '!');
			if (isNot)
				position++;

			if (smiles[position] == 'H') {
				position++;
				position += range.parse(position, 1, 1);
				long flags = 0;
				if (range.min <= 0 && range.max >= 0)
					flags |= Molecule.cAtomQFNot0Hydrogen;
				if (range.min <= 1 && range.max >= 1)
					flags |= Molecule.cAtomQFNot1Hydrogen;
				if (range.min <= 2 && range.max >= 2)
					flags |= Molecule.cAtomQFNot2Hydrogen;
				if (range.min <= 3 && range.max >= 3)
					flags |= Molecule.cAtomQFNot3Hydrogen;

				if (isNot) {
					atomQueryFeatures |= flags;
					explicitHydrogens = HYDROGEN_ANY;
				}
				else {
					if (range.isSingle()) {
						explicitHydrogens = range.min;
					}
					else {
						atomQueryFeatures |= (Molecule.cAtomQFHydrogen & ~flags);
						explicitHydrogens = HYDROGEN_ANY;
					}
				}
				continue;
			}

			if (smiles[position] == 'D'      // number of explicit neighbours (incl. explicit H)
			 || smiles[position] == 'd') {   // (RDKit extension) number of non-H-neighbours
				// we translate both to the number of non-H neighbours (for 'D' we assume no explicit H to be present)
				position++;
				position += range.parse(position, 1, 1);
				long flags = 0;
				if (range.min <= 0 && range.max >= 0)
					flags |= Molecule.cAtomQFNot0Neighbours;
				if (range.min <= 1 && range.max >= 1)
					flags |= Molecule.cAtomQFNot1Neighbour;
				if (range.min <= 2 && range.max >= 2)
					flags |= Molecule.cAtomQFNot2Neighbours;
				if (range.min <= 3 && range.max >= 3)
					flags |= Molecule.cAtomQFNot3Neighbours;
				if (range.min <= 4 && range.max >= 4)
					flags |= Molecule.cAtomQFNot4Neighbours;

				if (flags != 0) {
					if (isNot)
						atomQueryFeatures |= flags;
					else if ((atomQueryFeatures & Molecule.cAtomQFNeighbours) != 0)
						atomQueryFeatures &= ~flags;
					else {
						flags = flags ^ Molecule.cAtomQFNeighbours;
						atomQueryFeatures |= flags;
					}
				}
				continue;
			}

			if (smiles[position] == 'z' && mAllowCactvs) {   // electro-negative neighbour count (CACTVS,RDKit extension)
				position++;
				position += range.parse(position, 1, 4);
				long flags = 0;
				if (range.min <= 0 && range.max >= 0)
					flags |= Molecule.cAtomQFNot0ENeighbours;
				if (range.min <= 1 && range.max >= 1)
					flags |= Molecule.cAtomQFNot1ENeighbour;
				if (range.min <= 2 && range.max >= 2)
					flags |= Molecule.cAtomQFNot2ENeighbours;
				if (range.min <= 3 && range.max >= 3)
					flags |= Molecule.cAtomQFNot3ENeighbours;
				if (range.min <= 4 && range.max >= 4)
					flags |= Molecule.cAtomQFNot4ENeighbours;

				if (flags != 0) {
					if (isNot)
						atomQueryFeatures |= flags;
					else if ((atomQueryFeatures & Molecule.cAtomQFENeighbours) != 0)
						atomQueryFeatures &= ~flags;
					else {
						flags = flags ^ Molecule.cAtomQFENeighbours;
						atomQueryFeatures |= flags;
					}
				}
				continue;
			}

			if (smiles[position] == 'X') {   // neighbour count including implicit hydrogens
				position++;
				position += range.parse(position, 1, 1);
				byte[] valences = Molecule.cAtomValence[atomicNo];
				if (valences == null)
					continue;

				int valence = valences[0];

				// if we have a locally defined charge, we update the valance properly
				int localCharge = parseCharge(smiles, position, skipCount);
				if (skipCount[0] != 0) {
					if (Molecule.isAtomicNoElectronegative(atomicNo))
						valence += localCharge;
					else if (atomicNo == 6)
						valence -= Math.abs(localCharge);
					else
						valence -= localCharge;
				}

				long flags = 0;
				// we convert into pi-electron count using standard valence
				if (valence-range.min <= 0 && valence-range.max >= 0)
					flags |= Molecule.cAtomQFNot0PiElectrons;
				if (valence-range.min <= 1 && valence-range.max >= 1)
					flags |= Molecule.cAtomQFNot1PiElectron;
				if (valence-range.min <= 2 && valence-range.max >= 2)
					flags |= Molecule.cAtomQFNot2PiElectrons;

				if (flags != 0) {
					if (isNot)
						atomQueryFeatures |= flags;
					else if ((atomQueryFeatures & Molecule.cAtomQFPiElectrons) != 0)
						atomQueryFeatures &= ~flags;
					else {
						flags = flags ^ Molecule.cAtomQFPiElectrons;
						atomQueryFeatures |= flags;
					}
				}
				continue;
			}

			if (smiles[position] == 'A' || smiles[position] == 'a') {
				atomQueryFeatures |= (isNot ^ smiles[position] == 'A') ? Molecule.cAtomQFNotAromatic : Molecule.cAtomQFAromatic;
				position++;
				continue;
			}

			if (smiles[position] == 'R') {
				position++;
				position += range.parse(position, 1, 3);
				long flags = 0;
				if (range.min <= 0 && range.max >= 0)
					flags |= Molecule.cAtomQFNotChain;
				if (range.min <= 1 && range.max >= 1)
					flags |= Molecule.cAtomQFNot2RingBonds;
				if (range.min <= 2 && range.max >= 2)
					flags |= Molecule.cAtomQFNot3RingBonds;
				if (range.min <= 3 && range.max >= 3)
					flags |= Molecule.cAtomQFNot4RingBonds;
				if (range.max > 3)
					mParentParser.smartsWarning((isNot?"!R":"R")+range.max);

				if (flags != 0) {
					if (isNot)
						atomQueryFeatures |= flags;
					else if ((atomQueryFeatures & Molecule.cAtomQFRingState) != 0)
						atomQueryFeatures &= ~flags;
					else {
						flags = flags ^ Molecule.cAtomQFRingState;
						atomQueryFeatures |= flags;
					}
				}
				continue;
			}

			if (smiles[position] == 'r') {
				position++;
				position += range.parse(position, 1, 1);
				if (range.isDefault) {
					if (isNot)
						atomQueryFeatures |= Molecule.cBondQFRingState & ~Molecule.cAtomQFNotChain;
					else
						atomQueryFeatures |= Molecule.cAtomQFNotChain;
					continue;
				}

				int ringSize = range.min;

				if (range.isRange())
					mParentParser.smartsWarning((isNot ? "!r" : "r") + range.toString());

				if (!isNot && ringSize >= 3 && ringSize <= 7)
					atomQueryFeatures |= (ringSize << Molecule.cAtomQFSmallRingSizeShift);
				else if (!range.isRange())
					mParentParser.smartsWarning((isNot ? "!r" : "r") + ringSize);
				continue;
			}

			if (smiles[position] == 'v') {
				position++;
				position += range.parse(position, 1, 1);

				int valence = range.min;

				if (range.isRange())
					mParentParser.smartsWarning((isNot ? "!v" : "v") + range.toString());

				if (!isNot && valence <= 14)
					abnormalValence = valence;
				else if (!range.isRange())
					mParentParser.smartsWarning((isNot ? "!v" : "v") + valence);
				continue;
			}

			if (smiles[position] == '^') {  // RDKit hybridisation is translated into number of pi-electrons
				position++;

				int hybridization = smiles[position++] - '0';

				if (hybridization < 1 || hybridization > 3)
					throw new Exception("SmilesParser: Unsupported hybridization. Position:"+position);

				long piElectrons = (hybridization == 1) ? Molecule.cAtomQFNot2PiElectrons
								 : (hybridization == 2) ? Molecule.cAtomQFNot1PiElectron : Molecule.cAtomQFNot0PiElectrons;

				if (!isNot)
					piElectrons = Molecule.cAtomQFPiElectrons & ~piElectrons;

				atomQueryFeatures |= piElectrons;

				continue;
			}

			if (smiles[position] == '$') {  // recursive SMARTS
				if (!isNot)
					throw new Exception("SmilesParser: non-negated recursive SMARTS relating to preceding atom are not supported yet. Position:"+position);

				if (excludeGroupList == null)
					excludeGroupList = new ArrayList<>();

				position += parseRecursiveGroup(smiles, position, excludeGroupList);
				continue;
			}

			if (allowSmarts && (smiles[position] == ';' || smiles[position] == '&')) { // we interpret high and low precendence AND the same way
				smartsFeatureFound = true;
				position++;
				continue;
			}

			if (allowSmarts && smiles[position] == ',' && isRepeatedAllowedORFeature(smiles, position, skipCount)) {    // we allow OR-logic for some query options if they have the same type
				smartsFeatureFound = true;
				position += skipCount[0] + 1;
				continue;
			}

			if (allowSmarts && smiles[position] == ',' && (mMode & SmilesParser.MODE_ENUMERATE_SMARTS) != 0) {
				smartsFeatureFound = true;
				position += 1;
				break;
			}

			if (smiles[position] == ',')
				throw new Exception("SmilesParser: alternative atom definitions not supported. (Tip: enumerate SMARTS): '"+(char)smiles[position]+"', position:"+position);

			throw new Exception("SmilesParser: unexpected character inside brackets: '"+(char)smiles[position]+"', position:"+position);
		}

		return position;
	}

	private boolean parseAtomLabelInBrackets(byte[] smiles, int position, int endIndex, AtomLabelInfo info) throws Exception {
		info.mayBeAromatic = true;
		info.mayBeAliphatic = true;
		if (smiles[position] == '#') {
			position++;
			smartsFeatureFound = true;
			info.atomicNo = 0;
			info.labelLength = 1;
			while (position < endIndex
					&& Character.isDigit(smiles[position])) {
				info.atomicNo = 10 * info.atomicNo + smiles[position] - '0';
				info.labelLength++;
				position++;
			}
			if (info.atomicNo == 0 || info.atomicNo >= Molecule.cAtomLabel.length)
				throw new Exception("SmilesParser: Atomic number out of range. position:"+(position-1));
			return true;
		}

		if (smiles[position] >= 'A' && smiles[position] <= 'Z') {
			info.labelLength = (smiles[position+1] >= 'a' && smiles[position+1] <= 'z') ? 2 : 1;
			info.atomicNo = Molecule.getAtomicNoFromLabel(new String(smiles, position, info.labelLength, StandardCharsets.UTF_8));
			if (info.labelLength == 2 && info.atomicNo == 0) {
				info.labelLength = 1;
				info.atomicNo = Molecule.getAtomicNoFromLabel(new String(smiles, position, info.labelLength, StandardCharsets.UTF_8));
			}
			info.mayBeAromatic = false;
			if (info.atomicNo == 0)
				throw new Exception("SmilesParser: Unknown atom label. position:"+(position-1));
			return true;
		}

		if ((smiles[position] == 'A' && smiles[position+1] == 's')
		 || (smiles[position] == 'S' && smiles[position+1] == 'e')) {
			info.labelLength = 2;
			info.atomicNo = Molecule.getAtomicNoFromLabel(new String(smiles, position, info.labelLength, StandardCharsets.UTF_8));
			info.mayBeAliphatic = false;
			return true;
		}

		if (smiles[position] == 'c'
		 || smiles[position] == 'n'
		 || smiles[position] == 'o'
		 || smiles[position] == 'p'
		 || smiles[position] == 's') {
			info.labelLength = 1;
			info.atomicNo = Molecule.getAtomicNoFromLabel(new String(smiles, position, info.labelLength, StandardCharsets.UTF_8));
			info.mayBeAliphatic = false;
			return true;
		}

		return false;
	}

	/**
	 * @param smiles
	 * @param position position of potential first charge symbol '+' or '-'
	 * @param characterCount receives number of characters needed for charge encoding
	 * @return extracted charge; 0: no charge defined or explicit charge=0 - distinguish by characterCount
	 */
	private int parseCharge(byte[] smiles, int position, int[] characterCount) {
		characterCount[0] = 0;
		if (smiles[position] == '+' || smiles[position] == '-') {
			byte symbol = smiles[position];
			int charge = 1;
			characterCount[0]++;
			while (smiles[position+characterCount[0]] == symbol) {
				charge++;
				characterCount[0]++;
			}
			if (charge == 1 && Character.isDigit(smiles[position+1])) {
				charge = smiles[position+1] - '0';
				characterCount[0]++;
			}
			return symbol == '+' ? charge : -charge;
		}
		return 0;
	}

	public int addParsedAtom(StereoMolecule mol, char theChar, int position) throws Exception {
		int atom = mol.addAtom(atomicNo);	// this may be a hydrogen, if defined as [H]
		mol.setAtomCharge(atom, charge);
		mol.setAtomMapNo(atom, mapNo, false);
		mol.setAtomAbnormalValence(atom, abnormalValence);
		if (atomQueryFeatures != 0) {
			if ((atomQueryFeatures & Molecule.cAtomQFAromatic) != 0) {
				atomQueryFeatures &= ~Molecule.cAtomQFAromatic;
				mol.setAtomMarker(atom, true);
			}
			else {
				mol.setAtomMarker(atom, false);
			}
			mol.setAtomQueryFeature(atom, atomQueryFeatures, true);
		}
		if (atomList != null) {
			int[] list = new int[atomList.size()];
			for (int i=0; i<atomList.size(); i++)
				list[i] = atomList.get(i);
			mol.setAtomList(atom, list);
			atomList.removeAll();
		}
		else {  // mark aromatic atoms
			if (Character.isLowerCase(theChar)) {
				if (atomicNo != 5 && atomicNo != 6 && atomicNo != 7 && atomicNo != 8 && atomicNo != 15 && atomicNo != 16 && atomicNo != 33 && atomicNo != 34)
					throw new Exception("SmilesParser: atomicNo " + atomicNo + " must not be aromatic. Position:"+(position-1));

				mol.setAtomMarker(atom, true);
			}
			else {
				mol.setAtomMarker(atom, false);
			}
		}
		if (excludeGroupList != null) {
			// Recursive group(s) as further substructure restriction relating to previously parsed atom within same square brackets.
			for (StereoMolecule group:excludeGroupList) {
				group.setAtomicNo(0, 0);    // use first atom as attachment point and neglect for now potential query features or negative atom lists
				mol.addSubstituent(group, atom);
			}
		}

		// put explicitHydrogen into atomCustomLabel to keep atom-relation when hydrogens move to end of atom list in handleHydrogen()
		if (explicitHydrogens != HYDROGEN_ANY && atomicNo != 1) {	// no custom labels for hydrogen to get useful results in getHandleHydrogenMap()
			byte[] bytes = new byte[1];
			bytes[0] = (byte)explicitHydrogens;
			mol.setAtomCustomLabel(atom, bytes);
		}

		return atom;
	}

	private int parseRecursiveGroup(byte[] smiles, int dollarIndex, ArrayList<StereoMolecule> groupList) throws Exception {
		if (smiles[dollarIndex+1] != '(')
			throw new Exception("SmilesParser: '$' for recursive SMARTS must be followed by '('. position:"+dollarIndex);

		int openBrackets = 1;
		int endIndex = dollarIndex+2;
		while (endIndex < smiles.length && openBrackets > 0) {
			if (smiles[endIndex] == '(')
				openBrackets++;
			else if (smiles[endIndex] == ')')
				openBrackets--;
			endIndex++;
		}

		if (openBrackets > 0)
			throw new Exception("SmilesParser: Missing closing ')' for recursive SMARTS. '('-position:"+(dollarIndex+1));

		StereoMolecule group = new StereoMolecule(16, 16);
		group.setFragment(true);
		SmilesParser parser = new SmilesParser(mMode);
		parser.setEnumerationPositionList(mParentParser.getEnumerationPositionList());
		parser.parse(group, smiles, dollarIndex+2, endIndex-1);
		groupList.add(group);

		if (smiles[dollarIndex-1] == '!')
			for (int atom=0; atom<group.getAllAtoms(); atom++)
				group.setAtomQueryFeature(atom, Molecule.cAtomQFExcludeGroup, true);

		return endIndex - dollarIndex;
	}

	/**
	 * If two subsequent features are delimited by comma (OR-logic), then we allow these
	 * - if they have the same type (and atom label, if an atom label is preceding), e.g. 'NX' in NX3 and NX4+
	 * - if the feature supports the logic of adding query features to previously given ones (D,R,X,z)
	 * @param smiles
	 * @param commaPosition
	 * @param skipCount int[1] to hold the number of characters to skip for atom label (0 if there is no atom label)
	 * @return true, if comma (OR-logic) is an allowed delimiter here
	 */
	private boolean isRepeatedAllowedORFeature(byte[] smiles, int commaPosition, int[] skipCount) {
		if (commaPosition < 3)
			return false;

		int index1 = commaPosition - 1;
		if (smiles[index1] == '+' || smiles[index1] == '-')
			index1--;

		if (!Character.isDigit(smiles[index1]))
			return false;

		index1--;

		if (smiles[index1] != 'D'
				&& smiles[index1] != 'R'
				&& smiles[index1] != 'X'
				&& smiles[index1] != 'z')
			return false;

		skipCount[0] = 0;
		while (index1 > 0 && Character.isLetter(smiles[index1-1])) {
			index1--;
			skipCount[0]++;
		}

		int index2 = commaPosition + 1;
		while (Character.isLetter(smiles[index1])) {
			if (smiles.length <= index2 || smiles[index1] != smiles[index2])
				return false;
			index1++;
			index2++;
		}
		return true;
	}

	public boolean atomQueryFeaturesFound() {
		return atomQueryFeatures != 0L || atomList != null;
	}

	public ArrayList<StereoMolecule> getExcludeGroupList() {
		return excludeGroupList;
	}

	public StereoMolecule getRecursiveGroup() {
		return recursiveSmartsList == null ? null : recursiveSmartsList.get(0);
	}

	private static class AtomLabelInfo {
		boolean mayBeAromatic,mayBeAliphatic;
		int atomicNo,labelLength;

		public AtomLabelInfo() {
			atomicNo = -1;
		}
	}

	private static class SmilesRange {
		private final byte[] smiles;
		private int pos;
		public int min,max;
		public boolean isDefault;

		public SmilesRange(byte[] smiles) {
			this.smiles = smiles;
		}

		public int parse(int position, int defaultMin, int defaultMax) {
			isDefault = false;
			pos = position;

			if (Character.isDigit(smiles[position])) {
				int val = parseInt();
				min = max = val;

				// If we have the same query feature, comma delimited and with different number, then we extend the range...
				int firstLetter = position-1;
				while (firstLetter > 1 && Character.isLetterOrDigit(smiles[firstLetter-1]))
					firstLetter--;
				while (smiles[pos] == ',') {
					boolean lettersMatch = true;
					int letterCount = position-firstLetter;
					for (int i=0; i<letterCount; i++) {
						if (smiles[firstLetter+i] != smiles[pos+1+i]) {
							lettersMatch = false;
							break;
						}
					}
					if (!lettersMatch)
						break;

					pos += 1+letterCount;
					val = parseInt();
					if (min > val)
						min = val;
					else if (max < val)
						max = val;
				}

				return pos - position;
			}

			if (smiles[position] == '{'
					&& Character.isDigit(smiles[position+1])) {
				pos++;
				min = parseInt();
				if (smiles[pos++] != '-')
					return 0;   // unexpected
				if (!Character.isDigit(smiles[pos]))
					return 0;   // unexpected
				max = parseInt();
				if (smiles[pos++] != '}')
					return 0;   // unexpected
				return pos - position;
			}

			min = defaultMin;
			max = defaultMax;
			isDefault = true;
			return 0;
		}

		public boolean isSingle() {
			return max == min;
		}

		public boolean isRange() {
			return max > min;
		}

		public String toString() {
			return "{"+min+"-"+max+"}";
		}

		private int parseInt() {
			int num = smiles[pos++] - '0';
			if (Character.isDigit(smiles[pos]))
				num = 10 * num + (smiles[pos++] - '0');
			return num;
		}
	}
}
