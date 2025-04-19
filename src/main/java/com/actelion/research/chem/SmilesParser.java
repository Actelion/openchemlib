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

import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.util.ArrayUtils;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;


public class SmilesParser {
	protected static final int SMARTS_MODE_MASK = 3;
	public static final int SMARTS_MODE_IS_SMILES = 0;
	public static final int SMARTS_MODE_GUESS = 1;
	public static final int SMARTS_MODE_IS_SMARTS = 2;

	public static final int MODE_SKIP_COORDINATE_TEMPLATES = 4;
	public static final int MODE_MAKE_HYDROGEN_EXPLICIT = 8;
	public static final int MODE_NO_CACTUS_SYNTAX = 16;  // if not set, then some CACTVS SMARTS extensions will be recognized and translated as close as possible
	public static final int MODE_SINGLE_DOT_SEPARATOR = 32;  // CONSIDER single dots '.' (rather than '..') as moelcule separator when parsing reactions
	public static final int MODE_CREATE_SMARTS_WARNING = 64;
	public static final int MODE_ENUMERATE_SMARTS = 128;

	private static final int INITIAL_CONNECTIONS = 16;
	private static final int MAX_CONNECTIONS = 100; // largest allowed one in SMILES is 99
	private static final int BRACKET_LEVELS = 32;
	private static final int MAX_AROMATIC_RING_SIZE = 15;

	// Unspecified hydrogen count within brackets means :=0 for SMILES and no-H-restriction for SMARTS.
	// Therefore, we have to distinguish from explicit H0, which defined query feature for SMARTS.
	protected static final int HYDROGEN_IMPLICIT_ZERO = 9;

	protected StereoMolecule mMol;
	private boolean[] mIsAromaticBond;
	private int mMode,mSmartsMode,mAromaticAtoms,mAromaticBonds,mCoordinateMode;
	private long mRandomSeed;
	private final boolean mCreateSmartsWarnings,mMakeHydrogenExplicit,mSingleDotSeparator;
	private StringBuilder mSmartsWarningBuffer;
	private boolean mSmartsFeatureFound;
	private ArrayList<EnumerationPosition> mEnumerationPositionList;

	/**
	 * 
	 * mBond indices; tracks bonds created by connection numbers, 
	 * as these need to be sorted for allene chirality
	 * 
	 * BH 2025.01.02 
	 * 
	 * 
	 */
	
	private Map<String, int[]> mapConnectionBonds;
	
	/**
	 * mAtom indices; tracks atoms with numeric connections, which need
	 * special numbering for allenes
	 * 
	 * BH 2025.01.02
	 * 
	 */
	private BitSet bsAtomHasConnections;
	public int[] mol2SmilesMap;
	protected byte[] smiles; // mainly for debugging
	/**
	 * Creates a new SmilesParser that doesn't allow SMARTS features to be present in
	 * parsed strings. SMARTS features cause an exception. The fragment flag of created
	 * molecules is never set.
	 */
	public SmilesParser() {
		this(SMARTS_MODE_IS_SMILES);
		}

	/**
	 * Creates a new SmilesParser that may or may not allow SMARTS features to be present in
	 * parsed strings. If smartsMode is SMARTS_MODE_IS_SMILES, then any SMARTS features cause
	 * an exception. If smartsMode is SMARTS_MODE_IS_SMARTS, then the input string is considered
	 * a SMARTS, e.g. 'CC' is taken as fragment of two non-aromatic carbon atoms connected by a
	 * single bond and without any implicit hydrogen atoms. If smartsMode is SMARTS_MODE_IS_GUESS,
	 * then the molecule is considered a substructure if any SMARTS features are discovered.
	 * Depending on whether SMARTS features are found, created molecules have the fragment flag set
	 * or not set.
	 * @param mode one of SMARTS_MODE... and optionally other mode flags
	 */
	public SmilesParser(int mode) {
		mMode = mode & ~SMARTS_MODE_MASK;
		mSmartsMode = mode & SMARTS_MODE_MASK;
		mSingleDotSeparator = (mode & MODE_SINGLE_DOT_SEPARATOR) != 0;
		mCreateSmartsWarnings = (mode & MODE_CREATE_SMARTS_WARNING) != 0;
		mMakeHydrogenExplicit = ((mode & MODE_MAKE_HYDROGEN_EXPLICIT) != 0);
		mCoordinateMode = CoordinateInventor.MODE_DEFAULT;
		if ((mode & MODE_SKIP_COORDINATE_TEMPLATES) != 0)
			mCoordinateMode |= CoordinateInventor.MODE_SKIP_DEFAULT_TEMPLATES;
		if (mMakeHydrogenExplicit)
			mCoordinateMode &= ~CoordinateInventor.MODE_REMOVE_HYDROGEN;
		}

	/**
	 * Depending on the parse() parameters, the SmilesParser may or may not generate new atom coordinates
	 * after parsing the SMILES. In difficult cases the employed CoordinateInventor uses random decisions
	 * when optimizing colliding coordinates. In strained and bridged ring systems, generated coordinates
	 * may not correctly represent all E/Z-bond configurations.
	 * Calling this method with a seed != 0 causes the creation of reproducible atom coordinates.
	 * @param seed value different from 0 in order to always create the same reproducible atom coordinates
	 */
	public void setRandomSeed(long seed) {
		mRandomSeed = seed;
		}

	public StereoMolecule parseMolecule(String smiles) {
		return smiles == null ? null : parseMolecule(smiles.getBytes(StandardCharsets.UTF_8));
		}

	/**
	 * Convenience method to quickly obtain a StereoMolecule from a SMILES string.
	 * If you process many SMILES, then the parse() methods are preferred, because
	 * they avoid the steady instantiation new StereoMolecules.
	 * @param smiles
	 * @return
	 */
	public StereoMolecule parseMolecule(byte[] smiles) {
		StereoMolecule mol = new StereoMolecule();
		try {
			parse(mol, smiles);
			}
		catch (Exception e) {
			e.printStackTrace();
			return null;
			}
		return mol;
		}

	public static boolean isReactionSmiles(byte[] smiles) {
		return isReactionSmiles(smiles, null);
		}

	public static boolean isReactionSmiles(byte[] smiles, int[] catalystCountHolder) {
		int count = 0;
		int index = -1;

		while (count < 3) {
			index = ArrayUtils.indexOf(smiles, (byte)'>', index + 1);
			while (index>0 && smiles[index - 1] == (byte)'-')
				index = ArrayUtils.indexOf(smiles, (byte)'>', index + 1);

			if (index == -1)
				break;

			count++;

			if (catalystCountHolder != null && count == 1) {
				catalystCountHolder[0] = 0;
				if (index+1<smiles.length && smiles[index+1] != '>') {
					catalystCountHolder[0] = 1;
					for (int i=index+1; i<smiles.length && (smiles[i] != '>' || smiles[i-1] == '-'); i++)
						if (smiles[i] == '.' && smiles[i-1] != '.')
							catalystCountHolder[0]++;
					}
				}
			}

		return count == 2;
		}

	public Reaction parseReaction(String smiles) throws Exception {
		return smiles == null ? null : parseReaction(smiles.getBytes(StandardCharsets.UTF_8));
		}

	public Reaction parseReaction(byte[] smiles) throws Exception {
		int index1 = ArrayUtils.indexOf(smiles, (byte)'>');
		while (index1 > 0 && smiles[index1-1] == (byte)'-')
			index1 = ArrayUtils.indexOf(smiles, (byte)'>', index1+1);

		int index2 = (index1 == -1) ? -1 : ArrayUtils.indexOf(smiles, (byte)'>', index1+1);
		while (index2 > 0 && smiles[index2-1] == (byte)'-')
			index2 = ArrayUtils.indexOf(smiles, (byte)'>', index2+1);

		if (index2 == -1)
			throw new Exception("Missing one or both separators ('>').");
		if (ArrayUtils.indexOf(smiles, (byte)'>', index2+1) != -1)
			throw new Exception("Found more than 2 separators ('>').");

		Reaction rxn = new Reaction();

		int part = 0;
		int index = 0;
		int closingGroupBracketIndex = -1;
		while (index < smiles.length) {
			while (index<smiles.length && smiles[index] == '.')
				index++;

			if (smiles[index] == '(') {   // brackets may be used to group unconnected fragments into one molecule (in case of reactions)
				if (closingGroupBracketIndex != -1)
					throw new Exception("Second open group bracket found before closing first one.");

				index++;
				int level = 0;
				for (int i=index; i<smiles.length; i++) {
					if (smiles[i] == '(') {
						level++;
						}
					else if (smiles[i] == ')') {
						if (level-- == 0) {
							closingGroupBracketIndex = i;
							break;
							}
						}
					}
				}

			int end = index;
			while (end<smiles.length
				&& smiles[end] != '>'
				&& !(smiles[end] == '.' && ((mSingleDotSeparator && closingGroupBracketIndex==-1) || closingGroupBracketIndex==end-1 || end+1==smiles.length || smiles[end+1] == '.')))
				end++;

			int molend = end;
			if (closingGroupBracketIndex == end-1) {
				molend--;
				closingGroupBracketIndex = -1;
				}

			if (index != molend) {
				StereoMolecule mol = new StereoMolecule();
				parse(mol, smiles, index, molend);
				if (mSmartsMode == SMARTS_MODE_GUESS && mSmartsFeatureFound)
					return new SmilesParser(mMode | SMARTS_MODE_IS_SMARTS).parseReaction(smiles);

				if (part == 0)
					rxn.addReactant(mol);
				else if (part == 1)
					rxn.addCatalyst(mol);
				else
					rxn.addProduct(mol);
				}

			index = end;
			while (index < smiles.length && smiles[index] == '>') {
				index++;
				part++;
				}
			}

		return rxn;
		}

	protected ArrayList<EnumerationPosition> getEnumerationPositionList() {
		return mEnumerationPositionList;
	}

	protected void setEnumerationPositionList(ArrayList<EnumerationPosition> l) {
		mEnumerationPositionList = l;
	}

	public String[] enumerateSmarts(String smarts) throws Exception {
		mEnumerationPositionList = new ArrayList<>();
		mSmartsMode = SMARTS_MODE_IS_SMARTS;
		mMode |= MODE_ENUMERATE_SMARTS;

		ArrayList<String> smartsList = new ArrayList<>();
		smartsList.add(smarts);

		try {
			parse(new StereoMolecule(), smarts);
		}
		catch (Exception e) {
			System.out.println(e.getMessage());
		}

		EnumerationPosition[] options = mEnumerationPositionList.toArray(new EnumerationPosition[0]);
		Arrays.sort(options);

		for (EnumerationPosition option : options) {
			ArrayList<String> enumeration = new ArrayList<>();
			for (String s : smartsList)
				option.enumerate(this, s.getBytes(StandardCharsets.UTF_8), enumeration);

			smartsList = enumeration;
		}

		return smartsList.toArray(new String[0]);
	}

	/**
	 * If createSmartsWarning in the constructor was passed as true, then this method
	 * returns a list of all SMARTS features, which could not be interpreted in the most recently
	 * parsed SMILES/SMARTS pattern.
	 * @return
	 */
	public String getSmartsWarning() {
		return mSmartsWarningBuffer == null ? "" : "Unresolved SMARTS features:"+mSmartsWarningBuffer;
		}

	/**
	 * Parses the given smiles into the molecule, creates proper atom coordinates
	 * to reflect correct double bond geometries and translates tetrahedral and allene
	 * parities into up/down-bonds. SMARTS features are neglected unless
	 * setAllowSmartsFeatures(true) was called before parsing.
	 * @param mol
	 * @param smiles
	 * @throws Exception
	 */
	public void parse(StereoMolecule mol, String smiles) throws Exception {
		parse(mol, smiles.getBytes(StandardCharsets.UTF_8), true, true);
		}

	public void parse(StereoMolecule mol, byte[] smiles) throws Exception {
		parse(mol, smiles, true, true);
		}

	public void parse(StereoMolecule mol, byte[] smiles, int position, int endIndex) throws Exception {
		parse(mol, smiles, position, endIndex, true, true);
		}

	public void parse(StereoMolecule mol, byte[] smiles, boolean createCoordinates, boolean readStereoFeatures) throws Exception {
		parse(mol, smiles, 0, smiles.length, createCoordinates, readStereoFeatures);
		}

	public void parse(StereoMolecule mol, byte[] smiles, int position, int endIndex, boolean createCoordinates,
			boolean readStereoFeatures) throws Exception {
		this.smiles = smiles;
		mMol = mol;
		mMol.clear();
		mapConnectionBonds = new HashMap<>();
		bsAtomHasConnections = new BitSet();

		if (mSmartsWarningBuffer != null)
			mSmartsWarningBuffer.setLength(0);

		mAromaticAtoms = 0;
		mSmartsFeatureFound = false;
		boolean allowSmarts = (mSmartsMode != SMARTS_MODE_IS_SMILES);

		TreeMap<Integer, THParity> parityMap = null;

		int[] baseAtom = new int[BRACKET_LEVELS];
		baseAtom[0] = -1;

		int[] ringClosureAtom = new int[INITIAL_CONNECTIONS];
		int[] ringClosurePosition = new int[INITIAL_CONNECTIONS];
		int[] ringClosureBondType = new int[INITIAL_CONNECTIONS];
		int[] ringClosureBondQueryFeatures = new int[INITIAL_CONNECTIONS];
		for (int i = 0; i < INITIAL_CONNECTIONS; i++)
			ringClosureAtom[i] = -1;

		int atomMass = 0;
		int fromAtom = -1;
		boolean squareBracketOpen = false;
		boolean isDoubleDigit = false;
		boolean hasLeadingBracket = false;
		int bracketLevel = 0;
		int bondType = Molecule.cBondTypeSingle;
		int bondQueryFeatures = 0;

		while (smiles[position] <= 32)
			position++;

		while (position < endIndex) {
			char theChar = (char) smiles[position++];

			// if there is an atom symbol,
			if (Character.isLetter(theChar) || theChar == '*' || theChar == '?'
					|| (theChar == '!' && allowSmarts && squareBracketOpen)
					|| (theChar == '#' && allowSmarts && squareBracketOpen)
					|| (theChar == '$' && allowSmarts && squareBracketOpen)) {
				SmilesAtomParser atomParser = new SmilesAtomParser(this, mMode | mSmartsMode);

				if (!squareBracketOpen) {
					position = atomParser.parseAtomOutsideBrackets(smiles, position, endIndex, allowSmarts);
				} else if ((mMode & MODE_ENUMERATE_SMARTS) != 0) {
					EnumerationPosition ep = new EnumerationPosition(position - 1);
					position = atomParser.parseAtomInsideBrackets(smiles, position, endIndex, true, true);
					if (smiles[position - 1] != ']') { // we have multiple options and create an option list
						while (smiles[position - 1] != ']') {
							position = atomParser.parseAtomInsideBrackets(smiles, position + 1, endIndex, true, true);
							ep.increase();
						}
						mEnumerationPositionList.add(ep);
					}
				} else {
					position = atomParser.parseAtomInsideBrackets(smiles, position, endIndex, allowSmarts, false);
				}

				squareBracketOpen = false;

				if (atomParser.getRecursiveGroup() != null) {
					fromAtom = baseAtom[bracketLevel];

					baseAtom[bracketLevel] = mol.getAllAtoms();
					mol.addMolecule(atomParser.getRecursiveGroup());

					if (fromAtom != -1 && bondType != Molecule.cBondTypeDeleted) {
						int bond = mMol.addBond(fromAtom, fromAtom, bondType);
						if (bondQueryFeatures != 0) {
							mSmartsFeatureFound = true;
							mMol.setBondQueryFeature(bond, bondQueryFeatures, true);
							mMol.adaptBondTypeToQueryFeatures(bond);
						}
					}

					// Reset bond type and query features to default.
					bondType = Molecule.cBondTypeSingle;
					bondQueryFeatures = 0;

					continue;
				}

				///////////////////////////////////////////////////////////////////////////////
				// At this position the atom is determined and the square bracket is closed! //
				///////////////////////////////////////////////////////////////////////////////

				if (atomParser.atomicNo == -1 && theChar != '?')
					throw new Exception("SmilesParser: unknown element label found. Position:" + (position - 1));

				if (atomParser.atomQueryFeaturesFound())
					mSmartsFeatureFound = true;

				int atom = atomParser.addParsedAtom(mMol, theChar, position);
				if (mMol.isMarkedAtom(atom))
					mAromaticAtoms++;

				fromAtom = baseAtom[bracketLevel];
				if (fromAtom != -1 && bondType != Molecule.cBondTypeDeleted) {
					int bond = mMol.addBond(fromAtom, atom, bondType);
					if (bondQueryFeatures != 0) {
						mSmartsFeatureFound = true;
						mMol.setBondQueryFeature(bond, bondQueryFeatures, true);
						mMol.adaptBondTypeToQueryFeatures(bond);
					}
				}

				// Reset bond type and query features to default.
				bondType = Molecule.cBondTypeSingle;
				bondQueryFeatures = 0;

				baseAtom[bracketLevel] = atom;
				if (atomMass != 0) {
					mMol.setAtomMass(atom, atomMass);
					atomMass = 0;
				}

				if (readStereoFeatures) {
					THParity parity = (parityMap == null) ? null : parityMap.get(fromAtom);
					if (parity != null) // if previous atom is a stereo center
						parity.addNeighbor(atom, position, atomParser.atomicNo == 1 && mMol.getAtomMass(atom) == 0);

					if (atomParser.parityFound) { // if this atom is a stereo center
						if (parityMap == null)
							parityMap = new TreeMap<>();

						// using position as hydrogenPosition is close enough
						int hydrogenCount = (atomParser.explicitHydrogens == HYDROGEN_IMPLICIT_ZERO) ? 0
								: atomParser.explicitHydrogens;
						// BH this is great, but it does not solve the problem for allene neighbors
						parityMap.put(atom, new THParity(atom, position - 2, fromAtom, hydrogenCount, position - 1,
								atomParser.isClockwise));
					}
				}

				continue;
			}

			if (theChar == '.') {
				baseAtom[bracketLevel] = -1;
				bondType = Molecule.cBondTypeDeleted;
				continue;
			}

			if (isBondSymbol(theChar)) {
				if (squareBracketOpen)
					throw new Exception("SmilesParser: unexpected bond symbol inside square brackets: '" + theChar
							+ "', position:" + (position - 1));

				int excludedBonds = 0;
				while (isBondSymbol(theChar)) {
					if (theChar == '!') {
						theChar = (char) smiles[position++];
						if (theChar == '@')
							bondQueryFeatures |= Molecule.cBondQFNotRing;
						else if ((theChar == '-' && smiles[position] == '>')
								|| (theChar == '<' && smiles[position] == '-')) {
							excludedBonds |= Molecule.cBondTypeMetalLigand;
							position++;
						} else if (theChar == '-')
							excludedBonds |= Molecule.cBondTypeSingle;
						else if (theChar == '=')
							excludedBonds |= Molecule.cBondTypeDouble;
						else if (theChar == '#')
							excludedBonds |= Molecule.cBondTypeTriple;
						else if (theChar == '$')
							excludedBonds |= Molecule.cBondTypeQuadruple;
						else if (theChar == ':')
							excludedBonds |= Molecule.cBondTypeDelocalized;
						else
							throw new Exception("SmilesParser: bond symbol '" + theChar
									+ "' not allowed after '!'. Position:" + (position - 1));
					} else {
						if (theChar == '@')
							bondQueryFeatures |= Molecule.cBondQFRing;
						else if (theChar == '=')
							bondType = Molecule.cBondTypeDouble;
						else if (theChar == '#')
							bondType = Molecule.cBondTypeTriple;
						else if (theChar == '$')
							bondType = Molecule.cBondTypeQuadruple;
						else if (theChar == ':')
							bondType = Molecule.cBondTypeDelocalized;
						else if (theChar == '~')
							bondQueryFeatures |= Molecule.cBondTypeSingle | Molecule.cBondTypeDouble
									| Molecule.cBondTypeTriple | Molecule.cBondTypeDelocalized;
						else if (theChar == '/') {
							if (readStereoFeatures)
								bondType = Molecule.cBondTypeUp; // encode slash temporarily in bondType
						} else if (theChar == '\\') {
							if (readStereoFeatures)
								bondType = Molecule.cBondTypeDown; // encode slash temporarily in bondType
						}

						// Smiles extention 'dative bond'
						else if ((theChar == '-' && smiles[position] == '>')
								|| (theChar == '<' && smiles[position] == '-')) {
							bondType = Molecule.cBondTypeMetalLigand;
							position++;
						}

						if (smiles[position] == ',') {
							bondQueryFeatures |= bondSymbolToQueryFeature(
									bondType == Molecule.cBondTypeMetalLigand ? '>' : theChar);
							while (smiles[position] == ',') {
								if ((smiles[position + 1] == '<' && smiles[position + 2] == '-')
										|| (smiles[position + 1] == '-' && smiles[position + 2] == '>')) {
									bondQueryFeatures |= bondSymbolToQueryFeature('>');
									position += 3;
								} else {
									bondQueryFeatures |= bondSymbolToQueryFeature((char) smiles[position + 1]);
									position += 2;
								}
							}
						}
					}

					if (smiles[position] == ';') {
						position++;
						theChar = (char) smiles[position++];
						continue;
					}

					if (excludedBonds != 0)
						bondQueryFeatures |= Molecule.cBondQFBondTypes & ~excludedBonds;

					break;
				}

				continue;
			}

			if (theChar <= ' ') { // we stop reading at whitespace
				position = endIndex;
				continue;
			}

			if (Character.isDigit(theChar)) {
				int number = theChar - '0';
				if (squareBracketOpen) {
					while (position < endIndex && Character.isDigit(smiles[position])) {
						number = 10 * number + smiles[position] - '0';
						position++;
					}
					atomMass = number;
				} else {
					int bondTypePosition = isDoubleDigit ? position - 3 : position - 2;
					boolean hasBondType = (smiles[bondTypePosition] == '-' || smiles[bondTypePosition] == '/'
							|| smiles[bondTypePosition] == '\\' || smiles[bondTypePosition] == '='
							|| smiles[bondTypePosition] == '#' || smiles[bondTypePosition] == '$'
							|| smiles[bondTypePosition] == ':' || smiles[bondTypePosition] == '>'
							|| smiles[bondTypePosition] == '~');
					if (isDoubleDigit && position < endIndex && Character.isDigit(smiles[position])) {
						number = 10 * number + smiles[position] - '0';
						isDoubleDigit = false;
						position++;
					}
					if (number >= ringClosureAtom.length) {
						if (number >= MAX_CONNECTIONS)
							throw new Exception("SmilesParser: ringClosureAtom number out of range: " + number);

						int oldSize = ringClosureAtom.length;
						int newSize = ringClosureAtom.length;
						while (newSize <= number)
							newSize = Math.min(MAX_CONNECTIONS, newSize + INITIAL_CONNECTIONS);

						ringClosureAtom = Arrays.copyOf(ringClosureAtom, newSize);
						ringClosurePosition = Arrays.copyOf(ringClosurePosition, newSize);
						ringClosureBondType = Arrays.copyOf(ringClosureBondType, newSize);
						ringClosureBondQueryFeatures = Arrays.copyOf(ringClosureBondQueryFeatures, newSize);
						for (int i = oldSize; i < newSize; i++)
							ringClosureAtom[i] = -1;
					}
					if (ringClosureAtom[number] == -1) {
						ringClosureAtom[number] = baseAtom[bracketLevel];
						ringClosurePosition[number] = position - 1;
						ringClosureBondType[number] = hasBondType ? bondType : -1;
						ringClosureBondQueryFeatures[number] = hasBondType ? bondQueryFeatures : 0;
					} else {
						if (ringClosureAtom[number] == baseAtom[bracketLevel])
							throw new Exception("SmilesParser: ring closure to same atom");

						if (readStereoFeatures && parityMap != null) {
							THParity parity = parityMap.get(ringClosureAtom[number]);
							if (parity != null)
								parity.addNeighbor(baseAtom[bracketLevel], ringClosurePosition[number], false);
							parity = parityMap.get(baseAtom[bracketLevel]);
							if (parity != null)
								parity.addNeighbor(ringClosureAtom[number], position - 1, false);
						}

						if (ringClosureBondType[number] != -1)
							bondType = ringClosureBondType[number];
						else if (bondType == Molecule.cBondTypeUp) // interpretation inverts, if we have the slash bond
																	// at the second closure digit rather than at the
																	// first
							bondType = Molecule.cBondTypeDown;
						else if (bondType == Molecule.cBondTypeDown)
							bondType = Molecule.cBondTypeUp;
						// ringClosureAtom is the parent atom, i.e. the baseAtom of the first occurrence
						// of the closure digit
						int a1 = ringClosureAtom[number];
						int a2 = baseAtom[bracketLevel];
						int bond = addConnection(a1, a2, bondType, ringClosurePosition[number], position - 1);
						if (ringClosureBondQueryFeatures[number] != 0)
							bondQueryFeatures = ringClosureBondQueryFeatures[number];
						if (bondQueryFeatures != 0) {
							mSmartsFeatureFound = true;
							mMol.setBondQueryFeature(bond, ringClosureBondQueryFeatures[number], true);
							mMol.adaptBondTypeToQueryFeatures(bond);
						}
						ringClosureAtom[number] = -1; // for number re-usage
					}
					bondType = Molecule.cBondTypeSingle;
					bondQueryFeatures = 0;
				}
				continue;
			}

			if (theChar == '+') {
				throw new Exception("SmilesParser: '+' found outside brackets. Position:" + (position - 1));
			}

			if (theChar == '(') {
				if (baseAtom[bracketLevel] == -1) {
					// Leading '(' are superfluous and not good style, but we allow and ignore them
					// including their closing counterparts
					hasLeadingBracket = true;
					continue;
				}
				bracketLevel++;
				if (baseAtom.length == bracketLevel)
					baseAtom = Arrays.copyOf(baseAtom, baseAtom.length + BRACKET_LEVELS);
				baseAtom[bracketLevel] = baseAtom[bracketLevel - 1];
				continue;
			}

			if (theChar == ')') {
				if (bracketLevel == 0) {
					if (!hasLeadingBracket)
						throw new Exception(
								"SmilesParser: Closing ')' without opening counterpart. Position:" + (position - 1));
					baseAtom[0] = -1;
					hasLeadingBracket = false; // we allow for a new leading '(', e.g. after '.'
					continue;
				}
				bracketLevel--;
				continue;
			}

			if (theChar == '[') {
				squareBracketOpen = true;
				continue;
			}

			if (theChar == ']') {
				throw new Exception("SmilesParser: closing bracket at unexpected position:" + (position - 1));
			}

			if (theChar == '%') {
				isDoubleDigit = true;
				continue;
			}

			/*
			 * if (theChar == '.') { if (bracketLevel != 0) throw new
			 * Exception("SmilesParser: '.' found within brackets"); baseAtom[0] = -1; //
			 * for (int i=0; i<ringClosureAtom.length; i++) we allow ringClosures between
			 * fragments separated by '.' // ringClosureAtom[i] = -1; continue; }
			 */

			throw new Exception("SmilesParser: unexpected character outside brackets: '" + theChar + "' position:"
					+ (position - 1));
		}

		// Check for unsatisfied open bonds
		if (bondType != Molecule.cBondTypeSingle)
			throw new Exception("SmilesParser: dangling open bond");
		for (int rca : ringClosureAtom)
			if (rca != -1)
				throw new Exception("SmilesParser: dangling ring closure.");

		int[] handleHydrogenAtomMap = mMol.getHandleHydrogenMap();

		// If the number of explicitly defined hydrogens conflicts with the occupied and
		// default valence,
		// then try to change radical state to compensate. If that is impossible, then
		// set an abnormal valence.
		mMol.setHydrogenProtection(true); // We may have a fragment. Therefore, prevent conversion of explicit H into a
											// query feature.
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
			if (mMol.getAtomCustomLabel(atom) != null) { // if we have the exact number of hydrogens
				int explicitHydrogen = mMol.getAtomCustomLabelBytes(atom)[0];

				if (mSmartsFeatureFound || mSmartsMode == SMARTS_MODE_IS_SMARTS) {
					if (explicitHydrogen != HYDROGEN_IMPLICIT_ZERO) { // in SMARTS there is not implicit zero hydrogen!
						if (mMakeHydrogenExplicit) {
							for (int i = 0; i < explicitHydrogen; i++)
								mMol.addBond(atom, mMol.addAtom(1), 1);
						} else {
							if (explicitHydrogen == 0)
								mMol.setAtomQueryFeature(atom, Molecule.cAtomQFHydrogen & ~Molecule.cAtomQFNot0Hydrogen,
										true);
							if (explicitHydrogen == 1)
								mMol.setAtomQueryFeature(atom, Molecule.cAtomQFHydrogen & ~Molecule.cAtomQFNot1Hydrogen,
										true);
							if (explicitHydrogen == 2)
								mMol.setAtomQueryFeature(atom, Molecule.cAtomQFHydrogen & ~Molecule.cAtomQFNot2Hydrogen,
										true);
							if (explicitHydrogen == 3)
								mMol.setAtomQueryFeature(atom, Molecule.cAtomQFHydrogen & ~Molecule.cAtomQFNot3Hydrogen,
										true);
						}
					}
				} else {
					if (explicitHydrogen == HYDROGEN_IMPLICIT_ZERO)
						explicitHydrogen = 0;

					if (!mMol.isMetalAtom(atom) && (!mMol.isMarkedAtom(atom)
							|| (mMol.getAtomicNo(atom) == 6 && mMol.getAtomCharge(atom) == 0))) {
						// We don't correct aromatic non-carbon atoms, because for these the number of
						// explicit hydrogens encodes whether a pi-bond needs to be placed at the atom
						// when resolving aromaticity. Same applies for charged carbon.
						byte[] valences = Molecule.getAllowedValences(mMol.getAtomicNo(atom));
						boolean compatibleValenceFound = false;
						int usedValence = mMol.getOccupiedValence(atom);
						usedValence -= mMol.getElectronValenceCorrection(atom, usedValence);
						usedValence += explicitHydrogen;
						if (mMol.isMarkedAtom(atom)) // later we will use a valence for the pi-bond
							usedValence++;
						for (byte valence : valences) {
							if (usedValence <= valence) {
								compatibleValenceFound = true;
								if (valence == usedValence + 2)
									mMol.setAtomRadical(atom, Molecule.cAtomRadicalStateT);
								else if (valence == usedValence + 1)
									mMol.setAtomRadical(atom, Molecule.cAtomRadicalStateD);
								else if (valence != usedValence || valence != valences[0])
									mMol.setAtomAbnormalValence(atom, usedValence);
								break;
							}
						}
						if (!compatibleValenceFound)
							mMol.setAtomAbnormalValence(atom, usedValence);
					}

					if (mMakeHydrogenExplicit || !mMol.supportsImplicitHydrogen(atom))
						for (int i = 0; i < explicitHydrogen; i++)
							mMol.addBond(atom, mMol.addAtom(1), 1);
				}
			} else if (!mMakeHydrogenExplicit && (mSmartsFeatureFound || mSmartsMode == SMARTS_MODE_IS_SMARTS)) {
				// if we don't have a hydrogen count on the atom, but we have explicit hydrogen
				// atoms
				// and if we decode a SMARTS, then we convert explicit hydrogens into an 'at
				// least n hydrogen'
				int explicitHydrogen = mMol.getExplicitHydrogens(atom);
				if (explicitHydrogen >= 1)
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot0Hydrogen, true);
				if (explicitHydrogen >= 2)
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot1Hydrogen, true);
				if (explicitHydrogen >= 3)
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot2Hydrogen, true);
				if (explicitHydrogen >= 4)
					mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot3Hydrogen, true);
			}
		}

		if (!mMakeHydrogenExplicit && (mSmartsFeatureFound || mSmartsMode == SMARTS_MODE_IS_SMARTS))
			mMol.removeExplicitHydrogens(false);

		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);

		correctValenceExceededNitrogen(); // convert pyridine oxides and nitro into polar structures with valid nitrogen
											// valences

		locateAromaticDoubleBonds(allowSmarts, mSmartsFeatureFound);

		mMol.removeAtomCustomLabels();
		mMol.setHydrogenProtection(false);

		if (readStereoFeatures) {
			assignKnownEZBondParities();

			if (parityMap != null) {
				for (THParity parity : parityMap.values())
					mMol.setAtomParity(handleHydrogenAtomMap[parity.mCentralAtom],
							parity.calculateParity(handleHydrogenAtomMap), false);

				mMol.setParitiesValid(0);
			}
		}

		// defines unknown EZ parities as such, i.e. prevent coordinate generation to
		// create implicit EZ-parities
		mMol.setParitiesValid(0);

		if (createCoordinates) {
			CoordinateInventor inventor = new CoordinateInventor(mCoordinateMode);
			if (mRandomSeed != 0)
				inventor.setRandomSeed(mRandomSeed);
			inventor.invent(mMol);

			if (readStereoFeatures)
				mMol.setUnknownParitiesToExplicitlyUnknown();
		}

		if (mSmartsFeatureFound || mSmartsMode == SMARTS_MODE_IS_SMARTS) {
			mMol.setFragment(true);
			mMol.validateAtomQueryFeatures();
			mMol.validateBondQueryFeatures();
		}
	}

	private int addConnection(int a1, int a2, int bondType, int pos1, int pos2) {
		// BH 2025.01.02 for allene parities
		int bond = mMol.addBond(a1,  a2, bondType);
		mapConnectionBonds.put(a1 + "_" + a2, new int[] {pos1, pos2});
		mapConnectionBonds.put(a2 + "_" + a1, new int[] {pos1, pos2});
		bsAtomHasConnections.set(a1);
		bsAtomHasConnections.set(a2);
		//System.out.println("bond " + bond + " connects " + a1 + " and " + a2 + " pos" + pos1 + " " + pos2);
		return bond;
	}

	/**
	 * @return true if the previously parsed SMILES contained a SMARTS feature and was not parsed with SMARTS_MODE_IS_SMILES
	 */
	public boolean isSmarts() {
		return mSmartsFeatureFound;
	}

	private boolean isBondSymbol(char theChar) {
		return theChar == '-'
			|| theChar == '='
			|| theChar == '#'
			|| theChar == '$'
			|| theChar == ':'
			|| theChar == '/'
			|| theChar == '\\'
			|| theChar == '<'
			|| theChar == '~'
			|| theChar == '!'
			|| theChar == '@';
		}

	private int bondSymbolToQueryFeature(char symbol) {
		return symbol == '=' ? Molecule.cBondTypeDouble
			 : symbol == '#' ? Molecule.cBondTypeTriple
			 : symbol == '$' ? Molecule.cBondTypeQuadruple
			 : symbol == ':' ? Molecule.cBondTypeDelocalized
			 : symbol == '>' ? Molecule.cBondTypeMetalLigand
			 : symbol == '~' ? Molecule.cBondQFBondTypes : Molecule.cBondTypeSingle;
		}

	protected void smartsWarning(String feature) {
		if (mCreateSmartsWarnings) {
			if (mSmartsWarningBuffer == null)
				mSmartsWarningBuffer = new StringBuilder();

			mSmartsWarningBuffer.append(" ");
			mSmartsWarningBuffer.append(feature);
			}
		}

	private void locateAromaticDoubleBonds(boolean allowSmartsFeatures, boolean smartsFeatureFound) throws Exception {
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
		mIsAromaticBond = new boolean[mMol.getBonds()];
		mAromaticBonds = 0;

		// all explicitly defined aromatic bonds are taken
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondType(bond) == Molecule.cBondTypeDelocalized) {
				mMol.setBondType(bond, Molecule.cBondTypeSingle);
				mIsAromaticBond[bond] = true;
				mAromaticBonds++;
				}
			}

		boolean[] isAromaticRingAtom = new boolean[mMol.getAtoms()];

		// assume all bonds of small rings to be aromatic if the ring consists of aromatic atoms only
		RingCollection ringSet = new RingCollection(mMol, RingCollection.MODE_SMALL_AND_LARGE_RINGS);
		boolean[] isAromaticRing = new boolean[ringSet.getSize()];
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			int[] ringAtom = ringSet.getRingAtoms(ring);
			isAromaticRing[ring] = true;
			for (int j : ringAtom) {
				if (!mMol.isMarkedAtom(j)) {
					isAromaticRing[ring] = false;
					break;
				}
			}
			if (isAromaticRing[ring]) {
				for (int j : ringAtom) isAromaticRingAtom[j] = true;

				int[] ringBond = ringSet.getRingBonds(ring);
				for (int j : ringBond) {
					if (!mIsAromaticBond[j]) {
						mIsAromaticBond[j] = true;
						mAromaticBonds++;
					}
				}
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

		// If both atoms of a bond are marked as aromatic and
		// if none of the two atoms is a member of a fully aromatic ring,
		// then assume the bond to be an aromatic one.
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (!mIsAromaticBond[bond]) {
				int atom1 = mMol.getBondAtom(0, bond);
				int atom2 = mMol.getBondAtom(1, bond);
				if (!isAromaticRingAtom[atom1]
				 && !isAromaticRingAtom[atom2]
				 && mMol.isMarkedAtom(atom1)
				 && mMol.isMarkedAtom(atom2)) {
					mIsAromaticBond[bond] = true;
					mAromaticBonds++;
					}
				}
			}

		mMol.ensureHelperArrays(Molecule.cHelperRings);	// to accommodate for the structure changes

		if (mSmartsMode == SMARTS_MODE_IS_SMARTS
		 || (mSmartsMode == SMARTS_MODE_GUESS && mSmartsFeatureFound))
			protectExplicitAromaticBonds();

		// Since Smiles don't have aromaticity information about bonds, we assume that all
		// bonds of a ring are aromatic if all of its atoms are aromatic. This is not always true
		// (e.g. in fbc@@@LdbbbbbRJvcEBMIpTqrAD@@@@@@@@), which may lead to wrong resolution of
		// conjugated double bonds leaving unpaired single aromatic atoms.
		// We cache the (untrustworthy) isAromaticBond array to later find paths between single
		// aromatic atoms.
		boolean[] isAromaticBond = new boolean[mMol.getBonds()];
		if (mMol.getBonds()>=0) System.arraycopy(mIsAromaticBond, 0, isAromaticBond, 0, mMol.getBonds());

			// Some Smiles contain 'aromatic' rings with atoms not being compatible
			// with a PI-bond. These include: tertiary non-charged nitrogen, [nH],
			// sulfur, non-charged oxygen, charged carbon, etc...
			// All these atoms and attached bonds are marked as handled to avoid
			// attached bonds to be promoted (changed to double bond) later.
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			if (isAromaticRing[ring]) {
				int[] ringAtom = ringSet.getRingAtoms(ring);
				for (int k : ringAtom) {
					if (!qualifiesForPi(k)) {
						if (mMol.isMarkedAtom(k)) {
							mMol.setAtomMarker(k, false);// mark: atom aromaticity handled
							mAromaticAtoms--;
						}
						for (int j = 0; j<mMol.getConnAtoms(k); j++) {
							int connBond = mMol.getConnBond(k, j);
							if (mIsAromaticBond[connBond]) {
								mIsAromaticBond[connBond] = false;
								mAromaticBonds--;
							}
						}
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

		/* Some SMILES still contain kekulized aromatic rings with lowercase hetero atoms, e.g. C1=CC=C[se]1
		// If we recognize those rings to be aromatic, we remove the aromaticity marker from all ring atoms.
		if (mAromaticAtoms != 0) {
			RingCollection daylightTypeRingSet = new RingCollection(mMol, RingCollection.MODE_SMALL_RINGS_AND_AROMATICITY | RingCollection.MODE_INCLUDE_TAUTOMERIC_BONDS);
			for (int ring=0; ring<daylightTypeRingSet.getSize(); ring++) {
				if (daylightTypeRingSet.isAromatic(ring)) {
					int[] ringAtom = daylightTypeRingSet.getRingAtoms(ring);
					for (int atom:ringAtom)
						if (mMol.isMarkedAtom(atom)) {
							mMol.setAtomMarker(atom, false);
							mAromaticAtoms--;
						}
					}
				}
			} taken out, because OpenChemLib should not start interpreting invalid SMILES; TLS 20-Oct-2021 */

		while (mAromaticAtoms >= 2)
			if (!connectConjugatedRadicalPairs(isAromaticBond))
				break;

		if (allowSmartsFeatures) {
			if (mAromaticAtoms != 0) {
				for (int atom=0; atom<mMol.getAtoms(); atom++) {
					if (mMol.isMarkedAtom(atom)) {
						mMol.setAtomMarker(atom, false);
						mMol.setAtomQueryFeature(atom, Molecule.cAtomQFAromatic, true);
						mAromaticAtoms--;
						smartsFeatureFound = true;
						}
					}
				}
			if (mAromaticBonds != 0) {
				for (int bond=0; bond<mMol.getBonds(); bond++) {
					if (mIsAromaticBond[bond]) {
						mIsAromaticBond[bond] = false;
						mMol.setBondType(bond, Molecule.cBondTypeDelocalized);
						mAromaticBonds--;
						smartsFeatureFound = true;
						}
					}
				}
			}
		else {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mMol.isMarkedAtom(atom) && mMol.getImplicitHydrogens(atom) != 0) {
					mMol.setAtomMarker(atom, false);
					mMol.setAtomRadical(atom, Molecule.cAtomRadicalStateD);
					mAromaticAtoms--;
					}
				}
			}

		if ((mSmartsMode == SMARTS_MODE_IS_SMILES)
		 || (mSmartsMode == SMARTS_MODE_GUESS && !smartsFeatureFound)) {
			if (mAromaticAtoms != 0)
				throw new Exception("Assignment of aromatic double bonds failed");
			if (mAromaticBonds != 0)
				throw new Exception("Assignment of aromatic double bonds failed");
			}
		}


	private boolean connectConjugatedRadicalPairs(boolean[] isAromaticBond) {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.isMarkedAtom(atom)) {
				int[] graphLevel = new int[mMol.getAtoms()];
				int[] graphAtom = new int[mMol.getAtoms()];
				int[] graphParent = new int[mMol.getAtoms()];

				graphAtom[0] = atom;
				graphLevel[atom] = 1;
				graphParent[atom] = -1;
				int current = 0;
				int highest = 0;
				while (current <= highest) {
					int bondOrder = ((graphLevel[graphAtom[current]] & 1) == 1) ? 1 : 2;
					for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
						int bond = mMol.getConnBond(graphAtom[current], i);
						if (mMol.getBondOrder(bond) == bondOrder && isAromaticBond[bond]) {
							int candidate = mMol.getConnAtom(graphAtom[current], i);
							if (graphLevel[candidate] == 0) {
								if (bondOrder == 1 && mMol.isMarkedAtom(candidate)) {
									int parent = graphAtom[current];
									while (parent != -1) {
										mMol.setBondType(mMol.getBond(candidate,  parent), bondOrder == 1 ?
												Molecule.cBondTypeDouble : Molecule.cBondTypeSingle);
										bondOrder = 3 - bondOrder;
										candidate = parent;
										parent = graphParent[parent];
										}

									mMol.setAtomMarker(atom, false);
									mMol.setAtomMarker(candidate, false);
									mAromaticAtoms -= 2;
									return true;
									}

								graphAtom[++highest] = candidate;
								graphParent[candidate] = graphAtom[current];
								graphLevel[candidate] = graphLevel[graphAtom[current]]+1;
								}
							}
						}
					current++;
					}
				}
			}
		return false;
		}

	private void addLargeAromaticRing(int bond) {
		int[] graphLevel = new int[mMol.getAtoms()];
		int[] graphAtom = new int[mMol.getAtoms()];
		int[] graphBond = new int[mMol.getAtoms()];
		int[] graphParent = new int[mMol.getAtoms()];

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
						for (int j=0; j<=highest; j++) {
							if (!mIsAromaticBond[graphBond[i]]) {
								mIsAromaticBond[graphBond[i]] = true;
								mAromaticBonds++;
								}
							}
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
		if (!RingCollection.qualifiesAsAromaticAtomicNo(mMol.getAtomicNo(atom)))
			return false;

		if (mMol.getAtomicNo(atom) == 6) {
			if (!mMol.isMarkedAtom(atom))	// already member of another aromatic ring
				return false;

			if (mMol.getAtomCharge(atom) > 0)
				return false;
			}

		int explicitHydrogens = (mMol.getAtomCustomLabel(atom) == null || mMol.getAtomCustomLabelBytes(atom)[0] == HYDROGEN_IMPLICIT_ZERO) ?
								0 : mMol.getAtomCustomLabelBytes(atom)[0];
		int freeValence = mMol.getFreeValence(atom) - explicitHydrogens;
		if (freeValence < 1)
			return false;

		if (mMol.getAtomicNo(atom) == 16
		 || mMol.getAtomicNo(atom) == 34
		 || mMol.getAtomicNo(atom) == 52) {
			if (mMol.getConnAtoms(atom) == 2 && mMol.getAtomCharge(atom) <= 0)
				return false;
			return freeValence != 2;	// e.g. -S(=O)- correction to account for tetravalent S,Se
			}

		return true;
		}


	private void promoteBond(int bond) {
		if (mMol.getBondType(bond) == Molecule.cBondTypeSingle)
			mMol.setBondType(bond, Molecule.cBondTypeDouble);

		for (int i=0; i<2; i++) {
			int bondAtom = mMol.getBondAtom(i, bond);
			if (mMol.isMarkedAtom(bondAtom)) {
				mMol.setAtomMarker(bondAtom, false);
				mAromaticAtoms--;
				}
			for (int j=0; j<mMol.getConnAtoms(bondAtom); j++) {
				int connBond = mMol.getConnBond(bondAtom, j);
				if (mIsAromaticBond[connBond]) {
					mIsAromaticBond[connBond] = false;
					mAromaticBonds--;
					}
				}
			}
		}

	/**
	 * In general, we kekulize the bonds of any aromatic ring even if the input is SMARTS.
	 * For non-ring chains of lower case atoms we also assume that alternating double-single bonds
	 * are the desired outcome. If the input is SMARTS, however, an aromatic bond may be meant
	 * as query feature. For now, we only consider separated aromatic bonds as intentional
	 * query feature. We might consider doing that for non-ring chains as well.
	 */
	private void protectExplicitAromaticBonds() {
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mIsAromaticBond[bond]) {
				boolean isSingleAromaticBond = true;
				for (int i=0; i<2 && isSingleAromaticBond; i++) {
					int bondAtom = mMol.getBondAtom(i, bond);
					for (int j=0; j<mMol.getConnAtoms(bondAtom) && isSingleAromaticBond; j++) {
						if (bond != mMol.getConnBond(bondAtom, j)
						 && mIsAromaticBond[mMol.getConnBond(bondAtom, j)])
							isSingleAromaticBond = false;
					}
				}

				if (isSingleAromaticBond) {
					mMol.setBondType(bond, Molecule.cBondTypeDelocalized);
					mAromaticBonds--;
					for (int i=0; i<2; i++) {
						int bondAtom = mMol.getBondAtom(i, bond);
						if (mMol.isMarkedAtom(bondAtom)) {
							mMol.setAtomMarker(bondAtom, false);
							mAromaticAtoms--;
						}
					}
				}
			}
		}
	}



	/**
	 * Locate bonds that have no aromatic neighbour bond on at least one side.
	 * Convert these terminal or single aromatic bonds to a double bond
	 * and their aromatic neighbour bonds to a single bond.
	 * Do this until no terminal or single aromatic ring bonds remains.
 	 */
	private void promoteObviousBonds() {
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
						mMol.setAtomAbnormalValence(atom, -1);
						break;
						}
					}
				}
			}
		}

	private boolean assignKnownEZBondParities() {
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
							if (refAtom[i] == -1
							 && (mMol.getBondType(connBond) == Molecule.cBondTypeUp
							  || mMol.getBondType(connBond) == Molecule.cBondTypeDown)) {
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
					// if both bonds connected to the double bond atoms have the same slash direction, we have Z
					// (assuming that the parent atom (i.e. bondAtom[0]) in both cases is the double bond atom)
					boolean isZ = mMol.getBondType(refBond[0]) == mMol.getBondType(refBond[1]);

					// We need to correct, because slash or backslash refer to the double bonded
					// atom and not to the double bond itself as explained in opensmiles.org:
					//     F/C=C/F and C(\F)=C/F are both E
					// bondAtom[0] is always the parent in graph to bondAtom[1]. We use this to correct:
					for (int i=0; i<2; i++)
						if (refAtom[i] == mMol.getBondAtom(0, refBond[i]))
							isZ = !isZ;

					// E/Z configuration in the StereoMolecule refer to those neighbors with
					// lower atom index. Thus, we adapt for this:
					for (int i=0; i<2; i++)
						if (otherAtom[i] != -1
						 && otherAtom[i] < refAtom[i])
							isZ = !isZ;

					mMol.setBondParity(bond, isZ ? Molecule.cBondParityZor2
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

	private class EnumerationPosition implements Comparable<EnumerationPosition> {
		int mPosition,mCount;

		/**
		 * @param position position of first option in original smarts
		 */
		public EnumerationPosition(int position) {
			mPosition = position;
			mCount = 1;
			}

		public void increase() {
			mCount++;
			}

		public void enumerate(SmilesParser parser, byte[] smarts, ArrayList<String> enumeration) throws Exception {
			ArrayList<String> optionList = new ArrayList<>();

			int start = mPosition;
			SmilesAtomParser atomParser = new SmilesAtomParser(parser, mMode | mSmartsMode);
			int end = atomParser.parseAtomInsideBrackets(smarts, start+1, smarts.length, true, true)-1;
			if (smarts[end] != ']') {  // we have multiple options and create an option list
				optionList.add(new String(smarts, start, end-start));
				while (smarts[end] != ']') {
					start = end+1;
					end = atomParser.parseAtomInsideBrackets(smarts, start+1, smarts.length, true, true)-1;
					optionList.add(new String(smarts, start, end-start));
				}
			}

			for (String option : optionList)
				enumeration.add(new String(smarts, 0, mPosition) + option + new String(smarts, end, smarts.length-end));
			}

		@Override
		public int compareTo(EnumerationPosition o) {
			return Integer.compare(o.mPosition, mPosition);
		}
	}

	private class THParity {
		private static final int PSEUDO_ATOM_HYDROGEN = Integer.MAX_VALUE - 1;
		private static final int PSEUDO_ATOM_LONE_PAIR = Integer.MAX_VALUE;

		private class ParityNeighbour {
			int mAtom,mPosition;

			public ParityNeighbour(int atom, int position) {
				mAtom = atom;
				mPosition = position;
			}
			
			public String toString() {
				// makes debugging much easier!
				return "[" + (mAtom == PSEUDO_ATOM_HYDROGEN ? "h" 
						: mAtom == PSEUDO_ATOM_LONE_PAIR ? "lp" 
						: mMol.getAtomLabel(mAtom))
						 + mPosition + "]";
			}
		}


		int mCentralAtom,mCentralAtomPosition;
		boolean mIsClockwise,mError;
		ArrayList<ParityNeighbour> mNeighbourList;

		/**
		 * Instantiates a new parity object during smiles traversal.
		 * @param centralAtom index of atom processed
		 * @param centralAtomPosition position in SMILES of central atom
		 * @param fromAtom index of parent atom of centralAtom (-1 if centralAtom is first atom in smiles)
		 * @param explicitHydrogen Daylight syntax: hydrogen atoms defined within square bracket of other atom
		 * @param hydrogenPosition position in SMILES of central atom
		 * @param isClockwise true if central atom is marked with @@ rather than @
		 */
		public THParity(int centralAtom, int centralAtomPosition, int fromAtom, int explicitHydrogen, int hydrogenPosition, boolean isClockwise) {
			if (explicitHydrogen != 0 && explicitHydrogen != 1) {
				mError = true;
				}
			else {
				mCentralAtom = centralAtom;
				mCentralAtomPosition = centralAtomPosition;
				mIsClockwise = isClockwise;
				mNeighbourList = new ArrayList<>();

				// If we have a fromAtom, an explicit hydrogen, or a lone pair,
				// then add it as a normal neighbour.
				if (fromAtom != -1)
					addNeighbor(fromAtom, centralAtomPosition-1, false);

				if (fromAtom != -1 && explicitHydrogen == 1)
					addNeighbor(PSEUDO_ATOM_HYDROGEN, centralAtomPosition+1, false);
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
		 * @param atom
		 * @param position
		 */
		public void addNeighbor(int atom, int position, boolean unused) {
			if (!mError) {
				if (mNeighbourList.size() == 4) {
					mError = true;
					return;
				}

				mNeighbourList.add(new ParityNeighbour(atom, position));
			}
		}

		public int calculateParity(int[] handleHydrogenAtomMap) {
			if (mError)
				return Molecule.cAtomParityUnknown;

			// We need to translate smiles-parse-time atom indexes to those that the
			// molecule
			// uses after calling handleHydrogens, which is called from
			// ensureHelperArrays().
			for (ParityNeighbour neighbour : mNeighbourList)
				if (neighbour.mAtom != PSEUDO_ATOM_HYDROGEN && neighbour.mAtom != PSEUDO_ATOM_LONE_PAIR)
					neighbour.mAtom = handleHydrogenAtomMap[neighbour.mAtom];

			boolean isInverse = false;
			switch (mNeighbourList.size()) {
			case 2:
				isInverse = isInverseOrderAllene(handleHydrogenAtomMap);
				break;
			case 3:
				// All hydrogens atoms within SMILES all stereo centers all hydrogens must be
				// explicit (as explicit atoms or as H count in square brackets).
				// Therefore, three neighbour atoms is a rare situation, e.g. CC[S@](=O)C or
				// frozen out CC[N@H]C
				// In these cases we add the electron pair as pseudo neighbour
				mNeighbourList.add(new ParityNeighbour(PSEUDO_ATOM_LONE_PAIR, mCentralAtomPosition));
				// Note that this case also covers initial [C@H] as in alanine [C@H](N)(C)C(=O)O.
				// A lone pair will suffice in place of pseudo-hydrogen for these purposes.
				//$FALL-THROUGH$
			case 4:
				isInverse = isInverseOrderTH();
				break;
			default:
				return Molecule.cAtomParityUnknown;
			}
			//System.out.println(mNeighbourList);
			return (mIsClockwise ^ isInverse) ? Molecule.cAtomParity1 : Molecule.cAtomParity2;
			/*
			 * System.out.println();
			 * System.out.println("central:"+mCentralAtom+(mIsClockwise?" @@":" @")+" from:"
			 * +((mFromAtom == -1)?"none":Integer.toString(mFromAtom))+" with "
			 * +mImplicitHydrogen+" hydrogens");
			 * System.out.print("neighbors: "+mNeighborAtom[0]+"("+mNeighborPosition[0]+(
			 * mNeighborIsHydrogen[0]?",H":",non-H")+")"); for (int i=1; i<mNeighborCount;
			 * i++) System.out.print(", "+mNeighborAtom[i]+"("+mNeighborPosition[i]+(
			 * mNeighborIsHydrogen[i]?",H":",non-H")+")"); System.out.println();
			 * System.out.println("parity:"+parity);
			 */
		}

		/**
		 * Determine the atom index order parity for an allene. 
		 * 
		 * @return true if the SMILES atom ordering does not match the Molecule's atom ordering.
		 * @author hansonr
		 */
		private boolean isInverseOrderAllene(int[] smiles2MolMap) {
			boolean inversion = false;
			if (mol2SmilesMap == null) {
				mol2SmilesMap = new int[smiles2MolMap.length];
				for (int i = mol2SmilesMap.length; --i >= 0;) {
					mol2SmilesMap[smiles2MolMap[i]] = i; 
				}
			}
			inversion = checkAlleneInversion(smiles2MolMap, 0, inversion);
			inversion = checkAlleneInversion(smiles2MolMap, 1, inversion);
			return inversion;
		}


		/**
		 * Determine for a given end, whether there is a switching of the order of
		 * attached groups. The two important considerations we need to make are (a) if
		 * there are bracketed explicit H atoms that need to be considered and (b) if
		 * there are "ring" connection numbers. If no such complications exist, the
		 * order in the Molecule will be the same as in the SMILES.
		 * 
		 * @param index     0 (early end) or 1 (late end)
		 * @param inversion
		 * @return inversion after possible further inversion
		 * @author hansonr
		 */
		private boolean checkAlleneInversion(int[] smiles2MolMap, int index, boolean inversion) {
			final StereoMolecule mMol = SmilesParser.this.mMol; // avoiding multiple synthetic calls
			ParityNeighbour a1 = mNeighbourList.get(index);
			int smilesAtom = a1.mAtom;
			int centerAtom = smiles2MolMap[this.mCentralAtom];
			boolean hasImplicitH = mMol.getImplicitHydrogens(smilesAtom) > 0 && mMol.getConnAtoms(smilesAtom) == mMol.getAllConnAtoms(smilesAtom);
			boolean hasConnections = bsAtomHasConnections.get(mol2SmilesMap[smilesAtom]);
			//System.out.println(hasImplicitH + " " + hasConnections + " ca=" + mMol.getConnAtoms(smilesAtom) + " aca="
			//		+ mMol.getAllConnAtoms(smilesAtom) + " ih=" + mMol.getImplicitHydrogens(smilesAtom));

//			for (int i = 0, n = mMol.getAllConnAtoms(smilesAtom); i < n; i++) {
//				int b = mMol.getConnAtom(smilesAtom, i);
//				if (b == centerAtom)
//					continue;
//				System.out.println("m" + index + "_" + smilesAtom + "-" + b + " " + mMol.getAtomicNo(b));
//			}

			if (!hasImplicitH && !hasConnections) {
				// ignore FC or CF
				// but OC([H])= is special for some reason
				int b = mMol.getConnAtom(smilesAtom, 2); // H atom
				if (mMol.getAtomicNo(b) == 1) {
				//	System.out.println("check H " + new String(smiles));
					// BH
					// four tested cases:
					// [H]C(O)=[C@@]=CF  NO inversion!??
					// OC([H])=[C@]=CF   YES inversion!??
					// OC=[C@]=C([H])F   YES inversion OK
					// OC(F)=[C@]=C(O)[H] NO inversion OK
					//
					// Why we need to invert OC([H]) and not [H]C(O) makes no sense to me at all

					b = mol2SmilesMap[b];
					if (index == 0 ? b > smilesAtom 
							: mol2SmilesMap[mMol.getConnAtom(smilesAtom, 1)] > b) {
						inversion = !inversion;
					}
				}
				return inversion;
			}

			if (hasImplicitH) {
				// don't worry about leading F[CH]
				// just reverse trailing =[CH](F) (index==0), and [CH](F)= (index=1)
				// The H will not be in the molecule.
				// for index=0, it will be connected atom 0
				// for index=1, it will be connected atom 1
				int a2 = mMol.getConnAtom(smilesAtom, index);
				int b = mol2SmilesMap[a2];
				if (b > smilesAtom) {
					inversion = !inversion;
//					System.out.println("inverting CX");
					if (mMol.getAtomicNo(a2) == 1) {
						// this makes no sense to me
						// OC([H]) has a fluke? requires inversion, even though it should not.
						inversion = !inversion;
					}
				}
			}
			if (hasConnections) {
				// looking for C12 or C1(X) or [CH]1
				// connections will be earlier or later in the molecule.
				// the bonds will be to the connected atoms, in unusable order

				int n = mMol.getAllConnAtoms(smilesAtom);
				/**
				 * the order of the atoms is by position, but that is not significant for
				 * connections
				 */
				int[] connAtoms = new int[n - 1];
				/**
				 * the order of bonds is irrelevant
				 */
				int[] connBonds = new int[n - 1];
				int nConnected = 0;
				for (int i = 0, p = 0; i < n; i++) {
					int b = mMol.getConnAtom(smilesAtom, i);
					if (b == centerAtom)
						continue;
					connAtoms[p] = b;
					int bond = mMol.getBond(smilesAtom, b);
					int[] positions = mapConnectionBonds.get(mol2SmilesMap[smilesAtom] + "_" + mol2SmilesMap[b]);
					connBonds[p++] = (positions == null ? -1 : bond);
					if (positions != null)
						nConnected++;
				}
//				System.out.println("atoms " + Arrays.toString(connAtoms));
//				System.out.println("bonds " + Arrays.toString(connBonds));

				if (hasImplicitH) {
					// [CH]1
					if (connAtoms[0] < smilesAtom) {
						// earlier atom, but later than H here
						inversion = !inversion;
					}
				} else {
					// c1 or C12 or C1(X)
					// this gets tricky. Note that OpenSmiles
					// specifies only C1(X) is valid, not C(X)1
					// branched_atom ::= atom ringbond* branch*
					// but here the reader effectively treats these as C(X)1

					// So the question is, where is that atom?

					switch (nConnected) {
					case 1:
						if (mMol.getImplicitHydrogens(smilesAtom) > 0) {
							// C1= no problem; the connected atom is an H
							break;
						}
						// XC1 or C1(X), but I don't know if that is always the first connection
						boolean isFirst = (connBonds[0] >= 0);
						int connAtom = connAtoms[isFirst ? 0 : 1];
						int otherAtom = connAtoms[isFirst ? 1 : 0];
						if (otherAtom > smilesAtom) {
							// C1(X)
							if (connAtom > smilesAtom) {
								if (mMol.getAtomicNo(otherAtom) != 1) {
									// if the other atom is H, it will be already reversed
									// by this point.
									inversion = !inversion;
								}
							}
						} else {
							// XC1
							if (connAtom < smilesAtom)
								inversion = !inversion;
						}
						break;
					case 2:
						// trickier, because we have a lot of options with C12= or C21=
						// and the exact position of the numbers 1 and 2 are not known.
						int b1 = getOtherAtom(connBonds[0], smilesAtom);
						int b2 = getOtherAtom(connBonds[1], smilesAtom);
						// positions are the byte positions of the two ends of each connection
						int[] positions1 = mapConnectionBonds.get(smilesAtom + "_" + mol2SmilesMap[b1]);
						int[] positions2 = mapConnectionBonds.get(smilesAtom + "_" + mol2SmilesMap[b2]);
//						System.out.println(Arrays.toString(positions1));
//						System.out.println(Arrays.toString(positions2));
						int p1, p2;
						int c = mCentralAtomPosition;
						if (index == 0) {
							// first terminus C21=[C@]=.... with 2 and 1 pointing somewhere
							//
							// Two possibilities in each case.
							//
							// ....1....C12=[C@]=.....
							// ..p[0]..p[1]...........
							//
							// .........C12=[C@]=..1..
							// ........p[0]......p[1].
							//
							// If one end is late, then we know we want the other end,
							// and if both are early, we want the later one.
							p1 = (positions1[1] < c ? positions1[1] : positions1[0]);
							p2 = (positions2[1] < c ? positions2[1] : positions2[0]);
						} else {
							// second terminus
							// ....1...=[C@]=C12............
							// ..p[0].......p[1]............

							// ........=[C@]=C12........1...
							// .............p[0]......p[1]..

							p1 = (positions1[0] < c ? positions1[1] : positions1[0]);
							p2 = (positions2[0] < c ? positions2[1] : positions2[0]);
						}
						// inversion if this does not match the order of Molecule atoms
						if ((p1 < p2) != (b1 < b2)) {
							inversion = !inversion;
						}
						break;
					}
				}
			}
//			System.out.println("=" + index + " " + inversion + " " + new String(smiles));
			return inversion;
		}

		private int getOtherAtom(int b, int a) {
			return (mMol.mBondAtom[0][b] == a ? mMol.mBondAtom[1][b] : mMol.mBondAtom[0][b]);
		}

		private boolean isInverseOrderTH() {
			boolean inversion = false;
			for (int i=1; i<mNeighbourList.size(); i++) {
				for (int j=0; j<i; j++) {
					if (mNeighbourList.get(j).mAtom > mNeighbourList.get(i).mAtom)
						inversion = !inversion;
					if (mNeighbourList.get(j).mPosition > mNeighbourList.get(i).mPosition)
						inversion = !inversion;
				}
			}
			return inversion;
		}
	}

	private static void testStereo() {
		final String[][] data = { { "F/C=C/I", "F/C=C/I" },
								  { "F/C=C\\I", "F/C=C\\I" },
								  { "C(=C/I)/F", "F/C=C\\I" },
								  { "[H]C(/F)=C/I", "F/C=C\\I" },
								  { "C(=C\\1)/I.F1", "F/C=C/I" },
								  { "C(=C1)/I.F/1", "F/C=C/I" },
								  { "C(=C\\F)/1.I1", "F/C=C/I" },
								  { "C(=C\\F)1.I\\1", "F/C=C/I" },
								  { "C\\1=C/I.F1", "F/C=C/I" },
								  { "C1=C/I.F/1", "F/C=C/I" },
								  { "C(=C\\1)/2.F1.I2", "F/C=C/I" },
								  { "C/2=C\\1.F1.I2", "F/C=C/I" },
								  { "C/1=C/C=C/F.I1", "F/C=C/C=C\\I" },
								  { "C1=C/C=C/F.I\\1", "F/C=C/C=C\\I" },
								  { "C(/I)=C/C=C/1.F1", "F/C=C/C=C\\I" },
								  { "C(/I)=C/C=C1.F\\1", "F/C=C/C=C\\I" },

								  { "[C@](Cl)(F)(I)1.Br1", "F[C@](Cl)(Br)I" },
								  { "Br[C@](Cl)(I)1.F1", "F[C@](Cl)(Br)I" },
								  { "[C@H](F)(I)1.Br1", "F[C@H](Br)I" },
								  { "Br[C@@H](F)1.I1", "F[C@H](Br)I" },

								  { "C[S@@](CC)=O", "CC[S@](C)=O" },
								  { "[S@](=O)(C)CC", "CC[S@](C)=O" },

								  { "F1.OC=[C@]=C1", "OC=[C@]=CF" },
								  { "OC=[C@]=C1F.[H]1", "OC=[C@]=CF" },
								  { "[H]C(O)=[C@@]=CF", "OC=[C@]=CF" },
								  { "C(O)=[C@@]=CF", "OC=[C@]=CF" },
								  { "OC=[C@@]=C(F)[H]", "OC=[C@]=CF" },
								  { "CC(F)=[C@@]=CO", "CC(F)=[C@@]=CO" },
								  { "OC=[C@]=C(C)F", "CC(F)=[C@@]=CO" },
								  { "OC=[C@]=C(C)F", "CC(F)=[C@@]=CO" },
								  { "CC(F)=[C@@]=CO", "CC(F)=[C@@]=CO" },
								  { "CC(F)=[C@]=C(O)[H]", "CC(F)=[C@@]=CO" },
								  { "CC(F)=[C@]=C(O)Cl", "CC(F)=[C@]=C(O)Cl" },
								  { "ClC(O)=[C@]=C(F)C", "CC(F)=[C@]=C(O)Cl" },
								  { "OC(Cl)=[C@]=C(C)F", "CC(F)=[C@]=C(O)Cl" },
								  { "C1(Cl)=[C@]=C(C)F.O1", "CC(F)=[C@]=C(O)Cl" },
								  { "C(O)(Cl)=[C@]=C(C)F", "CC(F)=[C@]=C(O)Cl" },
								  { "[C@](=C(C)(F))=C(O)Cl", "CC(F)=[C@]=C(O)Cl" },
		};
		StereoMolecule mol = new StereoMolecule();
		for (String[] test:data) {
			try {
				new SmilesParser().parse(mol, test[0]);
				String smiles = new IsomericSmilesCreator(mol).getSmiles();
				System.out.print("IN:"+test[0]+" OUT:"+smiles);
				if (!test[1].equals(smiles))
					System.out.println(" EXPECTED: "+test[1]+" ERROR!");
				else
					System.out.println(" OK");
				}
			catch (Exception e) {
				e.printStackTrace();
				}
			}
		}

	public static void main(String[] args) {
		testStereo();

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
								  
								  { "OC=[C@]=CF", "allene-1", "gJQHBIAIVVb`@" },
								  { "OC([H])=[C@]=CF", "allene-1", "gJQHBIAIVVb`@" },
								  { "OC=[C@]=C([H])F", "allene-1", "gJQHBIAIVVb`@" },
								  
								  { "F1.OC=[C@]=C1", "allene-1", "gJQHBIAIVVb`@" },
								  { "OC=[C@]=C1F.[H]1", "allene-1", "gJQHBIAIVVb`@" },
								  { "[H]C(O)=[C@@]=CF", "allene-1", "gJQHBIAIVVb`@" },
								  { "C(O)=[C@@]=CF", "allene-1", "gJQHBIAIVVb`@" },
								  { "OC=[C@@]=C(F)[H]", "allene-1", "gJQHBIAIVVb`@" },
								  { "CC(F)=[C@@]=CO", "allene-2", "gGQHJIAIgfZJ@" },
								  { "OC=[C@]=C(C)F", "allene-2", "gGQHJIAIgfZJ@" },
								  { "OC=[C@]=C(C)F", "allene-2", "gGQHJIAIgfZJ@" },
								  { "CC(F)=[C@@]=CO", "allene-2", "gGQHJIAIgfZJ@" },
								  { "CC(F)=[C@]=C(O)[H]", "allene-2", "gGQHJIAIgfZJ@" },
								  { "CC(F)=[C@]=C(O)Cl", "allene-3", "gNqDDHbrBS[TmSH@" },
								  { "ClC(O)=[C@]=C(F)C", "allene-3", "gNqDDHbrBS[TmSH@" },
								  { "OC(Cl)=[C@]=C(C)F", "allene-3", "gNqDDHbrBS[TmSH@" },
								  { "C1(Cl)=[C@]=C(C)F.O1", "allene-3", "gNqDDHbrBS[TmSH@" },
								  { "C(O)(Cl)=[C@]=C(C)F", "allene-3", "gNqDDHbrBS[TmSH@" },
								  { "[C@](=C(C)(F))=C(O)Cl", "allene-3", "gNqDDHbrBS[TmSH@" },

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
