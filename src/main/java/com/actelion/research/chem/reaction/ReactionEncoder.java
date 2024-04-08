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

import com.actelion.research.chem.*;
import com.actelion.research.util.ArrayUtils;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

public class ReactionEncoder
{
	public static final char MOLECULE_DELIMITER = ' ';
	public static final char PRODUCT_IDENTIFIER = '!';
	public static final char CATALYST_DELIMITER = '+';	// character must not collide with idcode or coordinate encodings
	public static final char OBJECT_DELIMITER = '#';

	public static final String MOLECULE_DELIMITER_STRING = " ";
	public static final String OBJECT_DELIMITER_STRING = "#";

	public static final int INCLUDE_MAPPING = 1;
	public static final int INCLUDE_COORDS = 2;
	public static final int INCLUDE_DRAWING_OBJECTS = 4;
	public static final int INCLUDE_CATALYSTS = 8;

	public static final int INCLUDE_ALL = 15;
	public static final int INCLUDE_RXN_CODE_ONLY = 0;
	public static final int INCLUDE_DEFAULT = INCLUDE_MAPPING | INCLUDE_COORDS;

	public static final int RETAIN_REACTANT_AND_PRODUCT_ORDER = 16;


    private ReactionEncoder()
    {}

	/**
	 * Creates a String containing a canonical reaction code by
	 * creating idcodes of every reactant and product and
	 * concatenating them in lexically sorted order. This creates
	 * a canonical reaction code. The drawback is, however, that
	 * the original order of reactants and products may be changed.
	 * If mapping information is available this will be encoded
	 * in a 2nd string. Otherwise, this will be an empty string.
	 * Coordinates, if available, will be encoded in a 3rd string.
	 * If there are drawing objects assigned to this reaction
	 * then these are encoded in a 4th string.
	 * If the reaction contains catalysts, they are encoded as 5th string.
	 *
	 * @return String[5] with reaction code, mapping, coordinates, drawing objects, catalysts
	 */
	public static String[] encode(Reaction reaction, boolean keepAbsoluteCoordinates) {
		return encode(reaction, keepAbsoluteCoordinates, true);
	}

	/**
	 * Creates a canonical or non-canonical String containing a reaction
	 * code by creating idcodes of every reactant and product and
	 * concatenating them in original or canonical order.
	 * If mapping information is available this will be encoded
	 * in a 2nd string. Otherwise, this will be null.
	 * Coordinates, if available, will be encoded in a 3rd string.
	 * If there are drawing objects assigned to this reaction
	 * then these are encoded in a 4th string.
	 * If the reaction contains catalysts, they are encoded as 5th string.
	 *
	 * @param reaction
	 * @param keepAbsoluteCoordinates
	 * @param sortByIDCode whether to sort reactant and product idcodes to produce a canonical reaction code
	 * @return String[5] with reaction code, mapping, coordinates, drawing objects, catalysts
	 */
	public static String[] encode(Reaction reaction, boolean keepAbsoluteCoordinates, boolean sortByIDCode) {
		if (reaction == null
			|| reaction.getReactants() == 0
			|| reaction.getProducts() == 0) {
			return null;
		}

		String[] idcode = new String[reaction.getMolecules()];
		String[] mapping = new String[reaction.getMolecules()];
		String[] coords = new String[reaction.getMolecules()];

		for (int i = 0; i < reaction.getMolecules(); i++) {
			StereoMolecule mol = reaction.getMolecule(i);

			// reactants may not use cAtomQFRxnParityHint
			if (mol.isFragment() && i < reaction.getReactants())
				for (int atom=0; atom<mol.getAllAtoms(); atom++)
					mol.setAtomQueryFeature(atom, Molecule.cAtomQFRxnParityHint, false);

			Canonizer canonizer = new Canonizer(mol);
			idcode[i] = canonizer.getIDCode();
			if (idcode[i] == null) {
				return null;
			}

			mapping[i] = canonizer.getEncodedMapping();
			coords[i] = canonizer.getEncodedCoordinates(keepAbsoluteCoordinates);
		}

		StringBuilder idcodeSequence = new StringBuilder();
		StringBuilder coordsSequence = new StringBuilder();
		StringBuilder mappingSequence = new StringBuilder();

		for (int i = 0; i < reaction.getReactants(); i++) {
			int index = i;
			if (sortByIDCode) {
				String maxString = "";
				index = -1;
				for (int j = 0; j < reaction.getReactants(); j++) {
					if (maxString.compareTo(idcode[j]) < 0) {
						maxString = idcode[j];
						index = j;
					}
				}
			}
			if (i > 0) {
				idcodeSequence.append(MOLECULE_DELIMITER);
				mappingSequence.append(MOLECULE_DELIMITER);
				coordsSequence.append(MOLECULE_DELIMITER);
			}
			idcodeSequence.append(idcode[index]);
			mappingSequence.append(mapping[index]);
			coordsSequence.append(coords[index]);
			idcode[index] = "";
		}

		idcodeSequence.append(PRODUCT_IDENTIFIER);
		mappingSequence.append(MOLECULE_DELIMITER);
		coordsSequence.append(MOLECULE_DELIMITER);

		for (int i = reaction.getReactants(); i < reaction.getMolecules(); i++) {
			int index = i;
			if (sortByIDCode) {
				String maxString = "";
				index = -1;
				for (int j = reaction.getReactants(); j < reaction.getMolecules(); j++) {
					if (maxString.compareTo(idcode[j]) < 0) {
						maxString = idcode[j];
						index = j;
					}
				}
			}
			if (i > reaction.getReactants()) {
				idcodeSequence.append(MOLECULE_DELIMITER);
				mappingSequence.append(MOLECULE_DELIMITER);
				coordsSequence.append(MOLECULE_DELIMITER);
			}
			idcodeSequence.append(idcode[index]);
			mappingSequence.append(mapping[index]);
			coordsSequence.append(coords[index]);
			idcode[index] = "";
		}

		String[] result = new String[5];
		result[0] = idcodeSequence.toString();
		if (mappingSequence.length() > reaction.getMolecules() - 1)   // delimiters only
		{
			result[1] = mappingSequence.toString();
		}
		if (coordsSequence.length() > reaction.getMolecules() - 1)   // delimiters only
		{
			result[2] = coordsSequence.toString();
		}
		if (reaction.getDrawingObjects() != null) {
			result[3] = reaction.getDrawingObjects().toString();
		}
		if (reaction.getCatalysts() != 0) {
			result[4] = encodeCatalysts(reaction, keepAbsoluteCoordinates);
		}

		return result;
	}

	private static String encodeCatalysts(Reaction reaction, boolean keepAbsoluteCoordinates) {
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<reaction.getCatalysts(); i++) {
			if (sb.length() != 0)
				sb.append(CATALYST_DELIMITER);
			Canonizer canonizer = new Canonizer(reaction.getCatalyst(i));
			sb.append(canonizer.getIDCode());
			if (keepAbsoluteCoordinates) {
				sb.append(" ");
				sb.append(canonizer.getEncodedCoordinates(true));
			}
		}
		return sb.toString();
	}

	/**
	 * Creates a String containing a reaction code by creating idcodes of every reactant and product and
	 * concatenating them in original (if mode includes RETAIN_REACTANT_AND_PRODUCT_ORDER) or in
	 * lexical order. In the latter case this string is a canonical reaction encoding.
	 * If mapping information is available this will be encoded in a 2nd string.
	 * Coordinates, if available, will be encoded in a 3rd string.
	 * If there are drawing objects assigned to this reaction then these are encoded in a 4th string.
	 *
	 * @return One String with reaction code, coordinates, mapping, drawing objects as defined by mode.
	 */
	public static String encode(Reaction reaction, boolean keepAbsoluteCoordinates, int mode) {
		String[] result = encode(reaction, keepAbsoluteCoordinates, (mode & RETAIN_REACTANT_AND_PRODUCT_ORDER) == 0);
		if (result == null) {
			return null;
		}

		StringBuffer buf = new StringBuffer(result[0]);

		if (mode != 0) {
			buf.append(OBJECT_DELIMITER);
			if ((mode & INCLUDE_MAPPING) != 0
				&& result.length > 1
				&& result[1] != null) {
				buf.append(result[1]);
			}
		}

		mode &= ~INCLUDE_MAPPING;
		if (mode != 0) {
			buf.append(OBJECT_DELIMITER);
			if ((mode & INCLUDE_COORDS) != 0
				&& result.length > 2
				&& result[2] != null) {
				buf.append(result[2]);
			}
		}

		mode &= ~INCLUDE_COORDS;
		if (mode != 0) {
			buf.append(OBJECT_DELIMITER);
			if ((mode & INCLUDE_DRAWING_OBJECTS) != 0
				&& result.length > 3
				&& result[3] != null) {
				buf.append(result[3]);
			}
		}

		mode &= ~INCLUDE_DRAWING_OBJECTS;
		if (mode != 0) {
			buf.append(OBJECT_DELIMITER);
			if ((mode & INCLUDE_CATALYSTS) != 0
				&& result.length > 4
				&& result[4] != null) {
				buf.append(result[4]);
			}
		}

		return buf.toString();
	}

	/**
	 * Creates a Reaction object by interpreting a reaction code,
	 * mapping, coordinates and drawing objects that were earlier created
	 * by this class.
	 * If rxnCoords are relative or null, and if ensureCoordinates==true
	 * then all reactants and products are placed automatically along a
	 * horizontal line.
	 *
	 * @return Reaction
	 */
	public static Reaction decode(String rxnCode, String rxnMapping, String rxnCoords,
								  String rxnObjects, String rxnCatalysts, boolean ensureCoordinates, Reaction rxn) {
		if (rxnCode == null || rxnCode.length() == 0) {
			return null;
		}

		boolean isProduct = false;
		int idcodeIndex = 0;
		int mappingIndex = 0;
		int coordsIndex = 0;

		int productIndex = rxnCode.indexOf(PRODUCT_IDENTIFIER);
		if (productIndex == -1) {
			return null;
		}

		if (rxn == null)
			rxn = new Reaction();
		else
			rxn.clear();

		while (idcodeIndex != -1) {
			if (idcodeIndex > productIndex) {
				isProduct = true;
			}

			int delimiterIndex = rxnCode.indexOf(MOLECULE_DELIMITER, idcodeIndex);
			if (!isProduct
				&& (delimiterIndex > productIndex || delimiterIndex == -1)) {
				delimiterIndex = productIndex;
			}

			String idcode = null;
			if (delimiterIndex == -1) {
				idcode = rxnCode.substring(idcodeIndex);
				idcodeIndex = -1;
			} else {
				idcode = rxnCode.substring(idcodeIndex, delimiterIndex);
				idcodeIndex = delimiterIndex + 1;
			}

			String mapping = null;
			if (rxnMapping != null && rxnMapping.length() != 0) {
				delimiterIndex = rxnMapping.indexOf(MOLECULE_DELIMITER, mappingIndex);
				if (delimiterIndex == -1) {
					mapping = rxnMapping.substring(mappingIndex);
				} else {
					mapping = rxnMapping.substring(mappingIndex, delimiterIndex);
					mappingIndex = delimiterIndex + 1;
				}
			}

			String coords = null;
			if (rxnCoords != null && rxnCoords.length() != 0) {
				delimiterIndex = rxnCoords.indexOf(MOLECULE_DELIMITER, coordsIndex);
				if (delimiterIndex == -1) {
					coords = rxnCoords.substring(coordsIndex);
				} else {
					coords = rxnCoords.substring(coordsIndex, delimiterIndex);
					coordsIndex = delimiterIndex + 1;
				}
			}

			IDCodeParser parser = new IDCodeParser(ensureCoordinates);
			StereoMolecule mol = parser.getCompactMolecule(idcode, coords);

			if (mapping != null) {
				parser.parseMapping(mapping.getBytes(StandardCharsets.UTF_8));
			}

			if (isProduct) {
				rxn.addProduct(mol);
			} else {
				rxn.addReactant(mol);
			}
		}

		if (rxnObjects != null && rxnObjects.length() != 0) {
			rxn.setDrawingObjects(new DrawingObjectList(rxnObjects));
		}

		if (rxnCatalysts != null && rxnCatalysts.length() != 0) {
			IDCodeParser parser = new IDCodeParser(ensureCoordinates);
			int index1 = 0;
			int index2 = rxnCatalysts.indexOf(CATALYST_DELIMITER);
			while (index2 != -1) {
				rxn.addCatalyst(parser.getCompactMolecule(rxnCatalysts.substring(index1, index2)));
				index1 = index2+1;
				index2 = rxnCatalysts.indexOf(CATALYST_DELIMITER, index1);
			}
			rxn.addCatalyst(parser.getCompactMolecule(rxnCatalysts.substring(index1)));
		}

		return rxn;
	}

	/**
	 * Creates a Reaction object by interpreting a reaction code,
	 * mapping, coordinates and drawing objects that were earlier created
	 * by this class.
	 * If rxnCoords are relative or null, and if ensureCoordinates==true
	 * then all reactants and products are placed automatically along a
	 * horizontal line.
	 *
	 * @return Reaction
	 */
	public static Reaction decode(byte[] rxnCode, byte[] rxnMapping, byte[] rxnCoords,
								  String rxnObjects, byte[] rxnCatalysts, boolean ensureCoordinates) {
		if (rxnCode == null || rxnCode.length == 0) {
			return null;
		}

		boolean isProduct = false;
		int idcodeIndex = 0;
		int mappingIndex = 0;
		int coordsIndex = 0;

		int productIndex = indexOf(rxnCode, PRODUCT_IDENTIFIER);
		if (productIndex == -1)
			return null;

		Reaction rxn = new Reaction();

		while (idcodeIndex != -1) {
			if (idcodeIndex > productIndex)
				isProduct = true;

			int delimiterIndex = indexOf(rxnCode, MOLECULE_DELIMITER, idcodeIndex);
			if (!isProduct && (delimiterIndex > productIndex || delimiterIndex == -1))
				delimiterIndex = productIndex;

			int idcodeStart = idcodeIndex;
			idcodeIndex = (delimiterIndex == -1) ? -1 : delimiterIndex + 1;

			int mappingStart = -1;
			if (rxnMapping != null && mappingIndex < rxnMapping.length) {
				mappingStart = (rxnMapping[mappingIndex] == MOLECULE_DELIMITER) ? -1 : mappingIndex;
				delimiterIndex = indexOf(rxnMapping, MOLECULE_DELIMITER, mappingIndex);
				if (delimiterIndex != -1)
					mappingIndex = delimiterIndex + 1;
			}

			int coordsStart = -1;
			if (rxnCoords != null && rxnCoords.length != 0) {
				coordsStart = coordsIndex;
				delimiterIndex = indexOf(rxnCoords, MOLECULE_DELIMITER, coordsIndex);
				if (delimiterIndex != -1)
					coordsIndex = delimiterIndex + 1;
			}

			IDCodeParser parser = new IDCodeParser(ensureCoordinates);
			parser.neglectSpaceDelimitedCoordinates();
			StereoMolecule mol = parser.getCompactMolecule(rxnCode, rxnCoords, idcodeStart, coordsStart);

			if (mappingStart != -1)
				parser.parseMapping(rxnMapping, mappingStart);

			if (isProduct)
				rxn.addProduct(mol);
			else
				rxn.addReactant(mol);
		}

		if (rxnObjects != null && rxnObjects.length() != 0) {
			rxn.setDrawingObjects(new DrawingObjectList(rxnObjects));
		}

		if (rxnCatalysts != null && rxnCatalysts.length != 0) {
			IDCodeParser parser = new IDCodeParser(ensureCoordinates);
			int index1 = 0;
			int index2 = indexOf(rxnCatalysts, CATALYST_DELIMITER);
			while (index2 != -1) {
				rxn.addCatalyst(parser.getCompactMolecule(rxnCatalysts, index1));
				index1 = index2+1;
				index2 = indexOf(rxnCatalysts, CATALYST_DELIMITER, index1);
			}
			rxn.addCatalyst(parser.getCompactMolecule(rxnCatalysts, index1));
		}

		return rxn;
	}

	private static int indexOf(byte[] bytes, char ch) {
		for (int i=0; i<bytes.length; i++)
			if (bytes[i] == ch)
				return i;

		return -1;
		}

	private static int indexOf(byte[] bytes, char ch, int start) {
		for (int i=start; i<bytes.length; i++)
			if (bytes[i] == ch)
				return i;

		return -1;
	}

	public static Reaction decode(String s, boolean ensureCoordinates) {
		return decode(s, ensureCoordinates, null);
		}

	/**
	 * Creates a Reaction object by interpreting a reaction code,
	 * mapping, coordinates and drawing objects that were earlier created
	 * by this class and are passed OBJECT_DELIMITER-delimited within
	 * one string.
	 * If rxnCoords are relative or null, and if ensureCoordinates==true
	 * then all reactants and products are placed automatically along a
	 * horizontal line.
	 *
	 * @return Reaction
	 */
	public static Reaction decode(String s, boolean ensureCoordinates, Reaction rxn) {
		if (s == null)
			return null;

		String rxnCode = s;
		String rxnMapping = null;
		String rxnCoords = null;
		String rxnObjects = null;
		String rxnCatalysts = null;
		int index1 = s.indexOf(OBJECT_DELIMITER);
		if (index1 == -1) {
			rxnCode = s;
		} else {
			rxnCode = s.substring(0, index1);
			int index2 = s.indexOf(OBJECT_DELIMITER, index1 + 1);
			if (index2 == -1) {
				rxnMapping = s.substring(index1 + 1);
			} else {
				rxnMapping = s.substring(index1 + 1, index2);
				int index3 = s.indexOf(OBJECT_DELIMITER, index2 + 1);
				if (index3 == -1) {
					rxnCoords = s.substring(index2 + 1);
				} else {
					rxnCoords = s.substring(index2 + 1, index3);
					int index4 = s.indexOf(OBJECT_DELIMITER, index3 + 1);
					if (index4 == -1) {
						rxnObjects = s.substring(index3 + 1);
					} else {
						rxnObjects = s.substring(index3 + 1, index4);
						rxnCatalysts = s.substring(index4 + 1);
					}
				}
			}
		}

		return decode(rxnCode, rxnMapping, rxnCoords, rxnObjects, rxnCatalysts, ensureCoordinates, rxn);
	}


	/**
	 * Creates a Reaction object by interpreting a reaction string encoded by this class.
	 * Include options define whether mapping, coordinates, catalysts, and drawing objects
	 # are included in the reaction object.
	 * @param s
	 * @param includeOptions
	 * @return Reaction
	 */
	public static Reaction decode(String s, int includeOptions, Reaction rxn) {
		if (s == null)
			return null;

		String rxnCode = s;
		String rxnMapping = null;
		String rxnCoords = null;
		String rxnObjects = null;
		String rxnCatalysts = null;
		int index1 = s.indexOf(OBJECT_DELIMITER);
		if (index1 == -1) {
			rxnCode = s;
		} else {
			rxnCode = s.substring(0, index1);
			int index2 = s.indexOf(OBJECT_DELIMITER, index1 + 1);
			if (index2 == -1) {
				rxnMapping = s.substring(index1 + 1);
			} else {
				rxnMapping = s.substring(index1 + 1, index2);
				int index3 = s.indexOf(OBJECT_DELIMITER, index2 + 1);
				if (index3 == -1) {
					rxnCoords = s.substring(index2 + 1);
				} else {
					rxnCoords = s.substring(index2 + 1, index3);
					int index4 = s.indexOf(OBJECT_DELIMITER, index3 + 1);
					if (index4 == -1) {
						rxnObjects = s.substring(index3 + 1);
					} else {
						rxnObjects = s.substring(index3 + 1, index4);
						rxnCatalysts = s.substring(index4 + 1);
					}
				}
			}
		}

		return decode(rxnCode,
			(includeOptions & INCLUDE_MAPPING) != 0 ? rxnMapping : null,
			(includeOptions & INCLUDE_COORDS) != 0 ? rxnCoords : null,
			(includeOptions & INCLUDE_DRAWING_OBJECTS) != 0 ? rxnObjects : null,
			(includeOptions & INCLUDE_CATALYSTS) != 0 ? rxnCatalysts : null,
			false, rxn);
	}

	/**
	 * Generates an array of all reactants and/or products of the encoded reaction string as bytes.
	 * If the string includes atom coordinates or if they are explicitly, these are used.
	 * At least one of includeReactants and includeProducts must be true.
	 * @param s encoded reaction
	 * @param includeReactants
	 * @param includeProducts
	 * @param includeMapping
	 * @return null (if reactants or products are missing) or StereoMolecule array with at least one molecule
	 */
	public static StereoMolecule[] decodeMolecules(String s, boolean includeCoords, boolean includeMapping, boolean includeReactants, boolean includeProducts) {
		if (s == null)
			return null;

		byte[] rxnCode = null;
		byte[] rxnMapping = null;
		byte[] rxnCoords = null;
		int index1 = s.indexOf(OBJECT_DELIMITER);
		if (index1 == -1) {
			rxnCode = s.getBytes(StandardCharsets.UTF_8);
		} else {
			rxnCode = s.substring(0, index1).getBytes(StandardCharsets.UTF_8);
			if (includeMapping || includeCoords) {
				int index2 = s.indexOf(OBJECT_DELIMITER, index1 + 1);
				if (index2 == -1) {
					if (includeMapping)
						rxnMapping = s.substring(index1 + 1).getBytes(StandardCharsets.UTF_8);
				} else {
					if (includeMapping)
						rxnMapping = s.substring(index1 + 1, index2).getBytes(StandardCharsets.UTF_8);
					if (includeCoords) {
						int index3 = s.indexOf(OBJECT_DELIMITER, index2 + 1);
						if (index3 == -1) {
							rxnCoords = s.substring(index2 + 1).getBytes(StandardCharsets.UTF_8);
						} else {
							rxnCoords = s.substring(index2 + 1, index3).getBytes(StandardCharsets.UTF_8);
							}
						}
					}
				}
			}

		return decodeMolecules(rxnCode, rxnCoords, rxnMapping, includeReactants, includeProducts);
		}

		/**
		 * Generates an array of all reactants and/or products of the encoded reaction string as bytes.
		 * If the string includes atom coordinates or if they are explicitly, these are used.
		 * At least one of includeReactants and includeProducts must be true.
		 * @param rxnBytes may contain atom coordinates
		 * @param coords may be null
		 * @param mapping may be null
		 * @param includeReactants
		 * @param includeProducts
		 * @return null (if reactants or products are missing) or StereoMolecule array with at least one molecule
		 */
	public static StereoMolecule[] decodeMolecules(byte[] rxnBytes, byte[] coords, byte[] mapping, boolean includeReactants, boolean includeProducts) {
		if (rxnBytes == null || rxnBytes.length == 0)
			return null;

		int reactantEnd = ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.PRODUCT_IDENTIFIER);
		if (reactantEnd <= 0)
			return null;

		int productIndex = reactantEnd + 1;
		int productEnd = ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.OBJECT_DELIMITER, productIndex);
		if (productEnd == -1)
			productEnd = rxnBytes.length;

		if (productIndex == productEnd)
			return null;

		int coordsIndex = 0;
		if (coords == null) {
			coordsIndex = 1+ArrayUtils.indexOf(rxnBytes, (byte) ReactionEncoder.OBJECT_DELIMITER, productEnd + 1);
			if (coordsIndex != 0)
				coords = rxnBytes;
		}
		if (coords != null) {
			if (!includeReactants) {
				int reactantCount = 1;
				for (int i=0; i<reactantEnd; i++)
					if (rxnBytes[i] == ReactionEncoder.MOLECULE_DELIMITER)
						reactantCount++;
				for (int i=0; reactantCount != 0 && i<coords.length; i++) {
					if (coords[i] == ReactionEncoder.MOLECULE_DELIMITER) {
						coordsIndex = i + 1;
						reactantCount--;
					}
				}
			}
		}

		int mappingIndex = 0;
		if (mapping != null) {
			if (!includeReactants) {
				int reactantCount = 1;
				for (int i=0; i<reactantEnd; i++)
					if (rxnBytes[i] == ReactionEncoder.MOLECULE_DELIMITER)
						reactantCount++;
				for (int i=0; reactantCount != 0 && i<mapping.length; i++) {
					if (mapping[i] == ReactionEncoder.MOLECULE_DELIMITER) {
						mappingIndex = i + 1;
						reactantCount--;
					}
				}
			}
		}

		ArrayList<StereoMolecule> moleculeList = new ArrayList<>();
		if (includeReactants) {
			int reactantIndex = 0;
			do {
				IDCodeParser parser = new IDCodeParser();
				parser.neglectSpaceDelimitedCoordinates();
				StereoMolecule reactant = parser.getCompactMolecule(rxnBytes, coords, reactantIndex, coordsIndex);
				if (reactant.getAllAtoms() != 0)
					moleculeList.add(reactant);

				reactantIndex = 1+ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.MOLECULE_DELIMITER, reactantIndex);
				if (coords != null)
					coordsIndex = 1+ArrayUtils.indexOf(coords, (byte)ReactionEncoder.MOLECULE_DELIMITER, coordsIndex);
				if (mapping != null) {
					parser.parseMapping(mapping, mappingIndex);
					mappingIndex = 1+ArrayUtils.indexOf(mapping, (byte)ReactionEncoder.MOLECULE_DELIMITER, mappingIndex);
				}
			} while (reactantIndex != 0 && reactantIndex < reactantEnd);
		}
		if (includeProducts) {
			do {
				IDCodeParser parser = new IDCodeParser();
				parser.neglectSpaceDelimitedCoordinates();
				StereoMolecule product = parser.getCompactMolecule(rxnBytes, coords, productIndex, coordsIndex);
				if (product.getAllAtoms() != 0)
					moleculeList.add(product);

				productIndex = 1+ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.MOLECULE_DELIMITER, productIndex);
				if (coords != null)
					coordsIndex = 1+ArrayUtils.indexOf(coords, (byte)ReactionEncoder.MOLECULE_DELIMITER, coordsIndex);
				if (mapping != null) {
					parser.parseMapping(mapping, mappingIndex);
					mappingIndex = 1+ArrayUtils.indexOf(mapping, (byte)ReactionEncoder.MOLECULE_DELIMITER, mappingIndex);
				}
			} while (productIndex != 0 && productIndex < productEnd);
		}

		return moleculeList.size() == 0 ? null : moleculeList.toArray(new StereoMolecule[0]);
	}

	/**
	 * Generates an array of all catalysts of the encoded reaction string as bytes.
	 * If the string includes atom coordinates, these are used.
	 * @param rxnBytes
	 * @return null or StereoMolecule array with at least one molecule
	 */
	public static StereoMolecule[] decodeCatalysts(byte[] rxnBytes) {
		if (rxnBytes == null || rxnBytes.length == 0)
			return null;

		int index = 0;
		for (int i=0; i<4; i++) {
			index = 1 + ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.OBJECT_DELIMITER, index);
			if (index == 0)
				return null;
		}

		if (index == rxnBytes.length)
			return null;

		ArrayList<StereoMolecule> catalystList = new ArrayList<StereoMolecule>();
		while (index != 0 && index < rxnBytes.length) {
			int nextIndex = 1+ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.CATALYST_DELIMITER, index);
			int coordsIndex = 1+ArrayUtils.indexOf(rxnBytes, (byte)' ', index);

			StereoMolecule catalyst = (coordsIndex != 0 && (nextIndex == 0 || nextIndex > coordsIndex)) ?
					  new IDCodeParser().getCompactMolecule(rxnBytes, rxnBytes, index, coordsIndex)
					: new IDCodeParser().getCompactMolecule(rxnBytes, null, index, -1);
			if (catalyst.getAllAtoms() != 0)
				catalystList.add(catalyst);

			index = nextIndex;
		}

		return catalystList.size() == 0 ? null : catalystList.toArray(new StereoMolecule[0]);
	}

	/**
	 * Generates an array of all reactants and/or products of the encoded reaction string as bytes.
	 * If the string includes atom coordinates or if they are explicitly, these are used.
	 * At least one of includeReactants and includeProducts must be true.
	 * @param rxnBytes may contain atom coordinates
	 * @return null (if reactants or products are missing) or StereoMolecule array with at least one molecule
	 */
	public static byte[][] getMoleculeIDCodes(byte[] rxnBytes, boolean includeReactants, boolean includeProducts) {
		if (rxnBytes == null || rxnBytes.length == 0)
			return null;

		int reactantEnd = ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.PRODUCT_IDENTIFIER);
		if (reactantEnd <= 0)
			return null;

		int productIndex = reactantEnd + 1;
		int productEnd = ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.OBJECT_DELIMITER, productIndex);
		if (productEnd == -1)
			productEnd = rxnBytes.length;

		if (productIndex == productEnd)
			return null;

		ArrayList<byte[]> moleculeList = new ArrayList<>();
		if (includeReactants) {
			int reactantIndex = 0;
			while (reactantIndex < reactantEnd) {
				int index2 = ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.MOLECULE_DELIMITER, reactantIndex);
				if (index2 == -1)
					index2 = reactantEnd;
				moleculeList.add(Arrays.copyOfRange(rxnBytes, reactantIndex, index2));
				reactantIndex = 1 + index2;
			}
		}
		if (includeProducts) {
			while (productIndex < productEnd) {
				int index2 = ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.MOLECULE_DELIMITER, productIndex);
				if (index2 == -1)
					index2 = productEnd;
				moleculeList.add(Arrays.copyOfRange(rxnBytes, productIndex, index2));
				productIndex = 1 + index2;
			}
		}

		return moleculeList.size() == 0 ? null : moleculeList.toArray(new byte[0][]);
	}
}
