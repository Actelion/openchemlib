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

package com.actelion.research.chem.reaction;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.DrawingObjectList;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.ArrayUtils;

import java.util.ArrayList;

public class ReactionEncoder
{
	public static final char MOLECULE_DELIMITER = ' ';
	public static final char PRODUCT_IDENTIFIER = '!';
	public static final char CATALYST_DELIMITER = '!';
	public static final char OBJECT_DELIMITER = '#';

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
	 * Creates a String containing a unique reaction code by
	 * creating idcodes of every reactant and product and
	 * concatenating them in lexical order.
	 * If mapping information is available this will be encoded
	 * in a 2nd string. Otherwise this will be null.
	 * Coordinates, if available, will be encoded in a 3rd string.
	 * If there are drawing objects assigned to this reaction
	 * then these are encoded in a 4th string.
	 *
	 * @return String[4] with reaction code, coordinates, mapping, drawing objects
	 */
	public static String[] encode(Reaction reaction, boolean keepAbsoluteCoordinates) {
		return encode(reaction, keepAbsoluteCoordinates, true);
	}

	/**
	 * Creates a non-unique String containing a reaction code by
	 * creating idcodes of every reactant and product and
	 * concatenating them in original order.
	 * If mapping information is available this will be encoded
	 * in a 2nd string. Otherwise this will be null.
	 * Coordinates, if available, will be encoded in a 3rd string.
	 * If there are drawing objects assigned to this reaction
	 * then these are encoded in a 4th string.
	 * If the reaction contains catalysts, they are encoded as 5th string.
	 *
	 * @param reaction
	 * @param keepAbsoluteCoordinates
	 * @param sortByIDCode
	 * @return String[5] with reaction code, coordinates, mapping, drawing objects, catalysts
	 */
	private static String[] encode(Reaction reaction, boolean keepAbsoluteCoordinates, boolean sortByIDCode) {
		if (reaction == null
			|| reaction.getReactants() == 0
			|| reaction.getProducts() == 0) {
			return null;
		}

		String[] idcode = new String[reaction.getMolecules()];
		String[] mapping = new String[reaction.getMolecules()];
		String[] coords = new String[reaction.getMolecules()];

		for (int i = 0; i < reaction.getMolecules(); i++) {
			Canonizer canonizer = new Canonizer(reaction.getMolecule(i));
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
								  String rxnObjects, String rxnCatalysts, boolean ensureCoordinates) {
		if (rxnCode == null || rxnCode.length() == 0) {
			return null;
		}

		boolean isProduct = false;
		int idcodeIndex = 0;
		int mappingIndex = 0;
		int coordsIndex = 0;
		boolean reactionLayoutRequired = false;

		int productIndex = rxnCode.indexOf(PRODUCT_IDENTIFIER);
		if (productIndex == -1) {
			return null;
		}

		Reaction rxn = new Reaction();
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

			if (!reactionLayoutRequired && (coords == null || !parser.coordinatesAreAbsolute(coords)))
				reactionLayoutRequired = true;

			if (mapping != null) {
				parser.parseMapping(mapping.getBytes());
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

		rxn.setReactionLayoutRequired(reactionLayoutRequired);

		return rxn;
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
	public static Reaction decode(String s, boolean ensureCoordinates) {
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

		return decode(rxnCode, rxnMapping, rxnCoords, rxnObjects, rxnCatalysts, ensureCoordinates);
	}


	/**
	 * Creates a Reaction object by interpreting a reaction string encoded by this class.
	 * Include options define whether mapping, coordinates, catalysts, and drawing objects
	 # are included in the reaction object.
	 * @param s
	 * @param includeOptions
	 * @return Reaction
	 */
	public static Reaction decode(String s, int includeOptions) {
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
			false);
	}

	/**
	 * Generates an array of all products of the encoded reaction string as bytes.
	 * If the string includes atom coordinates, these are used.
	 * @param rxnBytes
	 * @return null or StereoMolecule array with at least one molecule
	 */
	public static StereoMolecule[] decodeProducts(byte[] rxnBytes) {
		if (rxnBytes == null || rxnBytes.length == 0)
			return null;

		int productIndex = 1+ ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.PRODUCT_IDENTIFIER);
		if (productIndex == 0)
			return null;

		int productEnd = ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.OBJECT_DELIMITER, productIndex);
		if (productEnd == -1)
			productEnd = rxnBytes.length;

		if (productIndex == productEnd)
			return null;

		byte[] coords = null;
		int coordsIndex = 1+ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.OBJECT_DELIMITER, productEnd+1);
		if (coordsIndex != 0) {
			int reactantIndex = 0;
			while (reactantIndex < productIndex) {	// advance coordinate index one step for every reactant
				reactantIndex = 1+ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.MOLECULE_DELIMITER, reactantIndex);
				coordsIndex = 1+ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.MOLECULE_DELIMITER, coordsIndex);
			}
			coords = rxnBytes;
		}


		ArrayList<StereoMolecule> productList = new ArrayList<StereoMolecule>();
		while (productIndex != -1 && productIndex < productEnd) {
			StereoMolecule product = new IDCodeParser().getCompactMolecule(rxnBytes, coords, productIndex, coordsIndex);
			if (product.getAllAtoms() != 0)
				productList.add(product);

			productIndex = 1+ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.MOLECULE_DELIMITER, productIndex);
			coordsIndex = 1+ArrayUtils.indexOf(rxnBytes, (byte)ReactionEncoder.MOLECULE_DELIMITER, coordsIndex);
		}

		return productList.size() == 0 ? null : productList.toArray(new StereoMolecule[0]);
	}
}
