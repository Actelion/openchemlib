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

package com.actelion.research.chem.reaction;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.DrawingObjectList;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;

public class ReactionEncoder
{
	public static final char MOLECULE_DELIMITER = ' ';
	public static final char PRODUCT_IDENTIFIER = '!';
	public static final char OBJECT_DELIMITER = '#';

	public static final int INCLUDE_MAPPING = 1;
	public static final int INCLUDE_COORDS = 2;
	public static final int INCLUDE_DRAWING_OBJECTS = 4;
	public static final int RETURN_RXN_CODE_ONLY = 0;
	public static final int RETURN_DEFAULT = INCLUDE_MAPPING | INCLUDE_COORDS;


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
	public static String[] encode(Reaction reaction, boolean keepAbsoluteCoordinates)
	{
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

		StringBuffer idcodeSequence = new StringBuffer();
		StringBuffer coordsSequence = new StringBuffer();
		StringBuffer mappingSequence = new StringBuffer();

		for (int i = 0; i < reaction.getReactants(); i++) {
			String maxString = "";
			int maxIndex = -1;
			for (int j = 0; j < reaction.getReactants(); j++) {
				if (maxString.compareTo(idcode[j]) < 0) {
					maxString = idcode[j];
					maxIndex = j;
				}
			}
			if (i > 0) {
				idcodeSequence.append(MOLECULE_DELIMITER);
				mappingSequence.append(MOLECULE_DELIMITER);
				coordsSequence.append(MOLECULE_DELIMITER);
			}
			idcodeSequence.append(idcode[maxIndex]);
			mappingSequence.append(mapping[maxIndex]);
			coordsSequence.append(coords[maxIndex]);
			idcode[maxIndex] = "";
		}

		idcodeSequence.append(PRODUCT_IDENTIFIER);
		mappingSequence.append(MOLECULE_DELIMITER);
		coordsSequence.append(MOLECULE_DELIMITER);

		for (int i = reaction.getReactants(); i < reaction.getMolecules(); i++) {
			String maxString = "";
			int maxIndex = -1;
			for (int j = reaction.getReactants(); j < reaction.getMolecules(); j++) {
				if (maxString.compareTo(idcode[j]) < 0) {
					maxString = idcode[j];
					maxIndex = j;
				}
			}
			if (i > reaction.getReactants()) {
				idcodeSequence.append(MOLECULE_DELIMITER);
				mappingSequence.append(MOLECULE_DELIMITER);
				coordsSequence.append(MOLECULE_DELIMITER);
			}
			idcodeSequence.append(idcode[maxIndex]);
			mappingSequence.append(mapping[maxIndex]);
			coordsSequence.append(coords[maxIndex]);
			idcode[maxIndex] = "";
		}

		String[] result = new String[4];
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

		return result;
	}

	/**
	 * Creates a String containing a unique reaction code by
	 * creating idcodes of every reactant and product and
	 * concatenating them in lexical order.
	 * If mapping information is available this will be encoded
	 * in a 2nd string.
	 * Coordinates, if available, will be encoded in a 3rd string.
	 * If there are drawing objects assigned to this reaction
	 * then these are encoded in a 4th string.
	 *
	 * @return One String with reaction code, coordinates, mapping, drawing objects
	 * as defined by whatToReturn.
	 */
	public static String encode(Reaction reaction, boolean keepAbsoluteCoordinates,
								int whatToReturn)
	{
		String[] result = encode(reaction, keepAbsoluteCoordinates);
		if (result == null) {
			return null;
		}

		StringBuffer buf = new StringBuffer(result[0]);
//		System.out.println("Buffer: 1:" + buf);
		if (whatToReturn != 0) {
			buf.append(OBJECT_DELIMITER);
			if ((whatToReturn & INCLUDE_MAPPING) != 0
				&& result.length > 1
				&& result[1] != null) {
				buf.append(result[1]);
			}
		}
//		System.out.println("Buffer: 2:" + buf);
		whatToReturn &= ~INCLUDE_MAPPING;
		if (whatToReturn != 0) {
			buf.append(OBJECT_DELIMITER);
			if ((whatToReturn & INCLUDE_COORDS) != 0
				&& result.length > 2
				&& result[2] != null) {
				buf.append(result[2]);
			}
		}
//		System.out.println("Buffer: 3:" + buf);
		whatToReturn &= ~INCLUDE_COORDS;
		if (whatToReturn != 0) {
			buf.append(OBJECT_DELIMITER);
			if ((whatToReturn & INCLUDE_DRAWING_OBJECTS) != 0
				&& result.length > 3
				&& result[3] != null) {
				buf.append(result[3]);
			}
		}

//		System.out.println("Buffer: 4:" + buf);
		return buf.toString();
	}

	/**
	 * Creates a Reaction object by interpreting a reaction code,
	 * mapping, coordinates and drawing objects that were earlier created
	 * by this class.
	 * If rxnCoords are relative or null, and if ensureCoordinates==true
	 * then all reactants and products are placed automatically along a
	 * horizontal line. In this case providing a valid Graphics ensure a
	 * more accurate molecule positioning.
	 *
	 * @return Reaction
	 */
	public static Reaction decode(String rxnCode, String rxnMapping, String rxnCoords,
								  String rxnObjects, boolean ensureCoordinates)
	{
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
	 * horizontal line. In this case providing a valid Graphics ensure a
	 * more accurate molecule positioning.
	 *
	 * @return Reaction
	 */
	public static Reaction decode(String s, boolean ensureCoordinates)
	{
		String rxnCode = s;
		String rxnMapping = null;
		String rxnCoords = null;
		String rxnObjects = null;
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
					rxnObjects = s.substring(index3 + 1);
				}
			}
		}

		return decode(rxnCode, rxnMapping, rxnCoords, rxnObjects, ensureCoordinates);
	}


	public static Reaction decode(String s, int type)
	{
		String rxnCode = s;
		String rxnMapping = null;
		String rxnCoords = null;
		String rxnObjects = null;
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
					rxnObjects = s.substring(index3 + 1);
				}
			}
		}

		return decode(rxnCode,
			(type & INCLUDE_MAPPING) == INCLUDE_MAPPING ? rxnMapping : null,
			(type & INCLUDE_COORDS) == INCLUDE_COORDS ? rxnCoords : null,
			(type & INCLUDE_DRAWING_OBJECTS) == INCLUDE_DRAWING_OBJECTS ? rxnObjects : null,
			false);
	}
}
