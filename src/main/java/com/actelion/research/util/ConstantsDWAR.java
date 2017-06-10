package com.actelion.research.util;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.SimpleDateFormat;
import java.util.Locale;

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
public class ConstantsDWAR {

	// This date format is recognized by the DataWarrior
	public static final SimpleDateFormat SIMPLE_DATE_FORMAT = new SimpleDateFormat("MMM dd, yyyy");

	// NumberFormat is not htread save.
	public static final String PATTERN_NF_DWAR4 ="0.0000";

	public static final String PATTERN_NF_DWAR = "0.000############";

	//
//	public static final DecimalFormat NF_DWAR4 = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.US));
//
//	public static final DecimalFormat NF_DWAR = new DecimalFormat("0.0", new DecimalFormatSymbols(Locale.US));



	// Separator for values in DataWarrior file.
	public static final String SEP_VALUE = "; ";
	
	// Line separator for multiple idcodes in a DataWarrior field
	public static final String SEP_LINE = "<NL>";
	

	
	// Constants from Molecule.cAtomLabel
	public static final int ATOM_LABEL_R1 = 142;
	public static final int ATOM_LABEL_R2 = 143;
	public static final int ATOM_LABEL_R3 = 144;
	
	private static int indexAtomLabel = 129;
	
	public static final int ATOM_LABEL_R4 = indexAtomLabel++;
	public static final int ATOM_LABEL_R5 = indexAtomLabel++;
	public static final int ATOM_LABEL_R6 = indexAtomLabel++;
	public static final int ATOM_LABEL_R7 = indexAtomLabel++;
	public static final int ATOM_LABEL_R8 = indexAtomLabel++;
	public static final int ATOM_LABEL_R9 = indexAtomLabel++;
	public static final int ATOM_LABEL_R10 = indexAtomLabel++;
	public static final int ATOM_LABEL_R11 = indexAtomLabel++;
	public static final int ATOM_LABEL_R12 = indexAtomLabel++;
	public static final int ATOM_LABEL_R13 = indexAtomLabel++;
	public static final int ATOM_LABEL_R14 = indexAtomLabel++;
	public static final int ATOM_LABEL_R15 = indexAtomLabel++;
	public static final int ATOM_LABEL_R16 = indexAtomLabel++;

	
	public static final int [] arrAtomLabelRGroups = {
		ATOM_LABEL_R1, 
		ATOM_LABEL_R2, 
		ATOM_LABEL_R3, 
		ATOM_LABEL_R4,
		ATOM_LABEL_R5,
		ATOM_LABEL_R6,
		ATOM_LABEL_R7,
		ATOM_LABEL_R8,
		ATOM_LABEL_R9,
		ATOM_LABEL_R10,
		ATOM_LABEL_R11,
		ATOM_LABEL_R12,
		ATOM_LABEL_R13,
		ATOM_LABEL_R14,
		ATOM_LABEL_R15,
		ATOM_LABEL_R16};

	/**
	 * Replaced by TAG_IDCODE2
	 */
	@Deprecated
	public static final String TAG_IDCODE = "idcode";
	
	public static final String TAG_IDCODE2 = "Structure";

	public static final String TAG_COOR2 = "idcoordinates2D";

	public static final String TAG_COOR3D = "idcoordinates3D";

	public static final String TAG_COOR = "idcoordinates";

	public static final String TAG_ID_QUERY = "Id Query";
	
	public static final String TAG_QUERY_IDENTIFIER = "QueryIdentifier";
	
	public static final String TAG_LIGAND_EFFICIENCY = "LE";

	public static final String TAG_STRUCTURE_INCLUDE = "Include";

	public static final String TAG_STRUCTURE_EXCLUDE = "Exclude";
	
	public static final String TAG_BONDS = "Bonds";
	
	public static final String TAG_ATOMS = "Atoms";
	
	public static final String TAG_NAME = "Name";
	
	public static final String TAG_RECORD_NO = "Record No";
	
	public static final String TAG_ACTNO = "Actelion No";

	public static final String TAG_SOURCE = "Source";




	public static final String [] TAG_NAMES = {TAG_NAME, TAG_ACTNO, TAG_RECORD_NO};

	public static final String ODE_EXTENSION = ".ode";

	public static final String DWAR_EXTENSION = ".dwar";

	public static final String REGEX_FILE_EXTENSION = "(.*\\"+ODE_EXTENSION+")|(.*\\"+DWAR_EXTENSION+")";

	public static final String SDF_EXTENSION = ".sdf";

	
	// Tags for reaction encoding
	// see ReactionEncoder
	// with reaction code, coordinates, mapping, drawing objects
	public static final String TAG_REACTION_CODE = "ReactionCode";
	public static final String TAG_REACTION_COORD = "ReactionCoordinates";
	public static final String TAG_REACTION_MAPPING = "ReactionMapping";
	public static final String TAG_REACTION_DRAW_OBJ = "ReactionDrawingObjects";
	
	public static final String CHARSET_ENCODING = "UTF8";


	

}
