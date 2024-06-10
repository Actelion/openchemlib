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

package com.actelion.research.util;

import java.text.SimpleDateFormat;

/**
 * 
 * ConstantsDWAR
 * @author Modest von Korff
 * Jan 15, 2013 MvK Start implementation
 */
public class ConstantsDWAR {

	// This date format is recognized by the DataWarrior
	public static final SimpleDateFormat SIMPLE_DATE_FORMAT = new SimpleDateFormat("MMM dd, yyyy");

	// NumberFormat is not thread save.
	public static final String PATTERN_NF_DWAR4 ="0.0000";

	public static final String PATTERN_NF_DWAR = "0.000############";

	// Program name for the class Platform
	public static final String PLATFORM_NAME_DATAWARRIOR = "datawarrior";

	public static final String PACKAGE_DATAWARRIOR_IDORSIA = "com.actelion.research.datawarrior.DataWarriorActelionLinux";


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

	// From SAR analysis
	public static final String TAG_IDCODE_SCAFFOLD = "Scaffold";

	public static final String TAG_COORDINATES_SCAFFOLD = "scaffoldAtomCoords";

	public static final String TAG_IDCODE2 = "Structure";

	public static final String TAG_COOR2 = "idcoordinates2D";

	public static final String IDCODE_EMPTY = "dH";



	/**
	 * Can be one or multiple sets of 3D coordinates.
	 */
	public static final String TAG_COOR3D = "idcoordinates3D";

	public static final String TAG_COOR3D_MMFF = "3D-Structure (low-energy random, mmff94s+)";

	/**
	 * Contains the idcode and the 3d coordinates.
	 */
	public static final String TAG_CONFORMERSET = "conformerSet";

	public static final String TAG_COOR = "idcoordinates";

	public static final String TAG_ID_QUERY = "Id Query";
	
	public static final String TAG_QUERY_IDENTIFIER = "QueryIdentifier";
	
	public static final String TAG_LIGAND_EFFICIENCY = "LE";

	public static final String TAG_STRUCTURE_INCLUDE = "Include";

	public static final String TAG_STRUCTURE_EXCLUDE = "Exclude";
	
	public static final String TAG_BONDS = "Bonds";
	
	public static final String TAG_ATOMS = "Atoms";

	public static final String TAG_MW = "Total Molweight";

	public static final String TAG_NAME = "Name";
	
	public static final String TAG_RECORD_NO = "Record No";

	@Deprecated // use TAG_CMP_ID instead!
	public static final String TAG_ACTNO = "Idorsia No";

	public static final String TAG_CMP_ID = "CompoundIdentifier";

	public static final String TAG_SOURCE = "Source";

	public static final String ATTR_YES = "Y";
	public static final String ATTR_NO = "N";



	public static final String [] TAG_NAMES = {TAG_NAME, TAG_ACTNO, TAG_RECORD_NO};

	public static final String ODE_EXTENSION = ".ode";

	// DataWarrior file
	public static final String DWAR_EXTENSION = ".dwar";

	public static final String REGEX_FILE_EXTENSION = "(.*\\"+ODE_EXTENSION+")|(.*\\"+DWAR_EXTENSION+")";

	// DataWarrior query file
	public static final String DWAQ_EXTENSION = ".dwaq";

	public static final String SDF_EXTENSION = ".sdf";

	
	// Tags for reaction encoding
	// see ReactionEncoder
	// with reaction code, coordinates, mapping, drawing objects
	public static final String TAG_REACTION_CODE = "ReactionCode";
	public static final String TAG_REACTION_COORD = "ReactionCoordinates";
	public static final String TAG_REACTION_MAPPING = "ReactionMapping";
	public static final String TAG_REACTION_DRAW_OBJ = "ReactionDrawingObjects";
	
	public static final String CHARSET_ENCODING = "UTF8"
			;
	public static final int MOLECULE_START_R_GROUPS = 129;


	

}
