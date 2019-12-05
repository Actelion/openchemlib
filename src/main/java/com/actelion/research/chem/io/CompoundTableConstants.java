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

package com.actelion.research.chem.io;

public interface CompoundTableConstants {
    String cColumnUnassignedItemText = "<Unassigned>";
    String cColumnUnassignedCode = "<none>";
    String cColumnNameRowList = "List '";

	// visible special columns
    String cColumnTypeIDCode = "idcode";
    String cColumnTypeRXNCode = "rxncode";

    String[] cParentSpecialColumnTypes = {
                                    cColumnTypeIDCode,
                                    cColumnTypeRXNCode };

        // non-parent special columns cannot be displayed
    String cColumnType2DCoordinates = "idcoordinates2D";
    String cColumnType3DCoordinates = "idcoordinates3D";
    String cColumnTypeAtomColorInfo = "atomColorInfo";
    String cColumnTypeReactionMapping = "atomMapping";
    String cColumnTypeReactionObjects = "reactionObjects";
        // in addition to these all DescriptorHandler.SHORT_NAMEs are valid column types

    String cReactionPartReaction = "reaction";  // this may only be used, if a molecule type descriptor is calculated from all merged reaction molecules
    String cReactionPartReactants = "reactants";
//    String cReactionPartPseudoCatalysts = "catalysts";
    String cReactionPartProducts = "products";
    String cReactionPartDelimiter = " of ";

    int cTextExclusionTypeContains = 1;
    int cTextExclusionTypeStartsWith = 2;
    int cTextExclusionTypeEquals = 3;
    int cTextExclusionTypeRegEx = 4;

    int cMaxTextCategoryCount = 65536;
    int cMaxDateOrDoubleCategoryCount = 16384;

    // summary mode for displaying values.
    int cSummaryModeNormal = 0;
    int cSummaryModeMean = 1;
    int cSummaryModeMedian = 2;
    int cSummaryModeMinimum = 3;
    int cSummaryModeMaximum = 4;
    int cSummaryModeSum = 5;
    String[] cSummaryModeText = { "All Values", "Mean Value", "Median Value", "Lowest Value", "Highest Value", "Sum" };
    String[] cSummaryModeCode = { "displayNormal","displayMean","displayMedian","displayMin","displayMax","displaySum" };

    int cDataTypeAutomatic = 0;
    int cDataTypeNumerical = 1;
    int cDataTypeDate = 2;
    int cDataTypeString = 3;
    String[] cDataTypeCode = {"auto", "num", "date", "text"};
    String[] cDataTypeText = {"Automatic", "Numerical", "Date", "Text"};

    // highlight mode for part-of-structure highlighting depending on current record similarity
    int cStructureHiliteModeFilter = 0;
    int cStructureHiliteModeCurrentRow = 1;
    int cStructureHiliteModeNone = 2;
    String[] cStructureHiliteModeText = { "Most Recent Filter", "Current Row Similarity", "No Highlighting" };
    String[] cStructureHiliteModeCode = { "hiliteFilter", "hiliteCurrent", "hiliteNone" };

    // highlight mode for part-of-reaction highlighting depending on current record similarity
    int cReactionHiliteModeReactionCenter = 0;
    int cReactionHiliteModeNone = 2;
    String[] cReactionHiliteModeText = { "Reaction Center", "No Highlighting" };
    String[] cReactionHiliteModeCode = { "hiliteReactionCenter", "hiliteNone" };

    String NEWLINE_REGEX = "\\r?\\n|\\r";	// regex for platform independent new line char(s)
	String NEWLINE_STRING = "<NL>";	// used in .dwar, .txt and .cvs files to indicate next line within a cell
	String TAB_STRING = "<TAB>";	// used in .dwar, .txt and .cvs files to indicate a tabulator within a cell
    String cEntrySeparator = "; ";
    byte[] cEntrySeparatorBytes = { ';', ' '};
    String cLineSeparator = "\n";
    byte cLineSeparatorByte = '\n'; // this must be equal to cLineSeparator
    String cRangeSeparation = " <= x < ";
    String cRangeNotAvailable = "<none>";
    String cDefaultDetailSeparator = "|#|";
    String cDetailIndexSeparator = ":";
    String cTextMultipleCategories = "<multiple categories>";

    String cColumnPropertySpecialType = "specialType";
    String cColumnPropertyParentColumn = "parent";
    String cColumnPropertyRelatedIdentifierColumn = "idColumn";
    String cColumnPropertyRelatedCatalystColumn = "catalystColumn";    // one could think of coupling solvent, condition, etc as well
    String[] cColumnRelationTypes = {cColumnPropertyRelatedIdentifierColumn, cColumnPropertyRelatedCatalystColumn};

    String cColumnPropertyGroupName = "groupName";
    String cColumnPropertyUseThumbNail = "useThumbNail";
    String cColumnPropertyImagePath = "imagePath";
    String cColumnPropertyIsFragment = "isFragment";    // specifies for structure & reaction columns, whether the fragment bit is set, when editing a new object
    String cColumnPropertyReactionPart = "reactionPart";
    String cColumnPropertyIsClusterNo = "isClusterNo";
    String cColumnPropertyDataMin = "dataMin";
    String cColumnPropertyDataMax = "dataMax";
    String cColumnPropertyCyclicDataMax = "cyclicDataMax";
    String cColumnPropertyDetailCount = "detailCount";
    String cColumnPropertyDetailName = "detailName";
    String cColumnPropertyDetailType = "detailType";
    String cColumnPropertyOrbitType = "orbitType";
    String cColumnPropertyDetailSource = "detailSource";
    String cColumnPropertyDetailSeparator = "detailSeparator";
    String cColumnPropertyDescriptorVersion = "version";
    String cColumnPropertyIsDisplayable = "isDisplayable";
    String cColumnPropertyBinBase = "binBase";
    String cColumnPropertyBinSize = "binSize";
    String cColumnPropertyBinIsLog = "binIsLog";
    String cColumnPropertyBinIsDate = "binIsDate";
    String cColumnPropertyLookupCount = "lookupCount";
    String cColumnPropertyLookupName = "lookupName";
    String cColumnPropertyLookupURL = "lookupURL";
    String cColumnPropertyLookupEncode = "lookupEncode";
    String cColumnPropertyLookupDetailURL = "lookupDetailURL";
    String cColumnPropertyLaunchCount = "launchCount";
    String cColumnPropertyLaunchName = "launchName";
    String cColumnPropertyLaunchCommand = "launchCommand";
    String cColumnPropertyLaunchOption = "launchOption";
    String cColumnPropertyLaunchDecoration = "launchDecoration";
	String cColumnPropertyLaunchAllowMultiple = "launchAllowMultiple";
    String cColumnPropertyOpenExternalName = "openExternalName";
    String cColumnPropertyOpenExternalPath = "openExternalPath";
    String cColumnPropertyReferencedColumn = "refColumn";
    String cColumnPropertyReferenceStrengthColumn = "refStrengthColumn";
    String cColumnPropertyReferenceType = "refType";
    String cColumnPropertyReferenceTypeRedundant = "redundant";	// a connection is always referenced on both records
    String cColumnPropertyReferenceTypeTopDown = "topdown";	// a connection is only referenced from top record
    String cColumnPropertyFormula = "formula";
    String cColumnPropertySuperposeMolecule = "superposeMol";	// idcode+coords to be displayed in every cell
    String cColumnPropertySuperpose = "superpose";  // cSuperposeValueReferenceRow or null
    String cColumnPropertySuperposeAlign = "align";  // cSuperposeAlignValueShape or null
    String cColumnPropertyCommentUploadStatus = "commentUploadStatus";
    String cColumnPropertyCommentDepartment = "commentDepartment";

    String cSuperposeValueReferenceRow = "refRow";  // "reference" or null
    String cSuperposeAlignValueShape = "shape";  // "reference" or null

    String cNativeFileHeaderStart = "<datawarrior-fileinfo>";
    String cNativeFileHeaderEnd = "</datawarrior-fileinfo>";
    String cNativeFileVersion = "version";
    String cNativeFileRowCount = "rowcount";
    String cNativeFileCreated = "created";

    String cColumnPropertyStart = "<column properties>";
    String cColumnPropertyEnd = "</column properties>";
    String cColumnName = "columnName";
    String cColumnProperty = "columnProperty";

    String cHitlistDataStart = "<hitlist data>";
    String cHitlistDataEnd = "</hitlist data>";
    String cHitlistName = "hitlistName";
    String cHitlistData = "hitlistData";

    String cDetailDataStart = "<detail data>";
    String cDetailDataEnd = "</detail data>";
    String cDetailID = "detailID";

    String cTemplateTagName = "datawarrior properties";
    String cPropertiesStart = "<"+cTemplateTagName+">";
    String cPropertiesEnd ="</"+cTemplateTagName+">";

    String cViewConfigTagName = "view configuration";

    String cDataDependentPropertiesStart = "<data dependent properties type=\"";
    String cDataDependentPropertiesEnd = "</data dependent properties>";

    String cExtensionNameFileExplanation = "explanation";
    String cExtensionNameMacroList = "macroList";

    String cFileExplanationStart = "<datawarrior "+cExtensionNameFileExplanation+">";
    String cFileExplanationEnd = "</datawarrior "+cExtensionNameFileExplanation+">";

    String cMacroListStart = "<datawarrior "+cExtensionNameMacroList+">";
    String cMacroListEnd = "</datawarrior "+cExtensionNameMacroList+">";

    String cAutoStartMacro = "autoStartMacro";
	}
