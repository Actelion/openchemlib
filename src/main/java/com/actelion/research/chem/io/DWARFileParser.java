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

package com.actelion.research.chem.io;

import com.actelion.research.chem.*;
import com.actelion.research.chem.descriptor.DescriptorConstants;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.descriptor.DescriptorHandlerStandard2DFactory;
import com.actelion.research.chem.descriptor.DescriptorHelper;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;
import com.actelion.research.io.BOMSkipper;
import com.actelion.research.util.BinaryDecoder;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;

public class DWARFileParser extends CompoundFileParser implements DescriptorConstants,CompoundTableConstants {

    public static final int MODE_COORDINATES_PREFER_2D = 1;
    public static final int MODE_COORDINATES_PREFER_3D = 2;
    public static final int MODE_COORDINATES_REQUIRE_2D = 3;
    public static final int MODE_COORDINATES_REQUIRE_3D = 4;
    private static final int MODE_COORDINATE_MASK = 7;
    public static final int MODE_BUFFER_HEAD_AND_TAIL = 8;
    public static final int MODE_EXTRACT_DETAILS = 16;

	private String[]		mFieldName;
	private String[]		mFieldData;
	private String			mLine,mCoordinate3DColumnName;
	private int[]			mFieldIndex;
    private int             mRecordCount,mMode;
	private int				mIDCodeColumn,mCoordinateColumn,mCoordinate2DColumn,mCoordinate3DColumn,
							mMoleculeNameColumn,mFragFpColumn;
    private TreeMap<String,Properties> mColumnPropertyMap;
	private TreeMap<String,SpecialField> mSpecialFieldMap;
	private TreeMap<String,Integer> mDescriptorColumnMap;
	private ArrayList<String> mHeadOrTailLineList;
	private HashMap<String,byte[]> mDetails;

    /**
     * Constructs a DWARFileParser from a file name with coordinate mode MODE_COORDINATES_PREFER_2D.
     * @param fileName
     */
	public DWARFileParser(String fileName) {
        try {
            mReader = new BufferedReader(new InputStreamReader(new FileInputStream(fileName), StandardCharsets.UTF_8));
			BOMSkipper.skip(mReader);
            mMode = MODE_COORDINATES_PREFER_2D;
            init();
            }
        catch (IOException e) {
            mReader = null;
            }
		}

    /**
     * Constructs a DWARFileParser from a File with coordinate mode MODE_COORDINATES_PREFER_2D.
     * @param file
     */
	public DWARFileParser(File file) {
        try {
            mReader = new BufferedReader(new InputStreamReader(new FileInputStream(file), StandardCharsets.UTF_8));
			BOMSkipper.skip(mReader);
            mMode = MODE_COORDINATES_PREFER_2D;
            init();
            }
        catch (IOException e) {
            mReader = null;
            }
		}

    /**
     * Constructs a DWARFileParser from a Reader with coordinate mode MODE_COORDINATES_PREFER_2D.
     * @param reader
     */
	public DWARFileParser(Reader reader) {
        try {
			mReader = (reader instanceof BufferedReader) ? (BufferedReader)reader : new BufferedReader(reader);
            mMode = MODE_COORDINATES_PREFER_2D;
            init();
            }
        catch (IOException e) {
            mReader = null;
            }
		}

    /**
     * Constructs a DWARFileParser from a file name with the specified coordinate mode.
     * @param fileName
     * @param mode one of 4 MODE_COORDINATE... modes
     */
    public DWARFileParser(String fileName, int mode) {
        try {
            mReader = new BufferedReader(new InputStreamReader(new FileInputStream(fileName), StandardCharsets.UTF_8));
			BOMSkipper.skip(mReader);
            mMode = mode;
            init();
            }
        catch (IOException e) {
            mReader = null;
            }
        }

    /**
     * Constructs a DWARFileParser from a File with the specified coordinate mode.
     * @param file
     * @param mode one of 4 MODE_COORDINATE... modes
     */
    public DWARFileParser(File file, int mode) {
        try {
            mReader = new BufferedReader(new InputStreamReader(new FileInputStream(file), StandardCharsets.UTF_8));
			BOMSkipper.skip(mReader);
            mMode = mode;
            init();
            }
        catch (IOException e) {
            mReader = null;
            }
        }

    /**
     * Constructs a DWARFileParser from a Reader with the specified coordinate mode.
     * @param reader
     * @param mode one of 4 MODE_COORDINATE... modes
     */
    public DWARFileParser(Reader reader, int mode) {
        try {
			mReader = (reader instanceof BufferedReader) ? (BufferedReader)reader : new BufferedReader(reader);
            mMode = mode;
            init();
            }
        catch (IOException e) {
            mReader = null;
            }
        }

    private String readHeadOrTailLine() throws IOException {
    	String line = mReader.readLine();
    	if ((mMode & MODE_BUFFER_HEAD_AND_TAIL) != 0 && line != null)
    		mHeadOrTailLineList.add(line);
    	return line;
    	}

    private void init() throws IOException {
	    setDescriptorHandlerFactory(DescriptorHandlerStandard2DFactory.getFactory());

    	if ((mMode & MODE_BUFFER_HEAD_AND_TAIL) != 0)
    		mHeadOrTailLineList = new ArrayList<>();

    	int coordinateMode = mMode & MODE_COORDINATE_MASK;
	
		String line = readHeadOrTailLine();
        if (line == null
         || !line.equals(cNativeFileHeaderStart))
            throw new IOException("no header found");

        mRecordCount = -1;
        line = readHeadOrTailLine();
        while (line != null
            && !line.equals(cNativeFileHeaderEnd)) {

            if (line.startsWith("<"+cNativeFileVersion)) {
                String version = extractValue(line);
                if (!version.startsWith("3.")
                 && !version.isEmpty())
                    throw new IOException("unsupported .dwar file version");
                }

            else if (line.startsWith("<"+cNativeFileRowCount)) {
                try {
                    mRecordCount = Integer.parseInt(extractValue(line));
                    }
                catch (NumberFormatException e) {}
                }
  
            line = readHeadOrTailLine();
            }

        line = readHeadOrTailLine();

        while (line != null
         && (line.equals(cFileExplanationStart)
          || line.equals(cMacroListStart))) {
            line = readHeadOrTailLine();
            while (line != null
                && !line.equals(cFileExplanationEnd)
                && !line.equals(cMacroListEnd))
                line = readHeadOrTailLine();
            line = readHeadOrTailLine();
        	}

        mColumnPropertyMap = new TreeMap<>();

        if (line != null
         && line.equals(cColumnPropertyStart)) {
            line = readHeadOrTailLine();
            String columnName = null;
            while (line != null
                && !line.equals(cColumnPropertyEnd)) {

                if (line.startsWith("<"+cColumnName)) {
                	columnName = extractValue(line);
                	mColumnPropertyMap.put(columnName, new Properties());
                    }

                else if (line.startsWith("<"+cColumnProperty)) {
                    String[] property = extractValue(line).split("\\t");
                    if(property.length==1) {
                    	mColumnPropertyMap.get(columnName).setProperty(property[0],"");
                    }
                    else {
                    	mColumnPropertyMap.get(columnName).setProperty(property[0], property[1]);
                    }
                }

                line = readHeadOrTailLine();
                }

            line = readHeadOrTailLine();
            }

        mSpecialFieldMap = new TreeMap<>();	// only take those columns that have a special type
        for (String columnName:mColumnPropertyMap.keySet()) {
        	Properties properties = mColumnPropertyMap.get(columnName);
        	String specialType = properties.getProperty(cColumnPropertySpecialType);
        	if (specialType != null)
        		mSpecialFieldMap.put(columnName, new SpecialField(
        				columnName,
        				specialType,
        				properties.getProperty(cColumnPropertyParentColumn),
        				properties.getProperty(cColumnPropertyRelatedIdentifierColumn),
        				properties.getProperty(cColumnPropertyDescriptorVersion)
        				));
        	}
        
        ArrayList<String> columnNameList = new ArrayList<String>();
        ArrayList<Integer> columnIndexList = new ArrayList<Integer>();

		if (line == null)
            throw new IOException("unexpected end of file");

		int fromIndex = 0;
		int toIndex = 0;
		int sourceColumn = 0;
		do {
			String columnName;
			toIndex = line.indexOf('\t', fromIndex);
			if (toIndex == -1) {
				columnName = line.substring(fromIndex);
				}
			else {
				columnName = line.substring(fromIndex, toIndex);
				fromIndex = toIndex+1;
				}

            if (mSpecialFieldMap.containsKey(columnName)) {
                mSpecialFieldMap.get(columnName).fieldIndex = sourceColumn;
                }
            else {
    			columnNameList.add(columnName);
    			columnIndexList.add(sourceColumn);
                }

			sourceColumn++;
			} while (toIndex != -1);

		mFieldName = new String[columnNameList.size()];
		mFieldIndex = new int[columnNameList.size()];
		for (int i=0; i<columnNameList.size(); i++) {
			mFieldName[i] = columnNameList.get(i);
			mFieldIndex[i] = columnIndexList.get(i);
			}

		mFieldData = new String[sourceColumn];

        mIDCodeColumn = -1;
        mCoordinateColumn = -1;
		mCoordinate2DColumn = -1;
		mCoordinate3DColumn = -1;
        mMoleculeNameColumn = -1;
        mFragFpColumn = -1;
        SpecialField idcodeColumn = mSpecialFieldMap.get("Structure");
        if (idcodeColumn == null
         || !idcodeColumn.type.equals(cColumnTypeIDCode)) {
            for (SpecialField specialColumn:mSpecialFieldMap.values()) {
                if (specialColumn.type.equals(cColumnTypeIDCode)) {
                    if (idcodeColumn == null
                     || idcodeColumn.fieldIndex > specialColumn.fieldIndex)
                        idcodeColumn = specialColumn;
                    }
                }
            }
        if (idcodeColumn != null) {
            if (idcodeColumn.idColumn != null) {
                for (int i=0; i<mFieldName.length; i++) {
                    if (idcodeColumn.idColumn.equals(mFieldName[i])) {
                        mMoleculeNameColumn = mFieldIndex[i];
                        break;
                        }
                    }
                }
            mIDCodeColumn = idcodeColumn.fieldIndex;
            for (SpecialField specialColumn:mSpecialFieldMap.values()) {
                if (idcodeColumn.name.equals(specialColumn.parent)) {
                    if (DESCRIPTOR_FFP512.shortName.equals(specialColumn.type)
                     && DescriptorHandlerLongFFP512.VERSION.equals(specialColumn.version))
                        mFragFpColumn = specialColumn.fieldIndex;
					else if (cColumnType2DCoordinates.equals(specialColumn.type))
						mCoordinate2DColumn = specialColumn.fieldIndex;
					else if (cColumnType3DCoordinates.equals(specialColumn.type)) {
						mCoordinate3DColumn = specialColumn.fieldIndex;
						mCoordinate3DColumnName = specialColumn.name;
						}

                    if (DescriptorHelper.isDescriptorShortName(specialColumn.type)
					 && DescriptorHelper.getDescriptorInfo(specialColumn.type).version.equals(specialColumn.version)) {
                        if (mDescriptorColumnMap == null)
                            mDescriptorColumnMap = new TreeMap<String,Integer>();
                        mDescriptorColumnMap.put(specialColumn.type, specialColumn.fieldIndex);
                        }
                    }
                }

			if (mCoordinate2DColumn != -1
			 && (coordinateMode == MODE_COORDINATES_REQUIRE_2D
			  || coordinateMode == MODE_COORDINATES_PREFER_2D
			  || (coordinateMode == MODE_COORDINATES_PREFER_3D && mCoordinate3DColumn == -1)))
				mCoordinateColumn = mCoordinate2DColumn;

			if (mCoordinate3DColumn != -1
			 && (coordinateMode == MODE_COORDINATES_REQUIRE_3D
			  || coordinateMode == MODE_COORDINATES_PREFER_3D
			  || (coordinateMode == MODE_COORDINATES_PREFER_2D && mCoordinate2DColumn == -1)))
				mCoordinateColumn = mCoordinate3DColumn;
            }
	    }

    /**
     * If you don't read any records after calling this method,
     * don't forget to call close() to close the underlying file.
     * @return whether the file contains chemical structures
     */
    public boolean hasStructures() {
    	return (mIDCodeColumn != -1);
    	}

    /**
     * @return whether the file contains chemical structures with explicit atom coordinates
     */
    public boolean hasStructureCoordinates() {
    	return (mCoordinateColumn != -1);
    	}

	/**
	 * @return whether the file contains chemical structures with explicit atom coordinates
	 */
	public boolean hasStructureCoordinates2D() {
		return (mCoordinate2DColumn != -1);
		}

	/**
	 * @return whether the file contains chemical structures with explicit atom coordinates
	 */
	public boolean hasStructureCoordinates3D() {
		return (mCoordinate3DColumn != -1);
		}

	public String getStructureCoordinates3DColumnName() {
		return mCoordinate3DColumnName;
		}

	public String[] getFieldNames() {
		return mFieldName;
		}

	/**
	 * @param columnName
	 * @return field index for special fields, e.g. to be used for getSpecialFieldData()
	 */
	public int getSpecialFieldIndex(String columnName) {
		for (SpecialField sf:mSpecialFieldMap.values())
			if (columnName.equals(sf.name))
				return sf.fieldIndex;

		return -1;
		}

	/**
	 * @param parentColumnName
	 * @param childType
	 * @return field index for special fields, e.g. to be used for getSpecialFieldData()
	 */
	public int getChildFieldIndex(String parentColumnName, String childType) {
		for (SpecialField sf:mSpecialFieldMap.values())
			if (parentColumnName.equals(sf.parent) && childType.equals(sf.type))
				return sf.fieldIndex;

		return -1;
		}

    public int getRowCount() {
        return mRecordCount;
        }

    /**
     * Provided that the mode contains MODE_BUFFER_HEAD_AND_TAIL, then this method
     * returns a list of all header/footer rows of the DWAR file. If this method is
     * called before all rows have been read, then the header lines including column
     * properties and the column title line are returned. If this method is
     * called after all rows have been read, then all lines after the data table, i.e. the
     * runtime properties, are returned.
     * @return
     */
    public ArrayList<String> getHeadOrTail() {
        return mHeadOrTailLineList;
        }

    /**
     * Provided that the mode contains MODE_EXTRACT_DETAILS, then this method
     * returns a map of all embedded detail objects of the DWAR file.
     * This method must not be called before all rows have been read.
     * @return
     */
    public HashMap<String,byte[]> getDetails() {
    	return mDetails;
    	}

    /**
     * Returns the entire line containing all row data
     * @return
     */
    public String getRow() {
        return mLine;
        }

	/**
	 * Returns the raw data in the following format:
	 * Columns are sorted according to the order of how they appear in DataWarrior:
	 *
	 * @param includeHeaderRow
	 * @return
	 */
	public String[][] getRawData(boolean includeHeaderRow, boolean structureAsSmiles) {

		String[] fn  = getFieldNames();
		TreeMap<String,SpecialField> sfn = getSpecialFieldMap();

		List<String> allFieldNames = new ArrayList<>();

		allFieldNames.addAll(Arrays.stream(fn).collect(Collectors.toList()));
		List<String> allFieldNamesForOutput = new ArrayList<>();
		allFieldNamesForOutput.addAll(allFieldNames);
		for(String sfi : sfn.keySet().stream().sorted(  (x,y) -> Integer.compare( sfn.get(x).fieldIndex , sfn.get(y).fieldIndex ) ).collect(Collectors.toList())) {
			if(sfn.get(sfi).type.equals(cColumnTypeIDCode)) {
				allFieldNames.add(sfn.get(sfi).fieldIndex, sfi);
				allFieldNamesForOutput.add(sfi);
			}
		}

		int nDataRows = getRowCount();
		int nOutputRows = nDataRows + (includeHeaderRow?1:0);

		String[][] rawData = new String[ nOutputRows ][ allFieldNames.size() ];

		if(includeHeaderRow) {
			for(int zi=0;zi<allFieldNamesForOutput.size();zi++) {
				String fi = allFieldNamesForOutput.get(zi);
				rawData[0][zi] = fi;
			}
		}

		next();

		for(int zi=0;zi<nDataRows;zi++) {
			int zOutput = zi + (includeHeaderRow?1:0);

			for(int zj=0;zj<allFieldNamesForOutput.size();zj++) {
				// determine if structure column or not
				String fi = allFieldNamesForOutput.get(zj); //allFieldNames.get(zj);//fn[zj];
				String rawData_i = "";

				if(sfn.containsKey(fi)) {
					rawData_i = getSpecialFieldData(getSpecialFieldIndex(fi));
					// special column..
					if( sfn.get(fi).type.equals(cColumnTypeIDCode) ) {
						if(structureAsSmiles) {
							try{
								StereoMolecule mi = new StereoMolecule();
								IDCodeParser icp = new IDCodeParser();
								icp.parse(mi,rawData_i);
								IsomericSmilesCreator isc = new IsomericSmilesCreator(mi);
								rawData_i = isc.getSmiles();
							}
							catch(Exception ex) {
								System.out.println("Exception for: "+rawData_i);
							}
						}
					}
				}
				else {
					// normal column:
					rawData_i = getFieldData( getFieldIndex( fi ) );
				}
				rawData[zOutput][zj] = rawData_i;
			}
			this.next();
		}
		return rawData;
	}

	public static void main(String args[]) {
		//DWARFileParser dwfp = new DWARFileParser("C:\\data\\ActelionFragmentLibrary_smallFragments.dwar");
		//DWARFileParser dwfp = new DWARFileParser("C:\\data\\hitexpclusteringtool\\input\\Ox1_HTS_215_Hits.dwar");
		DWARFileParser dwfp = new DWARFileParser("C:\\dev\\pyocl\\data\\Wikipedia_Compounds_6.dwar");

		String[][] rawData_a = dwfp.getRawData(true,false);
		dwfp = new DWARFileParser("C:\\dev\\pyocl\\data\\Wikipedia_Compounds_6.dwar");
		String[][] rawData_b = dwfp.getRawData(false,true);
		dwfp = new DWARFileParser("C:\\dev\\pyocl\\data\\Wikipedia_Compounds_6.dwar");
		String[][] rawData_c = dwfp.getRawData(true,true);
		System.out.println("DATA:\n");
		for(int zi=0;zi<rawData_a.length;zi++) {
			for(int zj=0;zj<rawData_a[0].length;zj++) {
				System.out.print(rawData_a[zi][zj]+" ;; ");
			}
			System.out.println("");
		}
		for(int zi=0;zi<rawData_b.length;zi++) {
			for(int zj=0;zj<rawData_b[0].length;zj++) {
				System.out.print(rawData_b[zi][zj]+" ;; ");
			}
			System.out.println();
		}
		for(int zi=0;zi<rawData_c.length;zi++) {
			for(int zj=0;zj<rawData_c[0].length;zj++) {
				System.out.print(rawData_c[zi][zj]+" ;; ");
			}
			System.out.println();
		}
	}


	protected boolean advanceToNext() {
		if (mReader == null)
			return false;

		mLine = null;
		try {
			mLine = mReader.readLine();
			if (mLine == null
			 || mLine.equals(cPropertiesStart)
			 || mLine.equals(cHitlistDataStart)
			 || mLine.equals(cDetailDataStart)
			 || mLine.startsWith(cDataDependentPropertiesStart)) {
				if ((mMode & MODE_BUFFER_HEAD_AND_TAIL) != 0) {
					mHeadOrTailLineList.clear();
					mHeadOrTailLineList.add(mLine);
					}
				while (mLine != null) {
					if (mLine.equals(cDetailDataStart)
					 && (mMode & MODE_EXTRACT_DETAILS) != 0)
						extractDetails();
					mLine = readHeadOrTailLine();
					}
			    mReader.close();
				return false;
				}
			}
		catch (IOException e) {
			return false;
			}

		int column = 0;
		int index1 = 0;
		int index2 = mLine.indexOf('\t');
		while (index2 != -1) {
			if (column<mFieldData.length)
				mFieldData[column] = mLine.substring(index1, index2).replace(NEWLINE_STRING, cLineSeparator).replace(TAB_STRING, "\t");
			column++;
			index1 = index2+1;
			index2 = mLine.indexOf('\t', index1);
			}
		if (column<mFieldData.length)
			mFieldData[column] = mLine.substring(index1).replace(NEWLINE_STRING, cLineSeparator).replace(TAB_STRING, "\t");

		return true;
		}

    private void extractDetails() {
		mDetails = new HashMap<String,byte[]>();
		try {
		    while (true) {
				String theLine = readHeadOrTailLine();
				if (theLine == null
				 || theLine.equals(cDetailDataEnd)) {
					break;
					}

				if (theLine.startsWith("<"+cDetailID)) {
					String detailID = extractValue(theLine);
					BinaryDecoder decoder = new BinaryDecoder(mReader);
					int size = decoder.initialize(8);
					byte[] detailData = new byte[size];
					for (int i=0; i<size; i++)
						detailData[i] = (byte)decoder.read();
					mDetails.put(detailID, detailData);
					}
				}
			}
		catch (Exception e) {
			mDetails = null;
			}
    	}

    /**
     * @return the row content of the first column containing chemical structures
     */
	public String getIDCode() {
        if (mIDCodeColumn == -1)
            return null;
		String s = mFieldData[mIDCodeColumn];
		return (s == null || s.isEmpty()) ? null : s;
		}

	/**
	 * This returns encoded atom coordinates according to the defined mode.
	 * If the compound file does not contain atom coordinates, then null is returned.
	 * If mode is one of MODE_COORDINATES_REQUIRE... and the required coordinate dimensionality
	 * (2D or 3D) is not available then null is returned.
	 * If mode is one of MODE_COORDINATES_PREFER... and the preferred coordinate dimensionality
	 * (2D or 3D) is not available then coordinates in another dimensionality are returned.
	 */
	public String getCoordinates() {
        if (mCoordinateColumn == -1)
            return null;
		String s = mFieldData[mCoordinateColumn];
		return (s == null || s.isEmpty()) ? null : s;
		}

	public String getCoordinates2D() {
		if (mCoordinate2DColumn == -1)
			return null;
		String s = mFieldData[mCoordinate2DColumn];
		return (s == null || s.isEmpty()) ? null : s;
		}

	public String getCoordinates3D() {
		if (mCoordinate3DColumn == -1)
			return null;
		String s = mFieldData[mCoordinate3DColumn];
		return (s == null || s.isEmpty()) ? null : s;
		}

	public Reaction getReaction(boolean ensureCoordinates) {
		String rxnCode = null;
		String rxnColumn = null;
		for (SpecialField sf : mSpecialFieldMap.values()) {
			if (CompoundTableConstants.cColumnTypeRXNCode.equals(sf.type)) {
				rxnColumn = sf.name;
				rxnCode = mFieldData[sf.fieldIndex];
				break;
			}
		}

		if (rxnCode == null)
			return null;

		String rxnMapping = null;
		String rxnCoords = null;
		String rxnCatalysts = null;
		for (SpecialField sf : mSpecialFieldMap.values()) {
			if (CompoundTableConstants.cColumnTypeReactionMapping.equals(sf.type) && rxnColumn.equals(sf.parent))
				rxnMapping = mFieldData[sf.fieldIndex];
			else if (CompoundTableConstants.cColumnType2DCoordinates.equals(sf.type) && rxnColumn.equals(sf.parent))
				rxnCoords = mFieldData[sf.fieldIndex];
			else if (CompoundTableConstants.cColumnTypeIDCode.equals(sf.type)
					&& sf.name.equals(getColumnProperties(rxnColumn).getProperty(cColumnPropertyRelatedCatalystColumn)))
				rxnCatalysts = mFieldData[sf.fieldIndex];
		}

		return rxnMapping == null && rxnCoords == null && rxnCatalysts == null ? ReactionEncoder.decode(rxnCode, ensureCoordinates)
				: ReactionEncoder.decode(rxnCode, rxnMapping, rxnCoords, null, rxnCatalysts, ensureCoordinates, null);
	}

	public String getMoleculeName() {
        if (mMoleculeNameColumn == -1)
            return null;
        String s = mFieldData[mMoleculeNameColumn];
        return (s == null || s.isEmpty()) ? null : s;
        }

    public Object getDescriptor(String shortName) {
        Integer column = (mDescriptorColumnMap == null) ? null : mDescriptorColumnMap.get(shortName);
        String s = (column == null) ? null : mFieldData[column.intValue()];
        return (s == null || s.isEmpty()) ? super.getDescriptor(shortName)
             : (getDescriptorHandlerFactory() == null) ? null
             : getDescriptorHandlerFactory().getDefaultDescriptorHandler(shortName).decode(s);
        }

    /**
     * @return the String encoded FragFp descriptor of the first column containing chemical structures
     */
	public String getIndex() {
        if (mFragFpColumn == -1)
            return null;
		String s = mFieldData[mFragFpColumn];
		return (s == null || s.isEmpty()) ? null : s;
		}

	public String getFieldData(int no) {
		return mFieldData[mFieldIndex[no]];
		}

	/**
	 * Returns a columnName->SpecialField map of all non-alphanumerical columns.
	 * SpecialField.type is one of the types defined in CompoundTableConstants:
	 * cColumnTypeIDCode,cColumnTypeRXNCode,cColumnType2DCoordinates,cColumnType3DCoordinates,
	 * cColumnTypeAtomColorInfo, and descriptor shortNames;
	 * @return special fields
	 */
	public TreeMap<String,SpecialField> getSpecialFieldMap() {
	    return mSpecialFieldMap;
	    }

	/**
	 * @param fieldIndex is available from special-field-TreeMap by getSpecialFieldMap().get(columnName).fieldIndex
	 * @return String encoded data content of special field, e.g. idcode
	 */
    public String getSpecialFieldData(int fieldIndex) {
        return mFieldData[fieldIndex];
        }

    /**
     * Returns the original column properties of any source column by column name.
     * @param columnName
     * @return
     */
    public Properties getColumnProperties(String columnName) {
    	return mColumnPropertyMap.get(columnName);
    	}

    private String extractValue(String theLine) {
        int index1 = theLine.indexOf("=\"") + 2;
        int index2 = theLine.indexOf("\"", index1);
        return theLine.substring(index1, index2);
        }

	public class SpecialField {
	    public String name;
	    public String type;
	    public String parent;
	    public String idColumn;
	    public String version;
	    public int fieldIndex;

	    public SpecialField(String name, String type, String parent, String idColumn, String version) {
	        this.name = name;
	        this.type = type;
	        this.parent = parent;
	        this.idColumn = idColumn;
	        this.version = version;
	        this.fieldIndex = -1;
	        }
	    }
    }
