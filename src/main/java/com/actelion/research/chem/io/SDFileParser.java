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

import com.actelion.research.chem.MolfileParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.UniqueStringList;
import com.actelion.research.io.BOMSkipper;

import java.io.*;

public class SDFileParser extends CompoundFileParser {
    private static final int DEFAULT_RECORDS_TO_INSPECT = 10240;
    private static final String[] cIDFieldNames = { "Idorsia No", "Actelion No", "ID", "IDNUMBER", "COMPOUND_ID", "NAME", "COMPND" };
	public static final String cNewLineString = "\n";

	private StringBuilder		mMolfileBuffer,mDataBuffer;
	private StereoMolecule		mMol;
	private String[]			mFieldName;
	private String[]			mFieldData;
	private int					mNoOfRecords,mIDFieldIndex;

	public SDFileParser(String fileName) {
		this(fileName, null);
		}

	public SDFileParser(String fileName, String[] fieldName) {
	    mNoOfRecords = -1;
	    
		mFieldName = fieldName;
		
		try {
			mReader = new BufferedReader(new InputStreamReader(new FileInputStream(fileName), "UTF-8"));
			BOMSkipper.skip(mReader);
		} catch (IOException e) {
			mReader = null;
		}
		
		
		init();
		}


	public SDFileParser(File file) {
		this(file, null);
		}


	public SDFileParser(File file, String[] fieldName) {
        mNoOfRecords = -1;
		mFieldName = fieldName;
		try {
    		mReader = new BufferedReader(new InputStreamReader(new FileInputStream(file), "UTF-8"));
			BOMSkipper.skip(mReader);
		} catch (IOException e) {
			mReader = null;
		}
		
		init();
	}


	public SDFileParser(Reader reader) {
		this(reader, null);
	}


	public SDFileParser(Reader reader, String[] fieldName) {
        mNoOfRecords = -1;
		mFieldName = fieldName;
		mReader = (reader instanceof BufferedReader) ? (BufferedReader)reader : new BufferedReader(reader);
		
		init();
		}


	private void init() {
		mMolfileBuffer = new StringBuilder(10240);
		mDataBuffer = new StringBuilder(10240);
		}
	
	private void extractAllFieldNames(int recordsToInspect) {
	    int records = 0;
//		TreeSet<String> fieldNameList = new TreeSet<String>(); Changed to keep the original order of field names. TLS 6Jan16
		UniqueStringList fieldNameList = new UniqueStringList();

		while (records < recordsToInspect) {
			String line;
			try {
				line = mReader.readLine();
				}
			catch (IOException e) {
				if (records < recordsToInspect)
					mNoOfRecords = records;
				break;
				}

			if (line == null) {
				if (records < recordsToInspect)
					mNoOfRecords = records;
				break;
				}

			if (line.startsWith("$$$$"))
				records++;

			if (line.startsWith(">")) {
				String fieldName = extractFieldName(line);
				if (fieldName != null)
					fieldNameList.addString(fieldName);
				}
			}

		try {
			mReader.close();
		    }
		catch (IOException e) {}

		mFieldName = fieldNameList.toArray();
		}


	/**
	 * Only accurate if getFieldNames() or getFieldNames(int) was called earlier
	 * and if the number of records of the SD-file is smaller than the number
	 * of records that were examined within the the getFieldNames() method.
	 * If not all records of the file were seen, then -1 is returned.
	 * For getRowCount() to reliably return the record count call getFieldNames(Integer.MAX_VALUE) first.
	 * @return number of rows or -1
	 */
	public int getRowCount() {
		return mNoOfRecords;
		}


	protected boolean advanceToNext() {
		if (mReader == null)
			return false;

// removed 13.8.2012 TLS; no need to read molfile in order to advance to the next record
//		if (mMolfileBuffer.length() != 0)
//			return true;

		mMolfileBuffer.setLength(0);
		mDataBuffer.setLength(0);
		
    	mMol = null;

		boolean molfileComplete = false;
		int fieldIndex = -1;
		String fieldName = null;
		String line;
		mFieldData = (mFieldName == null) ? null : new String[mFieldName.length];
		mIDFieldIndex = -1;

		do {
			try {
				line = mReader.readLine();
				if (line == null) {
	    			mMolfileBuffer.setLength(0);
		    		mReader.close();
			    	return false;
				    }
				}
			catch (IOException e) {
				mMolfileBuffer.setLength(0);
				return false;
				}

			if (!molfileComplete) {
				if (line.startsWith(">")) {	// to handle sd-record with molfiles without 'M  END'
					molfileComplete = true;
					mMolfileBuffer.append("M  END");
		    		mMolfileBuffer.append('\n');
		    		mDataBuffer.append(line);
		    		mDataBuffer.append('\n');
					}
				else {
					mMolfileBuffer.append(line);
		    		mMolfileBuffer.append('\n');
			    	if (line.startsWith("M  END"))
						molfileComplete = true;
					continue;
					}
				}
			else {
	    		mDataBuffer.append(line);
	    		mDataBuffer.append('\n');
				}

			if (mFieldName != null) {
				if (line.length() == 0) {
					fieldIndex = -1;
					}
				else if (fieldIndex == -1) {
					fieldName = extractFieldName(line);
					if (fieldName != null) {
					    // find fieldIndex to given fieldName
						fieldIndex = -1;
						for (int field=0; field<mFieldName.length; field++) {
							if (fieldName.equals(mFieldName[field])) {
								fieldIndex = field;
								break;
								}
							}

						// check whether field qualifies as compound identifier
						if (mIDFieldIndex == -1) {
                            for (String idName:cIDFieldNames) {
                                if (fieldName.equals(idName)) {
                                    mIDFieldIndex = fieldIndex;
                                    break;
                                    }
                                }
                            
                            }
						}
					}
				else {
					if (mFieldData[fieldIndex] == null) {
						mFieldData[fieldIndex] = line;
						}
					else {
						mFieldData[fieldIndex] = mFieldData[fieldIndex].concat(cNewLineString).concat(line);
						}
					}
				}
			} while (!line.startsWith("$$$$"));

		return true;
		}


	/**
	 * @return the molecule of the current record (null in case of parsing error)
	 */
	public StereoMolecule getMolecule() {
	    if (mMol != null)
	        return mMol;

	    mMol = new MolfileParser().getCompactMolecule(getNextMolFile());
	    if (mMol != null && (mMol.getName() == null || mMol.getName().length() == 0))
	        mMol.setName(getMoleculeName());
	    return mMol;
	    }


    public String getMoleculeName() {
        return (mIDFieldIndex != -1 && mFieldData != null) ?
                mFieldData[mIDFieldIndex] : (mMol != null) ? mMol.getName() : null;
        }


	/**
	 * Returns the molfile of the current record
	 * as one big String as it was read from the input file.
	 * Line endings are '\n'.
	 * @return 
	 */
    public String getNextMolFile() {
		String molfile = mMolfileBuffer.toString();
		return molfile;
		}


	/**
	 * Returns the field data of the current record
	 * as one big String as it was read from the input file.
	 * Line endings are '\n'.
	 * @return 
	 */
    public String getNextFieldData() {
		String fieldData = mDataBuffer.toString();
		return fieldData;
		}


	/**
	 * Returns a list of field names. If the field names were not passed
	 * to the constructor of SDFileParser, this method parses the file/reader
	 * to extract all field names and uses up this SDFileParser. In this case
	 * one needs to instantiate a new SDFileParser to sequentially iterate
	 * through the file/reader's records and supply the field name array to
	 * the constructor.
	 * @return array of field names
	 */
	public String[] getFieldNames() {
        if (mFieldName == null)
            extractAllFieldNames(DEFAULT_RECORDS_TO_INSPECT);

        return mFieldName;
	    }

    /**
     * Returns a list of field names. If the field names were not passed
     * to the constructor of SDFileParser, this method parses the file/reader
     * <recordsToInspect> records to extract all field names and uses up this
     * SDFileParser. In this case one needs to instantiate a new SDFileParser
     * to sequentially iterate through the file/reader's records and supply
     * the field name array to the constructor.
     * @return array of field names
     */
	public String[] getFieldNames(int recordsToInspect) {
        if (mFieldName == null)
            extractAllFieldNames(recordsToInspect);

        return mFieldName;
	    }

	/*	public boolean moreRecordsAvailable() {
		if (mFieldDataList == null || mFieldDataList.size() == 0)
			return false;

		mCurrentFieldData = mFieldDataList.get(0);
		mFieldDataList.remove(0);
		return true;
		}
*/

	public String getFieldData(int index) {
		if (mFieldData == null)
			return null;

		return mFieldData[index];
		}


	protected String extractFieldName(String line) {
		if (line.length() == 0
		 || line.charAt(0) != '>')
			return null;

		int index = 1;
		int openBracket = 0;
		int closeBracket = 0;
		while (index < line.length()) {
			if (line.charAt(index) == '<') {
				if (openBracket != 0)
					return null;
				openBracket = index;
				}
			else if (line.charAt(index) == '>') {
				if (closeBracket != 0)
					return null;
				closeBracket = index;
				}
			index++;
			}

		if (openBracket != 0 && openBracket < closeBracket)
			return line.substring(openBracket+1, closeBracket);

		// allow for MACCS-II field numbers, which have format DTn
		index = line.indexOf("DT", 1);
		if (index == -1)
			return null;

		int i = index+2;
		while (line.length()>i && Character.isDigit(line.charAt(i)))
			i++;
		
		return (i == index+2) ? null : line.substring(index, i);
		}
	}
