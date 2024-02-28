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
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.io.BOMSkipper;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.TreeMap;

/**
 * Quick and dirty RD-File parser that reads non-hierarchical RD-Files.
 */
public class RDFileParser {
    private static final int DEFAULT_RECORDS_TO_INSPECT = 1024;
	public static final String cNewLineString = "\n";

	private BufferedReader      mReader;
	private String[]			mFieldName;
	private String[]			mFieldData;
	private int					mNoOfRecords;
	private String              mLine,mIRegNo,mERegNo;
	private TreeMap<String,String> mDataMap;

	public RDFileParser(String fileName) {
	    mNoOfRecords = 0;
		try {
			mReader = new BufferedReader(new InputStreamReader(new FileInputStream(fileName), StandardCharsets.UTF_8));
			BOMSkipper.skip(mReader);
			readHeader();
		} catch (IOException e) {
			mReader = null;
		}
		}
		
	public RDFileParser(File file) {
        mNoOfRecords = 0;
		try {
    		mReader = new BufferedReader(new InputStreamReader(new FileInputStream(file), StandardCharsets.UTF_8));
			BOMSkipper.skip(mReader);
			readHeader();
		} catch (IOException e) {
			mReader = null;
		}
	}
	
	public RDFileParser(Reader reader) {
        mNoOfRecords = 0;
		mReader = (reader instanceof BufferedReader) ? (BufferedReader)reader : new BufferedReader(reader);
		try { readHeader(); } catch (IOException ioe) {}
		}
		
	private void readHeader() throws IOException {
		if (!"$RDFILE 1".equals(mReader.readLine())     // check for version 1
		 || mReader.readLine() == null                  // read and skip time stamp
		 || (mLine = mReader.readLine()) == null) {     // read first reaction tag
			mReader.close();
			mReader = null;
			return;
			}
		mDataMap = new TreeMap<>();
	}

	public boolean hasNext() {
		return mReader != null;
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

	public TreeMap<String,String> getFieldData() {
		return mDataMap;
		}

	public String getFieldData(String key) {
		return mDataMap.get(key);
		}

	public String getERegNo() {
		return mERegNo;
		}

	public String getIRegNo() {
		return mIRegNo;
		}

	/**
	 * RD-files may contains lists of molecules or lists of reactions.
	 * For RD-files containing reactions use this method.
	 * @return reaction or null, if the end of the file has been reached
	 */
	public Reaction getNextReaction() {
		if (mReader == null)
			return null;

		mIRegNo = null;
		mERegNo = null;
		Reaction rxn = null;
		mDataMap.clear();
		String key = null;

		while (mLine != null) {
			if (mLine.startsWith("$REREG ")) {
				mERegNo = mLine.substring(7).trim();
				}
			else if (mLine.startsWith("$RIREG ")) {
				mIRegNo = mLine.substring(7).trim();
				}
			else if (mLine.startsWith("$RFMT")) {
				rxn = new Reaction();
				try {
					if (!new RXNFileParser().parse(rxn, mReader)) {
						mReader = null;
						return null;
						}
					}
				catch (Exception e) {
					mReader = null;
					return null;
					}
				}
			else if (mLine.startsWith("$DTYPE ")) {
				key = mLine.substring(7).trim();
				}
			else if (mLine.startsWith("$DATUM ")) {
				if (key != null) {
					String value = mLine.substring(7).trim();
					if (value.length() != 0)
						mDataMap.put(key, value);
					key = null;
					}
				}

			try {
				mLine = mReader.readLine();
				if (mLine == null) {
					mReader.close();
					mReader = null;
					}
				}
			catch (IOException e) {
				return null;
				}

			if (isReactionNext())
				break;
			}

		mNoOfRecords++;
		return rxn;
		}

	public boolean isMoleculeNext() {
		return mLine != null && (mLine.startsWith("$MFMT") || mLine.startsWith("$MEGEG") || mLine.startsWith("$MIREG"));
	}

	public boolean isReactionNext() {
		return mLine != null && (mLine.startsWith("$RFMT") || mLine.startsWith("$REGEG") || mLine.startsWith("$RIREG"));
	}

	/**
	 * RD-files may contains lists of molecules or lists of reactions.
	 * For RD-files containing molecule use this method.
	 * @return molecule or null, if the end of the file has been reached
	 */
	public StereoMolecule getNextMolecule() {
		if (mReader == null)
			return null;

		mIRegNo = null;
		mERegNo = null;
		StereoMolecule mol = null;
		mDataMap.clear();
		String key = null;

		while (mLine != null) {
			if (mLine.startsWith("$MEREG ")) {
				mERegNo = mLine.substring(7).trim();
				}
			else if (mLine.startsWith("$MIREG ")) {
				mIRegNo = mLine.substring(7).trim();
				}
			else if (mLine.startsWith("$MFMT")) {
				mIRegNo = mLine.substring(7).trim();
				mol = new StereoMolecule();
				try {
					if (!new MolfileParser().parse(mol, mReader)) {
						mReader = null;
						return null;
						}
					}
				catch (Exception e) {
					mReader = null;
					return null;
					}
				}
			else if (mLine.startsWith("$DTYPE ")) {
				key = mLine.substring(7).trim();
				}
			else if (mLine.startsWith("$DATUM ")) {
				if (key != null) {
					String value = mLine.substring(7).trim();
					if (value.length() != 0)
						mDataMap.put(key, value);
					key = null;
					}
				}

			try {
				mLine = mReader.readLine();
				if (mLine == null) {
					mReader.close();
					mReader = null;
					}
				}
			catch (IOException e) {
				return null;
				}

			if (isMoleculeNext())
				break;
			}

		mNoOfRecords++;
		return mol;
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

	private void extractAllFieldNames(int recordsToInspect) {
		int records = 0;
		UniqueStringList fieldNameList = new UniqueStringList();

		while (records < recordsToInspect) {
			Object chem = isReactionNext() ? getNextReaction() :  getNextMolecule();
			if (chem == null)
				break;

			for (String key:mDataMap.keySet())
				fieldNameList.addString(key);

			records++;
			}

		mFieldName = fieldNameList.toArray();
		}
	}
