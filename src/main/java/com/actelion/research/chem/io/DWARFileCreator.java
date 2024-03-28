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

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.StereoMolecule;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;
import java.util.TreeMap;

import static com.actelion.research.chem.io.CompoundTableConstants.NEWLINE_STRING;
import static com.actelion.research.chem.io.CompoundTableConstants.TAB_STRING;

public class DWARFileCreator {
	private BufferedWriter mWriter;
	private DWARFileParser mMasterCopyParser;
	private String[] mRow;
	private ArrayList<String> mColumnTitleList;
	private TreeMap<Integer,Properties> mColumnPropertiesMap;

	/**
	 * Use a DWARFileCreator for writing native DataWarrior files without a CompoundTableModel.
	 * (if you have a populated CompoundTableModel, use the CompoundTableSaver instead).
	 * You may provide a freshly instantiated DWARFileParser as master copy that will provide header,
	 * explanation text, macros, column properties, and runtime properties for the new file.
	 * To use the DWARFileCreator you need to follow these steps:<br>
	 * - instantiate a new DWARFileCreator for a new output file<br>
	 * - optionally call setMasterCopy()<br>
	 * - define individual columns with addStructureColumn(), addDescriptorColumn(), and addAlphanumericalColumn()<br>
	 * - add custom column properties, if you need to with addColumnProperty()
	 * - call writeHeader() once to create the file and write file & table headers<br>
	 * - for every row call setRowStructure() and setRowValue() for cell values and then writeCurrentRow()<br>
	 * - if you didn't call setMasterCopy(), then optionally call writeTemplate() to add runtime properties<br>
	 * - call writeEnd() to close the file<br>
	 * @param writer make sure that it encodes as UTF-8
	 */
	public DWARFileCreator(BufferedWriter writer) {
		mWriter = writer;
		mColumnTitleList = new ArrayList<String>();
	}

	/**
	 * If the file to be created shall resemble another DataWarrior file regarding file
	 * explanation, macro content, columns names, column properties, runtime properties (template),
	 * then one may define a master copy with this method that serves as a blue print.
	 * @param parser DWARFileParser initialized with MODE_BUFFER_HEAD_AND_TAIL
	 */
	public void setMasterCopy(DWARFileParser parser) {
		mMasterCopyParser = parser;
		ArrayList<String> headAndTail = mMasterCopyParser.getHeadOrTail();
		for (String title:headAndTail.get(headAndTail.size()-1).split("\\t"))
			mColumnTitleList.add(title);
	}

	/**
	 * This adds a column to host canonical structure representations (idcodes).
	 * This call allocates the column and defines proper column properties.
	 * If you want DataWarrior to show structure in this column with the original
	 * atom coordinates, you need to also add a column for the encoded atom coordinates.
	 * Otherwise DataWarrior would generate 2D-coordinates for the structures atoms
	 * on the fly.
	 * @param name
	 * @param idColumnName null or column title of other column that hold a compound identifier
	 * @return new structure column index
	 */
	public int addStructureColumn(String name, String idColumnName) {
		int structureColumn = mColumnTitleList.size();
		String title = makeUnique(name);
		mColumnTitleList.add(title);
		addColumnProperty(structureColumn, "specialType", "idcode");
		if (idColumnName != null)
			addColumnProperty(structureColumn, "idColumn", idColumnName);
		return structureColumn;
	}

	/**
	 * Creates a new column to hold encoded 2-dimensional atom coordinates for the structures
	 * stored in the associated structure column. A structure column may not have more than one
	 * associated 2-dimensional atom coordinate column.
	 * This call allocates the column and defines proper column properties.
	 * @param structureColumn
	 * @return new coordinate column index
	 */
	public int add2DCoordinatesColumn(int structureColumn) {
		return addStructureChildColumn("idcoordinates2D", "idcoordinates2D", structureColumn);
	}

	/**
	 * Creates a new column to hold encoded 3-dimensional atom coordinates for the structures
	 * stored in the associated structure column. A structure column may have multiple
	 * associated 3-dimensional atom coordinate columns.
	 * This call allocates the column and defines proper column properties.
	 * @param name 3D-coordinate column names are used to distinguish multiple 3D-coordinate sets
	 * @param structureColumn
	 * @return new coordinate column index
	 */
	public int add3DCoordinatesColumn(String name, int structureColumn) {
		return addStructureChildColumn("idcoordinates3D", name, structureColumn);
	}

	/**
	 * Creates a new column to hold a chemical descriptor for the structures
	 * stored in the associated structure column. A structure column may have multiple
	 * associated descriptor columns, provided that these have distinct types.
	 * This call allocates the column and defines proper column properties.
	 * @param descriptorShortName name used to identify the descriptor type, e.g. 'FragFp'
	 * @param structureColumn
	 * @return new descriptor column index
	 */
	public int addDescriptorColumn(String descriptorShortName, String descriptorVersion, int structureColumn) {
		int column = addStructureChildColumn(descriptorShortName, descriptorShortName, structureColumn);
		addColumnProperty(column, "version", descriptorVersion);
		return column;
	}

	private int addStructureChildColumn(String specialType, String name, int structureColumn) {
		int column = mColumnTitleList.size();
		mColumnTitleList.add(makeUnique(name));
		addColumnProperty(column, "parent", mColumnTitleList.get(structureColumn));
		addColumnProperty(column, "specialType", specialType);
		return column;
	}

	/**
	 * Creates a new standard column to hold any alphanumerical content.
	 * @param name
	 * @return
	 */
	public int addAlphanumericalColumn(String name) {
		mColumnTitleList.add(makeUnique(name));
		return mColumnTitleList.size()-1;
	}

	/**
	 * This method may be used to define column properties, e.g. to associate the column
	 * with an URL for identifier lookups, to define a data range, or other special DataWarrior
	 * related features.
	 * @param column
	 * @param key
	 * @param value
	 */
	public void addColumnProperty(int column, String key, String value) {
		if (mColumnPropertiesMap == null)
			mColumnPropertiesMap = new TreeMap<Integer, Properties>();
		Properties cp = mColumnPropertiesMap.get(column);
		if (cp == null) {
			cp = new Properties();
			mColumnPropertiesMap.put(column, cp);
		}
		cp.setProperty(key, value);
	}

	/**
	 * This method, when called after calling setMasterCopy() or after completing
	 * defining columns, returns the list of defined column names in the correct order.
	 * Especially, when using a master copy, this method informs about the expected columns
	 * and their order to be used for adding the data.
	 * @return manually defined column titles or from master copy, if used
	 */
	public String[] getColumnTitles() {
		return mColumnTitleList.toArray(new String[0]);
		}

	/**
	 * This method, when called after calling setMasterCopy() or after completing
	 * defining columns, returns the column properties of the given column.
	 * @return manually defined column properties or from master copy, if used
	 */
	public Properties getColumnProperties(String columnTitle) {
		if (mMasterCopyParser != null) {
			return mMasterCopyParser.getColumnProperties(columnTitle);
			}
		else {
			for (int i=0; i<mColumnTitleList.size(); i++)
				if (mColumnTitleList.get(i).equals(columnTitle))
					return mColumnPropertiesMap.get(i);
			return null;
			}
		}

	/**
	 * Call this after defining columns and specifying column properties
	 * @param rowCount -1 if row count is not known
	 * @throws IOException
	 */
	public void writeHeader(int rowCount) throws IOException {
		if (mMasterCopyParser == null) {
			mWriter.write("<datawarrior-fileinfo>");
			mWriter.newLine();
			mWriter.write("<version=\"3.3\">");
			mWriter.newLine();
			if (rowCount > 0) {
				mWriter.write("<rowcount=\""+rowCount+"\">");
				mWriter.newLine();
			}
			mWriter.write("</datawarrior-fileinfo>");
			mWriter.newLine();
			writeColumnPropertiesAndTitles();
		}
		else {
			for (String line:mMasterCopyParser.getHeadOrTail()) {
				if (line.trim().matches("<rowcount=\"\\d+\">")) {
					if (rowCount < 0)
						continue;
					line = "<rowcount=\"" + rowCount + "\">";
					}
				if (line.trim().matches("<created=\"\\d+\">")) {
					line = "<created=\"" + System.currentTimeMillis() + "\">";
					}
				mWriter.write(line);
				mWriter.newLine();
			}
		}

		mRow = new String[mColumnTitleList.size()];
	}


	private void writeColumnPropertiesAndTitles() throws IOException {
		if (mColumnPropertiesMap != null) {
			mWriter.write("<column properties>");
			mWriter.newLine();
			for (int column:mColumnPropertiesMap.keySet()) {
				mWriter.write("<columnName=\""+mColumnTitleList.get(column)+"\">");
				mWriter.newLine();
				Properties cp = mColumnPropertiesMap.get(column);
				for (String key:cp.stringPropertyNames()) {
					String value = cp.getProperty(key);
					mWriter.write("<columnProperty=\""+key+"\t"+value+"\">");
					mWriter.newLine();
				}
			}
			mWriter.write("</column properties>");
			mWriter.newLine();
		}

		boolean isFirst = true;
		for (int i=0; i<mColumnTitleList.size(); i++) {
			if (isFirst)
				isFirst = false;
			else
				mWriter.write("\t");
			mWriter.write(mColumnTitleList.get(i));
		}
		mWriter.newLine();
	}

	/**
	 * Calculates the canonical structure representation as idcode from mol and puts
	 * it into the given idcodeColumn. If coordinateColumn != -1 then it also extracts
	 * the structures atom coordinates and puts them into the given coordinateColumn.
	 * It is assumed that the atom coordinates dimensionality (2D vs. 3D) matches that
	 * what was used to define the column.
	 * @param mol
	 * @param idcodeColumn
	 * @param coordsColumn
	 */
	public void setRowStructure(StereoMolecule mol, int idcodeColumn, int coordsColumn) {
		Canonizer canonizer = new Canonizer(mol);
		mRow[idcodeColumn] = canonizer.getIDCode();
		mRow[coordsColumn] = canonizer.getEncodedCoordinates();
	}

	/**
	 * Puts the given idcode into the given column. If you have defined a column for
	 * atom coordinates, you should also call setRowCoordinates() to store the corresponding
	 * coodinate encoding.
	 * @param idcode
	 * @param column
	 */
	public void setRowStructure(String idcode, int column) {
		mRow[column] = idcode;
	}

	/**
	 * Puts the given encoded atom coordinates into the given column.
	 * It is assumed that the atom coordinates dimensionality (2D vs. 3D) matches that
	 * what was used to define the column.
	 * @param coordinates
	 * @param column
	 */
	public void setRowCoordinates(String coordinates, int column) {
		mRow[column] = coordinates;
	}

	/**
	 * Puts the given value into the given column.
	 * @param value
	 * @param column
	 */
	public void setRowValue(String value, int column) {
		value = value.replaceAll("\\r?\\n|\\r", NEWLINE_STRING);
		value = value.replace("\t", TAB_STRING);
		mRow[column] = value;
		}

	/**
	 * Call this once per row after setting individual cell content with
	 * the respective setRowXXXX() methods.
	 * @throws IOException
	 */
	public void writeCurrentRow() throws IOException {
		boolean isFirst = true;
		for (int i=0; i<mRow.length; i++) {
			if (isFirst)
				isFirst = false;
			else
				mWriter.write("\t");
			if (mRow[i] != null) {
				mWriter.write(mRow[i]);
				mRow[i] = null;
			}
		}
		mWriter.newLine();
	}

	public void writeTemplate(ArrayList<String> properties) throws IOException {
		if (properties != null && properties.size() != 0) {
			mWriter.write(CompoundTableConstants.cPropertiesStart);
			mWriter.newLine();
			for (String propertyLine:properties) {
				mWriter.write(propertyLine);
				mWriter.newLine();
			}
			mWriter.write(CompoundTableConstants.cPropertiesStart);
			mWriter.newLine();
		}
	}

	public void writeEnd() throws IOException {
		if (mMasterCopyParser != null) {
			while (mMasterCopyParser.advanceToNext());
			for (String line:mMasterCopyParser.getHeadOrTail()) {
				mWriter.write(line);
				mWriter.newLine();
			}
		}

		mWriter.close();
	}

	/**
	 * Checks column name for uniqueness and modifies it if needed to be unique.
	 * @param name desired column name
	 * @return
	 */
	private String makeUnique(String name) {
		if (name == null || name.trim().length() == 0) {
			name = "Column 1";
		}
		else {
			name = name.trim().replaceAll("[\\x00-\\x1F]", "_");
		}

		while (columnNameExists(name)) {
			int index = name.lastIndexOf(' ');
			if (index == -1)
				name = name + " 2";
			else {
				try {
					int suffix = Integer.parseInt(name.substring(index+1));
					name = name.substring(0, index+1) + (suffix+1);
				}
				catch (NumberFormatException nfe) {
					name = name + " 2";
				}
			}
		}

		return name;
	}

	private boolean columnNameExists(String name) {
		for (String existing:mColumnTitleList)
			if (name.equalsIgnoreCase(existing))
				return true;
		return false;
	}
}
