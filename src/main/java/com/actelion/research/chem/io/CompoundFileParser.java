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

import java.io.BufferedReader;
import java.io.IOException;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorHandlerFactory;

public abstract class CompoundFileParser {
    private StereoMolecule mMol;
    private DescriptorHandlerFactory mDHFactory;
    private boolean mStructureUpToDate,mIDCodeUpToDate;
    private String mIDCode,mCoords;
	protected BufferedReader mReader;

    /**
     * Creates the proper parser for the given type of compound file (currently SD or DWAR).
     * @param fileName
     * @return parser or null, if the file doesn't exist or cannot be accessed
     */
	public static CompoundFileParser createParser(String fileName) {
        CompoundFileParser parser = null;
        int fileType = CompoundFileHelper.getFileType(fileName);
        if (fileType == CompoundFileHelper.cFileTypeDataWarrior)
            parser = new DWARFileParser(fileName);
        else if (fileType == CompoundFileHelper.cFileTypeSD)
            parser = new SDFileParser(fileName);

        return parser.mReader == null ? null : parser;
        }

    /**
     * Compiles all column names that contain alpha-numerical information.
     * Columns containing chemistry objects, coordinates or descriptors don't
     * appear in the list.
     * @return columns name array in the order of appearance
     */
    abstract public String[] getFieldNames();

    /**
     * Returns the cell content of the current row. Multi-line cell entries are
     * separated by a '\n' character.
     * @param column refers to alpha-numerical columns only, as getFieldNames()
     * @return
     */
    abstract public String getFieldData(int column);

    /**
     * Depending on data source returns the total row count or -1 if unknown
     * @return number of rows or -1
     */
    abstract public int getRowCount();

	/**
	 * Dont't call this method directly. Use next() instead.
	 * @return false if there is no next row
	 */
    abstract protected boolean advanceToNext();

    /**
     *
     * @return whether the file was found and open to accept next() calls
     */
    public boolean isOpen() {
        return mReader != null;
        }

    /**
     * Advances the row counter to the next row
     * @return false if there is no next row
     */
    public boolean next() {
        mStructureUpToDate = false;
        mIDCodeUpToDate = false;
        return advanceToNext();
        }

    /**
     * Closes the underlying reader. Call this, if you don't read all records of the file.
     * The reader is closed automatically after the last record has been read.
     */
    public final void close() {
    	if (mReader != null) {
    		try {
    			mReader.close();
    			}
    		catch (IOException ioe) {}
    		}
    	}

    /**
     * Either this method and getCoordinates() or getMolecule() must be overwritten!!!
     * @return idcode of first chemical structure column of the current row
     */
    public String getIDCode() {
        updateIDCodeAndCoords();
        return mIDCode;
        }

    /**
     * Either getIDCode and this method or getMolecule() must be overwritten!!!
     * @return idcoords of first chemical structure column of the current row
     */
    public String getCoordinates() {
        updateIDCodeAndCoords();
        return mCoords;
        }

    /**
     * @return name/id of (primary) chemical structure of the current row
     */
    abstract public String getMoleculeName();

    /**
     * If a requested descriptor is not available in a particuar
     * compound record, the parser can create one itself, provided its
     * DescriptorHandlerFactory knows the descriptor name. The default
     * DescriptorHandlerFactory is null, thus one needs to set one in
     * order to allow the parser to create descriptors.
     * @param factory
     */
    public void setDescriptorHandlerFactory(DescriptorHandlerFactory factory) {
        mDHFactory = factory;
        }

    /**
     * @return currently used DescriptorHandlerFactory
     */
    public DescriptorHandlerFactory getDescriptorHandlerFactory() {
        return mDHFactory;
        }

    /**
     * @param fieldName
     * @return index of the field with the given name, -1 if fieldName doesn't exist
     */
    public int getFieldIndex(String fieldName) {
        String[] name = getFieldNames();

        if (name != null)
            for (int i=0; i<name.length; i++)
                if (fieldName.equals(name[i]))
                    return i;

        return -1;
        }

    /**
     * If the file source contains encoded descriptors, then overwrite this method
     * to save the calculation time.
     * @param shortName
     * @return descriptor as int[] or whatever is the descriptors binary format
     */
    public Object getDescriptor(String shortName) {
        if (mDHFactory != null) {
            DescriptorHandler dh = mDHFactory.getDefaultDescriptorHandler(shortName);
            Object d = dh.createDescriptor(getMolecule());
            return dh.calculationFailed(d) ? null : d;
            }
        return null;
        }

    private void updateIDCodeAndCoords() {
        if (!mIDCodeUpToDate) {
            try {
                StereoMolecule mol = new StereoMolecule(getMolecule());
                mol.normalizeAmbiguousBonds();
                mol.canonizeCharge(true);
                Canonizer canonizer = new Canonizer(mol);
                mIDCode = canonizer.getIDCode();
                mCoords = canonizer.getEncodedCoordinates();
                }
            catch (Exception e) {
                mIDCode = null;
                mCoords = null;
                }
            mIDCodeUpToDate = true;
            }
        }

    /**
     * Either this method or getIDCode() and getCoordinates() must be overwritten!!!
     * @return the structure of the records (primary) molecule or null
     */
    public StereoMolecule getMolecule() {
        if (!mStructureUpToDate) {
            String idcode = getIDCode();
            String coords = getCoordinates();
            mMol = null;
            try {
	            mMol = new IDCodeParser(coords == null).getCompactMolecule(idcode, coords);
	            if (mMol != null)
	                mMol.setName(getMoleculeName());
            	}
            catch (Exception e) {}
            }
        mStructureUpToDate = true;
        return mMol;
        }
    }
