/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
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
     * @return the structure of the records (primary) molecule
     */
    public StereoMolecule getMolecule() {
        if (!mStructureUpToDate) {
            String idcode = getIDCode();
            String coords = getCoordinates();
            mMol = new IDCodeParser(coords == null).getCompactMolecule(idcode, coords);
            if (mMol != null)
                mMol.setName(getMoleculeName());
            }
        mStructureUpToDate = true;
        return mMol;
        }
    }
