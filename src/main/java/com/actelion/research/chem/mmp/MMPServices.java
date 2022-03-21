/*
 * Copyright (c) 2017
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
 * 3. Neither the name of the copyright holder nor the
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
 * @author Gregori Gerebtzoff
 */

package com.actelion.research.chem.mmp;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MMPServices {
	private Map<String, MMPReader> mmpReaders;
	
	public MMPServices() throws IOException {
		mmpReaders = new HashMap<String, MMPReader>();
	}
	
	/**
	 * Verifies if the data set exists
	 * @param datasetName Short name of the data set
	 * @return true/false if the data set exists
	 */
	private boolean verifyDatasetname(String datasetName) {
		return mmpReaders.containsKey(datasetName);
	}
	
	/**
	 * Reads a new MMP File
	 * @param br BufferedReader
	 * @param verbose Verbose
	 * @return short name of the data set
	 * @throws IOException
	 * @throws Exception
	 */
	public String readMMPFile(BufferedReader br, boolean verbose) throws IOException, Exception {
		MMPReader mmpReader = new MMPReader(br, verbose);
		String datasetName = mmpReader.getWhat("datasetName"); 
		mmpReaders.put(datasetName, mmpReader);
		return datasetName;
	}
	
	/**
	 * Gets the size of the chemical space for a specific data set
	 * @param datasetName Short name of the data set
	 * @param key idCode of the 'key' (constant part of the molecule)
	 * @return Size of the chemical space, -1 if the data set does not exist
	 */
	public int getChemicalSpaceSize(String datasetName, String key) {
		if (verifyDatasetname(datasetName)) {
			return getChemicalSpaceSize(datasetName, new String[]{key});
		}
		return -1;
	}
	
	/**
	 * Gets the size of the chemical space for a specific data set
	 * @param datasetName Short name of the data set
	 * @param keys Array of 'keys' idCodes (constant part of the molecule, one for single cut, two for double cuts)
	 * @return Size of the chemical space, -1 if the data set does not exist
	 */
	public int getChemicalSpaceSize(String datasetName, String[] keys) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getChemicalSpaceSize(keys);
		}
		return -1;
	}
	
	/**
	 * Gets the chemical space for a specific data set
	 * @param datasetName Short name of the data set
	 * @param keys Array of 'keys' idCodes (constant part of the molecule)
	 * @param value 'value' idCode (variable part of the molecule - not used yet). Can be null
	 * @param dataField Name of the data field for which data should be returned. Can be null
	 * @return a List of tab-delimited [idCodes, moleculeName, data] entries
	 */
	public List<String> getChemicalSpace(String datasetName, String[] keys, String value, String dataField) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getChemicalSpace(keys, value, dataField);
		}
		return null;
	}
	
	/**
	 * Generates the DWAR file of Matched Molecular Pairs for a specific data set and specific transformation 
	 * @param datasetName Short name of the data set
	 * @param idCode idCode of the seed molecule
	 * @param keys Array of 'keys' idCodes (constant part of the molecule)
	 * @param value1 seeded 'value' idCode (variable part of the molecule) 
	 * @param value2 target 'value' idCode (transformation)
	 * @param replacementSize Difference in number of heavy atoms between seed and target fragments
	 * @param properties List of data fields for which data should be returned
	 * @return a String containing the content of the whole DWAR file
	 */
	public String getMMPsDWAR(String datasetName, String idCode, String[] keys, String value1, String value2, int replacementSize, List<String> properties) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getMMPsDWAR(idCode, keys, value1, value2, replacementSize, properties);
		}
		return null;
	}
	
	/**
	 * Generates the DWAR file of the Chemical Space for a specific data set and a specific 'key'
	 * @param datasetName Short name of the data set
	 * @param idCode idCode of the seed molecule
	 * @param keys Array of 'keys' idCodes (constant part of the molecule)
	 * @param dataField Name of the data field for which data should be returned. Can be null.
	 * @return a String containing the content of the whole DWAR file
	 */
	public String getChemicalSpaceDWAR(String datasetName, String idCode, String[] keys, String dataField) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getChemicalSpaceDWAR(idCode, keys, dataField);
		}
		return null;
	}
	
	/**
	 * Gets the number of transformations for a specific data set, seed 'value' and deltas of heavy atoms
	 * @param datasetName Short name of the data set
	 * @param value1 idCode of the seed 'value' (variable part of the molecule)
	 * @param minAtoms minimal delta number of heavy atoms (compared to the seed fragment)
	 * @param maxAtoms maximal delta number of heavy atoms (compared to the seed fragment)
	 * @return the number of transformations, -1 if the data set does not exist
	 */
	public int getTransformationsSize(String datasetName, String value1, int minAtoms, int maxAtoms) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getTransformationsSize(value1, minAtoms, maxAtoms);
		}
		return -1;	
	}
	
	/**
	 * Returns a list of transformations
	 * @param datasetName Short name of the data set
	 * @param keys Array of 'keys' idCodes (constant part of the molecule)
	 * @param value1 seeded 'value' idCode (variable part of the molecule) 
	 * @param minAtoms minimal delta number of heavy atoms (compared to the seed fragment)
	 * @param maxAtoms maximal delta number of heavy atoms (compared to the seed fragment)
	 * @return List of String arrays ([seed, target, number of examples, transformed molecule exists])
	 */
	public List<String[]> getTransformationsTable(String datasetName, String[] keys, String value1, int minAtoms, int maxAtoms) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).transformationsListToTable(keys, value1, minAtoms, maxAtoms);
		}
		return null;
	}

	/**
	 * Generates the main JSON string for a seeded 'value'
	 * @param datasetName Short name of the data set
	 * @param idCode idCode of the whole seed molecule
	 * @param keys Array of 'keys' idCodes (constant part of the molecule)
	 * @param value1 seeded 'value' idCode (variable part of the molecule) 
	 * @param minAtoms minimal delta number of heavy atoms (compared to the seed fragment)
	 * @param maxAtoms maximal delta number of heavy atoms (compared to the seed fragment)
	 * @param sortBy SORT_BY_NUMBER_OF_EXAMPLES or SORT_BY_SIMILARITY
	 * @return a JSON string with all data
	 */
	public String getTransformationsJSON(String datasetName, String idCode, String[] keys, String value1, int minAtoms, int maxAtoms, String sortBy) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getTransformationsJSON(idCode, keys, value1, minAtoms, maxAtoms, sortBy);
		}
		return null;
	}
	
	/**
	 * Generates the DWAR file of the Transformations for a specific data set
	 * @param datasetName Short name of the data set
	 * @param idCode idCode of the whole seed molecule
	 * @param keys Array of 'keys' idCodes (constant part of the molecule)
	 * @param value1 seeded 'value' idCode (variable part of the molecule)
	 * @param minAtoms minimal delta number of heavy atoms (compared to the seed fragment)
	 * @param maxAtoms maximal delta number of heavy atoms (compared to the seed fragment)
	 * @param environmentSize Size of the local environment (0-5)
	 * @param properties List of data fields for which data should be returned
	 * @return a String containing the content of the whole DWAR file
	 */
	public String getTransformationsDWAR(String datasetName, String idCode, String[] keys, String value1, int minAtoms, int maxAtoms, int environmentSize, List<String> properties) { // idCode: whole molecule idCode
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getTransformationsDWAR(idCode, keys, value1, minAtoms, maxAtoms, environmentSize, properties);
		}
		return null;
	}
	
	/**
	 * Returns the idCode of a molecule from its name
	 * @param datasetName Short name of the data set
	 * @param molName Molecule name
	 * @return idCode string
	 */
	public String getIDCodeFromMolName(String datasetName, String molName) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getIDCodeFromMolName(molName);
		}
		return null;
	}
	
	/**
	 * Returns a list of available (numerical) data fields for a specific data set
	 * @param datasetName Short name of the data set
	 * @return List of field names
	 */
	public List<String> getDataFields(String datasetName) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getDataFields("fieldName");
		}
		return null;
	}
	
	/**
	 * Returns a list of long field names for each available numeric field
	 * @param datasetName Short name of the data set
	 * @return List of long field names (or short names if long names are not available)
	 */
	public List<String> getLongDataFields(String datasetName) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getDataFields("longFieldName");
		}
		return null;
	}
	
	/**
	 * Returns the list of categories for each available numeric field
	 * @param datasetName Short name of the data set
	 * @return List of categories, or 'other' if no categories are available
	 */
	public List<String> getCategoryNames(String datasetName) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getDataFields("categoryName");
		}
		return null;
	}
	
	/**
	 * Return the list of the 5% percentiles for each available numeric field
	 * @param datasetName Short name of the data set
	 * @return List of 5% percentiles
	 */
	public List<String> getPercentiles5(String datasetName) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getDataFields("percentile5");
		}
		return null;
	}
	
	/**
	 * Return the list of the 95% percentiles for each available numeric field
	 * @param datasetName Short name of the data set
	 * @return List of 95% percentiles
	 */	
	public List<String> getPercentiles95(String datasetName) {
		if (verifyDatasetname(datasetName)) {
			return mmpReaders.get(datasetName).getDataFields("percentile95");
		}
		return null;
	}
	
	/**
	 * Returns general informations about a specific data set
	 * @param datasetNames Ordered list of data set names; required to ensure that the order is identical to the one in the settings file
	 * @return Tab-delimited [short data set name, number of molecules, data generation date, one random molecule name]
	 */
	public String getDatasetInformations(ArrayList<String> datasetNames) {
		String retVal = "";
		for (String datasetName: datasetNames) {
			if (mmpReaders.containsKey(datasetName)) {
				String date = mmpReaders.get(datasetName).getWhat("date");
				String numberOfMolecules = mmpReaders.get(datasetName).getWhat("numberOfMolecules");
				String randomMoleculeName = mmpReaders.get(datasetName).getWhat("randomMoleculeName");
				if (retVal != "")
					retVal += "\n";
				retVal += datasetName + "\t" + numberOfMolecules + "\t" + date + "\t" + randomMoleculeName;
			}
		}
		return retVal;
	}
}