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

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.io.CompoundFileParser;
import com.actelion.research.chem.mmp.MMPFragmenter.MoleculeIndexID;

import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public class MMP {
	private static String VERSION = "1.1";									  // Version 1.1 has a modified Enumerator method (slightly faster reading)
	private static final String r1H = MMPFragmenter.createR1HMoleculeID();    // R1-H molecule IDCode (for Hydrogen replacements)
	private static final boolean TRANSFORM_UM_TO_PIC50 = true;                // transforms uM to pIC50 if a field name finishes with "_uM"
	private int maxValueAtoms = 0;                                            // size of the biggest fragment (value), for enumeration
	
	private MMPairs matchedMolecularPairs;                                    // container for all Matched Molecular Pairs (result of enumeration)
	private HashMap<Integer, HashMap<String, ArrayList<int[]>>> mMPIndex;     // valueAtoms -- keys - list of values ({valueIndex, moleculeIndex})
	private HashMap<String, ArrayList<MoleculeIndex>> wholeMoleculesIndex;    // IDCode - molIndex, molName, molData
	private List<List<String[]>> moleculesFragmentsID;                        // container for "clean" fragments (without R-groups) and MolName
	private MMPUniqueFragments mmpUniqueFragments;                            // Set of unique fragmentsID
	private MMPFragments mmpFragments;                                        // Container for molecule fragments (keys, value, molecule index)
	private int r1HIndex;                                                     // idCode of R-[H] molecule
	private String[] fieldNames;                                              // list of field names
	private boolean[] fieldNumerics;                                          // list of booleans (true if the field contains only numerical data)
	private float[] fieldPercentiles5;                                        // list of 5th percentiles
	private float[] fieldPercentiles95;                                       // list of 95th percentile
	private int moleculesRowCount;
	private String datasetName;
	
	static public class MoleculeIndex {
		int moleculeIndex;
		String moleculeName;
		String[] moleculeData;
		String moleculeIDCode;
		String moleculeIDCoord;
		
		public MoleculeIndex(int moleculeIndex, String moleculeIDCoord, String moleculeIDCode, String moleculeName, String[] moleculeData) {
			this.moleculeIndex = moleculeIndex;
			this.moleculeIDCoord = moleculeIDCoord;
			this.moleculeIDCode = moleculeIDCode;
			this.moleculeName = moleculeName;
			this.moleculeData = moleculeData;
		}
		
		public MoleculeIndex(int moleculeIndex, String moleculeName, String[] moleculeData) {
			this.moleculeIndex = moleculeIndex;
			this.moleculeName = moleculeName;
			this.moleculeData = moleculeData;
			this.moleculeIDCode = null;
			this.moleculeIDCoord = null;
		}
		
		public void setIDCode(String moleculeIDCode) {
			this.moleculeIDCode = moleculeIDCode;
		}
		
		public void setIDCoord(String moleculeIDCoord) {
			this.moleculeIDCoord = moleculeIDCoord;
		}
	}

	private static String getDateAndTime() {
		return new SimpleDateFormat("d-MMM-yyyy HH:mm:ss").format(Calendar.getInstance().getTime());
	}

	/**
	 * Helper function 
	 * @param str String to be analyzed
	 * @return true/false if the input String is numeric
	 */
	private static boolean isNumeric(String str) {
		if (str == null) {
			return true;
		}
		try {
			Double.parseDouble(str);
		}
		catch(NumberFormatException nfe) {  
			return false;  
		}  
		return true;  
	}
	
	/**
	 * Gets the mapped values corresponding to a defined size of the 'value' (variable part of the molecule)
	 * @param valueAtoms Number of heavy atoms of the 'value' (variable part of the molecule)
	 * @return Hashmap of keys - list of values ({valueIndex, moleculeIndex})
	 */
	private HashMap<String, ArrayList<int[]>> getIndex(int valueAtoms) {
		if (mMPIndex.containsKey(valueAtoms)) {
			return mMPIndex.get(valueAtoms);
		}
		return null;
	}
	
	/**
	 * Add new value to mMPIndex
	 * @param valueAtoms Number of heavy atoms of the 'value' (variable part of the molecule)
	 * @param keys Array of one (single cut) or two (double cut) strings for the 'key' (constant part(s) of the molecule) 
	 * @param values Array of two integers (valueIndex, moleculeIndex)
	 * @param isH true if the 'value' is a Hydrogen atom
	 * @return true/false if the new values were added
	 */
	private boolean addValues(int valueAtoms, String keys, int[] values, boolean isH) {
		boolean added = true;
		HashMap<String, ArrayList<int[]>> keysHash = null;
		ArrayList<int[]> valuesList = null;
		if (mMPIndex.containsKey(valueAtoms)) {
			keysHash = mMPIndex.get(valueAtoms);
		}
		else {
			keysHash = new HashMap<String, ArrayList<int[]>>();
		}
		if (keysHash.containsKey(keys)) {
			valuesList = keysHash.get(keys);
			if (valuesList == null)
				valuesList = new ArrayList<int[]>();
			if (isH == true && !valuesList.isEmpty()) {
				int[] lastItem = valuesList.get(valuesList.size() - 1);
				if (lastItem[0] != r1HIndex)
					valuesList.add(values);
				else
					added = false;
			}
			else {
				valuesList.add(values);
			}
		}
		else {
			valuesList = new ArrayList<int[]>();
			valuesList.add(values);
		}
		keysHash.put(keys, valuesList);
		mMPIndex.put(valueAtoms, keysHash);
		return added;
	}
	
	/**
	 * Generates a hash table of keys - list of values;<br>
	 * for double cuts, one key consists of a '\t'-separated string<br>
	 * of the two 'left' and 'right' fragmentsID.
	 * @param datasetName Name of the data set
	 * @param compoundFileParser Compound File Parser (SD Reader, database link, ...)
	 * @param verbose Verbose
	 * @throws IOException
	 */
	public MMP(String datasetName, CompoundFileParser compoundFileParser, boolean verbose) throws IOException {
		mMPIndex = new HashMap<Integer, HashMap<String, ArrayList<int[]>>>();
		wholeMoleculesIndex = new LinkedHashMap<String, ArrayList<MoleculeIndex>>();
		moleculesFragmentsID = new ArrayList<List<String[]>>();
		matchedMolecularPairs = new MMPairs();
		mmpUniqueFragments = new MMPUniqueFragments();
		r1HIndex = mmpUniqueFragments.addFragment(r1H, 0, null); // we force the number of atoms to 0 because [H] counts for 1
		mmpFragments = new MMPFragments();
		moleculesRowCount = 0;
		fieldNames = compoundFileParser.getFieldNames();
		fieldNumerics = new boolean[fieldNames.length];
		fieldPercentiles5 = new float[fieldNames.length];
		fieldPercentiles95 = new float[fieldNames.length];
		IDCodeParser idCodeParser = new IDCodeParser();
		boolean isSDFileParser = compoundFileParser.getClass().getName().contains("SDFileParser");
		@SuppressWarnings("unchecked")
		  ArrayList<Float>[] fieldDatas = (ArrayList<Float>[])new ArrayList[fieldNames.length];
		for (int i = 0; i < fieldNames.length; i++) {
			fieldDatas[i] = new ArrayList<Float>();
		}
		Arrays.fill(fieldNumerics, true);
		this.datasetName = datasetName;
		NumberFormat formatter = new DecimalFormat("#.##");
		int molCounter = 0;
		if (verbose) {
			if (compoundFileParser.getRowCount() != -1) {
				System.out.println(getDateAndTime() + ": fragmenting " + compoundFileParser.getRowCount() + " molecules...");
			}
			else {
				System.out.println(getDateAndTime() + ": fragmenting molecules...");
			}
		}
		while (compoundFileParser.next()) {
			StereoMolecule mol = compoundFileParser.getMolecule();
			mol.stripSmallFragments();
			String moleculeName = (compoundFileParser.getMoleculeName() == null) ? mol.getName() : compoundFileParser.getMoleculeName();
			String molID = compoundFileParser.getIDCode();
			String molIDCoord = compoundFileParser.getCoordinates();
			if (isSDFileParser) {
				// TODO: avoid re-parsing the IDCode to ensure that the saved molecule is the same as this one - problem occurs with SDF files...
				idCodeParser.parse(mol, molID);
			}
			String moleculeData[] = new String[fieldNames.length];
			for (int i=0; i<fieldNames.length; i++) {
				if (fieldNumerics[i] != false) {
					String fieldData = compoundFileParser.getFieldData(i);
					if (fieldData == null || fieldData.equals("N/A") || fieldData.equals("?")) {
						fieldData = null;
					}
					else if ((int) fieldData.charAt(0) == 65533) { // weird unicode question mark character
						fieldData = null;
					}
					moleculeData[i] = fieldData;
					if (!isNumeric(fieldData) && !isNumeric(fieldData.substring(1)) && !isNumeric(fieldData.substring(2))) { // for <, >, >= and <= symbols
						fieldNumerics[i] = false;
					}
					else if (fieldData != null && isNumeric(fieldData)) {
						if (TRANSFORM_UM_TO_PIC50 && fieldNames[i].endsWith("_uM")) {
							float data = round((float)-Math.log10(Float.parseFloat(fieldData)*1.0E-6), 3);
							fieldData = Float.toString(data);
							moleculeData[i] = fieldData;
							fieldDatas[i].add(data);
						}
						else {
							fieldDatas[i].add(Float.parseFloat(fieldData));
							moleculeData[i] = formatter.format(Float.parseFloat(fieldData));
						}
					}
					else if (fieldData != null && !fieldData.startsWith("<") && !fieldData.startsWith(">")) {
						fieldNumerics[i] = false;
					}
					else if (fieldData != null && (fieldData.startsWith("<=") || fieldData.startsWith(">="))) {
						if (TRANSFORM_UM_TO_PIC50 && fieldNames[i].endsWith("_uM")) {
							float data = round((float)-Math.log10(Float.parseFloat(fieldData.substring(2))*1.0E-6), 3);
							if (fieldData.startsWith(">=")) {
								moleculeData[i] = "<=" + Float.toString(data);
							}
							else {
								moleculeData[i] = ">=" + Float.toString(data);
							}
						}
						else {
							moleculeData[i] = fieldData.substring(0, 2) + formatter.format(Float.parseFloat(fieldData.substring(2)));
						}
					}
					else if (fieldData != null && (fieldData.startsWith("<") || fieldData.startsWith(">"))) {
						if (TRANSFORM_UM_TO_PIC50 && fieldNames[i].endsWith("_uM")) {
							float data = round((float)-Math.log10(Float.parseFloat(fieldData.substring(1))*1.0E-6), 3);
							if (fieldData.startsWith(">")) {
								moleculeData[i] = "<" + Float.toString(data);
							}
							else {
								moleculeData[i] = ">" + Float.toString(data);
							}
						}
						else {
							moleculeData[i] = fieldData.substring(0, 1) + formatter.format(Float.parseFloat(fieldData.substring(1)));
						}
					}
				}				
			}
			ArrayList<MoleculeIndex> moleculesIndex = new ArrayList<MoleculeIndex>();
			int molIndex = molCounter; 
			if (!wholeMoleculesIndex.containsKey(molID)) {
				MMPFragmenter mmp = new MMPFragmenter(mol);
				moleculesFragmentsID.add(mmp.getMoleculeFragmentsID());
				List<MoleculeIndexID> moleculeIndexesID = mmp.getMoleculeIndexesID(false);
				for (MoleculeIndexID moleculeIndexID: moleculeIndexesID) {
					String[] keysID = moleculeIndexID.getKeysID();
					String valueID = moleculeIndexID.getValueID();
					int valueAtoms = moleculeIndexID.getValueIDAtoms();
					int key1Index = mmpUniqueFragments.addFragment(keysID[0]);
					int valueIndex = mmpUniqueFragments.addFragment(valueID);
					moleculeIndexID.setValueIndex(valueIndex);
					if (keysID.length == 1) { // single cut
						addValues(valueAtoms, Integer.toString(key1Index) + "\t", new int[]{valueIndex, molCounter}, false);
						moleculeIndexID.setKeysIndex(new int[]{key1Index});
					}
					else { // double cut
						int key2Index = mmpUniqueFragments.addFragment(keysID[1]);
						addValues(valueAtoms, Integer.toString(key1Index) + "\t" + Integer.toString(key2Index), new int[]{valueIndex, molCounter}, false);
						moleculeIndexID.setKeysIndex(new int[]{key1Index, key2Index});
					}
					mmpFragments.addFragments(molCounter, moleculeIndexID);
					if (moleculeIndexID.getValueIDAtoms() > maxValueAtoms) {
						maxValueAtoms = moleculeIndexID.getValueIDAtoms();
					}
				}
				molCounter++;
			}
			else {
				moleculesIndex = wholeMoleculesIndex.get(molID);
				molIndex = moleculesIndex.get(0).moleculeIndex;
			}
			if (moleculesIndex.size() > 0) {
				moleculesIndex.add(new MoleculeIndex(molIndex, moleculeName, moleculeData));
			}
			else {
				moleculesIndex.add(new MoleculeIndex(molIndex, molIDCoord, molID, moleculeName, moleculeData));
			}
			wholeMoleculesIndex.put(molID, moleculesIndex);
			moleculesRowCount++;
			if (verbose) {
				if (moleculesRowCount % 1000 == 0) {
					System.out.println("# " + moleculesRowCount);
				}
				else if (moleculesRowCount % 100 == 0) {
					System.out.print("#");
				}
				else if (moleculesRowCount % 10 == 0) {
					System.out.print(".");
				}
			}
		}
		// Getting percentiles
		if (verbose) {
			System.out.println(" " + moleculesRowCount);
			System.out.println(getDateAndTime() + ": getting percentiles...");
		}
		for (int i=0; i<fieldNames.length; i++) {
			if (fieldNumerics[i] != false && fieldDatas[i].size() > 0) {
				Collections.sort(fieldDatas[i]);
				int index = (int)Math.floor(0.05 * fieldDatas[i].size()); // I use floor so that I don't have to correct for the array indexes starting at 0
				if (Math.round(0.05f * fieldDatas[i].size()) != 0.05f * fieldDatas[i].size()) {
					fieldPercentiles5[i] = fieldDatas[i].get(index);
				}
				else {
					fieldPercentiles5[i] = (fieldDatas[i].get(index) + fieldDatas[i].get(index+1)) / 2.0f;
				}
				index = (int)Math.floor(0.95 * fieldDatas[i].size());
				if (Math.round(0.95f * fieldDatas[i].size()) != 0.95f * fieldDatas[i].size()) {
					fieldPercentiles95[i] = fieldDatas[i].get(index);
				}
				else {
					fieldPercentiles95[i] = (fieldDatas[i].get(index) + fieldDatas[i].get(index+1)) / 2.0f;
				}
			}
//			else if (fieldNumerics[i] != false && fieldDatas[i].size() == 0) {
//				fieldNumerics[i] = false;
//			}
		}
		// Processing Hydrogens replacements
		if (verbose)
			System.out.println(getDateAndTime() + ": processing hydrogen replacements...");
		for (List<String[]> fragmentsID:moleculesFragmentsID) {
			for (String[] fragmentID:fragmentsID) { // {R-group, "clean"}
				if (wholeMoleculesIndex.containsKey(fragmentID[1])) {
					int moleculeIndex = wholeMoleculesIndex.get(fragmentID[1]).get(0).moleculeIndex;
					int fragmentIndex = mmpUniqueFragments.addFragment(fragmentID[0]);
					boolean added = addValues(0, Integer.toString(fragmentIndex) + "\t", new int[]{r1HIndex, moleculeIndex}, true);
					if (added) {
						// to do: replace -1 by the correct value (is it needed?)
						MoleculeIndexID moleculeIndexID = new MoleculeIndexID(new String[]{fragmentID[0]}, new int[]{fragmentIndex}, r1H, r1HIndex, null, 0, new int[]{-1}, new int[]{-1});
						mmpFragments.addFragments(moleculeIndex, moleculeIndexID); // String[] keysID, String valueID, int[] keysIDAtoms, int valueIDAtoms
					}
				}
			}
		}
		// Generating MMPs
		if (verbose)
			System.out.print(getDateAndTime() + ": generating MMPs");
		maxValueAtoms += 1;
		int counter = 0;
		int[][] combinations;
		if (VERSION == "1.0") {
			combinations = new int[maxValueAtoms*(maxValueAtoms+1)/2-1][2];
			for (int i=1; i<maxValueAtoms; i++) {
				for (int j=0; j<=i; j++) {
					combinations[counter][0] = i;
					combinations[counter][1] = j;
					counter++;
				}
			}
		}
		else {
			combinations = new int[maxValueAtoms*maxValueAtoms][2];
			for (int i=0; i<maxValueAtoms; i++) {
				for (int j=0; j<maxValueAtoms; j++) {
					combinations[counter][0] = i;
					combinations[counter][1] = j;
					counter++;
				}
			}
		}
		if (verbose)
			System.out.println(" (" + combinations.length + " combinations)...");
		MMPEnumerator mMPEnumerator = new MMPEnumerator();
		counter = 0;
		for (int[] combination:combinations) {
			if (combination[0] == combination[1] && combination[0] != 0) {
				mMPEnumerator = new MMPEnumerator(combination, getIndex(combination[0]), null, VERSION);
			}
			else if (combination[0] != combination[1]) {
				mMPEnumerator = new MMPEnumerator(combination, getIndex(combination[0]), getIndex(combination[1]), VERSION);
			}
			HashMap<String, List<String[]>> mMPs = mMPEnumerator.getMMPEnumeration();
//			matchedMolecularPairs.addMMPs(mMPs);
			if (mMPs != null && mMPs.size() > 0)
				matchedMolecularPairs.writeMMPEnumeration(mMPs);
			counter++;
			if (verbose) {
				if (counter % 1000 == 0) {
					System.out.println("# " + counter);
				}
				else if (counter % 100 == 0) {
					System.out.print("#");
				}
				else if (counter % 10 == 0) {
					System.out.print(".");
				}
			}
		}
		compoundFileParser.close();
		if (verbose)
			System.out.println(" " + counter + "\n" + getDateAndTime() + ": done.");
  	}
	
	/**
	 * Writes the Molecules block. A moleculeIndex column has been added<br>
	 * since molecules might be not unique &rarr; The index won't be unique.
	 * @param printWriter
	 */
	private void writeMolecules(PrintWriter printWriter) {
		String line = "moleculeIndex\tidcoordinates2D\tmolecule\tmoleculeName";
		printWriter.println("<molecules>");
		printWriter.println("<column properties>");
		printWriter.println("<columnName=\"moleculeIndex\">");
		printWriter.println("<columnName=\"idcoordinates2D\">");
		printWriter.println("<columnProperty=\"specialType	idcoordinates2D\">");
		printWriter.println("<columnProperty=\"parent	molecule\">");
		printWriter.println("<columnName=\"molecule\">");
		printWriter.println("<columnProperty=\"specialType	idcode\">");
		printWriter.println("<columnName=\"moleculeName\">");
		for (int i=0; i<fieldNames.length; i++) {
			if (fieldNumerics[i]) {
				if (TRANSFORM_UM_TO_PIC50 && fieldNames[i].endsWith("_uM")) {
					fieldNames[i] = "p" + fieldNames[i].substring(0, fieldNames[i].length()-3);
				}
				String columnName = fieldNames[i];
				String longName = fieldNames[i];
				String category = "other";
				String[] items = fieldNames[i].split("\t", -1);
				if (items.length > 1) {
					category = items[0];
					columnName = items[1];
					longName = items[1];
					if (items.length == 3) {
						longName = items[2];
					}
				}
				printWriter.println("<columnName=\"" + columnName + "\" percentile5=\"" + Math.floor(fieldPercentiles5[i]*10.0)/10.0 + "\" percentile95=\"" + Math.ceil(fieldPercentiles95[i]*10.0)/10.0 + "\" longName=\"" + longName + "\" category=\"" + category + "\">");
				line += "\t" + columnName;
			}
		}
		printWriter.println("</column properties>");
		printWriter.println(line);
		Iterator<String> it = wholeMoleculesIndex.keySet().iterator();
		while (it.hasNext()) {
			String molID = it.next();
			String molIDCoord = "";
			List<MoleculeIndex> moleculesIndex = wholeMoleculesIndex.get(molID);
			for (MoleculeIndex moleculeIndex: moleculesIndex) {
				if (molIDCoord == "" && moleculeIndex.moleculeIDCoord != null) {
					molIDCoord = moleculeIndex.moleculeIDCoord;
				}
				line = Integer.toString(moleculeIndex.moleculeIndex) + "\t" + molIDCoord + "\t" + molID + "\t" + moleculeIndex.moleculeName;
				for (int i=0; i<moleculeIndex.moleculeData.length; i++) {
					if (fieldNumerics[i]) {
						if (moleculeIndex.moleculeData[i] != null) {
							line += "\t" + moleculeIndex.moleculeData[i];
						}
						else {
							line += "\t"; // otherwise the word 'null' is written!
						}
					}
				}
				printWriter.println(line);
			}
		}
		printWriter.println("</molecules>");
	}
	
	/**
	 * Writes the header (general information) block and calls the writing of the different blocks
	 * @param printWriter
	 * @throws IOException
	 */
	public void writeMMPFile(PrintWriter printWriter) throws IOException {
		printWriter.println("<matchedmolecularpairs-fileinfo>");
		printWriter.println("<version=\"" + VERSION + "\">");
		DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy");
		Date date = new Date();
		printWriter.println("<date=\"" + dateFormat.format(date) + "\">");
		printWriter.println("<dataset=\"" + datasetName + "\">");
		printWriter.println("<moleculesrowcount=" + Integer.toString(moleculesRowCount) + ">");
		printWriter.println("<mmpuniquefragmentsrowcount=" + Integer.toString(mmpUniqueFragments.getUniqueFragmentsCount()) + ">");
		printWriter.println("<mmpfragmentsrowcount=" + Integer.toString(mmpFragments.getFragmentsCount()) + ">");
		printWriter.println("<mmprowcount=" + Integer.toString(matchedMolecularPairs.getMMPsCount()) + ">");
		printWriter.println("<keysminatoms=\"" + Integer.toString(MMPFragmenter.KEYS_MIN_ATOMS) + "\">");
		printWriter.println("</matchedmolecularpairs-fileinfo>");
		// Process fragments (adding FP) and save them
		writeMolecules(printWriter);
		mmpUniqueFragments.writeUniqueFragments(printWriter);
		mmpFragments.writeFragments(printWriter);
		matchedMolecularPairs.writeMMPs(printWriter);
		printWriter.close();
	}
	
	/**
	 * Helper function for rounding
	 * @param f Input value
	 * @param decimalPlace
	 * @return
	 */
	private static float round(float f, int decimalPlace) {
		BigDecimal bd = new BigDecimal(Float.toString(f));
		bd = bd.setScale(decimalPlace, BigDecimal.ROUND_HALF_UP);
	    return bd.floatValue();
	}
}
