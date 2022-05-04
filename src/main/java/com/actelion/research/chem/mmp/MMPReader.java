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

import com.actelion.research.chem.AbstractDepictor;
import com.actelion.research.chem.ExtendedDepictor;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.mmp.MMP.MoleculeIndex;
import com.actelion.research.chem.mmp.MMPUniqueFragments.MMPUniqueFragment;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.gui.swing.SwingDrawContext;
import com.actelion.research.util.Base64;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;
import java.util.*;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MMPReader {
	public static final String SORT_BY_SIMILARITY = "similarity";
	public static final String SORT_BY_NUMBER_OF_EXAMPLES = "results";
	
	private HashMap<String, ArrayList<MoleculeIndex>> wholeMoleculesIndex;   // Molecule idCode, [[molName, molDatas], ..]
	private List<String> molecules;                                          // Ordered List of molecules idCodes; molIndex is used in mmpFragmentsIndex
	private List<DataField> dataFields;                                      // Ordered List of numerical data fieds
	private List<String> uniqueFragmentsIndex;                               // Ordered List of unique fragments idCodes
	private MMPUniqueFragments mmpUniqueFragments;                           // List of unique fragments (index, number of atoms, fingerprints)
	private HashMap<String, List<int[]>> mmpFragmentsIndex;                  // keys (tab-delimited) - {valueFragmentIndex, molIndex}
	private HashMap<Integer, HashMap<Integer, List<int[]>>> mmpIndex;        // MMP container: <fragmentIndex of the first fragment ('value1'), <size of the second fragment, List of second fragments ('value2') and examples>  
	private String datasetName;
	private String date;
	private Integer keysMinAtoms;
	private String version;
	private MMPPropertyCalculator mPropertyCalculator;
	
	public class DataField {
		private String fieldName;
		private String longFieldName;
		private String categoryName;
		private String percentile5;
		private String percentile95;
		
		private DataField(String fieldName, String longFieldName, String categoryName, String percentile5, String percentile95) {
			this.fieldName = fieldName;
			this.categoryName = categoryName;
			this.longFieldName = longFieldName;
			this.percentile5 = percentile5;
			this.percentile95 = percentile95;
		}
	}
	
	private class MatchedMolecularPairExamples {
		ArrayList<MoleculeIndex> example1;
		ArrayList<MoleculeIndex> example2;
		int similarity;
		
		public MatchedMolecularPairExamples(ArrayList<MoleculeIndex> example1, ArrayList<MoleculeIndex> example2, int similarity) {
			this.example1 = example1;
			this.example2 = example2;
			this.similarity = similarity;
		}
	}
	
	public static class MatchedMolecularPair {
		int numberOfExamples;
		String value1;
		int value1Index;
		int value1Atoms;
		String value2;
		int value2Index;
		int value2Atoms;
		int similarity;
		Reaction transformation;
		String transformationString;
		List<MatchedMolecularPairExamples> mmpExamples;
		ArrayList<ArrayList<Double>> datas;
		ArrayList<Integer> similarities;
		int[][] n;
		int[][] increase;
		int[][] decrease;
		int[][] neutral;
		Double[][] average;
		Double[][] sd;
		Integer[][] numberOfIncrease;
		Integer[][] numberOfDecrease;
		Integer[][] numberOfNeutral;
		boolean targetExists;
		
		/**
		 * Creates a new Matched Molecular Pair
		 * @param value1 seed 'value' (variable part of the molecule)
		 * @param value1Index index of the seed 'value'
		 * @param value1Atoms number of heavy atoms of the seed 'value'
		 * @param value1FP fingerprints of the seed 'value'
		 * @param value2 target 'value' (variable part of the molecule)
		 * @param value2Index index of the target 'value'
		 * @param value2Atoms number of heavy atoms of the target 'value'
		 * @param value2FP fingerprints of the target 'value'
		 * @param mmpExamples List of examples
		 * @param numberOfFields number of numerical fields
		 * @param targetExists true/false if the target exists (the transformed seed molecule exists in the data set)
		 */
		public MatchedMolecularPair(String value1, int value1Index, int value1Atoms, String[] value1FP, String value2, int value2Index, int value2Atoms, String[] value2FP, List<MatchedMolecularPairExamples> mmpExamples, int numberOfFields, int targetExists) {
			this.numberOfExamples = mmpExamples.size();
			this.value1 = value1;
			this.value1Index = value1Index;
			this.value1Atoms = value1Atoms;
			this.value2 = value2;
			this.value2Index = value2Index;
			this.value2Atoms = value2Atoms;
			this.transformation = reactionFromTwoValues(value1, value2);
			this.transformationString = null;
			this.mmpExamples = mmpExamples;
			this.datas = new ArrayList<ArrayList<Double>>();
			this.similarities = new ArrayList<Integer>();
			this.n = new int[numberOfFields][6];
			this.average = new Double[numberOfFields][6];
			this.sd = new Double[numberOfFields][6];
			this.increase = new int[numberOfFields][6];
			this.decrease = new int[numberOfFields][6];
			this.neutral = new int[numberOfFields][6];
			this.targetExists = targetExists > -1 ? true : false;
			
			for (int i=0; i<numberOfFields; i++) {
				datas.add(new ArrayList<Double>());
			}
			for (int i=0; i<numberOfExamples; i++) {
				ArrayList<MoleculeIndex> example1 = mmpExamples.get(i).example1;
				ArrayList<MoleculeIndex> example2 = mmpExamples.get(i).example2;
				int similarity = mmpExamples.get(i).similarity;
				Double[] example1Data = averageData(example1, numberOfFields);
				Double[] example2Data = averageData(example2, numberOfFields);
				similarities.add(similarity);
				for (int j=0; j<numberOfFields; j++) {
					if (example1Data[j] != null && example2Data[j] != null) {
						for (int k=0; k<6; k++) {
							if (similarity == -1 || similarity >= k) {
								this.n[j][k] += 1;
							}
						}
						ArrayList<Double> data = datas.get(j);
						data.add(example2Data[j] - example1Data[j]);
						datas.set(j, data);
					}
				}
			}
			for (int i=0; i<numberOfFields; i++) {
				for (int j=0; j<6; j++) {
					this.average[i][j] = calcAverage(datas.get(i).subList(0,  this.n[i][j]));
					if (this.n[i][j] > 1) {
						this.sd[i][j] = calcStanDev(datas.get(i).subList(0,  this.n[i][j]));
					}
					else {
						this.sd[i][j] = null;
					}
					int[] tendency = calcTendency(datas.get(i).subList(0,  this.n[i][j]));
					this.increase[i][j] = tendency[0];
					this.decrease[i][j] = tendency[1];
					this.neutral[i][j] = tendency[2];
				}
			}
			this.similarity = 0;
			for (int i=value1FP.length-1; i>=0; i--) {
				if (value1FP[i].equals(value2FP[i])) {
					this.similarity = i + 1;
					break;
				}
			}
		}
		
		/**
		 * Generates the transformation string, to be displayed in DataWarrior
		 */
		private void calcTransformationString() {
			this.transformationString = idCodeFromTwoValues(value1, value2);
		}
		
		/**
		 * Calculates the average from a list of data
		 * @param data List of data
		 * @return average
		 */
		private static Double calcAverage(List<Double> data) {
			if (data.size() > 0) {
				double allData = 0.0;
				int numberOfData = 0;
				for (Double d: data) {
					numberOfData++;
					allData += d;
				}
				return allData / numberOfData;
			}
			return null;
		}
		
		/**
		 * Calculates the standard deviation from a list of data
		 * @param data List of data
		 * @return standard deviation
		 */
		private static double calcStanDev(List<Double> data) {
			return Math.pow(calcVariance(data), 0.5); 
		} 
		
		
		/**
		 * Calculates the variance from a list of data
		 * @param data List of data
		 * @return variance
		 */
		private static double calcVariance(List<Double> data) {
			int n = data.size();
			double total = 0;
			double sTotal = 0;
			double scalar = 1/(double)(n-1); 
			for (Double d: data) { 
				total += d; 
				sTotal += Math.pow(d, 2); 
			} 
			return (scalar*(sTotal - (Math.pow(total, 2)/n))); 
		}
		
		/**
		 * Calculates the tendency from a list of data: number of increase, decrease, and neutral
		 * @param data List of data
		 * @return array of [number of increase, number of decrease, number of neutral]
		 */
		private static int[] calcTendency(List<Double> data) {
			int[] retVal = new int[3]; // increase, decrease, neutral
			Arrays.fill(retVal, 0);
			for (Double d: data) {
				if (d >= 0.1) {
					retVal[0] += 1;
				}
				else if (d <= -0.1) {
					retVal[1] += 1;
				}
				else {
					retVal[2] += 1;
				}
			}
			return retVal;
		}
		
		/**
		 * Generates the Reaction object from two seed and target fragments 
		 * @param value1 seed 'value' (variable part of the molecule)
		 * @param value2 target 'value' (variable part of the molecule)
		 * @return Reaction object
		 */
		private Reaction reactionFromTwoValues(String value1, String value2) {
			StereoMolecule mol1 = new StereoMolecule();
			StereoMolecule mol2 = new StereoMolecule();
			IDCodeParser idCodeParser = new IDCodeParser();
			idCodeParser.parse(mol1, value1);
			idCodeParser.parse(mol2, value2);
			Reaction rxn = new Reaction(new StereoMolecule[]{mol1, mol2}, 1);
			return rxn;
		}
		
		/**
		 * Generates the idCode from the transformation of one seed 'value' to one target 'value', to be displayed in DataWarrior
		 * @param value1 seed 'value' (variable part of the molecule)
		 * @param value2 target 'value' (variable part of the molecule)
		 * @return idCode of the transformation
		 */
		private String idCodeFromTwoValues(String value1, String value2) {
			StereoMolecule mol1 = new StereoMolecule();
			StereoMolecule mol2 = new StereoMolecule();
			IDCodeParser idCodeParser = new IDCodeParser();
			idCodeParser.parse(mol1, value1);
			idCodeParser.parse(mol2, value2);
			Reaction rxn = new Reaction(new StereoMolecule[]{mol1, mol2}, 1);
			String[] rxnEncoder = ReactionEncoder.encode(rxn, false);
			return value1 + "!" + value2 + "##" + rxnEncoder[2];
//			return rxnEncoder[0] + "##" + rxnEncoder[2];
//			mol1.addSubstituent(mol2, -1);
//			return mol1.getIDCode();
		}
		
		/**
		 * In case of several compounds having the same structure, averages data between these compounds
		 * @param example List of data
		 * @param numberOfFields number of numerical fields
		 * @return Array of averages
		 */
		private Double[] averageData(ArrayList<MoleculeIndex> example, int numberOfFields) {
			int[] numberOfValues = new int[numberOfFields];
			Double[] values = new Double[numberOfFields];
			Arrays.fill(numberOfValues, 0);
			Arrays.fill(values, 0.0);
			for (int i=0; i<numberOfFields; i++) {
				for (MoleculeIndex moleculeIndex: example) {
					if (moleculeIndex.moleculeData != null && isNumeric(moleculeIndex.moleculeData[i])) {
						numberOfValues[i] += 1;
						values[i] += Double.parseDouble(moleculeIndex.moleculeData[i]);
					}
				}
			}
			for (int i=0; i<numberOfFields; i++) {
				if (numberOfValues[i] != 0) {
					values[i] /= numberOfValues[i];
				}
				else {
					values[i] = null;
				}
			}
			return values;
		}
		
		/**
		 * Helper function
		 * @param str input String
		 * @return true/false if the input string is numerical
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
	}

	public static final Comparator<MatchedMolecularPair> NUMBER_OF_EXAMPLES_SORT = new Comparator<MatchedMolecularPair>() {
		public int compare(MatchedMolecularPair matchedMolecularPair1, MatchedMolecularPair matchedMolecularPair2) {
			return matchedMolecularPair2.numberOfExamples - matchedMolecularPair1.numberOfExamples;
		}
	};

	public static final Comparator<MatchedMolecularPair> SIMILARITY_SORT = new Comparator<MatchedMolecularPair>() {
		public int compare(MatchedMolecularPair matchedMolecularPair1, MatchedMolecularPair matchedMolecularPair2) {
			return matchedMolecularPair2.similarity - matchedMolecularPair1.similarity;
		}
	};
	
	public static final Comparator<MatchedMolecularPairExamples> EXAMPLES_SIMILARITY_SORT = new Comparator<MatchedMolecularPairExamples>() {
		public int compare(MatchedMolecularPairExamples matchedMolecularPairExample1, MatchedMolecularPairExamples matchedMolecularPairExample2) {
			return matchedMolecularPairExample2.similarity - matchedMolecularPairExample1.similarity;
		}
	};
	
	public MMPReader(BufferedReader br, boolean verbose) throws IOException, Exception {
		mmpFragmentsIndex = new HashMap<String, List<int[]>>();
		mmpIndex = new HashMap<Integer, HashMap<Integer, List<int[]>>>();
		mmpUniqueFragments = new MMPUniqueFragments();
		mPropertyCalculator = new MMPPropertyCalculator();
		readMMPFile(br, verbose);
		br.close();
	}
	
	/**
	 * Reads the header block of a MMP file
	 * @param br
	 * @throws IOException
	 * @throws Exception
	 */
	private void readMMPFile(BufferedReader br, boolean verbose) throws IOException, Exception {
		HashMap<String, Integer> rowCounts = new HashMap<String, Integer>();
		int rowCountsCounter = 0;
		String strLine;
		Pattern pattern1 = Pattern.compile("<(.*?)=\"(.*?)\">");
		Pattern pattern2 = Pattern.compile("<(.*?rowcount)=([0-9]*?)>");
		while ((strLine = br.readLine()) != null && rowCountsCounter < 4) {
			Matcher matcher1 = pattern1.matcher(strLine);
			if (matcher1.find()) {
				if (matcher1.group(1).equals("dataset")) {
					datasetName = matcher1.group(2);
				}
				else if (matcher1.group(1).equals("date")) {
					date = matcher1.group(2);
				}
				else if (matcher1.group(1).equals("keysminatoms")) {
					keysMinAtoms = Integer.parseInt(matcher1.group(2));
				}
				else if (matcher1.group(1).equals("version")) {
					version = matcher1.group(2);
				}
			}
			else {
				Matcher matcher2 = pattern2.matcher(strLine);
				if (matcher2.find()) {
					rowCounts.put(matcher2.group(1), Integer.parseInt(matcher2.group(2)));
					rowCountsCounter++;
				}
			}
		}
		if (rowCountsCounter < 4) {
			throw new IOException("General: cannot find the four rowcount lines");
		}
		if (verbose) {
			System.out.println("The dataset contains " + rowCounts.get("moleculesrowcount") + " molecules, " + rowCounts.get("mmpuniquefragmentsrowcount") + " unique fragments, " + rowCounts.get("mmpfragmentsrowcount") + " molecules fragments combinations and " + rowCounts.get("mmprowcount") + " MMPs.");
			System.out.println("                  0        10        20        30        40        50        60        70        80        90       100");
		}
		readMolecules(br, rowCounts.get("moleculesrowcount"), verbose);
		readUniqueFragments(br, rowCounts.get("mmpuniquefragmentsrowcount"), verbose);
		readFragments(br, rowCounts.get("mmpfragmentsrowcount"), verbose);
		readMMPs(br, rowCounts.get("mmprowcount"), verbose);
	}
	
	/**
	 * Reads the Molecules block of a MMP file
	 * @param br
	 * @param rowCount Number of expected rows (from the header block)
	 * @throws IOException
	 * @throws Exception
	 */
	private void readMolecules(BufferedReader br, int rowCount, boolean verbose) throws IOException, Exception {
		wholeMoleculesIndex = new LinkedHashMap<String, ArrayList<MoleculeIndex>>();
		dataFields = new ArrayList<DataField>();
		molecules = new ArrayList<String>();
		keysMinAtoms = MMPFragmenter.KEYS_MIN_ATOMS;
		if (verbose)
			System.out.print("Molecules:        #");
		try {
			String strLine;
			boolean moleculesBlock = false;
			int linesToRead = 11;
			int entries = 0;
			int lastEntryIndex = -1;
			Pattern tagPattern = Pattern.compile("<(columnName=.*?)>");
		    Pattern attValue = Pattern.compile("(\\w+)=\"(.*?)\"");
			while ((strLine = br.readLine()) != null && entries < rowCount) {
				if (strLine.startsWith("<molecules>") || (moleculesBlock == true && linesToRead > 0)) {
					moleculesBlock = true;
					Matcher matcher = tagPattern.matcher(strLine);
					if (matcher.find()) {
						matcher = attValue.matcher(matcher.group(1));
						String fieldName = null;
						String longFieldName = null;
						String categoryName = null;
						String percentile5 = null;
						String percentile95 = null;
						while (matcher.find()) {
							if (matcher.group(1).equals("columnName")) {
								fieldName = matcher.group(2);
							}
							else if (matcher.group(1).equals("longName")) {
								longFieldName = matcher.group(2);
							}
							else if (matcher.group(1).equals("percentile5")) {
								percentile5 = matcher.group(2);
							}
							else if (matcher.group(1).equals("percentile95")) {
								percentile95 = matcher.group(2);
							}
							else if (matcher.group(1).equals("category")) {
								categoryName = matcher.group(2);
							}
						}
						if (fieldName != null && !fieldName.equals("moleculeIndex") && !fieldName.equals("idcoordinates2D") && !fieldName.equals("molecule") && !fieldName.equals("moleculeName")) {
							dataFields.add(new DataField(fieldName, longFieldName, categoryName, percentile5, percentile95));
							linesToRead++;
						}
					}
					linesToRead--;
				}
				else if (moleculesBlock == true && linesToRead == 0) {
					String[] items = strLine.split("\t", -1);
					if (items.length == dataFields.size() + 4) {
						items = strLine.split("\t", 5); // index, coordinates, idcode, name
						String[] data = items[4].split("\t", -1);
						int molIndex = Integer.parseInt(items[0]);
						MoleculeIndex moleculeIndex;
						ArrayList<MoleculeIndex> moleculesIndex = new ArrayList<MoleculeIndex>();
						if (molIndex == lastEntryIndex) {
							moleculeIndex = new MoleculeIndex(molIndex, items[3], data);
							moleculesIndex = wholeMoleculesIndex.get(items[2]);
						}
						else {
							moleculeIndex = new MoleculeIndex(molIndex, items[1], items[2], items[3], data);
							molecules.add(items[2]);
							if (molecules.size() != Integer.parseInt(items[0])+1) {
								System.out.println(molecules.size());
							}
						}
						moleculesIndex.add(moleculeIndex);
						wholeMoleculesIndex.put(items[2], moleculesIndex);
						printProgress(verbose, rowCount, entries);
						entries++;
						lastEntryIndex = molIndex;
					}
					else if (strLine.startsWith("</molecules>")) {
						throw new IOException("molecules: Bad number of entries");
					}
				}
			}
		}
		catch (IOException ioe) {
		}
		if (verbose)
			System.out.print("\n");
	}
	
	/**
	 * Reads the Unique Fragments block of a MMP file
	 * @param br
	 * @param rowCount Number of expected rows (from the header block)
	 */
	private void readUniqueFragments(BufferedReader br, int rowCount, boolean verbose) {
		uniqueFragmentsIndex = new ArrayList<String>(rowCount);
		if (verbose)
			System.out.print("Unique Fragments: #");
		try {
			String strLine;
			boolean mmpUniqueFragmentsBlock = false;
			int linesToRead = 17;
			int entries = 0;
			while ((strLine = br.readLine()) != null && entries < rowCount) {
				if (strLine.startsWith("<mmpUniqueFragments>") || (mmpUniqueFragmentsBlock == true && linesToRead > 0)) {
					mmpUniqueFragmentsBlock = true;
					linesToRead--;
				}
				else if (mmpUniqueFragmentsBlock == true && linesToRead == 0) {
					String[] items = strLine.split("\t", -1);
					if (items.length == 7) {
						uniqueFragmentsIndex.add(items[0]);
						mmpUniqueFragments.addFragment(items[0], Integer.parseInt(items[1]), new String[]{items[2], items[3], items[4], items[5], items[6]});
						printProgress(verbose, rowCount, entries);
						entries++;
					}
					else if (strLine.startsWith("</mmpUniqueFragments>")) {
						throw new IOException("mmpUniqueFragments: Bad number of entries");
					}
				}
			}
		}
		catch (IOException ioe) {
		}
		if (verbose)
			System.out.print("\n");
	}
	
	/**
	 * Reads the Fragments block of a MMP file
	 * @param br
	 * @param rowCount Number of expected rows (from the header block)
	 */	
	private void readFragments(BufferedReader br, int rowCount, boolean verbose)  {
		if (verbose)
			System.out.print("Fragments:        #");
		try {
			String strLine;
			boolean mmpFragmentsBlock = false;
			int linesToRead = 9;
			int entries = 0;
			while ((strLine = br.readLine()) != null && entries < rowCount) {
				if (strLine.startsWith("<mmpFragments>") || (mmpFragmentsBlock == true && linesToRead > 0)) {
					mmpFragmentsBlock = true;
					linesToRead--;
				}
				else if (mmpFragmentsBlock == true && linesToRead == 0) {
					String[] items = strLine.split("\t", -1);
					if (items.length == 5) {
						if (items[3].equals("1")) { // cutType
							addFragment(items[0], new int[]{Integer.parseInt(items[2]), Integer.parseInt(items[4])});
							// This is to index also {key-value} for smaller keys (used later to sort by similarity) that wouldn't be otherwise indexed
							if (mmpUniqueFragments.getFragmentAtoms(uniqueFragmentsIndex.get(Integer.parseInt(items[2]))) < keysMinAtoms) {
								addFragment(items[2], new int[]{Integer.parseInt(items[0]), Integer.parseInt(items[4])});
							}
						}
						else {
							addFragment(items[0] + "\t" + items[1], new int[]{Integer.parseInt(items[2]), Integer.parseInt(items[4])});
						}
						printProgress(verbose, rowCount, entries);
						entries++;
					}
				}
				else if (strLine.startsWith("</mmpFragments>")) {
					throw new IOException("mmpFragments: Bad number of entries");
				}
			}
		}
		catch (IOException ioe) {
		}
		if (verbose)
			System.out.print("\n");
	}
	
	/**
	 * Reads the Matched Molecular Pairs block of a MMP file
	 * @param br
	 * @param rowCount Number of expected rows (from the header block)
	 */
	private void readMMPs(BufferedReader br, int rowCount, boolean verbose) {
		if (verbose)
			System.out.print("MMPs:             #");
		try {
			String strLine;
			boolean mmpBlock = false;
			int linesToRead = 11;
			int entries = 0;
			int value1Atoms = -1;
			HashMap<Integer, HashMap<Integer, List<int[]>>> tempMMPIndex = null;
			while ((strLine = br.readLine()) != null && entries < rowCount)   {
				if (strLine.startsWith("<matchedMolecularPairs>") || (mmpBlock == true && linesToRead > 0)) {
					mmpBlock = true;
					linesToRead--;
				}
				else if (mmpBlock == true && linesToRead == 0) {
					String[] items = strLine.split("\t", -1);
					if (items.length == 7) {
						if (version == "1.1") {
							int val1Atoms = Integer.parseInt(items[1]);
							if (val1Atoms != value1Atoms) {
								if (tempMMPIndex != null) {
									mmpIndex.putAll(tempMMPIndex);
								}
								tempMMPIndex = new HashMap<Integer, HashMap<Integer, List<int[]>>>();
								value1Atoms = val1Atoms;
							}
							tempMMPIndex = addTempMMP(tempMMPIndex, Integer.parseInt(items[0]), Integer.parseInt(items[3]), Integer.parseInt(items[2]), items[6].split("\\|", -1));
						}
						else {
							addMMP(Integer.parseInt(items[0]), Integer.parseInt(items[3]), Integer.parseInt(items[2]), items[6].split("\\|", -1));
						}
						printProgress(verbose, rowCount, entries);
						entries++;
					}
				}
				else if (strLine.startsWith("</matchedMolecularPairs>")) {
					throw new IOException("matchedMolecularPairs: Bad number of entries");
				}
			}
			if (tempMMPIndex != null) {
				mmpIndex.putAll(tempMMPIndex);
			}
		}
		catch (IOException ioe) {
		}
		if (verbose)
			System.out.print("\n");
	}
	
	/**
	 * In verbose mode, print the process of reading the MMP file
	 * @param verbose
	 * @param rowCount Total number of rows
	 * @param entries Number of read entries
	 */
	private void printProgress(boolean verbose, int rowCount, int entries) {
		if (verbose) {
			double percentage = (entries+1) * 100.0 / rowCount;
			double one = 100.0 / rowCount;
			if (Math.floor(percentage) != Math.floor(percentage - one)) {
				if (Math.floor(percentage) % 10.0 == 0) {
					System.out.print("#");
				}
				else {
					System.out.print(".");
				}
			}
		}
	}

	/**
	 * Adds a new fragment
	 * @param keys tab-delimited keys (one for single cut, two for double cut)
	 * @param data [valueFragmentIndex, molIndex]
	 */
	private void addFragment(String keys, int[] data) {
		List<int[]> datas = new ArrayList<int[]>();
		if (mmpFragmentsIndex.containsKey(keys)) {
			datas = mmpFragmentsIndex.get(keys);
		}
		datas.add(data);
		mmpFragmentsIndex.put(keys, datas);
	}
	
	/**
	 * Adds a new Matched Molecular Pair (for version 1.1)
	 * @param tempMMPIndex temporary container for the MMPs
	 * @param value1FragmentIndex 'seed' fragment index
	 * @param value2Atoms 'target' number of heavy atoms
	 * @param value2
	 * @param examples List of examples; the first item is the target 'value' idCode
	 * @return the temporary container with the newly added MMP
	 */
	private HashMap<Integer, HashMap<Integer, List<int[]>>> addTempMMP(HashMap<Integer, HashMap<Integer, List<int[]>>> tempMMPIndex, int value1FragmentIndex, int value2Atoms, int value2, String[] examples) {
		int[] val2_examples = new int[examples.length * 2 + 1];
		val2_examples[0] = value2;
		int counter = 1;
		for (String value2AndExample: examples) {
			String[] items = value2AndExample.split(",");
			val2_examples[counter] = Integer.parseInt(items[0]);
			val2_examples[counter+1] = Integer.parseInt(items[1]);
			counter += 2;
		}
		HashMap<Integer, List<int[]>> values2 = new HashMap<Integer, List<int[]>>();
		List<int[]> values2OfSizeX = new ArrayList<int[]>();
		if (tempMMPIndex.containsKey(value1FragmentIndex)) {
			values2 = tempMMPIndex.get(value1FragmentIndex);
			if (values2.containsKey(value2Atoms)) {
				values2OfSizeX = values2.get(value2Atoms);
			}
		}
		values2OfSizeX.add(val2_examples);
		values2.put(value2Atoms, values2OfSizeX);
		tempMMPIndex.put(value1FragmentIndex, values2);
		return tempMMPIndex;
	}
	
	/**
	 * Adds a new Matched Molecular Pair (for version 1.0)
	 * @param value1FragmentIndex 'seed' fragment index
	 * @param value2Atoms 'target' number of heavy atoms
	 * @param value2
	 * @param examples List of examples; the first item is the target 'value' idCode
	 */
	private void addMMP(int value1FragmentIndex, int value2Atoms, int value2, String[] examples) {
		int[] val2_examples = new int[examples.length * 2 + 1];
		val2_examples[0] = value2;
		int counter = 1;
		for (String value2AndExample: examples) {
			String[] items = value2AndExample.split(",");
			val2_examples[counter] = Integer.parseInt(items[0]);
			val2_examples[counter+1] = Integer.parseInt(items[1]);
			counter += 2;
		}
		HashMap<Integer, List<int[]>> values2 = new HashMap<Integer, List<int[]>>();
		List<int[]> values2OfSizeX = new ArrayList<int[]>();
		if (mmpIndex.containsKey(value1FragmentIndex)) {
			values2 = mmpIndex.get(value1FragmentIndex);
			if (values2.containsKey(value2Atoms)) {
				values2OfSizeX = values2.get(value2Atoms);
			}
		}
		values2OfSizeX.add(val2_examples);
		values2.put(value2Atoms, values2OfSizeX);
		mmpIndex.put(value1FragmentIndex, values2);
	}
	
	/**
	 * Returns the fragment index from one fragment idCode
	 * @param fragment idCode of a fragment
	 * @return fragment index
	 */
	public Integer fragmentToFragmentIndex(String fragment) {
		if (uniqueFragmentsIndex.contains(fragment)) {
			return uniqueFragmentsIndex.indexOf(fragment);
		}
		return null;
	}

	/**
	 * Return the fragment indexes from one (single cut) or several fragment idCodes  
	 * @param fragments List of fragment idCodes
	 * @return List of fragment indexes
	 */
	public Integer[] fragmentToFragmentIndex(String[] fragments) {
		Integer retVal[] = new Integer[fragments.length];
		for (int i=0; i<fragments.length; i++) {
			retVal[i] = fragmentToFragmentIndex(fragments[i]);
		}
		return retVal;
	}
	
	/**
	 * Returns the number of heavy atoms of one fragment
	 * @param fragment idCode of a fragment
	 * @return number of heavy atoms
	 */
	public Integer fragmentToFragmentSize(String fragment) {
		return mmpUniqueFragments.getFragmentAtoms(fragment);
	}
	
	/**
	 * Returns the size of the chemical space corresponding to a specific 'key' (constant part of a molecule)
	 * @param key idCode of the 'key' (constant part of a molecule)
	 * @return the size of the chemical space
	 */
	public int getChemicalSpaceSize(String key) {
		return getChemicalSpaceSize(new String[]{key});
	}
	
	/**
	 * Returns the size of the chemical space corresponding to one constant part of a molecule
	 * @param keys array of one (single cut) or two (double cut) IDCodes of the 'keys' (constant part of the molecule)  
	 * @return the size of the chemical space.
	 */
	public int getChemicalSpaceSize(String[] keys) {
		int chemicalSpaceSize = 0;
		Integer[] keysIndex = fragmentToFragmentIndex(keys);
		String keysString = keysToKeysString(keysIndex);
		if (keysString != null && mmpFragmentsIndex.containsKey(keysString)) {
			List<int[]> chemicalSpace = mmpFragmentsIndex.get(keysString);
			Set<Integer> molList = new HashSet<Integer>();
			for (int[] chemSpace: chemicalSpace) {
				molList.add(chemSpace[1]);
			}
			// TODO: loop through molList and count "real" number of molecules (i.e. same structure, different names)
			chemicalSpaceSize = molList.size(); 
		}
		return chemicalSpaceSize;
	}
	
	/**
	 * Returns the list of idCode & molecules names representing the chemical space of one constant part of a molecule
	 * @param keys array of one (single cut) or two (double cut) IDCodes of the 'keys' (constant part of the molecule)
	 * @param value the variable part of the molecule; not used yet but might be used to identify the current compound
	 * @return an array of tab-delimited idCode, idCoord, molecule names [and data]
	 */
	public List<String> getChemicalSpace(String[] keys, String value, String dataField) {
		List<String> chemicalSpaceMolecules = new ArrayList<String>();
		List<int[]> chemicalSpace = new ArrayList<int[]>();
		Integer[] keysIndex = fragmentToFragmentIndex(keys);
		String keysString = keysToKeysString(keysIndex);
		int dataFieldIndex = -1;
		if (dataField != null) {
			for (int i=0; i<dataFields.size(); i++) {
				if (dataFields.get(i).fieldName.equals(dataField) || dataFields.get(i).longFieldName.equals(dataField)) {
					dataFieldIndex = i;
					break;
				}
			}
		}
		if (keysString != null && mmpFragmentsIndex.containsKey(keysString)) {
			chemicalSpace = mmpFragmentsIndex.get(keysString);
			Set<Integer> molList = new HashSet<Integer>();
			for (int[] chemSpace: chemicalSpace) {
				molList.add(chemSpace[1]);
			}
			for (Integer molIndex: molList) {				
				String idCode = molecules.get(molIndex);
				String idCoord = null;
				for (MoleculeIndex moleculeIndex: wholeMoleculesIndex.get(idCode)) {
					if (idCoord == null) {
						idCoord = moleculeIndex.moleculeIDCoord;
					}
					if (dataFieldIndex == -1) {
						chemicalSpaceMolecules.add(idCode + "\t" + idCoord + "\t" + moleculeIndex.moleculeName);
					}
					else {
						chemicalSpaceMolecules.add(idCode + "\t" + idCoord + "\t" + moleculeIndex.moleculeName + "\t" + moleculeIndex.moleculeData[dataFieldIndex]);
					}
				}
			}
		}
//		organize the array to put the current molecule on top
//		for (int i=0; i<chemicalSpace.size(); i++) {
//			if (chemicalSpace.get(i)[0] == value) {
//				Collections.swap(chemicalSpace, i, 0);
//			}
//		}
		return chemicalSpaceMolecules;
	}
	
	/**
	 * Generates the DWAR file for the chemical space for a specific 'key' and data field
	 * @param moleculeIDCode idCode of the seed molecule
	 * @param keys Array of one (single cut) or two (double cut) 'keys' idCodes (constant part of the molecule)
	 * @param dataField null or short name of a data field
	 * @return the whole string of the generated DWAR file
	 */
	public String getChemicalSpaceDWAR(String moleculeIDCode, String[] keys, String dataField) {
		int dataFieldIndex = -1;
		if (dataField != null) {
			for (int i=0; i<dataFields.size(); i++) {
				if (dataFields.get(i).fieldName.equals(dataField) || dataFields.get(i).longFieldName.equals(dataField)) {
					dataFieldIndex = i;
					dataField = dataFields.get(i).fieldName;
					break;
				}
			}
		}
		StringBuilder dWAR = new StringBuilder();
		dWAR.append("<datawarrior-fileinfo>\n");
		dWAR.append("<version=\"3.1\">\n");
		List<String> chemicalSpace = getChemicalSpace(keys, null, dataField);
		dWAR.append("<rowcount=\"" + chemicalSpace.size() + "\">\n");
		dWAR.append("</datawarrior-fileinfo>\n");
		dWAR.append("<column properties>\n");
		dWAR.append("<columnName=\"Structure\">\n");
		dWAR.append("<columnProperty=\"specialType	idcode\">\n");
		dWAR.append("<columnName=\"idcoordinates2D\">\n");
		dWAR.append("<columnProperty=\"specialType	idcoordinates2D\">\n");
		dWAR.append("<columnProperty=\"parent	Structure\">\n");
		dWAR.append("</column properties>\n");
		if (dataField != null) {
			dWAR.append("Structure\tidcoordinates2D\tActelion No\t" + dataField + "\n");
		}
		else {
			dWAR.append("Structure\tidcoordinates2D\tActelion No\n");
		}
		for (String chemSpace: chemicalSpace) {
			dWAR.append(chemSpace + "\n");
		}
		dWAR.append("<datawarrior properties>\n");
		if (dataFieldIndex != -1) {
			dWAR.append("<columnDescriptionCount=\"1\">\n");
			dWAR.append("<columnDescription_0=\"" + dataField + "\t" + dataFields.get(dataFieldIndex).longFieldName + "\">\n");
			dWAR.append("<colorMaxBackground__TableStructure=\"" + dataFields.get(dataFieldIndex).percentile95 + "\">\n");
			dWAR.append("<colorMinBackground__TableStructure=\"" + dataFields.get(dataFieldIndex).percentile5 + "\">\n");
			dWAR.append("<colorBackground__TableStructure_0=\"-13395457\">\n");
			dWAR.append("<colorBackground__TableStructure_1=\"-39322\">\n");
			dWAR.append("<colorColumnBackground__TableStructure=\"" + dataField + "\">\n");
			dWAR.append("<colorListModeBackground__TableStructure=\"straight\">\n");
		}
		dWAR.append("<columnWidth_Table_Structure=\"150\">\n");
		dWAR.append("<detailView=\"height[Data]=0.7;height[Structure]=0.3\">\n");
		dWAR.append("<filter0=\"#browser#	<disabled>\">\n");
		dWAR.append("<filter1=\"#structure#	Structure\">\n");
		
		dWAR.append("<filter2=\"#double#	" + dataField + "\">\n");
		dWAR.append("<mainViewCount=\"2\">\n");
		dWAR.append("<mainViewDockInfo0=\"root\">\n");
		dWAR.append("<mainViewDockInfo1=\"Table	right	0.5\">\n");
		dWAR.append("<mainViewName0=\"Table\">\n");
		dWAR.append("<mainViewName1=\"Structure\">\n");
		dWAR.append("<mainViewType0=\"tableView\">\n");
		dWAR.append("<mainViewType1=\"structureView\">\n");
		dWAR.append("<rightSplitting=\"0.6\">\n");
		dWAR.append("<rowHeight_Table=\"80\">\n");
		dWAR.append("<structureGridColumn_Structure=\"Structure\">\n");
		dWAR.append("<structureGridColumns_Structure=\"4\">\n");
		dWAR.append("</datawarrior properties>\n");
		return dWAR.toString();
	}
	
	/**
	 * Generates the DWAR file for the Matched Molecular Pairs corresponding to a specific seed 'value', target 'value', and list of data fields
	 * @param moleculeIDCode idCode of the seed molecule
	 * @param keys Array of one or two 'keys' idCodes (constant part of the molecule)
	 * @param value1 seed 'value' (variable part of the molecule)
	 * @param value2 target 'value' (variable part of the molecule)
	 * @param replacementSize size of the replacement (number of heavy atoms)
	 * @param properties List of data fields for which data should be retrieved
	 * @return the whole string of the generated DWAR file
	 */
	public String getMMPsDWAR(String moleculeIDCode, String[] keys, String value1, String value2, int replacementSize, List<String> properties) {
		List<MatchedMolecularPair> transformations = getTransformations(moleculeIDCode, keys, value1, replacementSize, replacementSize, null);
		NumberFormat formatter = new DecimalFormat("#.##");
		StringBuilder dWAR = new StringBuilder();
		dWAR.append("<datawarrior-fileinfo>\n");
		dWAR.append("<version=\"3.1\">\n");
		for (MatchedMolecularPair transformation: transformations) {
			if (transformation.value2.equals(value2)) {
				dWAR.append("<rowcount=\"" + transformation.numberOfExamples + "\">\n");
				break;
			}
		}
		dWAR.append("</datawarrior-fileinfo>\n");
		dWAR.append("<column properties>\n");
		dWAR.append("<columnName=\"Structure (1)\">\n");
		dWAR.append("<columnProperty=\"specialType	idcode\">\n");
		dWAR.append("<columnName=\"idcoordinates2D (1)\">\n");
		dWAR.append("<columnProperty=\"specialType	idcoordinates2D\">\n");
		dWAR.append("<columnProperty=\"parent	Structure (1)\">\n");
		dWAR.append("<columnName=\"Structure (2)\">\n");
		dWAR.append("<columnProperty=\"specialType	idcode\">\n");
		dWAR.append("<columnName=\"idcoordinates2D (2)\">\n");
		dWAR.append("<columnProperty=\"specialType	idcoordinates2D\">\n");
		dWAR.append("<columnProperty=\"parent	Structure (2)\">\n");
		dWAR.append("<columnName=\"FragFp 2\">\n");
		dWAR.append("<columnProperty=\"specialType	FragFp\">\n");
		dWAR.append("<columnProperty=\"parent	Structure (2)\">\n");
		dWAR.append("<columnProperty=\"version	1.2.1\">\n");
		dWAR.append("</column properties>\n");
		dWAR.append("FragFp 2	Structure (1)\tidcoordinates2D (1)\tStructure (2)\tidcoordinates2D (2)\tActelion No (1)\tActelion No (2)\tSimilarity");
		for (String property: properties) {
			dWAR.append("\t" + property + " (1)\t" + property + " (2)\t" + property + " (delta)");
		}
		dWAR.append("\n");
		for (MatchedMolecularPair transformation: transformations) {
			if (transformation.value2.equals(value2)) {
				// TODO: handling of multiple compounds with same structure (get(0), get(1)...)
				for (int j=0; j<transformation.numberOfExamples; j++) {
					MoleculeIndex example1 = transformation.mmpExamples.get(j).example1.get(0);
					MoleculeIndex example2 = transformation.mmpExamples.get(j).example2.get(0);
					dWAR.append("\t" + example1.moleculeIDCode + "\t" + example1.moleculeIDCoord + "\t" + example2.moleculeIDCode + "\t" + example2.moleculeIDCoord + "\t" + example1.moleculeName + "\t" + example2.moleculeName + "\t" + Integer.toString(transformation.similarities.get(j)));
					for (String property: properties) {
						int fieldIndex = -1;
						for (int i=0; i<dataFields.size(); i++) {
							if (dataFields.get(i).fieldName.equals(property) || dataFields.get(i).longFieldName.equals(property)) {
								fieldIndex = i;
								break;
							}
						}
						if (fieldIndex != -1) {
							String data1 = example1.moleculeData[fieldIndex];
							String data2 = "";
							if (example2.moleculeData != null) {
								data2 = example2.moleculeData[fieldIndex]; // virtual compound
							}
							String delta = "";
							if (!data1.equals("") && !data1.startsWith(">") && !data1.startsWith("<") && !data2.equals("") && !data2.startsWith(">") && !data2.startsWith("<")) {
								delta = formatter.format(Double.parseDouble(data2) - Double.parseDouble(data1));
							}
							dWAR.append("\t" + data1 + "\t" + data2 + "\t" + delta);
						}
					}
					dWAR.append("\n");
				}
			}
		}
		dWAR.append("<datawarrior properties>\n");
		if (properties.size() > 0) {
			dWAR.append("<axisColumn_2D View_0=\"" + properties.get(0) + " (1)\">\n");
			dWAR.append("<axisColumn_2D View_1=\"" + properties.get(0) + " (2)\">\n");
		}
		dWAR.append("<chartType_2D View=\"scatter\">\n");
		dWAR.append("<columnDescriptionCount=\"" + (properties.size()*3+1) + "\">\n");
		dWAR.append("<columnDescription_0=\"Similarity	Similarity of the local environment. 6: seed compound; 0-5: number of atoms identical to the seed compound; -1: key not found\">\n");
		int descriptionCounter = 1;
		for (String property: properties) {
			int fieldIndex = -1;
			for (int i=0; i<dataFields.size(); i++) {
				if (dataFields.get(i).fieldName.equals(property) || dataFields.get(i).longFieldName.equals(property)) {
					fieldIndex = i;
					break;
				}
			}
			if (fieldIndex != -1) {
				dWAR.append("<columnDescription_" + descriptionCounter + "=\"" + property + " (1)	First compound " + dataFields.get(fieldIndex).longFieldName + "\">\n");
				dWAR.append("<columnDescription_" + (descriptionCounter+1) + "=\"" + property + " (2)	Second compound " + dataFields.get(fieldIndex).longFieldName + "\">\n");
				dWAR.append("<columnDescription_" + (descriptionCounter+2) + "=\"" + property + " (delta)	Difference of " + dataFields.get(fieldIndex).longFieldName + "\">\n");
				descriptionCounter += 3;
			}
		}
		dWAR.append("<columnWidth_Table_Structure (1)=\"150\">\n");
		dWAR.append("<columnWidth_Table_Structure (2)=\"150\">\n");
		dWAR.append("<detailView=\"height[Data]=0.5;height[Structure (1)]=0.25;height[Structure (2)]=0.25\">\n");
		dWAR.append("<filter0=\"#browser#	#disabled#	Structure (1)\">\n");
		dWAR.append("<filter1=\"#structure#	Structure (2)\">\n");
		dWAR.append("<filter2=\"#category#	Similarity\">\n");
		int filterCounter = 3;
		for (String property: properties) {
			dWAR.append("<filter" + filterCounter + "=\"#double#	" + property + " (1)	#disabled#\">\n");
			dWAR.append("<filter" + (filterCounter+1) + "=\"#double#	" + property + " (2)	#disabled#\">\n");
			dWAR.append("<filter" + (filterCounter+2) + "=\"#double#	" + property + " (delta)	#disabled#\">\n");
			filterCounter += 3;
		}
		dWAR.append("<mainSplitting=\"0.8\">\n");
		dWAR.append("<mainView=\"Structure (2)\">\n");
		dWAR.append("<mainViewCount=\"4\">\n");
		dWAR.append("<mainViewDockInfo0=\"root\">\n");
		dWAR.append("<mainViewDockInfo1=\"Table	bottom	0.7\">\n");
		dWAR.append("<mainViewDockInfo2=\"Table	right	0.5\">\n");
		dWAR.append("<mainViewDockInfo3=\"Structure (1)	right	0.5\">\n");
		dWAR.append("<mainViewName0=\"Table\">\n");
		dWAR.append("<mainViewName1=\"2D View\">\n");
		dWAR.append("<mainViewName2=\"Structure (1)\">\n");
		dWAR.append("<mainViewName3=\"Structure (2)\">\n");
		dWAR.append("<mainViewType0=\"tableView\">\n");
		dWAR.append("<mainViewType1=\"2Dview\">\n");
		dWAR.append("<mainViewType2=\"structureView\">\n");
		dWAR.append("<mainViewType3=\"structureView\">\n");
		dWAR.append("<rightSplitting=\"0.6\">\n");
		dWAR.append("<rowHeight_Table=\"80\">\n");
		dWAR.append("<showNaNValues_2D View=\"true\">\n");
		dWAR.append("<structureGridColumn_Structure (1)=\"Structure (1)\">\n");
		dWAR.append("<structureGridColumn_Structure (2)=\"Structure (2)\">\n");
		dWAR.append("<structureGridColumns_Structure (1)=\"3\">\n");
		dWAR.append("<structureGridColumns_Structure (2)=\"3\">\n");
		dWAR.append("</datawarrior properties>\n");
		return dWAR.toString();
	}
	
	/**
	 * Returns the number of transformations for a defined variable part of a molecule
	 * @param value1 IDCode of the variable part of a molecule
	 * @param minAtoms minimal number of atoms in the replacement (relative to value1Atoms)
	 * @param maxAtoms maximal number of atoms in the replacement (relative to value1Atoms)
	 * @return an integer corresponding to the number of transformations
	 */
	public int getTransformationsSize(String value1, int minAtoms, int maxAtoms) {
		int mmpSize = 0;
		MMPUniqueFragment value1Fragment = mmpUniqueFragments.fragmentIDToFragment(value1);
		if (value1Fragment != null) {
			Integer value1Atoms = value1Fragment.getFragmentAtoms();
			Integer value1Index = value1Fragment.getFragmentIndex();
			if (value1Index != null && mmpIndex.containsKey(value1Index)) {
				HashMap<Integer, List<int[]>> mmps = mmpIndex.get(value1Index);
				for (int size=value1Atoms+minAtoms; size<=value1Atoms+maxAtoms; size++) {
					if (mmps.containsKey(size)) {
						mmpSize += mmps.get(size).size();
					}
				}
			}
		}
		return mmpSize;
	}
	
	/**
	 * Returns the list of transformations for a defined variable part of a molecule.
	 * @param moleculeIDCode IDCode of the input molecule
	 * @param keys IDCode of the constant part of the molecule
	 * @param value1 IDCode of the variable part of a molecule
	 * @param minAtoms minimal number of atoms in the replacement (relative to value1Atoms)
	 * @param maxAtoms maximal number of atoms in the replacement (relative to value1Atoms)
	 * @return a list of transformations currently sorted by decreasing number of examples for each transformation 
	 */
	private List<MatchedMolecularPair> getTransformations(String moleculeIDCode, String[] keys, String value1, int minAtoms, int maxAtoms, String sortBy) {
		List<MatchedMolecularPair> retVal = new ArrayList<MatchedMolecularPair>();
		Integer[] keyIndex = fragmentToFragmentIndex(keys);
		String keyIndexString = null;
		if (keyIndex.length == 1 && keyIndex[0] != null) {
			keyIndexString = Integer.toString(keyIndex[0]);
		}
		else if (keyIndex.length == 2 && keyIndex[0] != null && keyIndex[1] != null) {
			keyIndexString = Integer.toString(keyIndex[0]) + "\t" + Integer.toString(keyIndex[1]);
		}
		MMPUniqueFragment value1Fragment = mmpUniqueFragments.fragmentIDToFragment(value1);
		if (value1Fragment != null) {
			Integer value1Atoms = value1Fragment.getFragmentAtoms();
			Integer value1Index = value1Fragment.getFragmentIndex();
			// TODO: obtain 2 keysFragment for doublecuts and two keysFP, and use them in examplesToMolecules
			MMPUniqueFragment keysFragment = mmpUniqueFragments.fragmentIDToFragment(keys);
			String[] keysFP = null;
			if (keysFragment != null) {
				keysFP = keysFragment.getFragmentFP();	
			}
			String[] value1FP = value1Fragment.getFragmentFP();
			if (value1Index != null && mmpIndex.containsKey(value1Index)) {
				HashMap<Integer, List<int[]>> mmps = mmpIndex.get(value1Index);
				List<int[]> fragmentsIndex = null;
				if (keyIndexString != null && mmpFragmentsIndex.containsKey(keyIndexString)){
					fragmentsIndex = mmpFragmentsIndex.get(keyIndexString);
				}
				for (int size=value1Atoms+minAtoms; size<=value1Atoms+maxAtoms; size++) {
					if (mmps.containsKey(size)) {
						List<int[]> values2 = mmps.get(size);
						for (int[] value2_and_examples: values2) {
							int value2Index = value2_and_examples[0];
							String value2 = uniqueFragmentsIndex.get(value2Index);
//							int[] examples = Arrays.copyOfRange(value2_and_examples, 1, value2_and_examples.length);
							int targetExists = -1;
							if (fragmentsIndex != null) { // fragmentsIndex contains all [valueIndex, molIndex] for the seeded key 
								for (int[] fragmentIndex: fragmentsIndex) {
									if (fragmentIndex[0] == value2Index) {
										targetExists = fragmentIndex[1];
										break;
									}
								}
							}
							List<MatchedMolecularPairExamples> mmpExamples = examplesToMolecules(value2_and_examples, keysFP, value2Index, targetExists);
							if (targetExists == -1) { 
								ArrayList<MoleculeIndex> currents = wholeMoleculesIndex.get(moleculeIDCode);
								if (currents == null) {
									StereoMolecule virtualMol = new StereoMolecule();
									IDCodeParser idCodeParser = new IDCodeParser();
									idCodeParser.parse(virtualMol, moleculeIDCode);
									MoleculeIndex current = new MoleculeIndex(-1, "", moleculeIDCode, "", generateData(virtualMol));
//									current.setIDCode(moleculeIDCode);
									currents = new ArrayList<MoleculeIndex>();
									currents.add(current);
								}
								else {
									currents.get(0).setIDCode(moleculeIDCode);
								}
								StereoMolecule virtualMol = molFromKeyValue(keys, value2);
								MoleculeIndex virtual = new MoleculeIndex(-1, "", virtualMol.getIDCode(), "", generateData(virtualMol));
//								virtual.setIDCode(idCodeFromKeyValue(keys, value2));
								ArrayList<MoleculeIndex> virtuals = new ArrayList<MoleculeIndex>();
								virtuals.add(virtual);
								MatchedMolecularPairExamples matchedMolecularPairExamples = new MatchedMolecularPairExamples(currents, virtuals, 6);
								mmpExamples.add(0, matchedMolecularPairExamples);
							}
							Collections.sort(mmpExamples, EXAMPLES_SIMILARITY_SORT);
							MMPUniqueFragment value2Fragment = mmpUniqueFragments.fragmentIDToFragment(value2);
							String[] value2FP = value2Fragment.getFragmentFP();
							retVal.add(new MatchedMolecularPair(value1, value1Index, value1Atoms, value1FP, value2, value2Index, size, value2FP, mmpExamples, dataFields.size(), targetExists));
						}
					}
				}
			}
		}
		if (sortBy == null || sortBy.equals(SORT_BY_NUMBER_OF_EXAMPLES)) {
			Collections.sort(retVal, NUMBER_OF_EXAMPLES_SORT);
		}
		else if (sortBy.equals(SORT_BY_SIMILARITY)) {
			Collections.sort(retVal, SIMILARITY_SORT);
		}
		return retVal;
	}
	
	/**
	 * Generates the DWAR file for the Transformations corresponding to a specific seed 'value', number of atoms for the replacement, environment size, and list of data fields
	 * @param moleculeIDCode idCode of the seed molecule
	 * @param keys Array of 'keys' idCodes (constant part of the molecule)
	 * @param value1 seed 'value' idCode (variable part of the molecule)
	 * @param minAtoms minimal number of atoms in the replacement (relative to value1Atoms)
	 * @param maxAtoms maximal number of atoms in the replacement (relative to value1Atoms)
	 * @param environmentSize size of the environment (0-5)
	 * @param properties List of numerical data fields
	 * @return String of the whole DWAR file
	 */
	public String getTransformationsDWAR(String moleculeIDCode, String[] keys, String value1, int minAtoms, int maxAtoms, Integer environmentSize, List<String> properties) {
		if (environmentSize == null) {
			environmentSize = 0;
		}
		List<MatchedMolecularPair> transformations = getTransformations(moleculeIDCode, keys, value1, minAtoms, maxAtoms, null);
		NumberFormat formatter = new DecimalFormat("#.##");
		StringBuilder dWAR = new StringBuilder();
		dWAR.append("<datawarrior-fileinfo>\n");
		dWAR.append("<version=\"3.1\">\n");
		dWAR.append("<rowcount=\"" + transformations.size() + "\">\n");
		dWAR.append("</datawarrior-fileinfo>\n");
		dWAR.append("<column properties>\n");
		dWAR.append("<columnName=\"Transformation\">\n");
		dWAR.append("<columnProperty=\"specialType	rxncode\">\n");
		dWAR.append("<columnName=\"Product\">\n");
		dWAR.append("<columnProperty=\"specialType	idcode\">\n");
		dWAR.append("<columnName=\"Structure\">\n");
		dWAR.append("<columnProperty=\"specialType	idcode\">\n");
		dWAR.append("<columnName=\"idcoordinates2D\">\n");
		dWAR.append("<columnProperty=\"specialType	idcoordinates2D\">\n");
		dWAR.append("<columnProperty=\"parent	Structure\">\n");
		dWAR.append("</column properties>\n");
		dWAR.append("Transformation\tProduct\tDeltaAtoms\tStructure\tidcoordinates2D\tActelion No\tExists\tExamples");
		for (String property: properties) {
			dWAR.append("\t" + property + " Avg\t" + property + " SD\t" + property + " n");
		}
		dWAR.append("\n");
		for (MatchedMolecularPair transformation: transformations) {
			transformation.calcTransformationString();
			String buildingBlock = generateBuildingBlock(transformation.value2);
			dWAR.append(transformation.transformationString + "\t" + buildingBlock + "\t" + (transformation.value2Atoms - transformation.value1Atoms));
			if (!transformation.targetExists) {
				StereoMolecule virtualMol = molFromKeyValue(keys, transformation.value2);
				dWAR.append("\t" + virtualMol.getIDCode() + "\t\t");
			}
			else {
				dWAR.append("\t" + transformation.mmpExamples.get(0).example2.get(0).moleculeIDCode + "\t" + transformation.mmpExamples.get(0).example2.get(0).moleculeIDCoord + "\t" + transformation.mmpExamples.get(0).example1.get(0).moleculeName);
			}
			dWAR.append("\t" + transformation.targetExists + "\t" + transformation.numberOfExamples);
			for (String property: properties) {
				int fieldIndex = -1;
				for (int i=0; i<dataFields.size(); i++) {
					if (dataFields.get(i).fieldName.equals(property) || dataFields.get(i).longFieldName.equals(property)) {
						fieldIndex = i;
						break;
					}
				}
				if (fieldIndex != -1) {
					if (transformation.average[fieldIndex][environmentSize] != null) {
						dWAR.append("\t" + formatter.format(transformation.average[fieldIndex][environmentSize]));
					}
					else {
						dWAR.append("\t");
					}
					if (transformation.sd[fieldIndex][environmentSize] != null) {
						if (Double.isNaN(transformation.sd[fieldIndex][environmentSize])) {
							dWAR.append("\t0.0");
						}
						else {
							dWAR.append("\t" + formatter.format(transformation.sd[fieldIndex][environmentSize]));
						}
					}
					else {
						dWAR.append("\t");
					}
					dWAR.append("\t" + transformation.n[fieldIndex][environmentSize]);
				}
			}
			dWAR.append("\n");
		}
		dWAR.append("<datawarrior properties>\n");
		if (properties.size() > 0) {
			dWAR.append("<axisColumn_2D View_0=\"" + properties.get(0) + " Avg\">\n");
			if (properties.size() > 1) {
				dWAR.append("<axisColumn_2D View_1=\"" + properties.get(1) + " Avg\">\n");
			}
			else {
				dWAR.append("<axisColumn_2D View_1=\"" + properties.get(0) + " Avg\">\n");
			}
		}
		dWAR.append("<chartType_2D View=\"scatter\">\n");
		dWAR.append("<columnDescriptionCount=\"" + (properties.size()*3+3) + "\">\n");
		dWAR.append("<columnDescription_0=\"DeltaAtoms	Difference of number of heavy atoms between seed and target\">\n");
		dWAR.append("<columnDescription_1=\"Exists	Target structure exists in the dataset or not\">\n");
		dWAR.append("<columnDescription_2=\"Examples	Number of MMPs\">\n");
		int descriptionCounter = 3;
		for (String property: properties) {
			int fieldIndex = -1;
			for (int i=0; i<dataFields.size(); i++) {
				if (dataFields.get(i).fieldName.equals(property) || dataFields.get(i).longFieldName.equals(property)) {
					fieldIndex = i;
					break;
				}
			}
			if (fieldIndex != -1) {
				dWAR.append("<columnDescription_" + descriptionCounter + "=\"" + property + " Avg	Average of MMPs of " + dataFields.get(fieldIndex).longFieldName + "\">\n");
				dWAR.append("<columnDescription_" + (descriptionCounter+1) + "=\"" + property + " SD	SD of MMPs of " + dataFields.get(fieldIndex).longFieldName + "\">\n");
				dWAR.append("<columnDescription_" + (descriptionCounter+2) + "=\"" + property + " n	Number of MMPs of " + dataFields.get(fieldIndex).longFieldName + "\">\n");
				descriptionCounter += 3;
			}
		}
		dWAR.append("<columnWidth_Table_Structure=\"150\">\n");
		dWAR.append("<columnWidth_Table_Transformation=\"260\">\n");
		dWAR.append("<detailView=\"height[Data]=0.7;height[Structure]=0.3\">\n");
		dWAR.append("<filter0=\"#browser#	#disabled#	Transformation\">\n");
		dWAR.append("<filter1=\"#structure#	Structure\">\n");
		int filterCounter = 2;
		if (minAtoms != maxAtoms) {
			dWAR.append("<filter2=\"#double#	DeltaAtoms\">\n");
			filterCounter++;
		}
		dWAR.append("<filter" + filterCounter + "=\"#double#	Examples\">\n");
		filterCounter++;
		dWAR.append("<filter" + filterCounter + "=\"#category#	Exists\">\n");
		filterCounter++;
		for (String property: properties) {
			dWAR.append("<filter" + filterCounter + "=\"#double#	" + property + " Avg	#disabled#\">\n");
//			dWAR.append("<filter" + (filterCounter+1) + "=\"#double#	" + property + " n\">\n");
			filterCounter += 1;
		}
		dWAR.append("<mainSplitting=\"0.8\">\n");
		dWAR.append("<mainView=\"Structures\">\n");
		dWAR.append("<mainViewCount=\"3\">\n");
		dWAR.append("<mainViewDockInfo0=\"root\">\n");
		String tableName = "Transformations";
		if (environmentSize != 0) {
			tableName = "Transformations (environment size " + Integer.toString(environmentSize) + ")";
		}
		dWAR.append("<mainViewDockInfo1=\"" + tableName + "	bottom	0.7\">\n");
		dWAR.append("<mainViewDockInfo2=\"" + tableName + "	right	0.7\">\n");
		dWAR.append("<mainViewName0=\"" + tableName + "\">\n");
		dWAR.append("<mainViewName1=\"2D View\">\n");
		dWAR.append("<mainViewName2=\"Structures\">\n");
		dWAR.append("<mainViewType0=\"tableView\">\n");
		dWAR.append("<mainViewType1=\"2Dview\">\n");
		dWAR.append("<mainViewType2=\"structureView\">\n");
		dWAR.append("<rightSplitting=\"0.6\">\n");
		dWAR.append("<rowHeight_" + tableName + "=\"80\">\n");
		dWAR.append("<showNaNValues_2D View=\"true\">\n");
		dWAR.append("<structureGridColumn_Structures=\"Structure\">\n");
		dWAR.append("<structureGridColumns_Structures=\"3\">\n");
		dWAR.append("</datawarrior properties>\n");
		return dWAR.toString();
	}
	
	/**
	 * Returns a list of transformations
	 * @param keys Array of 'keys' idCodes (constant part of the molecule)
	 * @param value1 seed 'value' idCode (variable part of the molecule)
	 * @param minAtoms minimal number of atoms in the replacement (relative to value1Atoms)
	 * @param maxAtoms maximal number of atoms in the replacement (relative to value1Atoms)
	 * @return List of transformations
	 */
	public List<String[]> transformationsListToTable(String[] keys, String value1, int minAtoms, int maxAtoms) {
		List<MatchedMolecularPair> transformations = getTransformations(null, keys, value1, minAtoms, maxAtoms, null);
		List<String[]> retVal = new ArrayList<String[]>();
		for (MatchedMolecularPair transformation: transformations) {
			String isCurrent = transformation.targetExists ? "1" : "0";
			retVal.add(new String[]{transformation.value1, transformation.value2, Integer.toString(transformation.numberOfExamples), isCurrent});
		}
		return retVal;
	}
	
	/**
	 * Generates the main JSON file for Transformations and Matched Molecular Pairs
	 * @param moleculeIDCode idCode of the seed molecule
	 * @param keys Array of 'keys' idCodes (constant part of the molecule)
	 * @param value1 seed 'value' idCode (variable part of the molecule)
	 * @param minAtoms minimal number of atoms in the replacement (relative to value1Atoms)
	 * @param maxAtoms maximal number of atoms in the replacement (relative to value1Atoms)
	 * @param sortBy SORT_BY_SIMILARITY, SORT_BY_NUMBER_OF_EXAMPLES
	 * @return JSON string
	 */
	public String getTransformationsJSON(String moleculeIDCode, String[] keys, String value1, int minAtoms, int maxAtoms, String sortBy) {
		List<MatchedMolecularPair> transformations = getTransformations(moleculeIDCode, keys, value1, minAtoms, maxAtoms, sortBy);
		NumberFormat formatter = new DecimalFormat("#.##");
		StringBuilder jSON = new StringBuilder();
		jSON.append("{\"transformations\": [");
		int counter = 0;
		for (MatchedMolecularPair transformation: transformations) {
			if (counter > 0) {
				jSON.append(", ");
			}
			jSON.append( "\n\t{\"value1\": \"" + transformation.value1.replace("\\", "\\\\") + "\"");
			jSON.append(",\n\t \"value2\": \"" + transformation.value2.replace("\\", "\\\\") + "\"");
			jSON.append(",\n\t \"n\": " + transformation.numberOfExamples);
			jSON.append(",\n\t \"delta_atoms\": " + (transformation.value2Atoms - transformation.value1Atoms));
			jSON.append(",\n\t \"similarity\": " + transformation.similarity);
			String isCurrent = transformation.targetExists ? "true" : "false";
			jSON.append(",\n\t \"current\": " + isCurrent + "");
			try {
				jSON.append(",\n\t \"image\": \"" + getB64Image(getImage(transformation.value1, transformation.value2, 580, 266)) + "\"");
			}
			catch (Exception e) {
				// cannot generate image
			}
			// TODO: handling of multiple compounds with same structure (get(0), get(1)...)
			jSON.append(",\n\t \"compounds\": [");
			for (int j=0; j<transformation.numberOfExamples; j++) {
				if (j > 0) {
					jSON.append(", ");
				}
				jSON.append("[\"" + transformation.mmpExamples.get(j).example1.get(0).moleculeName + "\", \"" + transformation.mmpExamples.get(j).example2.get(0).moleculeName + "\"]");
			}
			jSON.append("]");
			jSON.append(",\n\t \"similarities\": [");
			for (int j=0; j<transformation.numberOfExamples; j++) {
				if (j > 0) {
					jSON.append(", ");
				}
				jSON.append(transformation.mmpExamples.get(j).similarity);
			}
			jSON.append("]");
			// TODO: handling of multiple compounds with same structure (get(0), get(1)...)
			jSON.append(",\n\t \"structures\": [");
			for (int j=0; j<transformation.numberOfExamples; j++) {
				if (j > 0) {
					jSON.append(", ");
				}
				jSON.append("[\"" + transformation.mmpExamples.get(j).example1.get(0).moleculeIDCode.replace("\\", "\\\\") + "\", \"" + transformation.mmpExamples.get(j).example2.get(0).moleculeIDCode.replace("\\", "\\\\") + "\"]");
			}
			jSON.append("]");
			jSON.append(",\n\t \"coordinates\": [");
			for (int j=0; j<transformation.numberOfExamples; j++) {
				if (j > 0) {
					jSON.append(", ");
				}
				jSON.append("[\"" + transformation.mmpExamples.get(j).example1.get(0).moleculeIDCoord.replace("\\", "\\\\") + "\", \"" + transformation.mmpExamples.get(j).example2.get(0).moleculeIDCoord.replace("\\", "\\\\") + "\"]");
			}
			jSON.append("]");
			jSON.append(",\n\t \"datas\": [");
			for (int i=0; i<transformation.n.length; i++) {
				if (i > 0) {
					jSON.append(", ");
				}
				jSON.append("{");
				for (int j=0; j<6; j++) {
					if (j > 0) {
						jSON.append(", ");
					}
					jSON.append( "\n\t\t\"similarity" + j + "\":");
					jSON.append( "\n\t\t{\"n\": " + transformation.n[i][j]);
					jSON.append(",\n\t\t \"increase\": " + transformation.increase[i][j]);
					jSON.append(",\n\t\t \"decrease\": " + transformation.decrease[i][j]);
					jSON.append(",\n\t\t \"neutral\": " + transformation.neutral[i][j]);
					if (transformation.average[i][j] != null) {
						jSON.append(",\n\t\t \"average\": " + formatter.format(transformation.average[i][j]));
					}
					else {
						jSON.append(",\n\t\t \"average\": " + transformation.average[i][j]);
					}
					if (transformation.sd[i][j] != null) {
						if (Double.isNaN(transformation.sd[i][j])) {
							jSON.append(",\n\t\t \"sd\": 0.0");
						}
						else {
							jSON.append(",\n\t\t \"sd\": " + formatter.format(transformation.sd[i][j]));
						}
					}
					else {
						jSON.append(",\n\t\t \"sd\": " + transformation.sd[i][j]);
					}
					jSON.append("}");
				}
				jSON.append(",\n\t\t \"data\": [");
				for (int j=0; j<transformation.numberOfExamples; j++) {
					ArrayList<MoleculeIndex> example1 = transformation.mmpExamples.get(j).example1;
					ArrayList<MoleculeIndex> example2 = transformation.mmpExamples.get(j).example2;
					if (j > 0) {
						jSON.append(",");
					}
					// TODO: handling of multiple compounds with same structure (get(0), get(1)...)
					String data1 = "";
					if (example1.get(0).moleculeData != null) {
						data1 = example1.get(0).moleculeData[i];
					}
					String data2 = "";
					if (example2.get(0).moleculeData != null) {
						data2 = example2.get(0).moleculeData[i]; // virtual compound
					}
					if (data1.equals("")) {
						data1 = "\"n.a.\"";
					}
					else if (data1.startsWith(">") || data1.startsWith("<")) {
						data1 = "\"" + data1 + "\"";
					}
					else {
						data1 = formatter.format(Double.parseDouble(data1));
					}
					if (data2.equals("")) {
						data2 = "\"n.a.\"";
					}
					else if (data2.startsWith(">") || data2.startsWith("<")) {
						data2 = "\"" + data2 + "\"";
					}
					else {
						data2 = formatter.format(Double.parseDouble(data2));
					}
					jSON.append("[" + data1  + ", " + data2 + "]");
				}
				jSON.append("]}");
			}
			jSON.append("]}");
			counter += 1;
		}
		jSON.append("]}");
	return jSON.toString();
	}
	
	/**
	 * For double cuts, the R-Group of the second fragment has to be correctly tagged 
	 * @param mol input fragment
	 * @return modified StereoMolecule
	 */
	private StereoMolecule changeR1ToR2(StereoMolecule mol) {
		for (int i=0; i<mol.getAtoms(); i++) {
			if (mol.getAtomicNo(i) == MMPFragmenter.FRAGMENT_ATOMIC_NO) {
				mol.setAtomicNo(i, MMPFragmenter.FRAGMENT_ATOMIC_NO+1);
				break;
			}
		}
		return mol;
	}
	
	/**
	 * Generates the idCode of the building block (with explicit hydrogens and R-group replaced by 'any' atom)
	 * @param value 'value' idCode (variable part of the molecule
	 * @return idCode of the building block
	 */
	private String generateBuildingBlock(String value) {
		StereoMolecule mol = new StereoMolecule();
		IDCodeParser idCodeParser = new IDCodeParser();
		idCodeParser.parse(mol, value);
		for (int i=0; i<mol.getAtoms(); i++) {
			if (mol.getAtomicNo(i) >= MMPFragmenter.FRAGMENT_ATOMIC_NO) {
				mol.setAtomQueryFeature(i, StereoMolecule.cAtomQFAny, true);
			}
			else {
				mol.setAtomQueryFeature(i, StereoMolecule.cAtomQFNoMoreNeighbours, true);
			}
		}
		return mol.getIDCode();
	}
	
	/**
	 * Generates a molecule from one 'value' and one or two 'keys'
	 * @param keys Array of 'keys' idCodes (constant part of the molecule)
	 * @param value 'value' idCode (variable part of the molecule)
	 * @return newly generated molecule
	 */
	private StereoMolecule molFromKeyValue(String[] keys, String value) {
		int [] firstAtoms = new int[2*keys.length];
		int [] rBonds = new int[2*keys.length];
		int [] rGroupAtoms = new int[2*keys.length];
		int [] secondAtoms = new int[2*keys.length];
		int [] rGroupsIndex = new int[2*keys.length];
		int [] rGroupCounters = new int[keys.length];
		int atomLabel;
		StereoMolecule mol1 = new StereoMolecule();
		StereoMolecule mol2 = new StereoMolecule();
		IDCodeParser idCodeParser = new IDCodeParser();
		idCodeParser.parse(mol1, value);
		idCodeParser.parse(mol2, keys[0]);
		mol1.addMolecule(mol2);
		if (keys.length == 2) {
			StereoMolecule mol3 = new StereoMolecule();
			idCodeParser.parse(mol3, keys[1]);
			mol3 = changeR1ToR2(mol3);
			mol1.addMolecule(mol3);
			rGroupCounters[1] = 2;
		}
		mol1.ensureHelperArrays(StereoMolecule.cHelperBitNeighbours);
		for (int bond=0; bond<mol1.getBonds(); bond++) {
			int atom1 = mol1.getBondAtom(0, bond);
			int atom2 = mol1.getBondAtom(1, bond);
			if (mol1.getAtomicNo(atom1) == 0 || mol1.getAtomicNo(atom1) >= 142 || mol1.getAtomicNo(atom2) == 0 || mol1.getAtomicNo(atom2) >= 142) {
				if (mol1.getAtomicNo(atom1) == 0 || mol1.getAtomicNo(atom1) >= 142) {
					atomLabel = mol1.getAtomicNo(atom1) - MMPFragmenter.FRAGMENT_ATOMIC_NO;
					rBonds[rGroupCounters[atomLabel]] = bond;
					firstAtoms[rGroupCounters[atomLabel]] = atom1;
					rGroupAtoms[rGroupCounters[atomLabel]] = 0;
					secondAtoms[rGroupCounters[atomLabel]] = atom2;
					rGroupsIndex[rGroupCounters[atomLabel]] = atomLabel;
				}
				else {
					atomLabel = mol1.getAtomicNo(atom2) - MMPFragmenter.FRAGMENT_ATOMIC_NO;
					rBonds[rGroupCounters[atomLabel]] = bond;
					firstAtoms[rGroupCounters[atomLabel]] = atom2;
					rGroupAtoms[rGroupCounters[atomLabel]] = 1;
					secondAtoms[rGroupCounters[atomLabel]] = atom1;
					rGroupsIndex[rGroupCounters[atomLabel]] = atomLabel;
				}
				rGroupCounters[atomLabel]++;
			}
		}
		mol1.setBondAtom(rGroupAtoms[0], rBonds[0], secondAtoms[1]);
		mol1.markBondForDeletion(rBonds[1]);
		mol1.markAtomForDeletion(firstAtoms[1]);
		mol1.markAtomForDeletion(firstAtoms[0]);
		if (keys.length == 2) {
			mol1.setBondAtom(rGroupAtoms[2], rBonds[2], secondAtoms[3]);
			mol1.markBondForDeletion(rBonds[3]);
			mol1.markAtomForDeletion(firstAtoms[3]);
			mol1.markAtomForDeletion(firstAtoms[2]);			
		}
		mol1.deleteMarkedAtomsAndBonds();
//		return mol1.getIDCode();
		return mol1;
	}
	
	/**
	 * Returns a list of MatchedMolecularPairExamples from a list of molecule names
	 * @param examples List of molecule indexes
	 * @param keys1FP fingerprints of the 'keys'
	 * @param value2Index index of the target 'value'
	 * @param targetExists true/false if the target molecule exists (after replacing the seed 'value' by the target 'value' in the seeded molecule)
	 * @return a list of MatchedMolecularPairExamples
	 */
	private List<MatchedMolecularPairExamples> examplesToMolecules(int[] examples, String[] keys1FP, int value2Index, int targetExists) {
		List<MatchedMolecularPairExamples> retVal = new ArrayList<MatchedMolecularPairExamples>();
		List<int[]> fragmentsIndex = null;
		if (mmpFragmentsIndex.containsKey(Integer.toString(value2Index))) {
			fragmentsIndex = mmpFragmentsIndex.get(Integer.toString(value2Index));
		}
		// pass keyIndex and use that instead of value2Index in case of too small value
		for (int i=1; i<examples.length; i+=2) { // examples[0] is the value2Index
			ArrayList<MoleculeIndex> example1 = molIndexToMolecule(examples[i]);
			ArrayList<MoleculeIndex> example2 = molIndexToMolecule(examples[i+1]);
			int similarity = -1;
			if (fragmentsIndex != null) { // fragmentsIndex contains all [valueIndex, molIndex] for the seeded key (in this case the value, to find the key) 
				for (int[] fragmentIndex: fragmentsIndex) {
					if (fragmentIndex[1] == examples[i+1]) {
						if (fragmentIndex[1] == targetExists) {
							similarity = 6;
						}
						else {
							similarity = 0;
							int keyIndex = fragmentIndex[0];
							String keyID = uniqueFragmentsIndex.get(keyIndex);
							MMPUniqueFragment keysFragment = mmpUniqueFragments.fragmentIDToFragment(keyID);
							String[] keys2FP = keysFragment.getFragmentFP();
							if (keys1FP != null) {
								for (int j=keys1FP.length-1; j>=0; j--) {
									if (keys1FP[j].equals(keys2FP[j])) {
										similarity = j + 1;
										break;
									}
								}
							}
						}
						break;
					}
				}
			}
			MatchedMolecularPairExamples matchedMolecularPairExamples = new MatchedMolecularPairExamples(example1, example2, similarity);
			retVal.add(matchedMolecularPairExamples);
		}
		return retVal;
	}
	
	/**
	 * Returns a list of MoleculeIndex objects from one molIndex
	 * @param molIndex
	 * @return
	 */
	private ArrayList<MoleculeIndex> molIndexToMolecule(int molIndex) {
		String moleculeIDCode = molecules.get(molIndex);
		ArrayList<MoleculeIndex> moleculesIndex = wholeMoleculesIndex.get(moleculeIDCode);
		// here I could put the IDCode just for the first entry, since I only read the first one later on...
		for (MoleculeIndex moleculeIndex: moleculesIndex) {
			moleculeIndex.setIDCode(moleculeIDCode);
		}
		return moleculesIndex;
	}
	
	/**
	 * Returns the idCode from one compound name, null if the compound is not found
	 * @param molName name of a compound
	 * @return null or idCode
	 */
	public String getIDCodeFromMolName(String molName) {
		for (Entry<String, ArrayList<MoleculeIndex>> cursor : wholeMoleculesIndex.entrySet()) {
			for (MoleculeIndex moleculeIndex: cursor.getValue()) {
				if (moleculeIndex.moleculeName.equals(molName)) {
					if (cursor.getValue().get(0).moleculeIDCoord != null) {
						return cursor.getKey() + "\t" + cursor.getValue().get(0).moleculeIDCoord;
					}
					return cursor.getKey();
				}
			}
		}
		return null;
//		Iterator<String> it = wholeMoleculesIndex.keySet().iterator();
//		while (it.hasNext()) {
//			String idCode = it.next();
//			List<MoleculeIndex> moleculesIndex = wholeMoleculesIndex.get(idCode);
//			for (MoleculeIndex moleculeIndex: moleculesIndex) {
//				if (moleculeIndex.moleculeName.equals(molName)) {
//					retVal = idCode;
//					return retVal;
//				}
//			}
//		}
//		return retVal;
	}

	private String keysToKeysString(Integer[] keys) {
		String keysString = null;
		if (keys.length == 1 && keys[0] != null) {
			keysString = Integer.toString(keys[0]);
		}
		else if (keys.length == 2 && keys[0] != null && keys[1] != null) {
			keysString = Integer.toString(keys[0]) + "\t" + Integer.toString(keys[1]);
		}
		return keysString;
	}
	
	/**
	 * Returns a list of fields
	 * @param what Requested data (fieldName, longFieldName, categoryName, percentile5, percentile95)
	 * @return List of fields
	 */
	public List<String> getDataFields(String what) {
		List<String> retVal = new ArrayList<String>();
		for (DataField dataField: dataFields) {
			if (what.equals("fieldName")) {
				retVal.add(dataField.fieldName);
			}
			else if (what.equals("longFieldName")) {
				if (dataField.longFieldName != null) {
					retVal.add(dataField.longFieldName);
				}
				else {
					retVal.add(dataField.fieldName);
				}
			}
			else if (what.equals("categoryName")) {
				if (dataField.categoryName != null) {
					retVal.add(dataField.categoryName);
				}
				else {
					retVal.add("Other");
				}
			}
			else if (what.equals("percentile5")) {
				retVal.add(dataField.percentile5);
			}
			else if (what.equals("percentile95")) {
				retVal.add(dataField.percentile95);
			}
		}
		return retVal;
	}
	
	/**
	 * Returns the requested general information of the data set
	 * @param what Requested information (datasetName, date, numberOfMolecules, randomMoleculeName)
	 * @return the requested information
	 */
	public String getWhat(String what) {
		if (what.equals("datasetName")) {
			return datasetName; 
		}
		else if (what.equals("date")) {
			return date;
		}
		else if (what.equals("numberOfMolecules")) {
			return Integer.toString(wholeMoleculesIndex.size());
		}
		else if (what.equals("randomMoleculeName")) {
			Random randomGenerator = new Random();
			if (molecules.size() > 0) {
				int index = randomGenerator.nextInt(molecules.size());
				ArrayList<MoleculeIndex> moleculesIndex = wholeMoleculesIndex.get(molecules.get(index));
				if (moleculesIndex.size() > 0) {
					return moleculesIndex.get(0).moleculeName;
				}
				return null;
			}
			return null;
		}
		return null;
	}
	
	/**
	 * Generates an image of the transformation from the seed 'value' to the target 'value'
	 * @param value1 seed 'value' idCode
	 * @param value2 target 'value' idCode
	 * @param width of the image
	 * @param height of the image
	 * @return BufferedImage object
	 */
	private static BufferedImage getImage(String value1, String value2, int width, int height) {
		BufferedImage bufferedImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		GenericRectangle viewRect = new GenericRectangle(0, 0, width, height);
		ExtendedDepictor extendedDepictor;
		StereoMolecule mol1 = new StereoMolecule();
		IDCodeParser idCodeParser = new IDCodeParser();
		idCodeParser.parse(mol1, value1);
		if (value2 == null) {
			extendedDepictor = new ExtendedDepictor(new StereoMolecule[]{mol1}, null);
		}
		else {
			StereoMolecule mol2 = new StereoMolecule();
			idCodeParser.parse(mol2, value2);
			Reaction rxn = new Reaction(new StereoMolecule[]{mol1, mol2}, 1);		
			extendedDepictor = new ExtendedDepictor(rxn, rxn.getDrawingObjects(), true);
		}
        Graphics2D graphics2D = (Graphics2D) bufferedImage.getGraphics();
		GenericDrawContext context = new SwingDrawContext(graphics2D);
        extendedDepictor.validateView(context, viewRect, AbstractDepictor.cModeInflateToMaxAVBL + 45);
        RenderingHints renderingHints = new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		graphics2D.addRenderingHints(renderingHints);
		extendedDepictor.paint(context);
        return bufferedImage;
	}
	
	/**
	 * Generates the 64 bits-encoded representation of an image
	 * @param bufferedImage
	 * @return 64 bits-encoded string
	 * @throws IOException
	 */
	private static String getB64Image(BufferedImage bufferedImage) throws IOException {
		ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
		OutputStream b64 = new Base64.OutputStream(outputStream);
		ImageIO.write(bufferedImage, "png", b64);
		return outputStream.toString("UTF-8");
	}
	
	private String[] generateData(StereoMolecule mol) {
		String[] datas = new String[dataFields.size()];
		Arrays.fill(datas, "");
		for (int i=0; i<dataFields.size(); i++) {
			if (dataFields.get(i).categoryName.equals("Calculated")) {
				String data = mPropertyCalculator.getCalculatedValue(dataFields.get(i).fieldName, mol);
				if (data != null) {
					datas[i] = data;
				}
			}
		}
		return datas;
	}
}
