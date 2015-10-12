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

package com.actelion.research.chem.prediction;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.SortedStringList;
import com.actelion.research.chem.StereoMolecule;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class ToxicityPredictor {
    public static final String[] RISK_NAME = { "unknown", "none", "low", "high" };
    public static final int cUnknownRisk = 0;
	public static final int cNoRisk = 1;
	public static final int cLowRisk = 2;
	public static final int cHighRisk = 3;

	public static final int cRiskTypes = 4;
	public static final int cRiskTypeMutagenic = 0;
	public static final int cRiskTypeTumorigenic = 1;
	public static final int cRiskTypeIrritant = 2;
	public static final int cRiskTypeReproductiveEffective = 3;

	public static final String[] cRiskNameA = {	"mutagenic",
												"tumorigenic",
												"irritant",
												"reproductive effective" };
	public static final String[] cRiskNameN = {	"Mutagenicity",
												"Tumorigenicity",
												"Irritating effects",
												"Reproductive effects" };

	private static ArrayList<String>[]	sHighRiskFragments;
	private static ArrayList<String>[]	sLowRiskFragments;
	private static SortedStringList[]	sRiskMolecules;
    private static boolean				sInitialized;

	@SuppressWarnings("unchecked")
	public ToxicityPredictor() {
		synchronized(ToxicityPredictor.class) {
			if (!sInitialized) {
				try {
					sRiskMolecules = new SortedStringList[cRiskTypes];
					sHighRiskFragments = new ArrayList[cRiskTypes];
					sLowRiskFragments = new ArrayList[cRiskTypes];
					sHighRiskFragments[0] = readIDCodeFile("/resources/toxpredictor/m1.txt");
					sHighRiskFragments[1] = readIDCodeFile("/resources/toxpredictor/t1.txt");
					sHighRiskFragments[2] = readIDCodeFile("/resources/toxpredictor/i1.txt");
					sHighRiskFragments[3] = readIDCodeFile("/resources/toxpredictor/r1.txt");
					sLowRiskFragments[0] = readIDCodeFile("/resources/toxpredictor/m2.txt");
					sLowRiskFragments[1] = readIDCodeFile("/resources/toxpredictor/t2.txt");
					sLowRiskFragments[2] = readIDCodeFile("/resources/toxpredictor/i2.txt");
					sLowRiskFragments[3] = readIDCodeFile("/resources/toxpredictor/r2.txt");
					sRiskMolecules[0] = readAndSortIDCodeFile("/resources/toxpredictor/m3.txt");
					sRiskMolecules[1] = readAndSortIDCodeFile("/resources/toxpredictor/t3.txt");
					sRiskMolecules[2] = readAndSortIDCodeFile("/resources/toxpredictor/i3.txt");
					sRiskMolecules[3] = readAndSortIDCodeFile("/resources/toxpredictor/r3.txt");
					sInitialized = true;
					}
				catch (Exception e) {
					System.out.println("Error ToxicityPredictor::initialize() " + e);
					}
				}
			}
		}


	public int assessRisk(StereoMolecule testMolecule, int riskType, ThreadMaster threadMaster) {
		if (!sInitialized)
			return cUnknownRisk;

		if (sRiskMolecules[riskType].contains(new Canonizer(testMolecule).getIDCode()))
			return cHighRisk;

		SSSearcher sss = new SSSearcher(SSSearcher.cMatchAtomCharge);
		StereoMolecule fragment = new StereoMolecule();
		for (int i=0; i<sHighRiskFragments[riskType].size(); i++) {
			if (threadMaster != null && threadMaster.threadMustDie())
				return cUnknownRisk;

			Thread.yield();

			new IDCodeParser(false).parse(fragment, sHighRiskFragments[riskType].get(i));
			sss.setMol(fragment, testMolecule);
			if (sss.isFragmentInMolecule())
				return cHighRisk;
			}

		for (int i=0; i<sLowRiskFragments[riskType].size(); i++) {
			if (threadMaster != null && threadMaster.threadMustDie())
				return cUnknownRisk;

			Thread.yield();

			new IDCodeParser(false).parse(fragment, sLowRiskFragments[riskType].get(i));
			sss.setMol(fragment, testMolecule);
			if (sss.isFragmentInMolecule())
				return cLowRisk;
			}

		return cNoRisk;
		}


	public ParameterizedStringList getDetail(StereoMolecule testMolecule, int riskType) {
		ParameterizedStringList theDetail = new ParameterizedStringList();

		if (!sInitialized) {
			theDetail.add("Toxicity predictor not properly initialized.",
								 ParameterizedStringList.cStringTypeText);
			return theDetail;
			}

		String idcode = new Canonizer(testMolecule).getIDCode();
		if (sRiskMolecules[riskType].contains(idcode)) {
			theDetail.add("This molecule is known to be "+cRiskNameA[riskType]+":",
								 ParameterizedStringList.cStringTypeText);
			theDetail.add(idcode, ParameterizedStringList.cStringTypeIDCode);
			return theDetail;
			}

		SSSearcher sss = new SSSearcher(SSSearcher.cMatchAtomCharge);
		boolean found = false;
		StereoMolecule fragment = new StereoMolecule();
		for (int i=0; i<sHighRiskFragments[riskType].size(); i++) {
			new IDCodeParser(false).parse(fragment, sHighRiskFragments[riskType].get(i));
			sss.setMol(fragment, testMolecule);
			if (sss.isFragmentInMolecule()) {
				if (!found)
					theDetail.add("High-risk fragments indicating "+cRiskNameN[riskType]+":",
										 ParameterizedStringList.cStringTypeText);

				found = true;
				theDetail.add(sHighRiskFragments[riskType].get(i), ParameterizedStringList.cStringTypeIDCode);
				}
			}

		found = false;
		for (int i=0; i<sLowRiskFragments[riskType].size(); i++) {
			new IDCodeParser(false).parse(fragment, sLowRiskFragments[riskType].get(i));
			sss.setMol(fragment, testMolecule);
			if (sss.isFragmentInMolecule()) {
				if (!found)
					theDetail.add("Medium-risk fragments indicating "+cRiskNameN[riskType]+":",
										 ParameterizedStringList.cStringTypeText);

				found = true;
				theDetail.add(sLowRiskFragments[riskType].get(i), ParameterizedStringList.cStringTypeIDCode);
				}
			}

		if (theDetail.getSize() == 0)
			theDetail.add("No indication for "+cRiskNameN[riskType]+" found.",
								 ParameterizedStringList.cStringTypeText);

		return theDetail;
		}


	private ArrayList<String> readIDCodeFile(String filename) throws Exception {
		BufferedReader theReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(filename)));
		ArrayList<String> fragments = new ArrayList<String>();
		while (true) {
			try {
				String idcode = theReader.readLine();
				if (idcode != null)
					fragments.add(idcode);
				else
					break;
				}
			catch (IOException e) {	break; }
			}
		theReader.close();

		return fragments;
		}


	private SortedStringList readAndSortIDCodeFile(String filename) throws Exception {
		BufferedReader theReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(filename)));
		SortedStringList molecules = new SortedStringList();
		while (true) {
			try {
				String idcode = theReader.readLine();
				if (idcode != null)
					molecules.addString(idcode);
				else
					break;
				}
			catch (IOException e) {	break; }
			}
		theReader.close();

		return molecules;
		}
	}


