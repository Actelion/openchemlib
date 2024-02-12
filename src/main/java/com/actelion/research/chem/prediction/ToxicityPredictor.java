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
import java.nio.charset.StandardCharsets;
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


	/**
	 * Before calculating any kind of property, make sure that the molecule's structure is standardized.
	 * Typically, molecules created by an IDCodeParser are standardized. Molecules generated from a
	 * SmilesParser or MolfileParser, or just drawn within an editor, should be standardized using the
	 * MoleculeStandardizer.
	 * @param mol
	 * @param riskType one of the four risk types cRiskType...
	 * @param threadMaster may be null
	 * @return toxicity risk class estimated from atom type specific increments
	 */
	public int assessRisk(StereoMolecule mol, int riskType, ThreadMaster threadMaster) {
		if (!sInitialized)
			return cUnknownRisk;

		if (sRiskMolecules[riskType].contains(new Canonizer(mol).getIDCode()))
			return cHighRisk;

		SSSearcher sss = new SSSearcher(SSSearcher.cMatchAtomCharge);
		StereoMolecule fragment = new StereoMolecule();
		for (int i=0; i<sHighRiskFragments[riskType].size(); i++) {
			if (threadMaster != null && threadMaster.threadMustDie())
				return cUnknownRisk;

			Thread.yield();

			new IDCodeParser(false).parse(fragment, sHighRiskFragments[riskType].get(i));
			sss.setMol(fragment, mol);
			if (sss.isFragmentInMolecule())
				return cHighRisk;
			}

		for (int i=0; i<sLowRiskFragments[riskType].size(); i++) {
			if (threadMaster != null && threadMaster.threadMustDie())
				return cUnknownRisk;

			Thread.yield();

			new IDCodeParser(false).parse(fragment, sLowRiskFragments[riskType].get(i));
			sss.setMol(fragment, mol);
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
		BufferedReader theReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(filename), StandardCharsets.UTF_8));
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
		BufferedReader theReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(filename), StandardCharsets.UTF_8));
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


