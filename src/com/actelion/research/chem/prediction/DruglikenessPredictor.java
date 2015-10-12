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
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;

public class DruglikenessPredictor {
	public static final double cDruglikenessUnknown = -999;

    private static boolean			sInitialized = false;
	private static IncrementTable	sIncrementTable;

	private ParameterizedStringList	mDetail;

	public DruglikenessPredictor() {
		synchronized(DruglikenessPredictor.class) {
			if (!sInitialized) {
		        try {
					sIncrementTable = new IncrementTable("/resources/druglikenessNoIndex.txt");
					sInitialized = true;
					}
				catch (Exception e) {
		            System.out.println("Unable to initialize DruglikenessPredictor");
					}
				}
			}
		}


	public double assessDruglikeness(StereoMolecule testMolecule, ThreadMaster threadMaster) {
		ParameterizedStringList detail = new ParameterizedStringList();

		if (!sInitialized) {
			detail.add("Druglikeness predictor not properly initialized.",
								ParameterizedStringList.cStringTypeText);
			return cDruglikenessUnknown;
			}

		detail.add("Found sub-structure fragments and their contributions:",
							ParameterizedStringList.cStringTypeText);
		detail.add("(yellow atoms carry at least one more substituent)",
							ParameterizedStringList.cStringTypeText);
		double nastyIncrementSum = 0.0;
		double incrementSum = 0.0;
		int fragmentCount = 0;
		SSSearcher sss = new SSSearcher(SSSearcher.cMatchAtomCharge);
		StereoMolecule fragment = new StereoMolecule();
		for (int i=0; i<sIncrementTable.getSize(); i++) {
			if (threadMaster != null && threadMaster.threadMustDie())
				return cDruglikenessUnknown;

			Thread.yield();

			new IDCodeParser(false).parse(fragment, sIncrementTable.getFragment(i));
			sss.setMol(fragment, testMolecule);
			if (sss.isFragmentInMolecule()) {
				double increment = sIncrementTable.getIncrement(i);
				if (increment < -1)
					nastyIncrementSum += increment;
				else {
					incrementSum += increment;
					fragmentCount++;
					}

				detail.add(sIncrementTable.getFragment(i),
								   ParameterizedStringList.cStringTypeIDCode);
				detail.add(Double.toString(increment),
								   ParameterizedStringList.cStringTypeDouble);
				}
			}

		if (fragmentCount == 0)
			return -1;

		double druglikeness = nastyIncrementSum + incrementSum / Math.sqrt(fragmentCount);

			// correct cut-off by also treating molecules
			// with more than 50 found fragments as more drug-like
		druglikeness = druglikeness + 0.0625 * (fragmentCount - 40);

		mDetail = detail;
		return druglikeness;
		}

	public String getDruglikenessString(StereoMolecule testMolecule) {
		if (!sInitialized)
			return "Druglikeness predictor not properly initialized.";

		double incrementSum = 0.0;
		int fragmentCount = 0;
		SSSearcher sss = new SSSearcher(SSSearcher.cMatchAtomCharge);
		StereoMolecule fragment = new StereoMolecule();
		for (int i=0; i<sIncrementTable.getSize(); i++) {
			new IDCodeParser(false).parse(fragment, sIncrementTable.getFragment(i));
			sss.setMol(fragment, testMolecule);
			if (sss.isFragmentInMolecule()) {
				incrementSum += sIncrementTable.getIncrement(i);
				fragmentCount++;
				}
			}

		double druglikeness = (fragmentCount == 0) ? -1 : incrementSum / Math.sqrt(fragmentCount);
		return druglikeness + "\t" + fragmentCount + "\t" + testMolecule.getAtoms();
		}

	/**
	 * If assessDruglikeness() was called multiple times in multiple threads, then
	 * getDetail() won't retrieve the expected detail.
	 * @return
	 */
	public ParameterizedStringList getDetail() {
		return mDetail;
		}
	}
