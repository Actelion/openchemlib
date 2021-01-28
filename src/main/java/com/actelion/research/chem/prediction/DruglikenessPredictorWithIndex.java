/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
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
import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.StereoMolecule;

public class DruglikenessPredictorWithIndex {
	public static final double cDruglikenessUnknown = -999;

	private static IncrementTableWithIndex	sIncrementTable;
    private static boolean					sInitialized = false;
	private static StereoMolecule[]			sFragmentList;

	public DruglikenessPredictorWithIndex() {
		synchronized(DruglikenessPredictorWithIndex.class) {
			if (!sInitialized) {
		        try {
					sIncrementTable = new IncrementTableWithIndex("/resources/druglikeness.txt");
					sFragmentList = new StereoMolecule[sIncrementTable.getSize()];
					for (int i=0; i<sIncrementTable.getSize(); i++)
						sFragmentList[i] = new IDCodeParser(false).getCompactMolecule(sIncrementTable.getFragment(i));
					sInitialized = true;
					}
				catch (Exception e) {
		            System.out.println("Unable to initialize DruglikenessPredictor");
					}
				}
			}
		}

	public double assessDruglikeness(StereoMolecule mol, long[] index, ThreadMaster threadMaster) {
		if (!sInitialized)
			return cDruglikenessUnknown;

		double nastyIncrementSum = 0.0;
		double incrementSum = 0.0;
		int fragmentCount = 0;
		SSSearcherWithIndex swi = new SSSearcherWithIndex(SSSearcher.cMatchAtomCharge);
		swi.setMolecule(mol, index);
		for (int i=0; i<sIncrementTable.getSize(); i++) {
			if (threadMaster != null && threadMaster.threadMustDie())
				return cDruglikenessUnknown;

			swi.setFragment(sFragmentList[i], sIncrementTable.getIndex(i));
			if (swi.isFragmentInMolecule()) {
				double increment = sIncrementTable.getIncrement(i);
				if (increment < -1)
					nastyIncrementSum += increment;
				else {
					incrementSum += increment;
					fragmentCount++;
					}
				}
			}

		if (fragmentCount == 0)
			return -1;

		double druglikeness = nastyIncrementSum + incrementSum / Math.sqrt(fragmentCount);

			// correct cut-off by also treating molecules
			// with more than 50 found fragments as more drug-like
		druglikeness = druglikeness + 0.0625 * (fragmentCount - 40);

		return druglikeness;
		}
	}
