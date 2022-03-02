/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
