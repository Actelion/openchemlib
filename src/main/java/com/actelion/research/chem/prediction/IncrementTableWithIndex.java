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

import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class IncrementTableWithIndex {
	private static final String cHeader = "<index version "+SSSearcherWithIndex.cIndexVersion+">";

    ArrayList<IncrementTableRecordWithIndex>	mRecords;

	protected IncrementTableWithIndex() {
		mRecords = new ArrayList<IncrementTableRecordWithIndex>();
		}


	protected IncrementTableWithIndex(String filename) throws Exception {
		BufferedReader theReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(filename)));

		String header = theReader.readLine();
		if (!header.equals(cHeader))
			throw new Exception("index version mismatch");

        DescriptorHandlerLongFFP512 descriptorHandler = new DescriptorHandlerLongFFP512();
		mRecords = new ArrayList<IncrementTableRecordWithIndex>();
		while (true) {
			String theLine = theReader.readLine();
			if (theLine == null)
				break;

			int firstTab = theLine.indexOf('\t');
			if (firstTab == -1)
				throw new Exception("line without TAB");
			int secondTab = theLine.indexOf('\t', firstTab+1);
			if (secondTab == -1)
				throw new Exception("line without second TAB");

            long[] index = (long[])descriptorHandler.decode(theLine.substring(0, firstTab));
			String idcode = theLine.substring(firstTab+1, secondTab);
			double increment = Double.valueOf(theLine.substring(secondTab+1)).doubleValue();

			mRecords.add(new IncrementTableRecordWithIndex(idcode, index, increment));
			}
		theReader.close();
		}


	protected void addElement(String idcode, long[] index, double increment) {
		mRecords.add(new IncrementTableRecordWithIndex(idcode, index, increment));
		}


	protected int getSize() {
		return mRecords.size();
		}


	protected String getFragment(int i) {
		return mRecords.get(i).mIDCode;
		}


	protected long[] getIndex(int i) {
		return mRecords.get(i).mIndex;
		}


	protected double getIncrement(int i) {
		return mRecords.get(i).mIncrement;
		}
	}


class IncrementTableRecordWithIndex {
	String				mIDCode;
	long[]				mIndex;
	double				mIncrement;

	protected IncrementTableRecordWithIndex(String idcode, long[] index, double increment) {
		mIDCode = idcode;
		mIndex = index;
		mIncrement = increment;
		}
	}