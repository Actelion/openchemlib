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

import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

public class IncrementTableWithIndex {
	private static final String cHeader = "<index version "+SSSearcherWithIndex.cIndexVersion+">";

    ArrayList<IncrementTableRecordWithIndex>	mRecords;

	protected IncrementTableWithIndex() {
		mRecords = new ArrayList<IncrementTableRecordWithIndex>();
		}


	protected IncrementTableWithIndex(String filename) throws Exception {
		BufferedReader theReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(filename), StandardCharsets.UTF_8));

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