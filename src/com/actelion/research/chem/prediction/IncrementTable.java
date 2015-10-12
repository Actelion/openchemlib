/*
 * @(#)IncrementTable.java
 *
 * Copyright 1997-2001 Actelion Ltd., Inc. All Rights Reserved.
 * 
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 * 
 * @author Thomas Sander
 */

package com.actelion.research.chem.prediction;

import java.io.*;
import java.util.ArrayList;

public class IncrementTable {
	ArrayList<IncrementTableRecord> mRecords;

	protected IncrementTable() {
		mRecords = new ArrayList<IncrementTableRecord>();
		}

	protected IncrementTable(String filename) throws Exception {
		BufferedReader theReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(filename)));
		mRecords = new ArrayList<IncrementTableRecord>();
		while (true) {
			String theLine = theReader.readLine();
			if (theLine == null)
				break;

			int tabPosition = theLine.indexOf('\t');
			if (tabPosition == -1)
				throw new Exception("line without TAB");

			String idcode = theLine.substring(0, tabPosition);
			double increment = Double.valueOf(theLine.substring(tabPosition+1)).doubleValue();

			mRecords.add(new IncrementTableRecord(idcode, increment));
			}
		theReader.close();
		}


	protected void addElement(String idcode, double increment) {
		mRecords.add(new IncrementTableRecord(idcode, increment));
		}


	protected int getSize() {
		return mRecords.size();
		}


	protected String getFragment(int i) {
		return mRecords.get(i).mIDCode;
		}


	protected double getIncrement(int i) {
		return mRecords.get(i).mIncrement;
		}
	}


class IncrementTableRecord {
	String	mIDCode;
	double	mIncrement;

	protected IncrementTableRecord(String idcode, double increment) {
		mIDCode = idcode;
		mIncrement = increment;
		}
	}