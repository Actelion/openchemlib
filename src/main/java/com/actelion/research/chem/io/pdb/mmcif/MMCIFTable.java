package com.actelion.research.chem.io.pdb.mmcif;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;

public class MMCIFTable {
	private String mName,mLine;
	private String[] mHeaderName;
	private TreeMap<String,Integer> mHeaderNameToIndexMap;

	public MMCIFTable(BufferedReader reader) throws IOException {
		mHeaderNameToIndexMap = new TreeMap<>();
		ArrayList<String> list = new ArrayList<>();
		mLine = reader.readLine();
		while (mLine.startsWith("_")) {
			mLine = mLine.trim();
			int index = mLine.indexOf('.');
			if (index == -1)
				throw new IOException("Missing '.' in table header line.");
			if (mName == null)
				mName = mLine.substring(0, index);
			else if (!mName.equals(mLine.substring(0, index)))
				throw new IOException("Inconsistent prefix in table header line.");
			if (mLine.indexOf('.', index+1) != -1)
				throw new IOException("Multiple '.' found in table header line.");
			String headerName = mLine.substring(index+1);
			mHeaderNameToIndexMap.put(headerName, list.size());
			list.add(headerName);
			mLine = reader.readLine();
		}
		mHeaderName = list.toArray(new String[0]);
	}

	public String getName() {
		return mName;
	}

	public String[] getHeaderNames() {
		return mHeaderName;
	}

	public int getIndex(String headerName) {
		return mHeaderNameToIndexMap.get(headerName);
	}

	public String[] parseRow(BufferedReader reader) throws IOException {
		if (mLine.startsWith("#"))
			return null;
		String[] row = mLine.split("\\s+");
		if (row.length != mHeaderName.length)
			throw new IOException("Inconsistent length of row entries ("+row.length+" and table headers ("+mHeaderName.length+").");
		mLine = reader.readLine();
		return row;
	}
}
