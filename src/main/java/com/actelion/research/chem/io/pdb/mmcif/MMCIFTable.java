package com.actelion.research.chem.io.pdb.mmcif;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;

public class MMCIFTable extends MMCIFElement {
	private String mName;
	private final String[] mHeaderName;
	private final TreeMap<String,Integer> mHeaderNameToIndexMap;

	public MMCIFTable(BufferedReader reader) throws IOException {
		mHeaderNameToIndexMap = new TreeMap<>();
		ArrayList<String> list = new ArrayList<>();
		String line = reader.readLine();
		while (line.startsWith("_")) {
			line = line.trim();
			int index = line.indexOf('.');
			if (index == -1)
				throw new IOException("Missing '.' in table header line.");
			if (mName == null)
				mName = line.substring(0, index);
			else if (!mName.equals(line.substring(0, index)))
				throw new IOException("Inconsistent prefix in table header line.");
			if (line.indexOf('.', index+1) != -1)
				throw new IOException("Multiple '.' found in table header line.");
			String headerName = line.substring(index+1);
			mHeaderNameToIndexMap.put(headerName, list.size());
			list.add(headerName);
			line = reader.readLine();
		}
		mHeaderName = list.toArray(new String[0]);
		super.init(line, reader);
	}

	public String getName() {
		return mName;
	}

	public String[] getHeaderNames() {
		return mHeaderName;
	}

	public int getIndex(String headerName) {
		Integer i = mHeaderNameToIndexMap.get(headerName);
		return i == null ? -1 : i;
	}

	public String[] parseRow(BufferedReader reader) throws IOException {
		if (isDone())
			return null;

		String[] row = new String[mHeaderName.length];
		for (int i=0; i<row.length; i++)
			row[i] = nextValue();

		if (row[0] == null)
			return null;

		if (row[row.length-1] == null)
			throw new IOException("Insufficient values for table '"+mName+"'.");

		return row;
	}

	public void skip(BufferedReader reader) throws IOException {
		while (!reader.readLine().trim().equals("#"));
	}
}
