package com.actelion.research.chem.io.pdb.mmcif;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.TreeMap;

public class MMCIFBlock extends MMCIFElement {
	private final String mName;
	private final TreeMap<String,String> mEntryMap;

	public MMCIFBlock(String line, BufferedReader reader) throws IOException {
		init(line, reader);
		int index1 = line.indexOf('.');
		mName = line.substring(1, index1);
		mEntryMap = new TreeMap<>();
		String fullKey = nextValue();
		while (!isDone()) {
			if (!fullKey.substring(1).startsWith(mName+"."))
				throw new IOException("Unexpected key in '"+mName+"' block.");

			String key = fullKey.substring(mName.length()+2);

			String value = nextValue();
			if (value == null)
				throw new IOException("Missing value for '"+key+"' in '"+mName+"' block.");

			mEntryMap.put(key, value);
			fullKey = nextValue();
		}
	}

	public boolean is(String name) {
		return mName.equals(name);
	}

	public String get(String key) {
		return mEntryMap.get(key);
	}
}
