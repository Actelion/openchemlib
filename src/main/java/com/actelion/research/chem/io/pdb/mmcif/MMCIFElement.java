package com.actelion.research.chem.io.pdb.mmcif;

import java.io.BufferedReader;
import java.io.IOException;

public abstract class MMCIFElement {
	private BufferedReader mReader;
	private String mLine;
	private int mLineIndex;
	private boolean mDone;

	public void init(String line, BufferedReader reader) {
		mReader = reader;
		mLine = line;
		mLineIndex = 0;
		mDone = false;
	}

	protected String nextValue() throws IOException {
		if (!mDone) {
			while (mLineIndex < mLine.length() && mLine.charAt(mLineIndex) == ' ')
				mLineIndex++;
			if (mLineIndex >= mLine.length())
				readNextLine();
		}

		if (mDone)
			return null;

		if (mLineIndex == 0 && mLine.startsWith(";")) {
			StringBuilder sb = new StringBuilder(mLine.substring(1));
			readNextLine();
			while (!mDone && !mLine.equals(";")) {
				sb.append(mLine);
				readNextLine();
			}
			readNextLine();
			return sb.toString();
		}

		if (mLine.charAt(mLineIndex) == '\'') {
			mLineIndex++;
			int startIndex = mLineIndex;
			while (mLineIndex < mLine.length() && mLine.charAt(mLineIndex) != '\'')
				mLineIndex++;
			mLineIndex++;
			return mLine.substring(startIndex, mLineIndex-1);
		}

		int startIndex = mLineIndex;
		while (mLineIndex < mLine.length() && mLine.charAt(mLineIndex) != ' ')
			mLineIndex++;
		return mLine.substring(startIndex, mLineIndex);
	}

	private void readNextLine() throws IOException {
		mLine = mReader.readLine().trim();
		mLineIndex = 0;
		if (mLine.equals("#"))
			mDone = true;
	}

	public boolean isDone() {
		return mDone;
	}
}
