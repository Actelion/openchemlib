package com.actelion.research.io;

import java.io.IOException;
import java.io.Reader;

public class BOMSkipper {
	/**
	 * Some text editors (e.g. Notepad on Windows) write a BOM as first character,
	 * when writing an UTF-8 encoded text file. BOMs are unicode characters encoded
	 * as 2,3, or 4-byte sequence specific for the unicode type (UTF-8, 16, ...) and
	 * little vs. big-endian byte order. Java Readers and InputStreams don't filter
	 * BOM out of the stream. Thus, we need to do it ourselfes.
	 * @param reader
	 */
	public static void skip(Reader reader) {
		try {
			reader.mark(1);
			char[] possibleBOM = new char[1];

			if (reader.read(possibleBOM) == 1 && possibleBOM[0] != '\ufeff')
				reader.reset();
			}
		catch (IOException ioe) {}
		}
	}