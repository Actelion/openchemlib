/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
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
* 3. Neither the name of the the copyright holder nor the
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
*/

package com.actelion.research.chem.io;

import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.MolfileV3Creator;
import com.actelion.research.chem.StereoMolecule;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

/**
 * A streaming writer for SD files (.sdf) with optional gzip compression (.sdf.gz).
 * Molecules are written one at a time, so memory usage stays constant regardless
 * of the number of molecules being written.
 *
 * <p>Usage example:</p>
 * <pre>
 * try (SDFileWriter writer = new SDFileWriter("output.sdf.gz")) {
 *     for (StereoMolecule mol : molecules) {
 *         writer.writeMolecule(mol);
 *     }
 * }
 * </pre>
 *
 * <p>To include SD properties:</p>
 * <pre>
 * try (SDFileWriter writer = new SDFileWriter("output.sdf")) {
 *     for (StereoMolecule mol : molecules) {
 *         Map&lt;String, String&gt; properties = new LinkedHashMap&lt;&gt;();
 *         properties.put("MOLECULAR_WEIGHT", "180.16");
 *         properties.put("COMPOUND_ID", "ABC-123");
 *         writer.writeMolecule(mol, properties);
 *     }
 * }
 * </pre>
 */
public class SDFileWriter implements AutoCloseable {
	private static final String NEWLINE = "\n";
	private static final String RECORD_SEPARATOR = "$$$$";

	private final BufferedWriter mWriter;
	private final boolean mUseV3000;
	private boolean mClosed;

	/**
	 * Creates an SDFileWriter for the given file path.
	 * If the filename ends with ".gz", output will be gzip-compressed.
	 * Uses V2000 molfile format by default.
	 * @param fileName path to the output file
	 * @throws IOException if the file cannot be opened for writing
	 */
	public SDFileWriter(String fileName) throws IOException {
		this(fileName, false);
	}

	/**
	 * Creates an SDFileWriter for the given file path.
	 * If the filename ends with ".gz", output will be gzip-compressed.
	 * @param fileName path to the output file
	 * @param useV3000 if true, uses V3000 molfile format; otherwise V2000
	 * @throws IOException if the file cannot be opened for writing
	 */
	public SDFileWriter(String fileName, boolean useV3000) throws IOException {
		this(createOutputStream(fileName), useV3000);
	}

	/**
	 * Creates an SDFileWriter for the given File.
	 * If the filename ends with ".gz", output will be gzip-compressed.
	 * Uses V2000 molfile format by default.
	 * @param file the output file
	 * @throws IOException if the file cannot be opened for writing
	 */
	public SDFileWriter(File file) throws IOException {
		this(file, false);
	}

	/**
	 * Creates an SDFileWriter for the given File.
	 * If the filename ends with ".gz", output will be gzip-compressed.
	 * @param file the output file
	 * @param useV3000 if true, uses V3000 molfile format; otherwise V2000
	 * @throws IOException if the file cannot be opened for writing
	 */
	public SDFileWriter(File file, boolean useV3000) throws IOException {
		this(createOutputStream(file), useV3000);
	}

	/**
	 * Creates an SDFileWriter that writes to the given Writer.
	 * The caller is responsible for any compression wrapping.
	 * Uses V2000 molfile format by default.
	 * @param writer the output writer
	 */
	public SDFileWriter(Writer writer) {
		this(writer, false);
	}

	/**
	 * Creates an SDFileWriter that writes to the given Writer.
	 * The caller is responsible for any compression wrapping.
	 * @param writer the output writer
	 * @param useV3000 if true, uses V3000 molfile format; otherwise V2000
	 */
	public SDFileWriter(Writer writer, boolean useV3000) {
		mWriter = (writer instanceof BufferedWriter) ? (BufferedWriter) writer : new BufferedWriter(writer);
		mUseV3000 = useV3000;
		mClosed = false;
	}

	/**
	 * Creates an SDFileWriter that writes to the given OutputStream.
	 * The caller is responsible for any compression wrapping.
	 * Uses V2000 molfile format by default.
	 * @param outputStream the output stream
	 */
	public SDFileWriter(OutputStream outputStream) {
		this(outputStream, false);
	}

	/**
	 * Creates an SDFileWriter that writes to the given OutputStream.
	 * The caller is responsible for any compression wrapping.
	 * @param outputStream the output stream
	 * @param useV3000 if true, uses V3000 molfile format; otherwise V2000
	 */
	public SDFileWriter(OutputStream outputStream, boolean useV3000) {
		this(new OutputStreamWriter(outputStream, StandardCharsets.UTF_8), useV3000);
	}

	/**
	 * Writes a molecule as an SD record without any properties.
	 * @param mol the molecule to write
	 * @throws IOException if an I/O error occurs
	 */
	public void writeMolecule(StereoMolecule mol) throws IOException {
		writeMolecule(mol, null);
	}

	/**
	 * Writes a molecule as an SD record with the given properties.
	 * @param mol the molecule to write
	 * @param properties SD field name-value pairs (may be null); use a LinkedHashMap to preserve order
	 * @throws IOException if an I/O error occurs
	 */
	public void writeMolecule(StereoMolecule mol, Map<String, String> properties) throws IOException {
		if (mClosed)
			throw new IOException("SDFileWriter is already closed");

		// Write the molfile block
		if (mUseV3000)
			new MolfileV3Creator(mol).writeMolfile(mWriter);
		else
			new MolfileCreator(mol).writeMolfile(mWriter);

		// Write SD properties
		if (properties != null) {
			for (Map.Entry<String, String> entry : properties.entrySet()) {
				mWriter.write(">  <" + entry.getKey() + ">");
				mWriter.write(NEWLINE);
				String value = entry.getValue();
				if (value != null) {
					mWriter.write(value);
					mWriter.write(NEWLINE);
				}
				mWriter.write(NEWLINE);
			}
		}

		// Write record separator
		mWriter.write(RECORD_SEPARATOR);
		mWriter.write(NEWLINE);
	}

	/**
	 * Flushes the underlying writer without closing it.
	 * @throws IOException if an I/O error occurs
	 */
	public void flush() throws IOException {
		if (!mClosed)
			mWriter.flush();
	}

	/**
	 * Closes the writer and releases resources.
	 * @throws IOException if an I/O error occurs
	 */
	@Override
	public void close() throws IOException {
		if (!mClosed) {
			mClosed = true;
			mWriter.close();
		}
	}

	private static OutputStream createOutputStream(String fileName) throws IOException {
		OutputStream os = new FileOutputStream(fileName);
		if (fileName.toLowerCase().endsWith(".gz"))
			os = new GZIPOutputStream(os);
		return os;
	}

	private static OutputStream createOutputStream(File file) throws IOException {
		OutputStream os = new FileOutputStream(file);
		if (file.getName().toLowerCase().endsWith(".gz"))
			os = new GZIPOutputStream(os);
		return os;
	}
}
