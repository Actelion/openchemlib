/*
 * Copyright (c) 2017
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
 * @author Gregori Gerebtzoff
 */

package com.actelion.research.chem.mmp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MMPairs {
	private File temp = File.createTempFile("matchedmolecularpairs", ".tmp");
	private PrintWriter mmpWriter = new PrintWriter(temp);
	private int numberOfLines = 0;
	
	MMPairs() throws IOException {
		mmpWriter.println("value1FragmentIndex\tvalue1Atoms\tvalue2FragmentIndex\tvalue2Atoms\tcutType\tnumberOfExamples\texamples");
	}
	
	/**
	 * Writes the Matched Molecular Pairs into a temporary file
	 * @param mMPs Matched Molecular Pairs to be written
	 * @throws IOException
	 */
	public void writeMMPEnumeration(HashMap<String, List<String[]>> mMPs) throws IOException {
		if (mMPs != null && mMPs.size() > 0) { 
			for (Map.Entry<String, List<String[]>> entry: mMPs.entrySet()) {
				String examples = "";
				for (String[] example: entry.getValue()) {
					if (examples == "") {
						examples = example[0] + "," + example[1];
					}
					else {
						examples += "|" + example[0] + "," + example[1];
					}
				}
				mmpWriter.println(entry.getKey() + "\t" + Integer.toString(entry.getValue().size()) + "\t" + examples);
				numberOfLines++;
			}
		}
	}
	
	/**
	 * Get the number of Matched Molecular Pairs
	 * @return
	 */
	public int getMMPsCount() {
		return numberOfLines;
	}
	
	/**
	 * Writes the Matched Molecular Pairs block
	 * @param printWriter
	 */
	public void writeMMPs(PrintWriter printWriter) throws IOException  {
		mmpWriter.close();
		printWriter.println("<matchedMolecularPairs>");
		printWriter.println("<column properties>");
		printWriter.println("<columnName=\"value1FragmentIndex\">");
		printWriter.println("<columnName=\"value1Atoms\">");
		printWriter.println("<columnName=\"value2FragmentIndex\">");
		printWriter.println("<columnName=\"value2Atoms\">");
		printWriter.println("<columnName=\"cutType\">");
		printWriter.println("<columnName=\"numberOfExamples\">");
		printWriter.println("<columnName=\"examples\">");
		printWriter.println("</column properties>");
		BufferedReader br = new BufferedReader(new FileReader(temp));
		String strLine;
		while ((strLine = br.readLine()) != null)   {
			printWriter.println(strLine);
		}
		br.close();
		temp.delete();
		printWriter.println("</matchedMolecularPairs>");
	}
}
