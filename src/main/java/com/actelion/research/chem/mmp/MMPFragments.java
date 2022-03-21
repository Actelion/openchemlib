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

import com.actelion.research.chem.mmp.MMPFragmenter.MoleculeIndexID;

import java.io.*;
import java.util.List;


public class MMPFragments {
	private File temp = File.createTempFile("fragments", ".tmp");
	private PrintWriter fragmentsWriter = new PrintWriter(temp);
	private int numberOfLines = 0;
	
	public MMPFragments() throws IOException {
		fragmentsWriter.println("key1FragmentIndex\tkey2FragmentIndex\tvalueFragmentIndex\tcutType\tmoleculeIndex");
	}
	
	/**
	 * Adds all molecule fragmentation variations
	 * @param moleculeIndex Index of the molecule
	 * @param moleculeIndexesID All fragmentation variations
	 */
	public void addFragments(int moleculeIndex, List<MoleculeIndexID> moleculeIndexesID) {
		for (MoleculeIndexID moleculeIndexID: moleculeIndexesID) {
			addFragments(moleculeIndex, moleculeIndexID);
		}
	}
	
	/**
	 * Adds a new molecule fragment
	 * @param moleculeIndex Index of the molecule
	 * @param moleculeIndexID One fragmentation variation of the molecule
	 */
	public void addFragments(int moleculeIndex, MoleculeIndexID moleculeIndexID) {
		int[] keysSize = moleculeIndexID.getKeysIDAtoms();
		int[] keysIndex = moleculeIndexID.getKeysIndex();
		if (keysSize.length == 1) { // single cut
			fragmentsWriter.println(Integer.toString(keysIndex[0]) + "\t\t" + Integer.toString(moleculeIndexID.getValueIndex()) + "\t1\t" + Integer.toString(moleculeIndex));
		}
		else { // double cut
			fragmentsWriter.println(Integer.toString(keysIndex[0]) + "\t" + Integer.toString(keysIndex[1]) + "\t" +  Integer.toString(moleculeIndexID.getValueIndex()) + "\t2\t" + Integer.toString(moleculeIndex));
		}
		numberOfLines++;
	}
	
	/**
	 * Writes the Molecules Fragments block
	 * @param printWriter
	 */
	public void writeFragments(PrintWriter printWriter) throws IOException  {
		fragmentsWriter.close();
		printWriter.println("<mmpFragments>");
		printWriter.println("<column properties>");
		printWriter.println("<columnName=\"key1FragmentIndex\">");
		printWriter.println("<columnName=\"key2FragmentIndex\">");
		printWriter.println("<columnName=\"valueFragmentIndex\">");
		printWriter.println("<columnName=\"cutType\">");
		printWriter.println("<columnName=\"moleculeIndex\">");
		printWriter.println("</column properties>");
		BufferedReader br = new BufferedReader(new FileReader(temp));
		String strLine;
		while ((strLine = br.readLine()) != null)   {
			printWriter.println(strLine);
		}
		br.close();
		temp.delete();
		printWriter.println("</mmpFragments>");
	}
	
	/**
	 * Get the number of fragments
	 * @return
	 */
	public int getFragmentsCount() {
		return numberOfLines;
	}	
}
