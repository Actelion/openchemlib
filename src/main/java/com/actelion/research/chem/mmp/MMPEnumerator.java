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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

public class MMPEnumerator {
	private HashMap<String, List<String[]>> mMPEnumeration;
	
	public HashMap<String, List<String[]>> getMMPEnumeration() {
		return mMPEnumeration;
	}
	
	/**
	 * Adds new Matched Molecular Pair to mMPEnumeration
	 * @param values Tab-delimited string of [pair1, pair2, pair1Size, pair2Size, cutType]
	 * @param data Tab-delimited string of [example1a, example1b, leftKey, null|rightKey]
	 */
	private void addMMP(String values, String[] data) {
		List<String[]> datas = new ArrayList<String[]>();
		if (mMPEnumeration.containsKey(values)) {
			datas = mMPEnumeration.get(values);
		}
		datas.add(data);
		mMPEnumeration.put(values, datas);
	}
	
	public MMPEnumerator() throws IOException { }
	
	/**
	 * Creates a new MMPEnumerator
	 * @param combination [A, B] array of two combination sizes (A=size of seed, B=size of replacement; sizes are number of heavy atoms)
	 * @param keysHash1 HashMap of the keys of size A 
	 * @param keysHash2 HashMap of the keys of size B
	 * @param version Version of the enumerator (1.0 or 1.1)
	 */
	public MMPEnumerator(int[] combination, HashMap<String, ArrayList<int[]>> keysHash1, HashMap<String, ArrayList<int[]>> keysHash2, String version) throws IOException {
		mMPEnumeration = new HashMap<String, List<String[]>>();
		if ((combination[0] == combination[1] && keysHash1 != null) || (combination[0] != combination[1] && keysHash1 != null && keysHash2 != null)) {			
			if (combination[0] != combination[1]) {
				Iterator<String> it = keysHash1.keySet().iterator();
				while (it.hasNext()) {
					String keysString = it.next();
					if (keysHash2.containsKey(keysString)) {
						ArrayList<int[]> valuesList1 = keysHash1.get(keysString);
						ArrayList<int[]> valuesList2 = keysHash2.get(keysString);
						String cutType = Integer.toString(keysString.split("\t").length);
						for (int[] values1:valuesList1) {
							for (int[] values2:valuesList2) {
								addMMP(Integer.toString(values1[0]) + "\t" + Integer.toString(combination[0]) + "\t" + Integer.toString(values2[0]) + "\t" + Integer.toString(combination[1]) + "\t" + cutType, new String[]{Integer.toString(values1[1]), Integer.toString(values2[1]), keysString});
								if (version == "1.0")
									addMMP(Integer.toString(values2[0]) + "\t" + Integer.toString(combination[1]) + "\t" + Integer.toString(values1[0]) + "\t" + Integer.toString(combination[0]) + "\t" + cutType, new String[]{Integer.toString(values2[1]), Integer.toString(values1[1]), keysString});
							}
						}
					}
				}
			}
			else {
				Iterator<String> it = keysHash1.keySet().iterator();
				while (it.hasNext()) {
					String keysString = it.next();
					ArrayList<int[]> valuesList = keysHash1.get(keysString);
					String numberOfCuts = Integer.toString(keysString.split("\t").length);
					if (valuesList.size() > 1) {
						for (int i=0; i<valuesList.size()-1; i++) {
							for (int j=i+1; j<valuesList.size(); j++) {
								int[] values1 = valuesList.get(i);
								int[] values2 = valuesList.get(j);
								if (values1[0] != values2[0]) { // needed to avoid replacing one fragment by the same
									addMMP(Integer.toString(values1[0]) + "\t" + Integer.toString(combination[0]) + "\t" + Integer.toString(values2[0]) + "\t" + Integer.toString(combination[1]) + "\t" + numberOfCuts, new String[]{Integer.toString(values1[1]), Integer.toString(values2[1]), keysString});
									addMMP(Integer.toString(values2[0]) + "\t" + Integer.toString(combination[1]) + "\t" + Integer.toString(values1[0]) + "\t" + Integer.toString(combination[0]) + "\t" + numberOfCuts, new String[]{Integer.toString(values2[1]), Integer.toString(values1[1]), keysString});
								}
							}
						}
					}
				}
			}
		}
	}
}
