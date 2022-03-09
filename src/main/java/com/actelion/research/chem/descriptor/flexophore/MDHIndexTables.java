/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * @author Modest v. Korff
 */

package com.actelion.research.chem.descriptor.flexophore;

import java.util.ArrayList;
import java.util.List;

public class MDHIndexTables {
	
	
	private static MDHIndexTables INSTANCE = null;
	
	private static final int MAX_SIZE = 50;
	
	private List<int [][]> liArrAtPairsBonds;
	
	private List<int [][]> liConnectionTable;
	
	private MDHIndexTables() {
		init();
	}
	
	public static MDHIndexTables getInstance(){
		if(INSTANCE==null)
			INSTANCE=new MDHIndexTables();
		
		return INSTANCE;
	}
	
	private void init(){
		
		liArrAtPairsBonds = new ArrayList<int[][]>();
		liArrAtPairsBonds.add(0, null);
		liArrAtPairsBonds.add(1, null);
		for (int i = 2; i < MAX_SIZE+1; i++) {
			liArrAtPairsBonds.add(i, createBondTable(i));
		}
		
		liConnectionTable = new ArrayList<int[][]>();
		liConnectionTable.add(0, null);
		liConnectionTable.add(1, null);
		for (int i = 2; i < MAX_SIZE+1; i++) {
			liConnectionTable.add(i, createConnectionTable(i));
		}
	}
	
	private int [][] createBondTable(int nodes){
		
		int bonds = ((nodes * nodes)-nodes) / 2;
		
		int [][] arrAtPairsBonds = new int [2][bonds];
		int cc=0;
		for (int i = 0; i < nodes; i++) {
			for (int j = i+1; j < nodes; j++) {
				arrAtPairsBonds[0][cc]=i;
				arrAtPairsBonds[1][cc]=j;
				cc++;
			}
		}
		
		return arrAtPairsBonds;
	}
	
	private int [][] createConnectionTable(int nodes){
		
		int [][] arrAtPairsBonds = liArrAtPairsBonds.get(nodes);
		
		int bonds = ((nodes * nodes)-nodes) / 2;
		int [][] arrConnBond = new int [nodes][nodes-1];
		int [] arrNumBonds = new int [nodes];
		for(int bnd=0; bnd< bonds; bnd++) {
			int at1 = arrAtPairsBonds[0][bnd];
			int at2 = arrAtPairsBonds[1][bnd];
			
			arrConnBond[at1][arrNumBonds[at1]++]=bnd;
			arrConnBond[at2][arrNumBonds[at2]++]=bnd;
			
		}
		return arrConnBond;
	}
	
	public int [][] getAtomPairsBondsTable(int nodes){
		return liArrAtPairsBonds.get(nodes);
	}
	
	public int [][] getConnectionTable(int nodes){
		return liConnectionTable.get(nodes);
	}

}
