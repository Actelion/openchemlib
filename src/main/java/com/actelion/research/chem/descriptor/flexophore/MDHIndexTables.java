/*
 * Copyright (c) 2020.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
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
