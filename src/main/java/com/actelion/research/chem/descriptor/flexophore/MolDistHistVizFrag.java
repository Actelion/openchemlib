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

import com.actelion.research.chem.Molecule3D;

public class MolDistHistVizFrag extends MolDistHistViz {

	
	private static final long serialVersionUID = 24112010;

	private int [] arrIndexParentNodes;
	
	public MolDistHistVizFrag() {
		// TODO Auto-generated constructor stub
	}

	public MolDistHistVizFrag(int length) {
		super(length);
		// TODO Auto-generated constructor stub
	}

	public MolDistHistVizFrag(int length, Molecule3D ff) {
		super(length, ff);
		// TODO Auto-generated constructor stub
	}

	public MolDistHistVizFrag(MolDistHistVizFrag mdhvf) {
		super(mdhvf);
		if(mdhvf.arrIndexParentNodes != null) {
			arrIndexParentNodes = new int [mdhvf.arrIndexParentNodes.length];
			System.arraycopy(mdhvf.arrIndexParentNodes, 0, arrIndexParentNodes, 0, mdhvf.arrIndexParentNodes.length);
		}
	}

	public MolDistHistVizFrag(MolDistHist mdh) {
		super(mdh);
		// TODO Auto-generated constructor stub
	}

	public int[] getArrIndexParentNodes() {
		return arrIndexParentNodes;
	}
	
	public MolDistHistVizFrag copy(){
		
		MolDistHistVizFrag mdhvf = new MolDistHistVizFrag(this);
		
		if(arrIndexParentNodes != null) {
			int [] arr = new int [arrIndexParentNodes.length];
			System.arraycopy(arrIndexParentNodes, 0, arr, 0, arrIndexParentNodes.length);
			mdhvf.setArrIndexParentNodes(arr);
		}
		
		return mdhvf;
	}

	public void setArrIndexParentNodes(int[] arrIndexParentNodes) {
		this.arrIndexParentNodes = arrIndexParentNodes;
	}

}
