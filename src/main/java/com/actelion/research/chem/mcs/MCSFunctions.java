package com.actelion.research.chem.mcs;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

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
public class MCSFunctions {

	private MCS mcs;

	/**
	 * 
	 */
	public MCSFunctions() {
		mcs = new MCS();

	}
	
	public StereoMolecule getRemainingStructure(StereoMolecule molProduct, StereoMolecule molEduct1){
		
		boolean [] arrBondProduct = new boolean [molProduct.getBonds()];
		
		boolean [] arrBondReactant2 = null;
		
		mcs.set(molProduct, molEduct1);
				
		mcs.getMCSBondArray(arrBondProduct, arrBondReactant2);
		
		StereoMolecule molRemainingSubstituent = new StereoMolecule(molProduct);
				
		molRemainingSubstituent.ensureHelperArrays(Molecule.cHelperRings);
				
		boolean [] arrBondProductInverse = new boolean [molProduct.getBonds()];
		
		for (int i = 0; i < arrBondProduct.length; i++) {
			
			if(!arrBondProduct[i]){
				arrBondProductInverse[i]=true;
			}
		}
		
		molProduct.copyMoleculeByBonds(molRemainingSubstituent, arrBondProductInverse, true, null);
		
		molRemainingSubstituent.ensureHelperArrays(Molecule.cHelperRings);
				
		return molRemainingSubstituent;
		
	}
	
	public StereoMolecule getMaximumCommonSubstructure(List<StereoMolecule> li) throws Exception {
		
		StereoMolecule molMCS = li.get(0);
		
		List<StereoMolecule> liMCS = new ArrayList<StereoMolecule>();
		
		liMCS.add(molMCS);
				
		for (int i = 1; i < li.size(); i++) {
			
			StereoMolecule mol = li.get(i);
			
			mcs.set(mol, molMCS);
			
			molMCS = mcs.getMCS();
			
			if(molMCS==null){
				System.err.println("No common substructure found. break!");
			}
		}
				
		return molMCS;
		
	}
	
	public double getScore(StereoMolecule mol, StereoMolecule frag){
		
		mcs.set(mol, frag);
		
		if(mcs.getMCS()==null){
			return 0;
		}
		
		return mcs.getScore();
		
	}



}
