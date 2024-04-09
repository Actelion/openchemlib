package com.actelion.research.chem.mcs;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.properties.complexity.IBitArray;
import com.actelion.research.chem.shredder.Fragment;

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
public class BondVector2IdCode {
	
	private StereoMolecule mol;
	
	private List<int []> liRingSets;
	
	public BondVector2IdCode(StereoMolecule mol) {
		this.mol = mol;
		
		liRingSets = new ArrayList<int []>();
		RingCollection rc = mol.getRingSet();
		int rings = rc.getSize();
		for (int i = 0; i < rings; i++) {
			int [] arrIndexBnd = rc.getRingBonds(i);
			liRingSets.add(arrIndexBnd);
		}
	}
	
	/**
	 * 
	 * @param fragDefByBonds
	 * @return true if the fragment contains parts of a ring.
	 */
	public boolean containsFragmentOpenRing(IBitArray fragDefByBonds){
		
		boolean openRing = false;
		for(int [] arrIndexBnd : liRingSets){
			int ccOverlap=0;
			for (int i = 0; i < arrIndexBnd.length; i++) {
				if(fragDefByBonds.isBitSet(arrIndexBnd[i])){
					ccOverlap++;
				}
			}
			if((ccOverlap > 0) && ccOverlap < arrIndexBnd.length) {
				openRing = true;
				break;
			}
		}
		
		return openRing;
	}
	public String getFragmentIdCode(IBitArray fragDefByBonds){
		
		StereoMolecule frag = convert(fragDefByBonds, false);
		Canonizer can = new Canonizer(frag);
		String idcode = can.getIDCode();
		return idcode; 
	}
	public Fragment getFragment(IBitArray fragDefByBonds){
				
		StereoMolecule frag = convert(fragDefByBonds, false);
		Canonizer can = new Canonizer(frag);
		String idcode = can.getIDCode();
		Fragment fragment = new Fragment(idcode);
		fragment.setMol(frag);
		fragment.setSize(frag.getBonds());
		return fragment; 
	}
	
	public Fragment getFragment(IBitArray fragDefByBonds, boolean addWildcards){
		
		StereoMolecule frag = convert(fragDefByBonds, addWildcards);
		Canonizer can = new Canonizer(frag);
		String idcode = can.getIDCode();
		Fragment fragment = new Fragment(idcode);
		fragment.setMol(frag);
		fragment.setSize(frag.getBonds());
		return fragment; 
	}
	
	private StereoMolecule convert(IBitArray fragDefByBonds, boolean addWildcards){
		
		int bonds = mol.getBonds();
		int atoms = mol.getAtoms();
		boolean [] arrBonds = new boolean [bonds];
		boolean [] arrAtoms = new boolean [atoms];
		int bondsFragment = 0;
		
		for (int i = 0; i < bonds; i++) {
			if(fragDefByBonds.isBitSet(i)){
				arrBonds[i] = true;
				bondsFragment++;
				arrAtoms[mol.getBondAtom(0, i)] = true;
				arrAtoms[mol.getBondAtom(1, i)] = true;
			}
		}
		
		int atomsFrag = 0;
		for (int i = 0; i < arrAtoms.length; i++) {
			if(arrAtoms[i]){
				atomsFrag++;
			}
		}

		StereoMolecule fragSubBonds = new StereoMolecule(atomsFrag, bondsFragment);
		int [] indexAtoms = mol.copyMoleculeByBonds(fragSubBonds, arrBonds, true, null);

		// Add ring and aromatic info.
		// Added 07.07.2014
		// 
		
		int indexAtomNew = 0;
		for (int i = 0; i < indexAtoms.length; i++) {
			if(indexAtoms[i]>-1) {
			if((mol.getAtomQueryFeatures(indexAtoms[i]) & Molecule.cAtomQFNotChain) > 0){
				fragSubBonds.setAtomQueryFeature(indexAtomNew, Molecule.cAtomQFNotChain, true);
			}
			if((mol.getAtomQueryFeatures(indexAtoms[i]) & Molecule.cAtomQFAromatic) > 0){
				fragSubBonds.setAtomQueryFeature(indexAtomNew, Molecule.cAtomQFAromatic, true);
			}
			indexAtomNew++;
			}
		}
				
		if(addWildcards) {
			boolean [] arrAtomCopied2Fragment = new boolean [mol.getAtoms()];
			for (int i = 0; i < indexAtoms.length; i++) {
				if(indexAtoms[i] > -1)
					arrAtomCopied2Fragment[i] = true;
			}
			
			for (int i = 0; i < indexAtoms.length; i++) {
				if(indexAtoms[i] > -1) {
					
					int atIndexOld = i;
					int nConnected = mol.getConnAtoms(atIndexOld);
					for (int j = 0; j < nConnected; j++) {
						int indexAtConn = mol.getConnAtom(atIndexOld, j);
						if(!arrAtomCopied2Fragment[indexAtConn]){
							int atWildCard = fragSubBonds.addAtom(0);
							int atIndexNew = indexAtoms[i];
							fragSubBonds.addBond(atIndexNew, atWildCard, Molecule.cBondTypeSingle);
							fragSubBonds.setAtomQueryFeature(atWildCard, Molecule.cAtomQFAny, true);
						}
					}
				}
			}
		}
		fragSubBonds.ensureHelperArrays(Molecule.cHelperRings);
		return fragSubBonds;
	}
	
	public String getFragmentIdCodeCarbonSkeleton(IBitArray fragDefByBonds){
		StereoMolecule frag = convert(fragDefByBonds, false);
		for (int i = 0; i < frag.getAtoms(); i++) {
			frag.setAtomicNo(i, 6);
		}
		Canonizer can = new Canonizer(frag);
		String idcode = can.getIDCode();
		return idcode; 
	}
}
