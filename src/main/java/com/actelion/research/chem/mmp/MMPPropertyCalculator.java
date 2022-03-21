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

import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.MolecularFlexibilityCalculator;
import com.actelion.research.chem.prediction.CLogPPredictor;
import com.actelion.research.chem.prediction.PolarSurfaceAreaPredictor;
import com.actelion.research.util.DoubleFormat;

import java.util.TreeSet;

public class MMPPropertyCalculator {
	String propertyName;
	String shortDisplayedPropertyName;
	String longDisplayedPropertyName;
	int calculatedPropertyAttributeIndex;
	String value;
	
	public MMPPropertyCalculator() {}
	
	/**
	 * @param propertyName Name of the calculated property
	 * @param calculatedPropertyAttributeIndex Column index of the calculated property
	 */
	private void CalculateCompound(String propertyName, StereoMolecule mol, int calculatedPropertyAttributeIndex) {
		this.propertyName = propertyName;
		this.calculatedPropertyAttributeIndex = calculatedPropertyAttributeIndex;
		this.value = null;
		if (propertyName.toLowerCase().equals("total_weight") || propertyName.toLowerCase().equals("mw")) {
			this.shortDisplayedPropertyName = "MW";
			this.longDisplayedPropertyName = "Molecular Weight";
			double totalWeight = new MolecularFormula(mol).getRelativeWeight();
			value = DoubleFormat.toString(totalWeight, 6, true);
		}
		else if (propertyName.toLowerCase().equals("logp")) {
			this.shortDisplayedPropertyName = "LogP";
			this.longDisplayedPropertyName = "LogP";
			CLogPPredictor predictor = new CLogPPredictor();
			float cLogP = predictor.assessCLogP(mol);
			value = Float.toString(cLogP);
		}
		else if (propertyName.toLowerCase().equals("acceptors")) {
			this.shortDisplayedPropertyName = "Acceptors";
			this.longDisplayedPropertyName = "Acceptors";
			int count = 0;
			for (int atom=0; atom<mol.getAllAtoms(); atom++)
				if (mol.getAtomicNo(atom) == 7 || mol.getAtomicNo(atom) == 8)
					count++;
			value = ""+count;
		}
		else if (propertyName.toLowerCase().equals("donors")) {
			this.shortDisplayedPropertyName = "Donors";
			this.longDisplayedPropertyName = "Donors";
			int count = 0;
			for (int atom=0; atom<mol.getAllAtoms(); atom++)
				if ((mol.getAtomicNo(atom) == 7 || mol.getAtomicNo(atom) == 8)
				 && mol.getAllHydrogens(atom) > 0)
					count++;
			value = ""+count;
		}
		else if (propertyName.toLowerCase().equals("psa")) {
			this.shortDisplayedPropertyName = "PSA";
			this.longDisplayedPropertyName = "Polar Surface Area";
			PolarSurfaceAreaPredictor predictor = new PolarSurfaceAreaPredictor();
			float psa = predictor.assessPSA(mol);
			value = Float.toString(psa);
		}
		else if (propertyName.toLowerCase().equals("shape")) {
			this.shortDisplayedPropertyName = "Shape";
			this.longDisplayedPropertyName = "Shape";
			value = DoubleFormat.toString(assessMolecularShape(mol));
		}
		else if (propertyName.toLowerCase().equals("flexibility")) {
			this.shortDisplayedPropertyName = "Flexibility";
			this.longDisplayedPropertyName = "Flexibility";
			MolecularFlexibilityCalculator predictor = new MolecularFlexibilityCalculator();
			float flexibility = predictor.calculateMolecularFlexibility(mol);
			value = Float.toString(flexibility);
		}
		else if (propertyName.toLowerCase().equals("complexity")) {
			this.shortDisplayedPropertyName = "Complexity";
			this.longDisplayedPropertyName = "Complexity";
			value = DoubleFormat.toString(assessMolecularComplexity(mol));
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("heavy_atoms")) {
			this.shortDisplayedPropertyName = "Heavy atoms";
			this.longDisplayedPropertyName = "Heavy atoms";
			value = ""+mol.getAtoms();
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("noncarbon_atoms") || propertyName.toLowerCase().equals("non-carbon atoms")) {
			this.shortDisplayedPropertyName = "Non-carbon atoms";
			this.longDisplayedPropertyName = "Non-carbon atoms";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperNeighbours);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (mol.getAtomicNo(atom) != 6)
					count++;
			value = ""+count;
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("metal_atoms")) {
			this.shortDisplayedPropertyName = "Metal atoms";
			this.longDisplayedPropertyName = "Metal atoms";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperNeighbours);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (mol.isMetalAtom(atom))
					count++;
			value = ""+count;
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("negative_atoms")) {
			this.shortDisplayedPropertyName = "Negative atoms";
			this.longDisplayedPropertyName = "Negative atoms";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperNeighbours);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (mol.isElectronegative(atom))
					count++;
			value = ""+count;
		}
		else if (propertyName.toLowerCase().equals("stereocenters")) {
			this.shortDisplayedPropertyName = "Stereocenters";
			this.longDisplayedPropertyName = "Stereocenters";
			value = ""+mol.getStereoCenterCount();
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("rotatable_bonds") || propertyName.toLowerCase().equals("rot. bonds")) {
			this.shortDisplayedPropertyName = "Rot. bonds";
			this.longDisplayedPropertyName = "Rotatable bonds";
			value = ""+mol.getRotatableBondCount();
		}
		else if (propertyName.toLowerCase().equals("rings")) {
			this.shortDisplayedPropertyName = "Rings";
			this.longDisplayedPropertyName = "Rings";
			mol.ensureHelperArrays(StereoMolecule.cHelperRings);
			value = ""+mol.getRingSet().getSize();
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("aromatic_rings") || propertyName.toLowerCase().equals("arom. rings")) {
			this.shortDisplayedPropertyName = "Arom. rings";
			this.longDisplayedPropertyName = "Aromatic rings";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperRings);
			RingCollection rc = mol.getRingSet();
			for (int i=0; i<rc.getSize(); i++)
				if (rc.isAromatic(i))
					count++;
			value = ""+count;
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("sp3_atoms")) {
			this.shortDisplayedPropertyName = "SP3 atoms";
			this.longDisplayedPropertyName = "SP3 atoms";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperRings);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if ((mol.getAtomicNo(atom) == 6 && mol.getAtomPi(atom) == 0)
				 || (mol.getAtomicNo(atom) == 7 && !mol.isFlatNitrogen(atom))
				 || (mol.getAtomicNo(atom) == 8 && mol.getAtomPi(atom) == 0 && !mol.isAromaticAtom(atom))
				 || (mol.getAtomicNo(atom) == 15)
				 || (mol.getAtomicNo(atom) == 16 && !mol.isAromaticAtom(atom)))
					count++;
			value = ""+count;
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("symmetric_atoms") || propertyName.toLowerCase().equals("symm. atoms")) {
			this.shortDisplayedPropertyName = "Symm. atoms";
			this.longDisplayedPropertyName = "Symmetric atoms";
			mol.ensureHelperArrays(StereoMolecule.cHelperSymmetrySimple);
			int maxRank = 0;
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (maxRank < mol.getSymmetryRank(atom))
					maxRank = mol.getSymmetryRank(atom);
			value = ""+(mol.getAtoms()-maxRank);
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("all_amides")) {
			this.shortDisplayedPropertyName = "All amides";
			this.longDisplayedPropertyName = "All amides";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperNeighbours);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (AtomFunctionAnalyzer.isAmide(mol, atom))
					count++;
			value = ""+count;
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("all_amines")) {
			this.shortDisplayedPropertyName = "All amines";
			this.longDisplayedPropertyName = "All amines";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperRings);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (AtomFunctionAnalyzer.isAmine(mol, atom))
					count++;
			value = ""+count;
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("alkyl_amines")) {
			this.shortDisplayedPropertyName = "Alkyl amines";
			this.longDisplayedPropertyName = "Alkyl amines";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperRings);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (AtomFunctionAnalyzer.isAlkylAmine(mol, atom))
					count++;
			value = ""+count;
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("aryl_amines")) {
			this.shortDisplayedPropertyName = "Aryl amines";
			this.longDisplayedPropertyName = "Aryl amines";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperRings);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (AtomFunctionAnalyzer.isArylAmine(mol, atom))
					count++;
			value = ""+count;
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("aromatic_nitrogens") || propertyName.toLowerCase().equals("arom. nitrogens")) {
			this.shortDisplayedPropertyName = "Arom. nitrogens";
			this.longDisplayedPropertyName = "Aromatic nitrogens";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperRings);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (mol.getAtomicNo(atom) == 7 && mol.isAromaticAtom(atom))
					count++;

			value = ""+count;
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("basic_nitrogens")) {
			this.shortDisplayedPropertyName = "Basic nitrogens";
			this.longDisplayedPropertyName = "Basic nitrogens";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperRings);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (AtomFunctionAnalyzer.isBasicNitrogen(mol, atom))
					count++;

			value = ""+count;
		}
		else if (propertyName.replace(" ",  "_").toLowerCase().equals("acidic_nitrogens")) {
			this.shortDisplayedPropertyName = "Acidic nitrogens";
			this.longDisplayedPropertyName = "Acidic nitrogens";
			int count = 0;
			mol.ensureHelperArrays(StereoMolecule.cHelperRings);
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (AtomFunctionAnalyzer.isAcidicOxygen(mol, atom))
					count++;

			value = ""+count;
		}
	}
	
	/**
	 * Returns the number of bonds of the shortest path between
	 * those two atoms with the largest topological distance.
	 * @param mol
	 * @return
	 */
	private double assessMolecularShape(StereoMolecule mol) {
		mol.ensureHelperArrays(StereoMolecule.cHelperRings);
		if (mol.getAtoms() == 0)
			return -1;
		if (mol.getBonds() == 0)
			return 0;
	
		int maxLength = 0;
		for (int atom=0; atom<mol.getAtoms(); atom++)
			if (mol.getConnAtoms(atom) == 1 || mol.isRingAtom(atom))
				maxLength = Math.max(maxLength, findHighestAtomDistance(mol, atom));
	
		return (double)(maxLength+1) / (double)mol.getAtoms();
		}
	
	/**
	 * Calculates the molecular complexity (warning: time consuming!)
	 * @param mol
	 * @return Molecular complexity
	 */
	private double assessMolecularComplexity(StereoMolecule mol) {
		final int MAX_BOND_COUNT = 7;
		int bondCount = Math.min(mol.getBonds()/2, MAX_BOND_COUNT);
	
		mol.ensureHelperArrays(StereoMolecule.cHelperRings);
	    StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());
	    TreeSet<String> fragmentSet = new TreeSet<String>();
	    int[] atomMap = new int[mol.getAllAtoms()];
	
	    boolean[][] bondsTouch = new boolean[mol.getBonds()][mol.getBonds()];
	    for (int atom=0; atom<mol.getAtoms(); atom++) {
	    	for (int i=1; i<mol.getConnAtoms(atom); i++) {
	        	for (int j=0; j<i; j++) {
	        		int bond1 = mol.getConnBond(atom, i);
	        		int bond2 = mol.getConnBond(atom, j);
	        		bondsTouch[bond1][bond2] = true;
	        		bondsTouch[bond2][bond1] = true;
	        		}
	    		}
	    	}
	
	    boolean[] bondIsMember = new boolean[mol.getBonds()];
	    int maxLevel = bondCount - 2;
	    int[] levelBond = new int[maxLevel+1];
		for (int rootBond=0; rootBond<mol.getBonds(); rootBond++) {
			bondIsMember[rootBond] = true;
			int level = 0;
			levelBond[0] = rootBond;
			while (true) {
				boolean levelBondFound = false;
				while (!levelBondFound && levelBond[level] < mol.getBonds()-1) {
					levelBond[level]++;
					if (!bondIsMember[levelBond[level]]) {
	    				for (int bond=rootBond; bond<mol.getBonds(); bond++) {
	    					if (bondIsMember[bond] && bondsTouch[bond][levelBond[level]]) {
	    						levelBondFound = true;
	    						break;
	    						}
	    					}
						}
					}
	
				if (levelBondFound) {
					bondIsMember[levelBond[level]] = true;
					if (level == maxLevel) {
					    mol.copyMoleculeByBonds(fragment, bondIsMember, true, atomMap);
					    fragmentSet.add(new Canonizer(fragment).getIDCode());
	    				bondIsMember[levelBond[level]] = false;
						}
					else {
						level++;
						levelBond[level] = rootBond;
						}
					}
				else {
					if (--level < 0)
						break;
					bondIsMember[levelBond[level]] = false;
					}
				}
			}
	
	    return Math.log(fragmentSet.size()) / bondCount;
	    }
	
	/**
	 * Calculates the topological distance to the topologically most remote atom.
	 * @param mol
	 * @param startAtom
	 * @return number of bonds from startAtom to remote atom
	 */
	private int findHighestAtomDistance(StereoMolecule mol, int startAtom) {
		int[] graphLevel = new int[mol.getAtoms()];
	    int[] graphAtom = new int[mol.getAtoms()];
	
	    graphAtom[0] = startAtom;
	    graphLevel[startAtom] = 1;
	
	    int current = 0;
	    int highest = 0;
	    while (current <= highest /* && graphLevel[current] <= maxLength */) {
	        int parent = graphAtom[current];
	        for (int i=0; i<mol.getConnAtoms(parent); i++) {
	            int candidate = mol.getConnAtom(parent, i);
	            if (graphLevel[candidate] == 0) {
	                graphAtom[++highest] = candidate;
	                graphLevel[candidate] = graphLevel[parent]+1;
	                }
	            }
	        current++;
	        }
	    return graphLevel[graphAtom[highest]] - 1;
		}
	
	public String getCalculatedValue(String propertyName, StereoMolecule mol) {
		CalculateCompound(propertyName, mol, -1);
		return value;
	}
	
	public String getCalculatedValue(String propertyName, StereoMolecule mol, int calculatedPropertyAttributeIndex) {
		CalculateCompound(propertyName, mol, calculatedPropertyAttributeIndex);
		return value;
	}
	
	public String getShortDisplayedPropertyName() {
		return shortDisplayedPropertyName;
	}

	public String getLongDisplayedPropertyName() {
		return longDisplayedPropertyName;
	}
	
	public int getCalculatedPropertyAttributeIndex() {
		return calculatedPropertyAttributeIndex;
	}
	
}


