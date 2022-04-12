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
 * @author Joel Wahl
 */

package com.actelion.research.chem.conf;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.pharmacophoretree.HungarianAlgorithm;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;

import java.util.*;

/**
 * employing RMSD calculation using a HungarianAlgorithm for the atom assignment to correct
 * for overestimated RMSD values obtained by conventional methods for molecules with symmetric substructures 
 * based on: dx.doi.org/10.1021/ci400534h
 *
 */

public class SymmetryCorrectedRMSDCalculator {
	private Conformer conf1;
	private Conformer conf2;

	
	public SymmetryCorrectedRMSDCalculator(Conformer conf1, Conformer conf2) {
		this.conf1 = conf1;
		this.conf2 = conf2;
	}
	
	public double calculate() {
		int[][] assignments = getAssignments();
		double rmsd = 0.0;
		for(int[] pair : assignments) { 
			Coordinates c1 = conf1.getCoordinates(pair[0]);
			Coordinates c2 = conf2.getCoordinates(pair[1]);
			rmsd+=c1.distanceSquared(c2);
		}
		rmsd/=assignments.length;
		rmsd = Math.sqrt(rmsd);
		
		return rmsd;
	}
	
	/**
	 * using a HungarianAlgorithm for the optimal assignment of atoms for the RMSD calculation
	 * Restriction: Only atoms having the same atom type can be assigned to each other
	 * Cost matrix: squared distance between the atoms
	 * @return
	 */
	
	private int[][] getAssignments() {
		List<int[]> assignments = new ArrayList<>();
		StereoMolecule mol1 = new StereoMolecule(conf1.getMolecule());
		StereoMolecule mol2 = new StereoMolecule(conf2.getMolecule());

		mol1.ensureHelperArrays(Molecule.cHelperParities);
		mol2.ensureHelperArrays(Molecule.cHelperParities);
		
		Integer[] atomTypes1 = getAtomTypes(mol1);
		Integer[] atomTypes2 = getAtomTypes(mol2);
	
		Set<Integer> uniqueAtomTypes1 = new HashSet<>();
		Collections.addAll(uniqueAtomTypes1, atomTypes1);
		Set<Integer> uniqueAtomTypes2 = new HashSet<>();
		Collections.addAll(uniqueAtomTypes2, atomTypes2);
		assert(uniqueAtomTypes1.equals(uniqueAtomTypes2));
		for(int at : uniqueAtomTypes1) { 
			List<Integer> occurences1 = new ArrayList<>();
			List<Integer> occurences2 = new ArrayList<>();
			for(int i=0;i<atomTypes1.length;i++) {
				int at1 = atomTypes1[i];
				if(at1==at)
					occurences1.add(i);
			}
			for(int i=0;i<atomTypes2.length;i++) {
				int at2 = atomTypes2[i];
				if(at2==at)
					occurences2.add(i);
			}
			assert(occurences1.size()==occurences2.size());
			double[][] costMatrix = new double[occurences1.size()][occurences2.size()];
			int counter1 = 0;
			for(int occ1 : occurences1) { 
				int counter2 = 0;
				for(int occ2 : occurences2) {
					double distSq = conf1.getCoordinates(occ1).distanceSquared(conf2.getCoordinates(occ2));
					costMatrix[counter1][counter2] = distSq;
					counter2++;
				}
				counter1++;
			}
			int[][] ass = HungarianAlgorithm.hgAlgorithm(costMatrix, "min");
			for(int[] pair : ass) {
				int p1 = occurences1.get(pair[0]);
				int p2 = occurences2.get(pair[1]);
				assignments.add(new int[] {p1,p2});
			}
 		}
		return assignments.toArray(new int[assignments.size()][2]);
		
	}
	
	private Integer[] getAtomTypes(StereoMolecule mol) {
		Integer[] atomTypes = new Integer[mol.getAtoms()];
		for(int a=0;a<mol.getAtoms();a++) {
			atomTypes[a] = InteractionAtomTypeCalculator.getAtomType(mol, a);
		}
		return atomTypes;
	}
	
	
	

}
