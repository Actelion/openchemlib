package com.actelion.research.chem.conf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.pharmacophoretree.HungarianAlgorithm;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;

/**
 * employing RMSD calculation using a HungarianAlgorithm for the atom assignment to correct
 * for overestimated RMSD values obtained by conventional methods for molecules with symmetric substructures 
 * based on: dx.doi.org/10.1021/ci400534h
 * @author joel
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
					costMatrix[counter2][counter1] = distSq;
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
