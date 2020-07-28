package com.actelion.research.chem.docking.scoring.idoscore;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.TorsionDB;

/**
 * knowledge-based term for the internal strain of a ligand, composed of the pair-potentials
 * (only atoms further away than 3 bonds from each other)
 * 
 * @author wahljo1
 *
 */

public class StrainTerm {
	
	private int[][] torsions;
	private List<int[]> atomPairs;
	private StereoMolecule mol;
	
	public StrainTerm(StereoMolecule mol) {
		this.mol = mol;
		atomPairs = new ArrayList<>();
		init();
	}
	
	public void init() {
		boolean[] isRotatableBond = new boolean[mol.getBonds()];
		TorsionDB.findRotatableBonds(mol,true, isRotatableBond);
		List<Integer> rotBonds = new ArrayList<Integer>();
		IntStream.range(0, isRotatableBond.length).forEach(e -> {
			if(isRotatableBond[e])
				rotBonds.add(e);
		});
	}
	
	
	

}
