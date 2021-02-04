package com.actelion.research.chem.docking;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.docking.scoring.idoscore.InteractionTerm;
import com.actelion.research.chem.interactionstatistics.InteractionDistanceStatistics;
import com.actelion.research.chem.interactionstatistics.SplineFunction;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;


public class ScoringTask {
	
	/**
	 * 
	 * @param receptor
	 * @param ligand
	 * @param receptorAtomTypes
	 * @param ligandAtomTypes
	 * @return
	 */
	public static double calcScore(StereoMolecule receptor, StereoMolecule ligand, int[] receptorAtomTypes, int[] ligandAtomTypes) {
		Set<Integer> receptorAtoms = new HashSet<Integer>(); 
		List<InteractionTerm> terms = new ArrayList<InteractionTerm>();
		MoleculeGrid molGrid = new MoleculeGrid(receptor);
		for(int l=0;l<ligand.getAtoms();l++) 
			receptorAtoms.addAll(molGrid.getNeighbours(ligand.getCoordinates(l), InteractionDistanceStatistics.CUTOFF_RADIUS));
		
		for(int p : receptorAtoms) {
			for(int l=0;l<ligand.getAtoms();l++) {
				Conformer recConf = new Conformer(receptor);
				Conformer ligConf = new Conformer(ligand);
				terms.add(InteractionTerm.create(recConf, ligConf, p,l, receptorAtomTypes, ligandAtomTypes));
			}
		}
		
		double score = 0.0;
		double[] gradient = new double[3*ligand.getAllAtoms()];
		for(InteractionTerm term : terms) {
			if(term!=null) {
			//System.out.println(score);
			score += term.getFGValue(gradient);
			}
		}
		
		return score;
		
		
	}
	/**
	 * calculate interaction of probe atom with receptor
	 * @param receptor
	 * @param probeAtomType
	 * @param c
	 * @param receptorAtomTypes
	 * @param ligandAtomTypes
	 * @param grid
	 * @return
	 */
	public static double calcScore(StereoMolecule receptor, int probeAtomType, Coordinates c, int[] receptorAtomTypes) {
		double score = 0.0;
		Set<Integer> receptorAtoms = new HashSet<Integer>(); 
		MoleculeGrid grid = new MoleculeGrid(receptor);
		receptorAtoms.addAll(grid.getNeighbours(c, InteractionDistanceStatistics.CUTOFF_RADIUS));
		for(int p : receptorAtoms) {
			SplineFunction f = InteractionDistanceStatistics.getInstance().getFunction(receptorAtomTypes[p], probeAtomType);
			double dist = c.distance(receptor.getCoordinates(p));
			score+=f.getFGValue(dist)[0];
		}
		
	
		
		return score;
		
		
	}
	

};