package com.actelion.research.chem.phesaflex;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.conf.TorsionRelevanceHelper;


/**
 * @author JW
 * Provides functionality to perform random dihedral angle perturbations on the 3D conformation of the molecule. 
 * Central bonds are perturbed by smaller values (given by torsion relevance), whereas terminal bonds can be perturbed by
 * 60 degrees
 *
 */
public class MetropolisMonteCarloHelper {
	
	private StereoMolecule mol;
	private static double MAX_ANGLE = 60.0/180.0*Math.PI;
	private static double MIN_ANGLE = 5.0/180.0*Math.PI;
	private static double TEMPERATURE = 0.0043; //move that reduces Tanimoto by 0.01 has 10% chance of acceptance
	private long seed = 1234L;
	private double rmax;
	private double rmin;
	private float[] torsionRelevance;
	private int[] rotatableBonds;
	private Random random;
	private int previousBond;
	private double previousAngle;
	private double slope;
	private BondRotationHelper bondRotationHelper;
	
	public MetropolisMonteCarloHelper(StereoMolecule mol) {
		this.mol = mol;
		
	}
	/**
	 * 
	 * @return: boolean, indicating if MMC Helper could successfully be initialized
	 */
	public boolean init() {
		boolean success = true;
		bondRotationHelper = new BondRotationHelper(mol);
		random = new Random(seed);
		boolean[] isRotatableBond = new boolean[mol.getBonds()];
		TorsionDB.findRotatableBonds(mol,true, isRotatableBond);
		List<Integer> rotBonds = new ArrayList<Integer>();
		IntStream.range(0, isRotatableBond.length).forEach(e -> {
			if(isRotatableBond[e])
				rotBonds.add(e);
		});
		rotatableBonds = rotBonds.stream().mapToInt(i->i).toArray();
		torsionRelevance = TorsionRelevanceHelper.getRelevance(mol, isRotatableBond);
		rmin = Float.MAX_VALUE;
		rmax = 0.0f;
		for(float relevance : torsionRelevance) {
			if (relevance<rmin)
				rmin = relevance;
			if(relevance>rmax)
				rmax = relevance;
		}
		slope = (MIN_ANGLE-MAX_ANGLE)/(rmax-rmin);
		if(rotatableBonds.length==0)
			success = false;
		return success;
	}
	
	public void step() {
		int prefactor = random.nextInt(2)<1 ? -1 : 1;
		previousBond = rotatableBonds[random.nextInt(rotatableBonds.length)];
		previousAngle = MAX_ANGLE+(torsionRelevance[previousBond]-rmin)*slope;
		previousAngle*=prefactor;
		bondRotationHelper.rotateSmallerSide(previousBond, previousAngle);
	}
	
	public void undoStep() {
		bondRotationHelper.rotateSmallerSide(previousBond, -previousAngle);
	}
	
	public boolean accept(double oldScore, double newScore) {
		boolean accept = false;
		if(newScore>oldScore)
			accept = true;
		else {
			double delta = -(newScore-oldScore);
			double p = Math.exp(-delta/TEMPERATURE);
			double rnd = random.nextDouble();
			if(rnd<p)
				accept = true;
		}
		return accept;
	}
	

}
