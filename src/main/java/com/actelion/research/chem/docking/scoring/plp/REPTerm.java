package com.actelion.research.chem.docking.scoring.plp;


import java.util.Map;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;

/**
 * repulsive term, see PLPTerm for reference
 *
 */

public class REPTerm implements PotentialEnergyTerm {
	
	public static final double A = 3.2;
	public static final double B = 5.0;
	public static final double B_sq = B*B;
	public static final double C = 0.05;
	public static final double D = 20;
	
	private int recAtom;
	private int ligAtom;

	private Conformer ligand;
	private Conformer receptor;
	
	
	private REPTerm(Conformer receptor, Conformer ligand, int recAtom, int ligAtom) {
		this.recAtom = recAtom;
		this.ligAtom = ligAtom;
		this.receptor = receptor;
		this.ligand = ligand;
		
	}
	
	public static REPTerm create(Conformer receptor, Conformer ligand, int recAtom, int ligAtom) {
		return new REPTerm(receptor,ligand,recAtom,ligAtom);
	}
		
	@Override
	public double getFGValue(double[] gradient) {
		final Coordinates ci = receptor.getCoordinates(recAtom);		
		final Coordinates ck = ligand.getCoordinates(ligAtom);				
		final Coordinates cr = ci.subC(ck);
		double r2 = cr.distSq();		
		Coordinates grad = new Coordinates();
		double energy = 0.0;
		if(r2>B_sq) {
			energy = 0; 
		} 
		else {
			double prefactor = 0.0;
			double r = Math.sqrt(r2);
			if(r<A) {
				prefactor = ((C-D)/A)*(1.0/r);
				grad = cr.scaleC(prefactor);
				energy = r*(C-D)/A + D;
			}
			else {
				prefactor = (-C/(B-A))*(1.0/r);
				grad = cr.scaleC(prefactor);
				energy = -C*(r-A)/(B-A) + C;
			}
		
		}
		gradient[3*ligAtom]-= grad.x;
		gradient[3*ligAtom+1]-= grad.y;
		gradient[3*ligAtom+2]-= grad.z;
	
		return energy;
	}
}


