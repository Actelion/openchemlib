package com.actelion.research.chem.docking.scoring.plp;

import java.util.HashMap;
import java.util.Map;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.interactionstatistics.SplineFunction;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;

/**
 * piecewise linear potential (PLP) as described in doi:10.1021/ci800298z
 * parameters from model 3 (M3)
 * @author wahljo1
 *
 */

public class PLPTerm implements PotentialEnergyTerm {
	
	public static final Map<String, Double> HBOND_TERM  = new HashMap<String, Double>();
	static {
		HBOND_TERM.put("A", 2.3);
		HBOND_TERM.put("B", 2.6);
		HBOND_TERM.put("C", 3.1);
		HBOND_TERM.put("D", 3.4);
		HBOND_TERM.put("E", -1.0);
		HBOND_TERM.put("F", 20.0);
	};
	
	public static final Map<String, Double> METAL_TERM  = new HashMap<String, Double>();
	static {
		METAL_TERM.put("A", 1.4);
		METAL_TERM.put("B", 2.2);
		METAL_TERM.put("C", 2.6);
		METAL_TERM.put("D", 2.8);
		METAL_TERM.put("E", -1.0);
		METAL_TERM.put("F", 20.0);
	};
	
	public static final Map<String, Double> BURIED_TERM  = new HashMap<String, Double>();
	static {
		BURIED_TERM.put("A", 3.4);
		BURIED_TERM.put("B", 3.6);
		BURIED_TERM.put("C", 4.5);
		BURIED_TERM.put("D", 5.5);
		BURIED_TERM.put("E", -0.1);
		BURIED_TERM.put("F", 20.0);
	};
	
	public static final Map<String, Double> NONPOLAR_TERM  = new HashMap<String, Double>();
	static {
		NONPOLAR_TERM.put("A", 3.4);
		NONPOLAR_TERM.put("B", 3.6);
		NONPOLAR_TERM.put("C", 4.5);
		NONPOLAR_TERM.put("D", 5.5);
		NONPOLAR_TERM.put("E", -0.4);
		NONPOLAR_TERM.put("F", 20.0);
	};
	
	private int recAtom;
	private int ligAtom;
	private double A;
	private double B;
	private double C;
	private double D;
	private double E;
	private double F;
	private double D_sq;
	private Conformer ligand;
	private Conformer receptor;

	
	private PLPTerm(Conformer receptor, Conformer ligand, int recAtom, int ligAtom, Map<String,Double> term) {
		A = term.get("A");
		B = term.get("B");
		C = term.get("C");
		D = term.get("D");
		E = term.get("E");
		F = term.get("F");
		D_sq = D*D;
		this.recAtom = recAtom;
		this.ligAtom = ligAtom;
		this.ligand = ligand;
		this.receptor = receptor;
		
	}
	
	public static PLPTerm create(Conformer receptor, Conformer ligand, int recAtom, int ligAtom, Map<String,Double> term) {
		return new PLPTerm(receptor, ligand, recAtom, ligAtom, term);
	}
	

	
	
	
	@Override
	public double getFGValue(double[] gradient) {
		final Coordinates ci = receptor.getCoordinates(recAtom);		
		final Coordinates ck = ligand.getCoordinates(ligAtom);				
		final Coordinates cr = ci.subC(ck);
		double r2 = cr.distSq();		
		Coordinates grad = new Coordinates();
		double energy = 0.0;
		if(r2>D_sq) {
			energy = 0; 
		} 
		else {
			double prefactor = 0.0;
			double r = Math.sqrt(r2);
			if(r<A) {
				prefactor = (-F/A)*(1.0/r);
				grad = cr.scaleC(prefactor);
				energy = (F*(A-r))/A;
			}
			else if(r<B) {
				prefactor = (E/(B-A))*(1.0/r);
				grad = cr.scaleC(prefactor);
				energy = (E*(r-A))/(B-A);
			}
			else if(r<C) {
				prefactor = 0.0;
				grad = cr.scaleC(prefactor);
				energy = E;
			}
			else if(r<=D) {
				prefactor = (-E/(D-C))*(1.0/r);
				grad = cr.scaleC(prefactor);
				energy = (E*(D-r))/(D-C);
			}
		}
		gradient[3*ligAtom]-= grad.x;
		gradient[3*ligAtom+1]-= grad.y;
		gradient[3*ligAtom+2]-= grad.z;
		return energy;
	}
}


