package com.actelion.research.chem.potentialenergy;

import com.actelion.research.chem.conf.Conformer;


public class PositionConstraint implements PotentialEnergyTerm {
	private double[] refPos;
	private double k;
	private double d;
	private int atom;
	private Conformer conf;
	
	
	public PositionConstraint (Conformer conf, int atom,  double k, double d) {
		refPos = new double[3];
		this.atom = atom;
		this.conf = conf;
		refPos[0] = conf.getX(atom);
		refPos[1] = conf.getY(atom);
		refPos[2] = conf.getZ(atom);
		
		this.k = k;
		this.d = d;
	}
	


	@Override
	public double getFGValue(double[] gradient) {
		double energy = 0.0;
		double[] pos = new double[] {conf.getX(atom), conf.getY(atom), conf.getZ(atom)};
		double dx = pos[0]-refPos[0];
		double dy = pos[1]-refPos[1];
		double dz = pos[2]-refPos[2];
		double dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
		double prefactor = 0.0;
		if(dist>d)
			prefactor = dist-d;
		else 
			prefactor = 0.0;
		energy+=0.5*k*prefactor*prefactor;
		gradient[atom] += prefactor*dx/Math.max(dist, 1E-08);
		gradient[atom+1] += prefactor*dy/Math.max(dist, 1E-08);
		gradient[atom+2] += prefactor*dz/Math.max(dist, 1E-08);
	
		return energy;
	}
		
	

}
