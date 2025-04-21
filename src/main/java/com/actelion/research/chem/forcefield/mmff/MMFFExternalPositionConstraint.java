package com.actelion.research.chem.forcefield.mmff;

import com.actelion.research.util.DoubleFormat;

public class MMFFExternalPositionConstraint implements EnergyTerm {
	private double[] refPos;
	private int constrainedAtom;
	private double k;
	private double d;
	
	
	public MMFFExternalPositionConstraint (int atom, double[] refPos, double k, double d) {
		this.refPos = refPos;
		this.constrainedAtom = atom;
		this.k = k;
		this.d = d;
	}

	@Override
	public double getEnergy(double[] pos) {
		return getEnergy(pos, null);
	}

	@Override
	public double getEnergy(double[] pos, StringBuilder detail) {
		double dx = pos[3*constrainedAtom]-refPos[0];
		double dy = pos[3*constrainedAtom+1]-refPos[1];
		double dz = pos[3*constrainedAtom+2]-refPos[2];
		double dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
		double prefactor = 0.0;
		if(dist>d)
			prefactor = dist-d;
		else 
			prefactor = 0.0;

		double energy = 0.5*k*prefactor*prefactor;

		if (detail != null)
			detail.append("extPosConstraint\t"+DoubleFormat.toString(dist)+"\t"+DoubleFormat.toString(d)+"\tall\t"+ DoubleFormat.toString(energy)+"\n");

		return energy;
	}

	@Override
	public void getGradient(double[] pos, double[] grad) {
			int a = constrainedAtom;
			double dx = pos[3*a]-refPos[0];
			double dy = pos[3*a+1]-refPos[1];
			double dz = pos[3*a+2]-refPos[2];
			double dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
			double prefactor = 0.0;
			if(dist>d)
				prefactor = dist-d;
			else 
				prefactor = 0.0;
			grad[a] += prefactor*dx/Math.max(dist, 1E-08);
			grad[a+1] += prefactor*dy/Math.max(dist, 1E-08);
			grad[a+2] += prefactor*dz/Math.max(dist, 1E-08);
	}
		
	

}
