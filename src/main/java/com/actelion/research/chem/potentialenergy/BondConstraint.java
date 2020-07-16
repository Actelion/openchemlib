package com.actelion.research.chem.potentialenergy;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;

public class BondConstraint implements PotentialEnergyTerm {
	
	private static double FORCE_CONSTANT = 50.0; //kcal/mol
	private Conformer conf;
	private int[] bondAtoms;
	private double targetValue;

	
	public BondConstraint(Conformer conf,int[] bondAtoms, double targetDistance) {
		this.conf = conf;
		this.bondAtoms = bondAtoms;
		this.targetValue = targetDistance;
	}


	@Override
	public double getFGValue(double[] gradient) {
		int a1 = bondAtoms[0];
		int a2 = bondAtoms[1];
		
		Coordinates c1 = conf.getCoordinates(a1);
		Coordinates c2 = conf.getCoordinates(a2);
		
		Coordinates d = c1.subC(c2);
		double dist = d.dist();
		Coordinates dGrad;
		double prefactor = FORCE_CONSTANT*(dist-targetValue);
		dGrad = d.scaleC(prefactor).scaleC(1.0/Math.max(dist, 1.0e-8));
		gradient[3*a1]  += dGrad.x;
		gradient[3*a1+1]  += dGrad.y;
		gradient[3*a1+2]  += dGrad.z;
		gradient[3*a2] -= dGrad.x;
		gradient[3*a2+1] -= dGrad.y;
		gradient[3*a2+2] -= dGrad.z;
		double distTerm = (dist-targetValue);
		return 0.5*FORCE_CONSTANT*distTerm*distTerm;
		
		
		
	}

}
