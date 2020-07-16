package com.actelion.research.chem.potentialenergy;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.torsionstrain.StatisticalTorsionTerm;

public class TorsionConstraint implements PotentialEnergyTerm {
	
	private static double FORCE_CONSTANT = 50.0; //kcal/mol
	private Conformer conf;
	private int[] torsionAtoms;
	private double targetValueMin;
	private double targetValueMax;
	
	/**
	 * Term that forces a torsion angle with a harmonic potential to a range of acceptable values (range given by the width) 
	 *
	 * @param conf
	 * @param torsionAtoms 
	 * @param targetValue : in degrees (0-360)
	 * @param width in degrees(>0)
	 */
	
	public TorsionConstraint(Conformer conf, int[] torsionAtoms, double targetValue, double width ) {
		this.conf = conf;
		this.torsionAtoms = torsionAtoms;
		this.targetValueMin = 2*Math.PI*(targetValue-width)/360.0;
		this.targetValueMax = 2*Math.PI*(targetValue+width)/360.0;
		if(targetValueMin < 0.0)
			targetValueMin+=2*Math.PI;
		if(targetValueMax < 0.0)
			targetValueMax=2*Math.PI;
	}
	
	double computeDihedralTerm(double dihedral)  {
		  double dihedralTarget = dihedral;
		  if (!(dihedral > targetValueMin && dihedral < targetValueMax)
		      && !(dihedral > targetValueMin && targetValueMin >  targetValueMax)
		      && !(dihedral < targetValueMax && targetValueMin > targetValueMax)) {
		    double minDihedralTarget = dihedral - targetValueMin;
		    double maxDihedralTarget = dihedral - targetValueMax;

		    dihedralTarget = (Math.abs(minDihedralTarget) < Math.abs(maxDihedralTarget)
		      ? targetValueMin : targetValueMax);
		  }
		  double dihedralTerm = dihedral - dihedralTarget;
		  return dihedralTerm;
		}
	
	
	
	
	@Override
	public double getFGValue(double[] gradient) {
		int a1 = torsionAtoms[0];
		int a2 = torsionAtoms[1];
		int a3 = torsionAtoms[2];
		int a4 = torsionAtoms[3];
		Coordinates c1 = conf.getCoordinates(a1);
		Coordinates c2 = conf.getCoordinates(a2);
		Coordinates c3 = conf.getCoordinates(a3);
		Coordinates c4 = conf.getCoordinates(a4);
		double dihedral = Coordinates.getDihedral(c1, c2, c3, c4);
		if(dihedral<0.0)
			dihedral+=2*Math.PI;
		double dihedralTerm = computeDihedralTerm(dihedral);
        double e = 0.5*dihedralTerm*dihedralTerm*FORCE_CONSTANT;
        double dEdPhi = FORCE_CONSTANT*dihedralTerm;
        
        StatisticalTorsionTerm.getCartesianTorsionGradient(torsionAtoms, conf,gradient, dEdPhi,new Coordinates[] {c1,c2,c3,c4},null);
        
        return e;
		
	}
	
	
	
	

}
