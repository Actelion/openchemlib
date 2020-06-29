package com.actelion.research.chem.potentialenergy;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;


public class AngleConstraint implements PotentialEnergyTerm {
	private static double FORCE_CONSTANT = 50.0; //kcal/mol
	private Conformer conf;
	private int[] angleAtoms;
	private double targetValue;

	
	public AngleConstraint(Conformer conf,int[] angleAtoms, double targetAngle) {
		this.conf = conf;
		this.angleAtoms = angleAtoms;
		this.targetValue = 2*Math.PI*targetAngle/360.0;
	}
	

	@Override
	public double getFGValue(double[] gradient) {
		int a1 = angleAtoms[0];
		int a2 = angleAtoms[1]; //middle atom
		int a3 = angleAtoms[2];
		
		Coordinates c1 = conf.getCoordinates(a1);
		Coordinates c2 = conf.getCoordinates(a2);
		Coordinates c3 = conf.getCoordinates(a3);
		
        Coordinates r0 = c1.subC(c2).unit();
        Coordinates r1 = c3.subC(c2).unit();

        double dist0 = c2.distance(c1);
        double dist1 = c3.distance(c2);

        double cosTheta = r0.cosAngle(r1);

        double sinThetaSq = 1.0 - cosTheta*cosTheta;
        double sinTheta = 1.0e-8;
        if (sinThetaSq > 0.0)
            sinTheta = Math.sqrt(sinThetaSq);

        double angleTerm = Math.acos(cosTheta) - targetValue;

        double dEdTheta = FORCE_CONSTANT*angleTerm;


        double dCos_dS[] = new double[]{
            1.0/dist0*(r1.x - cosTheta*r0.x),
            1.0/dist0*(r1.y - cosTheta*r0.y),
            1.0/dist0*(r1.z - cosTheta*r0.z),
            1.0/dist1*(r0.x - cosTheta*r1.x),
            1.0/dist1*(r0.y - cosTheta*r1.y),
            1.0/dist1*(r0.z - cosTheta*r1.z)
        };
        gradient[3*a1    ] += dEdTheta*dCos_dS[0]/(-sinTheta);
        gradient[3*a1 + 1] += dEdTheta*dCos_dS[1]/(-sinTheta);
        gradient[3*a1 + 2] += dEdTheta*dCos_dS[2]/(-sinTheta);

        gradient[3*a2    ] += dEdTheta*(-dCos_dS[0] - dCos_dS[3])/(-sinTheta);
        gradient[3*a2 + 1] += dEdTheta*(-dCos_dS[1] - dCos_dS[4])/(-sinTheta);
        gradient[3*a2 + 2] += dEdTheta*(-dCos_dS[2] - dCos_dS[5])/(-sinTheta);

        gradient[3*a3    ] += dEdTheta*dCos_dS[3]/(-sinTheta);
        gradient[3*a3 + 1] += dEdTheta*dCos_dS[4]/(-sinTheta);
        gradient[3*a3 + 2] += dEdTheta*dCos_dS[5]/(-sinTheta);
        
        return 0.5*FORCE_CONSTANT*angleTerm*angleTerm;
}

}
