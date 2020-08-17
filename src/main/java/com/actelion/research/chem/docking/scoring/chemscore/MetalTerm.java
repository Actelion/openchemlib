package com.actelion.research.chem.docking.scoring.chemscore;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;

public class MetalTerm implements PotentialEnergyTerm {
	
	private static final double D1 = 0.6;
	private static final double D2 = 0.8;

	
	private static final double PHI0 = Math.PI;
	private static final double PHI1 = 80*Math.PI/180.0;
	private static final double PHI2 = 90*Math.PI/180.0;

	private static final double ENERGY = -6.0;
	
	private Conformer ligand;
	private int acceptor;

	private Coordinates fitPoint;
	private int[] acceptorNeighbours;
	private double scale;

	
	
	private MetalTerm(Conformer ligand, int acceptor,  
			 int[] acceptorNeighbours, Coordinates fitPoint, double scale) {
		this.ligand = ligand;
		this.acceptor = acceptor;
		this.acceptor = acceptor;
		this.acceptorNeighbours = acceptorNeighbours;
		this.fitPoint = fitPoint;
		this.scale = scale;
		
	}
	
	public static MetalTerm create(Conformer ligand, int acceptor, int[] acceptorNeighbours,
			Coordinates fitPoint, double scale) {
		return new MetalTerm(ligand, acceptor, acceptorNeighbours,
				fitPoint, scale);
	}

	private double getDistTerm(double[] gradient) {
		System.out.println("check metal");
		Coordinates a;
		Coordinates grad = new Coordinates();
		double energy = 0.0;
		a = ligand.getCoordinates(acceptor);
		Coordinates r = a.subC(fitPoint);
		System.out.println(a);
		System.out.println(fitPoint);
		double d = r.dist();
		System.out.println("dist");
		System.out.println(d);
		if(d<D1) {
			energy = 1.0;
		}
		else if(d>D2) {
			energy = 0.0;
		}
		
		else {
			double prefactor = -1.0/(D2-D1)*(1.0/d);
			grad = r.scaleC(prefactor);

			gradient[3*acceptor]+= grad.x;
			gradient[3*acceptor+1]+= grad.y;
			gradient[3*acceptor+2]+= grad.z;

			energy = (D2-d)/(D2-D1);
		}
		
		return energy;
		
	}
	
	private double getAngleTerm(double[] gradient, int a1, int a2, double x0, double x1, double x2) {
		
		double energy = 0.0;
		Coordinates c1, c2, c3;
		c1 = ligand.getCoordinates(a1);

		c2 = ligand.getCoordinates(a2);

		c3 = fitPoint;

		
	    Coordinates r0 = c1.subC(c2).unit();
	    Coordinates r1 = c3.subC(c2).unit();

	    double dist0 = c2.distance(c1);
	    double dist1 = c3.distance(c2);

	    double cosTheta = r0.cosAngle(r1);

	    double angleTerm = Math.acos(cosTheta) - x0;

	    
	    if(angleTerm<x1)
	    	energy = 1.0;
	    else if(angleTerm>x2)
	    	energy = 0.0;
	    else {
	    	double prefactor = -1.0/(x2-x1); //derivative of energy term with respect to 
		    if(angleTerm<0) {
		    	prefactor = -prefactor;
		    	angleTerm = -angleTerm;
		    }
		    energy = (x2-angleTerm)/(x2-x1);

		    double sinThetaSq = 1.0 - cosTheta*cosTheta;
		    double sinTheta = 1.0e-8;
		    if (sinThetaSq > 0.0)
		        sinTheta = Math.sqrt(sinThetaSq);


		   double dCos_dS[] = new double[]{
				    1.0/dist0*(r1.x - cosTheta*r0.x),
		            1.0/dist0*(r1.y - cosTheta*r0.y),
		            1.0/dist0*(r1.z - cosTheta*r0.z),
		            1.0/dist1*(r0.x - cosTheta*r1.x),
		            1.0/dist1*(r0.y - cosTheta*r1.y),
		            1.0/dist1*(r0.z - cosTheta*r1.z)
		        };

		   gradient[3*a1    ] += prefactor*dCos_dS[0]/(-sinTheta);
		   gradient[3*a1 + 1] += prefactor*dCos_dS[1]/(-sinTheta);
		   gradient[3*a1 + 2] += prefactor*dCos_dS[2]/(-sinTheta);
		   
		   gradient[3*a2    ] += prefactor*(-dCos_dS[0] - dCos_dS[3])/(-sinTheta);
		   gradient[3*a2 + 1] += prefactor*(-dCos_dS[1] - dCos_dS[4])/(-sinTheta);
		   gradient[3*a2 + 2] += prefactor*(-dCos_dS[2] - dCos_dS[5])/(-sinTheta);

	    }
		return energy;
		
	}
	

	@Override
	public double getFGValue(double[] gradient) {
		List<Double> energies = new ArrayList<>();
		List<double[]> gradients = new ArrayList<double[]>();
		double energy = 0.0;
		double[] grad = new double[gradient.length];
		energy = getDistTerm(grad);
		double totEnergy = 0.0;
		if(energy!=0.0) {
			energies.add(energy);
			gradients.add(grad);

			for(int aa : acceptorNeighbours) {
				grad = new double[gradient.length];
				energy=getAngleTerm(gradient,aa,acceptor,PHI0,PHI1,PHI2);
				energies.add(energy);
				gradients.add(grad);
			}
		
			double[] totGrad = new double[gradient.length];
			System.out.println("engies");
			System.out.println(energies);
			totEnergy = scale*ENERGY;
			for(double eng : energies)
				totEnergy*=eng;
			//apply product rules
			for(int i=0;i<gradients.size();i++) {
				double[] g = gradients.get(i);
				for(int j=0;j<gradients.size();j++) {
					if(i==j)
						continue;
					double e = energies.get(j);
					double w = e*scale*ENERGY;
					for(int k=0;j<g.length;k++)
						g[k]*=w;
					
				}
				for(int l=0;l<totGrad.length;l++) {
					totGrad[l]+=g[l];
				}
			}
			for(int i=0;i<totGrad.length;i++) {
				gradient[i]+=totGrad[i];
			}
		}
		System.out.println("tot");
		System.out.println(totEnergy);
		return totEnergy;
	}
	
	
	

}
