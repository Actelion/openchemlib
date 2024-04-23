package com.actelion.research.chem.docking.scoring.chemscore;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;

public class SimpleMetalTerm implements PotentialEnergyTerm {
	
	private static final double D1 = 2.6;
	private static final double D2 = 3.0;
	private static final double CUTOFF_SQ = 9.0;

	
	private static final double PHI0 = Math.PI;
	private static final double PHI1 = 80*Math.PI/180.0;
	private static final double PHI2 = 90*Math.PI/180.0;

	private static final double ENERGY = -6.0;
	
	private Conformer receptor;
	private Conformer ligand;
	private int acceptor;
	private int metal;
	private int[] acceptorNeighbours;
	private double scale;

	
	
	private SimpleMetalTerm(Conformer receptor, Conformer ligand, int acceptor, int metal, 
			 int[] acceptorNeighbours, double scale) {
		this.receptor = receptor;
		this.ligand = ligand;
		this.acceptor = acceptor;
		this.metal = metal;
		this.acceptor = acceptor;
		this.acceptorNeighbours = acceptorNeighbours;
		this.scale = scale;
	}
	
	public static SimpleMetalTerm create(Conformer receptor, Conformer ligand, int acceptor, int metal, 
			 int[] acceptorNeighbours, double scale) {
		return new SimpleMetalTerm(receptor, ligand, acceptor,  metal, 
			 acceptorNeighbours, scale);
	}

	private double getDistTerm(double[] gradient) {
		Coordinates a;
		Coordinates m;
		Coordinates grad = new Coordinates();
		double energy = 0.0;
		a = ligand.getCoordinates(acceptor);
		m = receptor.getCoordinates(metal);
		Coordinates r = a.subC(m);
		double rSq = r.distSq();
		if(rSq<CUTOFF_SQ) {
	
			double d = Math.sqrt(rSq);

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
		}
		
		return energy;
		
	}
	
	private double getAngleTerm(double[] gradient, int a1, int a2, int a3, double x0, double x1, double x2) {
		
		double energy = 0.0;
		Coordinates c1, c2, c3;
		c1 = ligand.getCoordinates(a1);

		c2 = ligand.getCoordinates(a2);

		c3 = receptor.getCoordinates(a3);

		
	    Coordinates r0 = c1.subC(c2).unit();
	    Coordinates r1 = c3.subC(c2).unit();

	    double dist0 = c2.distance(c1);
	    double dist1 = c3.distance(c2);

	    double cosTheta = r0.cosAngle(r1);

	    double angleTerm = Math.acos(cosTheta) - x0;
	    boolean invert = false;
	    if(angleTerm<0) {
	    	angleTerm=-angleTerm;
	    	invert=true;
	    }
	    
	    if(angleTerm<x1)
	    	energy = 1.0;
	    else if(angleTerm>x2)
	    	energy = 0.0;
	    else {
	    	double prefactor = -1.0/(x2-x1); //derivative of energy term with respect to 
		    if(invert) {
		    	prefactor = -prefactor;
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
		if(energy!=0.0) {
			energies.add(energy);
			gradients.add(grad);

			for(int aa : acceptorNeighbours) {
				grad = new double[gradient.length];
				energy=getAngleTerm(gradient,aa,acceptor, metal,PHI0,PHI1,PHI2);
				energies.add(energy);
				gradients.add(grad);
			}
		}
		else 
			energies.add(energy);
		double[] totGrad = new double[gradient.length];
		double totEnergy = scale*ENERGY;
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
				for(int k=0;k<g.length;k++)
					g[k]*=w;
				
			}
			for(int l=0;l<totGrad.length;l++) {
				totGrad[l]+=g[l];
			}
		}
		for(int i=0;i<totGrad.length;i++) {
			gradient[i]+=totGrad[i];
		}
		return totEnergy;
	}

}
