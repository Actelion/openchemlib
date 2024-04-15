package com.actelion.research.chem.docking.scoring.chemscore;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.forcefield.mmff.Constants;
import com.actelion.research.chem.forcefield.mmff.Vector3;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;

public class HBTerm implements PotentialEnergyTerm {
	
	private static final double D0 = 1.85;
	private static final double D1 = 0.25;
	private static final double D2 = 0.65;
	private static final double CUTOFF_SQ = 6.25;

	
	private static final double PHI0 = Math.PI;
	private static final double PHI1 = 30*Math.PI/180.0;
	private static final double PHI2 = 80*Math.PI/180.0;

	
	private static final double PSI0 = Math.PI;
	private static final double PSI1 = 80*Math.PI/180.0;
	private static final double PSI2 = 100*Math.PI/180.0;
	
	private static final double ENERGY = -3.0;
	
	private Conformer receptor;
	private Conformer ligand;
	private int acceptor;
	private int donor;
	private int hydrogen;
	private int[] acceptorNeighbours;
	private boolean isLigAtomDonor;
	private boolean isLigAtomAcceptor;
	private double scale;
	
	
	private HBTerm(Conformer receptor, Conformer ligand, int acceptor, int donor, int hydrogen, 
			boolean isLigAtomAcceptor, boolean isLigAtomDonor, int[] acceptorNeighbours, double scale) {
		this.receptor = receptor;
		this.ligand = ligand;
		this.acceptor = acceptor;
		this.donor = donor;
		this.hydrogen = hydrogen;
		this.acceptorNeighbours = acceptorNeighbours;
		this.isLigAtomAcceptor = isLigAtomAcceptor;
		this.isLigAtomDonor = isLigAtomDonor;
		this.scale = scale;
		assert(isLigAtomDonor!=isLigAtomAcceptor); //donor and acceptor have to be from different molecules
	}
	
	
	public static HBTerm create(Conformer receptor, Conformer ligand, int acceptor, int donor, int hydrogen, 
			boolean isLigAtomAcceptor, boolean isLigAtomDonor, int[] acceptorNeighbours, double scale) {
		return new HBTerm(receptor, ligand, acceptor, donor, hydrogen, 
			isLigAtomAcceptor, isLigAtomDonor, acceptorNeighbours, scale);
	}
	
	
	/** d0 = 1.85;
	 *  d1 = 0.25;
	 *  d2 = 0.65
	 *  x = r-d0     r is the hydrogen-acceptor distance
	 *  if(|x|<d1 f = 1;
	 *  if(|x|>d2 f = 0;
	 *  else f(x) = (x2-x)/(x2-x1)
	 * @param gradient
	 * @return
	 */
	private double getDistTerm(double[] gradient) {
		Coordinates a;
		Coordinates h;
		Coordinates grad = new Coordinates();
		double energy = 0.0;
		if(isLigAtomAcceptor)
			a = ligand.getCoordinates(acceptor);
		else 
			a = receptor.getCoordinates(acceptor);
		if(isLigAtomDonor)
			h = ligand.getCoordinates(hydrogen);
		else 
			h = receptor.getCoordinates(hydrogen);
		Coordinates r = a.subC(h);
		double distSq = r.distSq();
		if(distSq<CUTOFF_SQ) {
			double d = Math.sqrt(distSq);
			boolean invert = false;
			double distTerm = d - D0;
			if(distTerm<0) {
				invert=true;
				distTerm = -distTerm;
			}
			if(distTerm<D1) {
				energy = 1.0;
			}
			else if(distTerm>D2) {
				energy = 0.0;
			}
			
			else {
				double prefactor = -1.0/(D2-D1)*(1.0/d);
				if(invert)
					prefactor*=-1;
				grad = r.scaleC(prefactor);
				grad.scale(scale);
				if(distTerm<0)
					grad.scaleC(-1.0);
				if(isLigAtomAcceptor) {//acceptor is ligand atom
					gradient[3*acceptor]+= grad.x;
					gradient[3*acceptor+1]+= grad.y;
					gradient[3*acceptor+2]+= grad.z;
				}
				else { // hydrogen is ligand atom	
					gradient[3*hydrogen]-= grad.x;
					gradient[3*hydrogen+1]-= grad.y;
					gradient[3*hydrogen+2]-= grad.z;
				}
				energy = (D2-distTerm)/(D2-D1);
			}
		}	
		return energy;
		
	}
	
	private double getAngleTerm(double[] gradient, int a1, int a2, int a3, boolean isLigAtomA1, boolean
			isLigAtomA2, boolean isLigAtomA3, double x0, double x1, double x2) {
		
		double energy = 0.0;
		Coordinates c1, c2, c3;
		if(isLigAtomA1)
			c1 = ligand.getCoordinates(a1);
		else
			c1 = receptor.getCoordinates(a1);
		if(isLigAtomA2)
			c2 = ligand.getCoordinates(a2);
		else
			c2 = receptor.getCoordinates(a2);
		if(isLigAtomA3)
			c3 = ligand.getCoordinates(a3);
		else
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
		   if(isLigAtomA1) {
		        gradient[3*a1    ] += prefactor*dCos_dS[0]/(-sinTheta);
		        gradient[3*a1 + 1] += prefactor*dCos_dS[1]/(-sinTheta);
		        gradient[3*a1 + 2] += prefactor*dCos_dS[2]/(-sinTheta);
		   }
		   if(isLigAtomA2) {
		        gradient[3*a2    ] += prefactor*(-dCos_dS[0] - dCos_dS[3])/(-sinTheta);
		        gradient[3*a2 + 1] += prefactor*(-dCos_dS[1] - dCos_dS[4])/(-sinTheta);
		        gradient[3*a2 + 2] += prefactor*(-dCos_dS[2] - dCos_dS[5])/(-sinTheta);
		   }
		   if(isLigAtomA3) {
		        gradient[3*a3    ] += prefactor*dCos_dS[3]/(-sinTheta);
		        gradient[3*a3 + 1] += prefactor*dCos_dS[4]/(-sinTheta);
		        gradient[3*a3 + 2] += prefactor*dCos_dS[5]/(-sinTheta);
		   }
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
			// first D-H-A angle
			grad = new double[gradient.length];
			energy= getAngleTerm(grad,donor,hydrogen,acceptor,
					isLigAtomDonor,isLigAtomDonor,isLigAtomAcceptor,
					PHI0,PHI1,PHI2);
			energies.add(energy);
			gradients.add(grad);
			for(int aa : acceptorNeighbours) {
				grad = new double[gradient.length];
				energy=getAngleTerm(grad,aa,acceptor,hydrogen,isLigAtomAcceptor,isLigAtomAcceptor,
						isLigAtomDonor,PSI0,PSI1,PSI2);
				energies.add(energy);
				gradients.add(grad);
			}
		}
		else 
			energies.add(0.0);
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
