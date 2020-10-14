/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/
package com.actelion.research.chem.conf.torsionstrain;

import java.text.DecimalFormat;
import java.util.Arrays;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.forcefield.mmff.Vector3;
import com.actelion.research.chem.interactionstatistics.SplineFunction;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;


/**
 * Represents a torsion potential as a function of the angle, derived from statistical torsion distributions from the COD/CSD
 * The dihedral angle is defined by a atom sequence of 4 atoms. If a core atom is sp3 hybridized and has two neighbours of the same symmetry rank,
 * then the torsion angle has to be defined by a virtual atom, occupying the position of the hypothetical third atom
 * 
 */
public class StatisticalTorsionTerm implements PotentialEnergyTerm  {


	private static final double EPS = 0.00001;
	public double rik2;
	private Conformer conf;
	private int atoms[];
	private int rearAtoms[][]; //to define virtual torsion atoms
	//Statistics
	private final SplineFunction f; 
	
	private StatisticalTorsionTerm(Conformer conf, int[] atoms, SplineFunction f) {
		this.conf = conf;
		this.f = f;			
		this.atoms = atoms;
		assessRearAtoms();
		
	}
	
	private void assessRearAtoms() {
		StereoMolecule mol = conf.getMolecule();
		rearAtoms = new int[2][];
		if(atoms[0]==-1) {
			rearAtoms[0] = new int[2];
			int index = 0;
			for(int i=0;i<mol.getConnAtoms(atoms[1]);i++) {
				if(mol.getConnAtom(atoms[1], i)==atoms[2])
					continue;
				else {
					rearAtoms[0][index] = mol.getConnAtom(atoms[1], i);
					index++;
				}
			}
		}
		
		if(atoms[3]==-1) {
			rearAtoms[1] = new int[2];
			int index = 0;
			for(int i=0;i<mol.getConnAtoms(atoms[2]);i++) {
				if(mol.getConnAtom(atoms[2], i)==atoms[1])
					continue;
				else {
					rearAtoms[1][index] = mol.getConnAtom(atoms[2], i);
					index++;
				}
			}
		}
		
	}
	
	public static StatisticalTorsionTerm create(Conformer conf, int[] atoms, String torsionID) {		
		SplineFunction f = StatisticalTorsionPotential.getInstance().getFunction(torsionID);
		if(f==null) {
			return null;
		}
		return new StatisticalTorsionTerm(conf,atoms, f);			
	}
	
	public static void getCartesianTorsionGradient(int[] atoms, Conformer conf, double[] gradient, double dEdPhi, Coordinates[] coords, int[][] rearAtoms) {
	
		Coordinates c1 = coords[0];
		Coordinates c2 = coords[1];
		Coordinates c3 = coords[2];
		Coordinates c4 = coords[3];
		
		int a1 = atoms[0];
		int a2 = atoms[1];
		int a3 = atoms[2];
		int a4 = atoms[3];
		
       Coordinates[] r = new Coordinates[]{
               c1.subC(c2),
               c3.subC(c2),
               c2.subC(c3),
               c4.subC(c3),
           };
       Coordinates[] t = new Coordinates[]{
               r[0].cross(r[1]),
               r[2].cross(r[3])
           };

           double[] d = new double[]{
               t[0].dist(),
               t[1].dist()
           };

           if (Math.abs(d[0]) < EPS|| Math.abs(d[1]) < EPS)
               return;

           t[0].unit();
           t[1].unit();

           double cosPhi = t[0].dot(t[1]);
           double sinPhiSq = 1.0 - cosPhi * cosPhi;
           double sinPhi = ((sinPhiSq > 0.0) ? Math.sqrt(sinPhiSq) : 0.0);


           double sinTerm = -dEdPhi * (Math.abs(sinPhi) < 0.00001
                   ? (1.0 / cosPhi) : (1.0 / sinPhi));

           double[] dCos_dT = new double[]{
               1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
               1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
               1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
               1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
               1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
               1.0 / d[1] * (t[0].z - cosPhi * t[1].z)
           };
           if(rearAtoms!=null) {
	           if(rearAtoms[0]!=null) { // chain rule, take into account derivate of virtual atom wrt rear atoms
	        	   double derX = sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
	        	   double derY = sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
	        	   double derZ = sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);
	        	   gradient[3*rearAtoms[0][0]] -= derX;
	        	   gradient[3*rearAtoms[0][0]+1] -= derY;
	        	   gradient[3*rearAtoms[0][0]+2] -= derZ;
	        	   gradient[3*rearAtoms[0][1]] -= derX;
	        	   gradient[3*rearAtoms[0][1]+1] -= derY;
	        	   gradient[3*rearAtoms[0][1]+2] -= derZ;
	           }
           
	           else {
	        	   gradient[3*a1+0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
	        	   gradient[3*a1+1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
	        	   gradient[3*a1+2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);
	           }
           }

           gradient[3*a2+0] += sinTerm * (dCos_dT[1] * (r[1].z - r[0].z)
                   + dCos_dT[2] * (r[0].y - r[1].y)
                   + dCos_dT[4] * (-r[3].z)
                   + dCos_dT[5] * (r[3].y));
           gradient[3*a2+1] += sinTerm * (dCos_dT[0] * (r[0].z - r[1].z)
                   + dCos_dT[2] * (r[1].x - r[0].x)
                   + dCos_dT[3] * (r[3].z)
                   + dCos_dT[5] * (-r[3].x));
           gradient[3*a2+2] += sinTerm * (dCos_dT[0] * (r[1].y - r[0].y)
                   + dCos_dT[1] * (r[0].x - r[1].x)
                   + dCos_dT[3] * (-r[3].y)
                   + dCos_dT[4] * (r[3].x));

           gradient[3*a3+0] += sinTerm * (dCos_dT[1] * (r[0].z)
                   + dCos_dT[2] * (-r[0].y)
                   + dCos_dT[4] * (r[3].z - r[2].z)
                   + dCos_dT[5] * (r[2].y - r[3].y));
           gradient[3*a3+1] += sinTerm * (dCos_dT[0] * (-r[0].z)
                   + dCos_dT[2] * (r[0].x)
                   + dCos_dT[3] * (r[2].z - r[3].z)
                   + dCos_dT[5] * (r[3].x - r[2].x));
           gradient[3*a3+2] += sinTerm * (dCos_dT[0] * (r[0].y)
                   + dCos_dT[1] * (-r[0].x)
                   + dCos_dT[3] * (r[3].y - r[2].y)
                   + dCos_dT[4] * (r[2].x - r[3].x));
           
           if(rearAtoms!=null) {
	           if(rearAtoms[1]!=null) {
	        	   double derX = sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
	        	   double derY = sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
	        	   double derZ = sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
	        	   gradient[3*rearAtoms[1][0]] -= derX;
	        	   gradient[3*rearAtoms[1][0]+1] -= derY;
	        	   gradient[3*rearAtoms[1][0]+2] -= derZ;
	        	   gradient[3*rearAtoms[1][1]] -= derX;
	        	   gradient[3*rearAtoms[1][1]+1] -= derY;
	        	   gradient[3*rearAtoms[1][1]+2] -= derZ;
	           }   
	           else {
	        	   gradient[3*a4+0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
	        	   gradient[3*a4+1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
	        	   gradient[3*a4+2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
	           }
           }
	}


	@Override
	public final double getFGValue(final double[] gradient) {
		for(int i=0;i<gradient.length;i++) {
			gradient[i] = 0.0;
		}
		Coordinates c1,c2,c3,c4;
		int a1 = atoms[0];
		int a2 = atoms[1];
		int a3 = atoms[2];
		int a4 = atoms[3];
		if(a1==-1) {
			Coordinates c = conf.getCoordinates(a2);
			Coordinates c1a = conf.getCoordinates(rearAtoms[0][0]);
			Coordinates c2a = conf.getCoordinates(rearAtoms[0][1]);
			Coordinates v1 = c1a.subC(c);
			Coordinates v2 = c2a.subC(c);
			c1 = v1.addC(v2);
			c1.scale(-1.0);
		}
		else {
			c1 = conf.getCoordinates(a1);
		}
		
		if(a4==-1) {
			Coordinates c = conf.getCoordinates(a3);
			Coordinates c1a = conf.getCoordinates(rearAtoms[1][0]);
			Coordinates c2a = conf.getCoordinates(rearAtoms[1][1]);
			Coordinates v1 = c1a.subC(c);
			Coordinates v2 = c2a.subC(c);
			c4 = v1.addC(v2);
			c4.scale(-1.0);
		}
		else {
			c4 = conf.getCoordinates(a4);
		}
		
		c2 = conf.getCoordinates(atoms[1]);
		c3 = conf.getCoordinates(atoms[2]);

		
		double dihedral = Coordinates.getDihedral(c1, c2, c3, c4);
		if(dihedral<0.0)
			dihedral+=2*Math.PI;
		double[] res = f.getFGValue(360.0*dihedral/2*Math.PI);
		double e = res[0];
		double dEdPhi = res[1];	
		getCartesianTorsionGradient(atoms, conf,gradient, dEdPhi,new Coordinates[] {c1,c2,c3,c4}, rearAtoms);
		
		return e;
	}

	

	

}
