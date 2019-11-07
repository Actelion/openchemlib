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
package com.actelion.research.chem.interactionstatistics;

import java.text.DecimalFormat;
import java.util.Arrays;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;


/**
 * ProteinLigandTerm is used to represent the energy between 2 atoms.
 * 
 * The function can be divider into 3 parts:
 * - Linear part for the very short distances <1A (to avoid extreme derivates)
 * - Statistics derived potential when the data is enough to have a dependable function  
 * - Lennard-Jones 8-4 VDW when stats are not available
 * 
 */
public class InteractionTerm  {

	//Taper to the null function close to cutoff distance
	private final static double CUTOFF = InteractionDistanceStatistics.CUTOFF_RADIUS - InteractionDistanceStatistics.BIN_SIZE;
	
	public double rik2;
	private double energy;
	private double factor;
	private StereoMolecule ligand;
	private StereoMolecule receptor;
	private int atoms[];	
	//Statistics
	private final DistanceDependentPairPotential f; 
	
	private InteractionTerm(StereoMolecule receptor, StereoMolecule ligand, int[] atoms, DistanceDependentPairPotential f, double factor) {
		this.receptor = receptor;
		this.ligand = ligand;
		this.f = f;			
		this.factor = factor;
		this.atoms = atoms;
	}
	
	public static InteractionTerm create(StereoMolecule receptor, StereoMolecule ligand, int p, int l, int[] receptorAtomTypes, int[] ligandAtomTypes) {		
		DistanceDependentPairPotential f = InteractionDistanceStatistics.getInstance().getFunction(receptorAtomTypes[p], ligandAtomTypes[l]);
		if(f==null) {
			return null;
		}
		return new InteractionTerm(receptor, ligand, new int[]{p, l}, f, 1.0);			
	}
	
	

	public final double getFGValue(final double[] gradient) {
		final Coordinates ci = receptor.getCoordinates(atoms[0]);		
		final Coordinates ck = ligand.getCoordinates(atoms[1]);				
		final Coordinates cr = ci.subC(ck);
		rik2 = cr.distSq();		
		double rik = 0.0;
		if(rik2>CUTOFF*CUTOFF) {
			energy = 0; 
		} else {
			double de=0;
			rik = Math.sqrt(rik2);

			double grad[] = f.getFGValue(rik);			
			energy = factor * grad[0];	 
			if(gradient!=null) de = factor * grad[1];				


			if(gradient!=null) {
			
				double deddt = (rik<=1? -10 : de) / rik;
				cr.scale(deddt);
				if(atoms[1]<gradient.length) {
					gradient[3*atoms[1]]-= cr.x;
					gradient[3*atoms[1]+1]-= cr.y;
					gradient[3*atoms[1]+2]-= cr.z;
				}
			}	
			

			
		}


		
		return energy;

	}
	

	

}
