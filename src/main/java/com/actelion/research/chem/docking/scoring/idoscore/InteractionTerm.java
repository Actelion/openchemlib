/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * 3. Neither the name of the copyright holder nor the
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
 * @author Joel Freyss
 */

package com.actelion.research.chem.docking.scoring.idoscore;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.interactionstatistics.InteractionDistanceStatistics;
import com.actelion.research.chem.interactionstatistics.SplineFunction;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;


/**
 * ProteinLigandTerm is used to represent the energy between 2 atoms.
 * 
 * The function can be divider into 3 parts:
 * - Linear part for the very short distances <1A (to avoid extreme derivates)
 * - Statistics derived potential when the data is enough to have a dependable function  
 * - Lennard-Jones 8-4 VDW when stats are not available
 * 
 */
public class InteractionTerm implements PotentialEnergyTerm {

	//Taper to the null function close to cutoff distance
	private final static double CUTOFF = InteractionDistanceStatistics.CUTOFF_RADIUS - InteractionDistanceStatistics.BIN_SIZE;
	private final static double CUTOFF_SQ = CUTOFF*CUTOFF;
	private double energy;
	private double factor;
	private Conformer ligand;
	private Conformer receptor;
	private int atoms[];	
	//Statistics
	private final SplineFunction f; 
	
	private InteractionTerm(Conformer receptor, Conformer ligand, int[] atoms, SplineFunction f, double factor) {
		this.receptor = receptor;
		this.ligand = ligand;
		this.f = f;			
		this.factor = factor;
		this.atoms = atoms;
	}
	


	public static InteractionTerm create(Conformer receptor, Conformer ligand, int p, int l, int[] receptorAtomTypes, int[] ligandAtomTypes) {		
		SplineFunction f = InteractionDistanceStatistics.getInstance().getFunction(receptorAtomTypes[p], ligandAtomTypes[l]);
		if(f==null) {
			return null;
		}
		return new InteractionTerm(receptor, ligand, new int[]{p, l}, f, 1.0);			
	}
	
	
	@Override
	public final double getFGValue(final double[] gradient) {
		final Coordinates ci = receptor.getCoordinates(atoms[0]);		
		final Coordinates ck = ligand.getCoordinates(atoms[1]);				
		final Coordinates cr = ci.subC(ck);
		double rik2 = cr.distSq();		
		double rik = 0.0;
		if(rik2>CUTOFF_SQ) {
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

				
				gradient[3*atoms[1]]-= cr.x;
				gradient[3*atoms[1]+1]-= cr.y;
				gradient[3*atoms[1]+2]-= cr.z;
			

				
			}	
			

			
		}


		
		return energy;

	}
	

	

}