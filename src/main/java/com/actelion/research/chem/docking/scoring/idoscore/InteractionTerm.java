/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Joel Freyss
 */
package com.actelion.research.chem.docking.scoring.idoscore;

import java.text.DecimalFormat;
import java.util.Arrays;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
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