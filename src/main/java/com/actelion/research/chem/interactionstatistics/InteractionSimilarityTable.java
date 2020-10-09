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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 *
 * distance gives the distance of the point of highest potential.
 * equivalence makes the difference, considering the distance and the potential
 * distance is absolute, equivalence is relative to another distance
 *
 * @author freyssj, modified by JW
 */
public class InteractionSimilarityTable {
	private static final double LOWER_BOUNDARY = 1.6; 
	private static final double UPPER_BOUNDARY = 4.5;
	private static final double WELL_DEPTH_CUTOFF = 2.0; 
	private static InteractionSimilarityTable instance = null; 
	private final InteractionDistanceStatistics stats;
	private final Map<Integer,Integer> keyToId;
	private final InteractionDescriptor[][] iDsToDescriptor;
	private final double[][] similarityTable;
	Set<Integer> atomKeys;
	

	/**
	 * 
	 * 
	 */
	public static class InteractionDescriptor {
		public int N;
		public double optimalDist;
		public double optimalStrength;
		
		public InteractionDescriptor(SplineFunction plf) {
			if(plf==null) return;
			N = plf.getTotalOccurences();
			try {
				//Find the minima of this function
				optimalStrength = 100;
				optimalDist = -1;
				for (double d = LOWER_BOUNDARY; d < UPPER_BOUNDARY; d+=.05) {
					double v = plf.getFGValue(d)[0];
					if(optimalStrength<0 && v>optimalStrength+0.5) break; //consider only first min
					if(v<optimalStrength) {
						optimalDist = d;
						optimalStrength = v;
					}
				}							
			} catch (Exception e) {
				e.printStackTrace();
			}			 
		}
		
		public double dist(InteractionDescriptor d2) {
			
			double s1 = optimalDist>0?optimalStrength:0;
			double m1 = optimalDist>0?optimalDist:UPPER_BOUNDARY;
			double s2 = d2.optimalDist>0?d2.optimalStrength:0;
			double m2 = d2.optimalDist>0?d2.optimalDist:UPPER_BOUNDARY;
			
			double dist1 = Math.abs(m1-m2)/(UPPER_BOUNDARY-LOWER_BOUNDARY);
			
			double dist2 = Math.abs(s2-s1)>WELL_DEPTH_CUTOFF ? 1.0 : Math.abs(s2-s1)/WELL_DEPTH_CUTOFF;
			
			return 0.5*dist1+0.5*dist2;
			 
		}
				
		@Override
		public String toString() {return new DecimalFormat("0.0").format(optimalDist)+":"+new DecimalFormat("0.00").format(optimalStrength)+" ["+N+"]";}

		
	}
		
	/**
	 * Private Constructor (singleton pattern)
	 * Creates a table of similarity between the interaction atom types
	 *  
	 *
	 */
	private InteractionSimilarityTable() {
		keyToId = new HashMap<Integer,Integer>();
		stats = InteractionDistanceStatistics.getInstance();
		
		//Prepare the proteinLigandIDs table
		atomKeys = stats.getAtomKeySet();
		int N = atomKeys.size();
		int index = 0;
		
		for(int key: atomKeys) {
			keyToId.putIfAbsent(key, index);
			index++;	
		}
		
		iDsToDescriptor = new InteractionDescriptor[N][N];
		
		for (int i : atomKeys) {
			for (int j : atomKeys) {
				SplineFunction plf = stats.getFunction(i, j);
				iDsToDescriptor[keyToId.get(i)][keyToId.get(j)] = new InteractionDescriptor(plf);
			}			
		}

		
		//Create the similarityTable table

		similarityTable = new double[N][N];
		for (int l1 = 0; l1 < N; l1++) {
			for (int l2 = l1; l2 < N; l2++) {
				double sum = 0;
				int total = 0;
				for (int i = 0; i < iDsToDescriptor.length; i++) {
					InteractionDescriptor id1 = iDsToDescriptor[i][l1];
					InteractionDescriptor id2 = iDsToDescriptor[i][l2];
					double coeff = id1.N+id2.N;
					sum   += id1.dist(id2) * coeff;
					total += coeff; 
				}
				similarityTable[l1][l2] = similarityTable[l2][l1] = total>0? sum/total: 5;		
			}			
		}

	}
	
	
	public static InteractionSimilarityTable getInstance() {
		if(instance==null) {
			synchronized(InteractionSimilarityTable.class) {
				if(instance==null) {
					instance = new InteractionSimilarityTable();
				}
			}
		}

		return instance;
	}
	

	public double getDistance(int type1, int type2) {
		int a = keyToId.get(InteractionDistanceStatistics.getInstance().getKey(type1));
		int b = keyToId.get(InteractionDistanceStatistics.getInstance().getKey(type2));
		return similarityTable[a][b];
	}
	
	/**
	 * D(LigandType_1, LigandType_2) = Sum( d( F(ProteinType_i, LigandType_1), F(ProteinType_i, LigandType_2)), i) 
	 * @param type1
	 * @param type2
	 * @return
	 */
	public double getDissimilarity(int type1, int type2) {
		int a = keyToId.get(InteractionDistanceStatistics.getInstance().getKey(type1));
		int b = keyToId.get(InteractionDistanceStatistics.getInstance().getKey(type2));

		double sum = 0;
		for (int i = 0; i < similarityTable.length; i++) {
			double diff = Math.abs( similarityTable[a][i] - similarityTable[b][i]);
			sum+=diff;
		}
		return sum/ similarityTable.length;
	}
	/**
	 * Compare similarity values of 2 types (across all lines)
	 * @param type
	 * @param maxDist
	 * @return
	 */
	public List<Integer> getEquivalentTypes(int type, double maxDist) {
		type = InteractionDistanceStatistics.getInstance().getKey(type);
		List<Integer> res = new ArrayList<Integer>();
		for (int type2 : atomKeys) {			
			if(getDissimilarity(type, type2)<maxDist) {
				res.add(type2);
			}
		}
		return res;		
	}
	
	

	
	
	
}
