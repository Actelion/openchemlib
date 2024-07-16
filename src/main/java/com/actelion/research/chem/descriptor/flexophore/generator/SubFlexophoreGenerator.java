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
 */

package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.calc.combinatorics.CombinationGenerator;
import com.actelion.research.chem.descriptor.flexophore.*;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * SubFlexophoreGenerator
 * @author Modest von Korff
 * Oct 12, 2012 MvK: Start implementation
 */
public class SubFlexophoreGenerator {


	private int minBinDistThresh;
	private int maxDistanceBinThresh;
	private int minNumDifferentInteractionTypes;

	private ViolatedConditionsCount violatedConditionsCount;

	private HashSet<Integer> hsInteractionType;

	/**
	 */
	public SubFlexophoreGenerator(SubFlexophoreGenerator a) {
		this.minBinDistThresh = a.minBinDistThresh;
		this.maxDistanceBinThresh = a.maxDistanceBinThresh;
		this.minNumDifferentInteractionTypes = a.minNumDifferentInteractionTypes;

		violatedConditionsCount = new ViolatedConditionsCount();

		hsInteractionType = new HashSet<>();
	}

	public SubFlexophoreGenerator(int minBinDistThresh, int maxDistanceBinThresh, int minNumDifferentInteractionTypes) {
		this.minBinDistThresh = minBinDistThresh;
		this.maxDistanceBinThresh = maxDistanceBinThresh;
		this.minNumDifferentInteractionTypes = minNumDifferentInteractionTypes;

		violatedConditionsCount = new ViolatedConditionsCount();

		hsInteractionType = new HashSet<>();
	}

	public List<MolDistHist> generateSubPharmacophoresCheckedRange(MolDistHist mdh, int size){

		List<int[]> liIndSub = CombinationGenerator.getAllOutOf(mdh.getNumPPNodes(), size);
		if(liIndSub == null)
			return null;

		List<MolDistHist> liMolDistHistSub = new ArrayList<>();

		for (int[] arrIndSub : liIndSub) {
			MolDistHist mdhSub = getSubFragmentCheckedRange(mdh, arrIndSub);
			if(mdhSub==null) {
				continue;
			}

			//
			// check minimum requirements
			//
			int nNodes = mdhSub.getNumPPNodes();
			hsInteractionType.clear();
			for (int k = 0; k < nNodes; k++) {
				PPNode node = mdhSub.getNode(k);
				int nIATypes = node.getInteractionTypeCount();
				for (int l = 0; l < nIATypes; l++) {
					int iaType = node.getInteractionType(l);
					hsInteractionType.add(iaType);
				}
			}

			if(hsInteractionType.size()<minNumDifferentInteractionTypes){
				violatedConditionsCount.ccViolatedMinDiffInteractionTypes++;
				continue;
			}

			liMolDistHistSub.add(mdhSub);

		}

		return liMolDistHistSub;
	}

	public MolDistHist getSubFragmentCheckedRange(MolDistHist mdh, int [] arrIndices){

		MolDistHist frag = new MolDistHist(arrIndices.length);

		for (int i = 0; i < arrIndices.length; i++) {
			PPNode node = new PPNode(mdh.getNode(arrIndices[i]));
			frag.addNode(node);
		}

		boolean minDistReached = false;
		boolean maxDistViolated = false;

		violated:
		for (int i = 0; i < arrIndices.length; i++) {
			for (int j = i+1; j < arrIndices.length; j++) {
				byte [] arrHist = mdh.getDistHist(arrIndices[i],arrIndices[j]);

				if(!minDistReached) {
					for (int length = minBinDistThresh; length < arrHist.length; length++) {
						if (arrHist[length] > 0) {
							minDistReached = true;
							break;
						}
					}
				}

				for (int k = arrHist.length - 1; k > maxDistanceBinThresh; k--) {
					if (arrHist[k] > 0) {
						maxDistViolated = true;
						break violated;
					}
				}
				frag.setDistHist(i,j,arrHist);
			}
		}

		if(!minDistReached || maxDistViolated){

			if(!minDistReached){
				violatedConditionsCount.ccMissedMinRange++;
			}
			if(maxDistViolated){
				violatedConditionsCount.ccViolatedMaxRange++;
			}

			return null;
		}
		frag.realize();
		return frag;
	}


	public ViolatedConditionsCount getViolatedConditionsCount() {
		return violatedConditionsCount;
	}


	public static List<MolDistHistViz> generateSubPharmacophores(MolDistHistViz mdhv, int minNumPPPoints, int maxNumPPPoints){
		
		List<HashSet<MolDistHistViz>> liHsMolDistHistViz = new ArrayList<HashSet<MolDistHistViz>>();
		
		for (int size = 0; size < maxNumPPPoints+1; size++) {
			liHsMolDistHistViz.add(new HashSet<MolDistHistViz>());
		}
					
		int maxNumPPPointsReal = Math.min(mdhv.getNumPPNodes(), maxNumPPPoints) + 1;
		
		for (int size = minNumPPPoints; size < maxNumPPPointsReal; size++) {
			List<MolDistHistViz> liMDHVSub = generateSubPharmacophores(mdhv, size);
			liHsMolDistHistViz.get(size).addAll(liMDHVSub);
		}
					
		List<MolDistHistViz> liMDHVSubAll = new ArrayList<MolDistHistViz>();
		
		for (HashSet<MolDistHistViz> hsMDHVSub : liHsMolDistHistViz) {
			liMDHVSubAll.addAll(hsMDHVSub);
		}
		
		
		return liMDHVSubAll;
	}







	public static List<MolDistHistViz> generateSubPharmacophores(List<MolDistHistViz> liMDHV, int minNumPPPoints, int maxNumPPPoints){
		
		List<HashSet<MolDistHistViz>> liHsMolDistHistViz = new ArrayList<HashSet<MolDistHistViz>>();
		
		for (int size = 0; size < maxNumPPPoints+1; size++) {
			liHsMolDistHistViz.add(new HashSet<MolDistHistViz>());
		}
			
		for (MolDistHistViz mdhv : liMDHV) {
			
			int maxNumPPPointsReal = Math.min(mdhv.getNumPPNodes(), maxNumPPPoints) + 1;
			
			for (int size = minNumPPPoints; size < maxNumPPPointsReal; size++) {
				List<MolDistHistViz> liMDHVSub = generateSubPharmacophores(mdhv, size);
				
				liHsMolDistHistViz.get(size).addAll(liMDHVSub);
				
			}
			
		}
		
		List<MolDistHistViz> liMDHVSubAll = new ArrayList<MolDistHistViz>();
		
		for (HashSet<MolDistHistViz> hsMDHVSub : liHsMolDistHistViz) {
			liMDHVSubAll.addAll(hsMDHVSub);
		}
		
		
		return liMDHVSubAll;
	}

	
	/**
	 * Generates a list of sub-fragments.
	 * @param mdhv
	 * @param size number of pharmacophore points
	 * @return unique list of features
	 */
	public static List<MolDistHistViz> generateSubPharmacophores(MolDistHistViz mdhv, int size){
		
		List<int[]> liIndSub = CombinationGenerator.getAllOutOf(mdhv.getNumPPNodes(), size);
		if(liIndSub == null)
			return null;
		
		HashSet<MolDistHistViz> hsMolDistHistViz = new HashSet<MolDistHistViz>();
		
		List<MolDistHistVizFrag> liMDH = getSubFragments(mdhv, liIndSub);
		for (MolDistHistViz mdhFrag : liMDH) {
			hsMolDistHistViz.add(mdhFrag);
		}
		
		List<MolDistHistViz> liMDHFeatures = new ArrayList<MolDistHistViz>(hsMolDistHistViz);
		
		return liMDHFeatures;
	}
	
	/**
	 * Inevitable pharmacophore points are considered (01.11.2012)
	 * @param mdh
	 * @param liIndices indices for pharmacophore points in <code>mdh</code>.
	 * @return
	 */
	public static List<MolDistHistVizFrag> getSubFragments(MolDistHistViz mdh, List<int []> liIndices){
		List<MolDistHistVizFrag> liFrags = new ArrayList<MolDistHistVizFrag>(liIndices.size());
		for (int [] arr : liIndices) {

			boolean inevitableAllInSub = true;
			if(mdh.getNumMandatoryPharmacophorePoints()>0){
				
				int nInevitableInSub = 0;
				for (int i = 0; i < arr.length; i++) {
					
					if(mdh.isMandatoryPharmacophorePoint(arr[i])){
						nInevitableInSub++;	
					}
				}
				
				if(nInevitableInSub < mdh.getNumMandatoryPharmacophorePoints()){
					inevitableAllInSub = false;
				}
			}
			
			if(inevitableAllInSub) {
			
				MolDistHistVizFrag frag = getSubFragment(mdh, arr);
				
				liFrags.add(frag);
			}
			
		}
		
		return liFrags;
	}


	/**
	 * @param mdh
	 * @param arrIndexNodes the original index is set in the index field of the new SubFragment.
	 * @return
	 */
	public static MolDistHistVizFrag getSubFragment(MolDistHistViz mdh, int [] arrIndexNodes){
		
		
		MolDistHistVizFrag frag = new MolDistHistVizFrag(arrIndexNodes.length, mdh.getMolecule());
		
		for (int i = 0; i < arrIndexNodes.length; i++) {
			PPNodeViz node = new PPNodeViz(mdh.getNode(arrIndexNodes[i]));
			frag.addNode(node);
		}
		
		for (int i = 0; i < arrIndexNodes.length; i++) {
			for (int j = i+1; j < arrIndexNodes.length; j++) {
				byte [] arrHist = mdh.getDistHist(arrIndexNodes[i],arrIndexNodes[j]);
				frag.setDistHist(i,j,arrHist);
			}
		}
		
		frag.setArrIndexParentNodes(arrIndexNodes);
		
		frag.realize();
			
		return frag;
	}

	public static List<MolDistHist> generateSubPharmacophores(MolDistHist mdh, int size){

		List<int[]> liIndSub = CombinationGenerator.getAllOutOf(mdh.getNumPPNodes(), size);
		if(liIndSub == null)
			return null;

		List<MolDistHist> liMolDistHistSub = new ArrayList<>();

		for (int[] arrIndSub : liIndSub) {
			MolDistHist mdhSub = getSubFragment(mdh, arrIndSub);
			liMolDistHistSub.add(mdhSub);
		}

		return liMolDistHistSub;
	}


	public static MolDistHist getSubFragment(MolDistHist mdh, int [] arrIndexNodes){

		MolDistHist frag = new MolDistHist(arrIndexNodes.length);

		for (int i = 0; i < arrIndexNodes.length; i++) {
			PPNode node = new PPNode(mdh.getNode(arrIndexNodes[i]));
			frag.addNode(node);
		}

		for (int i = 0; i < arrIndexNodes.length; i++) {
			for (int j = i+1; j < arrIndexNodes.length; j++) {
				byte [] arrHist = mdh.getDistHist(arrIndexNodes[i],arrIndexNodes[j]);
				frag.setDistHist(i,j,arrHist);
			}
		}

		frag.realize();

		return frag;
	}

	public static class ViolatedConditionsCount {

		private long ccMissedMinRange;
		private long ccViolatedMaxRange;
		private long ccViolatedMinDiffInteractionTypes;


		public ViolatedConditionsCount() {
			ccMissedMinRange=0;
			ccViolatedMaxRange=0;
			ccViolatedMinDiffInteractionTypes=0;
		}

		public void add(ViolatedConditionsCount v){
			ccMissedMinRange +=v.ccMissedMinRange;
			ccViolatedMaxRange +=v.ccViolatedMaxRange;
			ccViolatedMinDiffInteractionTypes +=v.ccViolatedMinDiffInteractionTypes;
		}

		public long getCcMissedMinRange() {
			return ccMissedMinRange;
		}

		public long getCcViolatedMaxRange() {
			return ccViolatedMaxRange;
		}

		public long getCcViolatedMinDiffInteractionTypes() {
			return ccViolatedMinDiffInteractionTypes;
		}
	}

}
