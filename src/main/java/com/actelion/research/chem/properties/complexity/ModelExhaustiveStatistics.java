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

package com.actelion.research.chem.properties.complexity;

import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.actelion.research.calc.Logarithm;
import com.actelion.research.util.Formatter;

public class ModelExhaustiveStatistics {
	
	public static final String TAG_ModelExhaustiveStatistics = "ModelExhaustiveStatistics";
	
	private static final String TAG_ATTR_NUM_BONDS_IN_FRAGMENT = "numBondsFragment"; 
	
	private static final String TAG_ATTR_NUM_FRAGMENTS = "numFragments"; 
	
	private static final String TAG_ATTR_NUM_UNIQUE = "numUniqueFrags"; 
	
	private static final String TAG_ATTR_RATIO_COVERED_BONDS = "ratioCoveredBonds"; 
	
	
	
	// Number of bonds in the fragment.
	private int numBondsInFrag;
	
	// Number of different fragments.
	private int nFragments;
	
	// Number of unique fragments.
	private int nUnique;
	
	
	private double ratioCoveredBonds;
		
	public ModelExhaustiveStatistics(int numBondsInFrag, int nFragments, int nUnique, double ratioCoveredBonds) {
				
		this.numBondsInFrag = numBondsInFrag;
		
		this.nFragments = nFragments;
		
		this.nUnique = nUnique;
		
		this.ratioCoveredBonds = ratioCoveredBonds;
	}
	
	public int getNumBondsInFragment() {
		return numBondsInFrag;
	}
	
	public int getFragments() {
		return nFragments;
	}
	
	public void setFragments(int nFragments) {
		this.nFragments = nFragments;
	}
	
	public int getUnique() {
		return nUnique;
	}
	
	public void setUnique(int nUnique) {
		this.nUnique = nUnique;
	}
	
	
	public double getRatioCoveredBonds() {
		return ratioCoveredBonds;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[# bonds in fragment=");
		sb.append(numBondsInFrag);
		sb.append(", fragments=");
		sb.append(nFragments);
		sb.append(", unique=");
		sb.append(nUnique);
		sb.append(", non overlapping=");
		sb.append(", ratio covered bonds=");
		sb.append(ratioCoveredBonds);
		sb.append("]");
		return sb.toString();
	}
	
	
	public static Comparator<ModelExhaustiveStatistics> getComparatorNumBonds() {
		
		Comparator<ModelExhaustiveStatistics> comp = new Comparator<ModelExhaustiveStatistics>(){
			 			
		public int compare(ModelExhaustiveStatistics o1, ModelExhaustiveStatistics o2) {
			
			
			if(o1.getNumBondsInFragment() > o2.getNumBondsInFragment()){
				return 1;
			}else if(o1.getNumBondsInFragment() < o2.getNumBondsInFragment()){
				return -1;
			}
			
			return 0;
		}};
		
		return comp;
	}
	
	public static String toString(ResultFragmentsStatistic resultFragmentsStatistic) {
		StringBuilder sb = new StringBuilder();
		
		List<ModelExhaustiveStatistics> liModelExhaustiveStatistics = resultFragmentsStatistic.getExhaustiveStatistics();

		Collections.sort(liModelExhaustiveStatistics, getComparatorNumBonds());
		
		sb.append("Num bonds");
		sb.append("\t");
		sb.append("Fragments");
		sb.append("\t");
		sb.append("Unique");
		sb.append("\t");
		sb.append("Sum isomorphs");
		sb.append("\t");
		sb.append("Ratio coverage by isomorph all");
		sb.append("\t");
		sb.append("Ratio coverage by isomorph max freq idcode");
		sb.append("\t");
		sb.append("ln(unique)");
		sb.append("\n");
		
		int sumFragments = 0;
		
		int sumUniqueFragments = 0;
		
		
		
		for (int i = 0; i < liModelExhaustiveStatistics.size(); i++) {
			ModelExhaustiveStatistics model = liModelExhaustiveStatistics.get(i);
								
			sumFragments += model.getFragments();
			
			sumUniqueFragments += model.getUnique();
			
			sb.append(model.getNumBondsInFragment());
			sb.append("\t");
			sb.append(model.getFragments());
			sb.append("\t");
			sb.append(model.getUnique());
			
			sb.append("\t");
			sb.append(Formatter.format3(model.getRatioCoveredBonds()));
			sb.append("\t");
			double yUnique = Logarithm.get(model.getUnique(), ObjectiveExhaustiveStatistics.BASE_LOG);
			sb.append(Formatter.format3(yUnique));
			sb.append("\n");
		}
		
		sb.append("Sum frags = " + sumFragments + "\tsum distinct " + sumUniqueFragments);
		sb.append("\n");

		int quarter = resultFragmentsStatistic.getBonds() / 4;
		
		sb.append("Bonds quarter " + quarter);
		sb.append("\n");
		sb.append("Bonds 1");
		sb.append("\t");
		sb.append("Bonds 2");
		sb.append("\t");
		sb.append("Slope");
		sb.append("\n");
		
		for (int i = 1; i < liModelExhaustiveStatistics.size(); i++) {
			
			ModelExhaustiveStatistics model0 = liModelExhaustiveStatistics.get(i-1);
			ModelExhaustiveStatistics model1 = liModelExhaustiveStatistics.get(i);
			
			double lnUnique0 = Math.log(model0.getUnique());
			
			double lnUnique1 = Math.log(model1.getUnique());
						
			double slope = lnUnique1 - lnUnique0;
						
			sb.append( model0.getNumBondsInFragment() + "\t" + model1.getNumBondsInFragment() + "\t" + Formatter.format3(slope));
						
			sb.append("\n");
		}

		
		return sb.toString();
	}
	
	public Element getXMLElement(Document doc) throws ParserConfigurationException, DOMException, IOException{
		
        Element nodeRoot  = doc.createElement(TAG_ModelExhaustiveStatistics);
        
        nodeRoot.setAttribute(TAG_ATTR_NUM_BONDS_IN_FRAGMENT, Integer.toString(numBondsInFrag));
        
        nodeRoot.setAttribute(TAG_ATTR_NUM_FRAGMENTS, Integer.toString(nFragments));
        
        nodeRoot.setAttribute(TAG_ATTR_NUM_UNIQUE, Integer.toString(nUnique));
        
        nodeRoot.setAttribute(TAG_ATTR_RATIO_COVERED_BONDS, Formatter.format4(ratioCoveredBonds));
                        
        
                
		return nodeRoot;
	}

	public static ModelExhaustiveStatistics readXMLElement(Element root) throws ParserConfigurationException, DOMException, IOException{
		
		int nBondsInFrag = Integer.parseInt(root.getAttribute(TAG_ATTR_NUM_BONDS_IN_FRAGMENT));
		
		int nFrag = Integer.parseInt(root.getAttribute(TAG_ATTR_NUM_FRAGMENTS));
				
		int nUnique = Integer.parseInt(root.getAttribute(TAG_ATTR_NUM_UNIQUE));
		
		double ratioCovered = Double.parseDouble(root.getAttribute(TAG_ATTR_RATIO_COVERED_BONDS));
		
		
		
		ModelExhaustiveStatistics modelExhaustiveStatistics = new ModelExhaustiveStatistics(nBondsInFrag, nFrag, nUnique, ratioCovered);
		
		
		return modelExhaustiveStatistics;
	}
	
	


}
