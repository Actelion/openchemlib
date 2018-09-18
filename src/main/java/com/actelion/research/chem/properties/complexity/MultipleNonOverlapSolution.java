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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import com.actelion.research.util.datamodel.IntVec;

public class MultipleNonOverlapSolution {
	
	public static final String TAG_MultipleNonOverlapSolution = "MultipleNonOverlapSolution";
	
	private static final String TAG_ATTR_NUM_EQUAL_IDCODES = "equalIdCodes"; 
	
	private static final String TAG_INDEX_SOLUTION_CHECKED = "IndexSolution";
	
	private static final String TAG_ATTR_INDEX_SOLUTION_CHECKED = "intVecSolution";
	
	private static final String TAG_CONTAINER = "Container";
	
	private static final String TAG_ATTR_CONTAINER = "intVecContainer";
	
	private static final String TAG_LIST_ISOMORPH = "ListNonOverlapIsomorph";
	
	private static final String TAG_ISOMORPH = "NonOverlapIsomorph";
	
	private static final String TAG_ATTR_ISOMORPH = "isomorph";
	
	// Bit list
	private IntVec arrIndexSolutionChecked;
	
	// Bonds from all solutions are collected in this instance.
	private FragmentDefinedByBondsIdCode container;
	
	// List with solutions
	private List<FragmentDefinedByBondsIdCode> liNonOverLappingIsomorphSubstruct;
	
	private int numEqualIdCodes;
	
	
	
	private MultipleNonOverlapSolution() {
		
	}
	
	/**
	 * 
	 * @param numEqualIdCodes
	 * @param livIdCodeInit
	 * @param indexSolution 
	 */
	public MultipleNonOverlapSolution(int numEqualIdCodes, FragmentDefinedByBondsIdCode livIdCodeInit, int indexSolution){
		
		this.numEqualIdCodes = numEqualIdCodes;
		
		arrIndexSolutionChecked = new IntVec(IntVec.getSizeForBits(numEqualIdCodes));
		
		arrIndexSolutionChecked.setBit(indexSolution);
		
		container = new FragmentDefinedByBondsIdCode(livIdCodeInit);
		
		container.addBits(livIdCodeInit.getBitArray());
		
		liNonOverLappingIsomorphSubstruct = new ArrayList<FragmentDefinedByBondsIdCode>();
		
		liNonOverLappingIsomorphSubstruct.add(livIdCodeInit);
	}
	
	public MultipleNonOverlapSolution(MultipleNonOverlapSolution solution){
		
		this.numEqualIdCodes = solution.numEqualIdCodes;

		arrIndexSolutionChecked = new IntVec(solution.arrIndexSolutionChecked);
		
		container = new FragmentDefinedByBondsIdCode(solution.container);
					
		liNonOverLappingIsomorphSubstruct = new ArrayList<FragmentDefinedByBondsIdCode>();
		
		liNonOverLappingIsomorphSubstruct.addAll(solution.liNonOverLappingIsomorphSubstruct);
	}
	
	boolean isSolutionIncluded(int indexSolutionIdCode){
		return arrIndexSolutionChecked.isBitSet(indexSolutionIdCode);
	}
	
	void addSolution(FragmentDefinedByBondsIdCode listWithIntVecIdCode, int indexSolution){
		
		arrIndexSolutionChecked.setBit(indexSolution);
		
		container.addBits(listWithIntVecIdCode.getBitArray());
			
		liNonOverLappingIsomorphSubstruct.add(listWithIntVecIdCode);
			
	}
	
	boolean isOverlap(FragmentDefinedByBondsIdCode listWithIntVecIdCode){
		
		if(container.isOverlappingBits(listWithIntVecIdCode.getBitArray())){
			
			return true;
		}
		
		return false;
	}
	
	public int hashCode() {
		return container.hashCode();
	}

	public boolean equals(Object obj) {
		
		if(!(obj instanceof MultipleNonOverlapSolution)){
			return false;
		}
		
		MultipleNonOverlapSolution m = (MultipleNonOverlapSolution)obj;
		
		if(!arrIndexSolutionChecked.equal(m.arrIndexSolutionChecked)){
			return false;
		}
		
		return true;
	}
	
	boolean isAllChecked(){
		
		int sum = 0;
		
		for (int i = 0; i < numEqualIdCodes; i++) {
			if(arrIndexSolutionChecked.isBitSet(i)){
				sum++;
			}
		}
		
		return (sum==numEqualIdCodes) ? true : false;
	}
	
	public IntVec getCheckerArray(){
		return arrIndexSolutionChecked;
	}
	
	/**
	 * @return the container
	 */
	public FragmentDefinedByBondsIdCode getContainer() {
		return container;
	}

	/**
	 * @return the liNonOverLappingIsomorphSubstruct
	 */
	public List<FragmentDefinedByBondsIdCode> getLiNonOverLappingIsomorphSubstruct() {
		return liNonOverLappingIsomorphSubstruct;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(arrIndexSolutionChecked.toStringBinaryDense());
		return sb.toString();
	}
	
	public Element getXMLElement(Document doc) throws ParserConfigurationException, DOMException, IOException{
				
        Element nodeRoot  = doc.createElement(TAG_MultipleNonOverlapSolution);
        
        nodeRoot.setAttribute(TAG_ATTR_NUM_EQUAL_IDCODES, Integer.toString(numEqualIdCodes));
        
        Element nodeIndexSolution = doc.createElement(TAG_INDEX_SOLUTION_CHECKED);

        nodeIndexSolution.setAttribute(TAG_ATTR_INDEX_SOLUTION_CHECKED, arrIndexSolutionChecked.write2String());
		
        Element nodeContainer  = doc.createElement(TAG_CONTAINER);

        nodeContainer.setAttribute(TAG_ATTR_CONTAINER, container.write2String());
                
        Element nodeListIsomorph  = doc.createElement(TAG_LIST_ISOMORPH);

        for (FragmentDefinedByBondsIdCode livIdCode : liNonOverLappingIsomorphSubstruct) {
			
            Element nodeIsomorph  = doc.createElement(TAG_ISOMORPH);
            
            nodeIsomorph.setAttribute(TAG_ATTR_ISOMORPH, livIdCode.write2String());

            nodeListIsomorph.appendChild(nodeIsomorph);
		}
                        
        nodeRoot.appendChild(nodeIndexSolution);
        nodeRoot.appendChild(nodeContainer);
        nodeRoot.appendChild(nodeListIsomorph);
        
        
		return nodeRoot;
	}

	public static MultipleNonOverlapSolution readXMLElement(Element root) throws ParserConfigurationException, DOMException, IOException{
		
		int numEqualIdCodes = Integer.parseInt(root.getAttribute(TAG_ATTR_NUM_EQUAL_IDCODES));
		
		
		Element nodeIndexSolution = null;
		Element nodeContainer = null;
		Element nodeListIsomorph = null;
		
		NodeList nl = root.getChildNodes();
		
		for (int i = 0; i < nl.getLength(); i++) {
			
			Node node = nl.item(i);
			
			String nodeName = node.getNodeName();
			
			if(MultipleNonOverlapSolution.TAG_INDEX_SOLUTION_CHECKED.equals(nodeName)){
				nodeIndexSolution = (Element)node;
			} else if(MultipleNonOverlapSolution.TAG_CONTAINER.equals(nodeName)){
				nodeContainer = (Element)node;
			} else if(MultipleNonOverlapSolution.TAG_LIST_ISOMORPH.equals(nodeName)){
				nodeListIsomorph = (Element)node;
			}
 
			
		}
		
		InputStream isIndexSolutionChecked = new ByteArrayInputStream(nodeIndexSolution.getAttribute(TAG_ATTR_INDEX_SOLUTION_CHECKED).getBytes());
		
		IntVec arrIndexSolutionChecked = IntVec.read(isIndexSolutionChecked);
		
		isIndexSolutionChecked.close();
		
		
		InputStream isContainer = new ByteArrayInputStream(nodeContainer.getAttribute(TAG_ATTR_CONTAINER).getBytes());

		FragmentDefinedByBondsIdCode container = FragmentDefinedByBondsIdCode.read(isContainer);
		
		isContainer.close();
		
		List<FragmentDefinedByBondsIdCode> liNonOverLappingIsomorphSubstruct = new ArrayList<FragmentDefinedByBondsIdCode>();
		
		
		NodeList nlIsomorph = nodeListIsomorph.getChildNodes();
		
		for (int i = 0; i < nlIsomorph.getLength(); i++) {
			
			Node node = nlIsomorph.item(i);
			
			String nodeName = node.getNodeName();
			
			if(TAG_ISOMORPH.equals(nodeName)){
				Element nodeIsomorph = (Element)node;
								
				InputStream isIsomorph = new ByteArrayInputStream(nodeIsomorph.getAttribute(TAG_ATTR_ISOMORPH).getBytes());

				FragmentDefinedByBondsIdCode livIdCodeIsomorph = FragmentDefinedByBondsIdCode.read(isIsomorph);
				
				isIsomorph.close();
				
				liNonOverLappingIsomorphSubstruct.add(livIdCodeIsomorph);
				
			}
 
			
		}
		
		MultipleNonOverlapSolution multipleNonOverlapSolution = new MultipleNonOverlapSolution();
		
		multipleNonOverlapSolution.numEqualIdCodes = numEqualIdCodes;
		
		multipleNonOverlapSolution.arrIndexSolutionChecked = arrIndexSolutionChecked;
		
		multipleNonOverlapSolution.container = container;
		
		multipleNonOverlapSolution.liNonOverLappingIsomorphSubstruct = liNonOverLappingIsomorphSubstruct;
		
		
		return multipleNonOverlapSolution;
	}
	
	
	
}