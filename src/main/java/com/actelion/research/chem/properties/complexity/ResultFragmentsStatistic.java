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
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;

public class ResultFragmentsStatistic {
	
	public static final String TAG_ResultFragmentsStatistic = "ResultFragmentsStatistic";
	
	private static final String TAG_ATTR_NAME = "Name";
	
	private static final String TAG_ATTR_IDCODE = "IdCode";
	
	private static final String TAG_ATTR_COORD = "Coordinates"; 
	
	private static final String TAG_LIST_ModelExhaustiveStatistics = "ListModelExhaustiveStatistics";
			
	private String name;
	
	private StereoMolecule mol;
	
	private List<ModelExhaustiveStatistics> liModelExhaustiveStatistics;
	
	/**
	 * @param liModelExhaustiveStatistics
	 */
	public ResultFragmentsStatistic(StereoMolecule mol,	List<ModelExhaustiveStatistics> liModelExhaustiveStatistics) {
		super();
		
		
		this.mol = mol;
		this.liModelExhaustiveStatistics = liModelExhaustiveStatistics;
	}


	/**
	 * @return the nAtomsMolecule
	 */
	public int getAtoms() {
		return mol.getAtoms();
	}


	/**
	 * @return the nBondsMolecule
	 */
	public int getBonds() {
		return mol.getBonds();
	}


	/**
	 * @return the liModelExhaustiveStatistics
	 */
	public List<ModelExhaustiveStatistics> getExhaustiveStatistics() {
		return liModelExhaustiveStatistics;
	}
	
	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}


	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}


	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		if(name!=null){
			sb.append("Name " + name);
			sb.append("\n");
		}
		
		sb.append("Atoms " + getAtoms() + ", bonds " + getBonds() + ".");
		
		sb.append("\n");
		
		sb.append(ModelExhaustiveStatistics.toString(this));
		
		
		return sb.toString();
	}
	
	public Element getXMLElement(Document doc) throws ParserConfigurationException, DOMException, IOException{
		
        Element nodeRoot  = doc.createElement(TAG_ResultFragmentsStatistic);
        
        nodeRoot.setAttribute(TAG_ATTR_NAME, name);
        
        Canonizer can = new Canonizer(mol);
        
        nodeRoot.setAttribute(TAG_ATTR_IDCODE, can.getIDCode());
        
        nodeRoot.setAttribute(TAG_ATTR_COORD, can.getEncodedCoordinates());
                                
        Element nodeListModelExhaustiveStatistics  = doc.createElement(TAG_LIST_ModelExhaustiveStatistics);

        for (ModelExhaustiveStatistics model : liModelExhaustiveStatistics) {
			
            Element nodeModelExhaustiveStatistics  = model.getXMLElement(doc);
            
            

            nodeListModelExhaustiveStatistics.appendChild(nodeModelExhaustiveStatistics);
		}
                        
        nodeRoot.appendChild(nodeListModelExhaustiveStatistics);
        
		return nodeRoot;
	}
	
	/**
	 * @return the mol
	 */
	public StereoMolecule getMol() {
		return mol;
	}


	public static ResultFragmentsStatistic readXMLElement(Element root) throws ParserConfigurationException, DOMException, IOException{
				
		String name = root.getAttribute(TAG_ATTR_NAME);
		
		String idcode = root.getAttribute(TAG_ATTR_IDCODE);
		
		String coord = root.getAttribute(TAG_ATTR_COORD);
		
		
		NodeList nl = root.getChildNodes();
		
		int len = nl.getLength();
		
		Element nodeListModelExhaustiveStatistics = null; 
		
		for (int i = 0; i < len; i++) {
			Node node = nl.item(i);
			
			String nodeName = node.getNodeName();
			
			if(TAG_LIST_ModelExhaustiveStatistics.equals(nodeName)){
				nodeListModelExhaustiveStatistics = (Element)node;
				break;
			}
			
		}
		
		
		NodeList nlModelExhaustiveStatistics = nodeListModelExhaustiveStatistics.getChildNodes();
		
		List<ModelExhaustiveStatistics> liModelExhaustiveStatistics = new ArrayList<ModelExhaustiveStatistics>();

		for (int i = 0; i < nlModelExhaustiveStatistics.getLength(); i++) {
			Node node = nlModelExhaustiveStatistics.item(i);
			
			String nodeName = node.getNodeName();
			
			if(ModelExhaustiveStatistics.TAG_ModelExhaustiveStatistics.equals(nodeName)){
				Element nodeModelExhaustiveStatistics = (Element)node;
				
				ModelExhaustiveStatistics modelExhaustiveStatistics = ModelExhaustiveStatistics.readXMLElement(nodeModelExhaustiveStatistics);
				
				liModelExhaustiveStatistics.add(modelExhaustiveStatistics);
				
			}
			
		}
		
		IDCodeParser parser = new IDCodeParser(false);
		
		StereoMolecule mol = parser.getCompactMolecule(idcode, coord);
		
		ResultFragmentsStatistic resultFragmentsStatistic = new ResultFragmentsStatistic(mol, liModelExhaustiveStatistics);
		
		resultFragmentsStatistic.setName(name);
		

		return resultFragmentsStatistic;
	}


}
