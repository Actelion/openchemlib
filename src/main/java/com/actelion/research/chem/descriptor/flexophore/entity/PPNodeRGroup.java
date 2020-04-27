package com.actelion.research.chem.descriptor.flexophore.entity;

import com.actelion.research.chem.descriptor.flexophore.PPNode;

/**
 * PPNodeStructure
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Feb 20, 2013 MvK Start implementation
 */
public class PPNodeRGroup extends PPNode {

	private String idcode;
	
	private String coordinates;
	
	public PPNodeRGroup() {
		super();
		
	}
	
	public PPNodeRGroup(PPNode node, String idcode) {
		super(node);
		this.idcode = idcode;
	}

	
	public String getIdCode() {
		return idcode;
	}
	
	public String getCoordinates() {
		return coordinates;
	}

	
	public void setIdCode(String idcode) {
		this.idcode = idcode;
	}
	
	public void set(String idcode, String coordinates) {
		this.idcode = idcode;
		this.coordinates = coordinates;
	}
	
	

}
