package com.actelion.research.chem.descriptor.flexophore.generator;

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.descriptor.flexophore.MolDistHistViz;
import com.actelion.research.chem.descriptor.flexophore.UnparametrizedAtomTypeException;
import com.actelion.research.util.datamodel.IDCodeCoord;

/**
 * FlexophoreCreateFunctions
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Jan 7, 2013 MvK Start implementation
 * Mar 01 2016 new Flexophore.
 */
public class FlexophoreCreateFunctions {

	public static MolDistHistViz create(String idcode) throws UnparametrizedAtomTypeException {
		
		IDCodeParser parser = new IDCodeParser();
						
		StereoMolecule mol = parser.getCompactMolecule(idcode);
		
		CoordinateInventor coordinateInventor = new CoordinateInventor();
		
		coordinateInventor.invent(mol);
		
		return createDescriptor(mol);
	}
	
	public static MolDistHistViz create(String idcode, String coord) throws UnparametrizedAtomTypeException {
		
		IDCodeParser parser = new IDCodeParser();
						
		StereoMolecule mol = parser.getCompactMolecule(idcode, coord);
		
		if(coord==null) {
			CoordinateInventor coordinateInventor = new CoordinateInventor();
			
			coordinateInventor.invent(mol);
		}
		
		return createDescriptor(mol);
	}
	
	public static MolDistHistViz create(IDCodeCoord idcodeCoord) throws UnparametrizedAtomTypeException {
				
		return create(idcodeCoord.getIdcode(), idcodeCoord.getCoordinates());
	}

	public static MolDistHistViz createDescriptor(StereoMolecule mol) {
		
    	StereoMolecule fragBiggest = mol;
    	
    	fragBiggest.stripSmallFragments();
    	
    	fragBiggest.ensureHelperArrays(StereoMolecule.cHelperCIP);
    	    	
    	MolDistHistViz mdh = null;
        try {
			mdh = CreatorMolDistHistViz.getInstance().create(mol);
		} catch (Exception e) {
			e.printStackTrace();
		}
        return mdh;
    }



}
