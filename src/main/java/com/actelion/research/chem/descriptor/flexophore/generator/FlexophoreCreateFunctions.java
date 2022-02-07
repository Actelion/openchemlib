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

import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.descriptor.flexophore.MolDistHistViz;
import com.actelion.research.chem.descriptor.flexophore.UnparametrizedAtomTypeException;
import com.actelion.research.util.datamodel.IDCodeCoord;

/**
 * FlexophoreCreateFunctions
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
