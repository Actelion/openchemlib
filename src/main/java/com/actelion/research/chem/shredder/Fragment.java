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

package com.actelion.research.chem.shredder;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.datamodel.IDCodeCoord;

/**
 * Fragment
 * Jan 18, 2013 MvK Start implementation
 */
public class Fragment extends IDCodeCoord {
	
    public static final String TAG_FREQUENCY_ONE_PER_MOLECULE = "FrequencyOnePerMolecule";
    
    public static final String TAG_FREQUENCY = "FrequencySumAll";
    
    public static final String TAG_RELATIVE_FREQUNCY = "Relative frequencySumAll";
    
    public static final String TAG_SIZE = "Size fragment"; 

    // This frequency term is only incremented ones per molecule.
	private int frequencyOnePerMol;
	
	// This frequency term counts all fragments. One molecule can have several identical fragments. 
	private int frequencySumAll;
	
	private int size;
	
	private StereoMolecule mol;
	

	public Fragment(String idcode) {
		super(idcode);
	}

	
	/**
	 * @param idcode
	 * @param coordinates
	 */
	public Fragment(String idcode, String coordinates) {
		super(idcode, coordinates);
	}



	/**
	 * @return the frequencyOnePermol
	 */
	public int getFrequencyOnePerMol() {
		return frequencyOnePerMol;
	}


	/**
	 * @param frequencyOnePerMol the frequencyOnePermol to set
	 */
	public void setFrequencyOnePerMol(int frequencyOnePerMol) {
		this.frequencyOnePerMol = frequencyOnePerMol;
	}


	/**
	 * @return the frequency
	 */
	public int getFrequencySumAll() {
		return frequencySumAll;
	}


	/**
	 * @param frequency the frequency to set
	 */
	public void setFrequencySumAll(int frequency) {
		this.frequencySumAll = frequency;
	}

	public void incrementFrequencySumAll() {
		this.frequencySumAll++;
	}
	
	public void incrementFrequencyOnePerMol() {
		this.frequencyOnePerMol++;
	}
	

	/**
	 * @return the size
	 */
	public int getSize() {
		return size;
	}


	/**
	 * @param size the size to set
	 */
	public void setSize(int size) {
		this.size = size;
	}
	
	public boolean equals(Object o) {
		
		boolean equals = true;
		
		if(!(o instanceof Fragment)){
			return false;
		}
		
		Fragment f = (Fragment) o;
		
		if(!idcode.equals(f.idcode)){
			equals = false;
		}
				
		return equals;
	}
	
	/**
	 * @return the mol
	 */
	public StereoMolecule getMol() {
		return mol;
	}


	/**
	 * @param mol the mol to set
	 */
	public void setMol(StereoMolecule mol) {
		this.mol = mol;
	}



}
