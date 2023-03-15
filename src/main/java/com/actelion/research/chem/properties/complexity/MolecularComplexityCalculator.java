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

import com.actelion.research.chem.StereoMolecule;

public class MolecularComplexityCalculator {

	public static final int TOTAL_CAPACITY = (int)(10*Math.pow(10,6));

	private ExhaustiveFragmentsStatistics exhaustiveFragmentsStatistics;
	
	public MolecularComplexityCalculator() {
				
		exhaustiveFragmentsStatistics = new ExhaustiveFragmentsStatistics(BitArray128.MAX_NUM_BITS, TOTAL_CAPACITY);
		
		exhaustiveFragmentsStatistics.setCollectFragmentIdCodes(false);
	}
	
	/**
	 * 
	 * @param mol
	 * @return NaN if complexity calculation failed or is not possible.
	 */
	public double calculate(StereoMolecule mol){
				
		int bondsNeededForComplexityCalc = ObjectiveExhaustiveStatistics.getNeededNumberOfBondsInFragment(mol);
		
		if(bondsNeededForComplexityCalc > BitArray128.MAX_NUM_BITS) {
			return Double.NaN;
		}
		
		ResultFragmentsStatistic resultFragmentsStatistic = exhaustiveFragmentsStatistics.create(mol, bondsNeededForComplexityCalc);
				
		SummaryFragments summaryFragments = new SummaryFragments(resultFragmentsStatistic); 
		
		return summaryFragments.getComplexityScore();

	}
	
	public void roundUp() throws Throwable{
		exhaustiveFragmentsStatistics.roundUp();
	}


}
