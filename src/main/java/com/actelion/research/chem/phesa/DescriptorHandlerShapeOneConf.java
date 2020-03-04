package com.actelion.research.chem.phesa;




/** 
 * @version: 1.0, February 2018
 * Author: J. Wahl
 * a molecular shape descriptor is a shape-ensemble (molecular shapes of generated conformers of a molecule)
 * it also contains information about the atom-connectivity,elements,... so for any molecular shape, the corresponding
 * molecule can be reconstructed 
 * there are different modi for calculating the shape similarity:
 * 0: only takes into account molecule shape for alignment and overlap calculation
 * 1: alignment solely based on shape, but overlap calculation incorporates pharmacophore overlaps
 * 2: both alignment and overlap calculation take a combined score of shape and pharmacophore overlap
 * 3: only the pharmacophore overlap is used for both aligment and overlap calculation -> this is the fastest method!
 * 19 April 2018: performance enhancement by using a cutoff for the calculation of atom-atom overlaps and preculated exp-values with linear interpolation
 * TODO: add Tversky index
*/


public class DescriptorHandlerShapeOneConf extends DescriptorHandlerShape {

	public DescriptorHandlerShapeOneConf() {
		super(true);
	}
	

	public DescriptorHandlerShapeOneConf(double ppWeight) {
		super(ppWeight);
	}
	
	public DescriptorHandlerShapeOneConf(int maxConfs,double ppWeight) {
		super(true,maxConfs,ppWeight);
	}
	
	@Override
	public DescriptorHandlerShape getThreadSafeCopy() {

		DescriptorHandlerShape dhs = new DescriptorHandlerShapeOneConf();
		dhs.ppWeight = ppWeight;
		dhs.flexible = flexible;
		dhs.maxConfs = maxConfs;
		return dhs;
	}

}