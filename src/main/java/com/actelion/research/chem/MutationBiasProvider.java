package com.actelion.research.chem;

public interface MutationBiasProvider {
	void setBiasReference(StereoMolecule referenceMolecule);
	double getBiasFactor(StereoMolecule mutatedMolecule);
	}
