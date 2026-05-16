package org.openmolecules.chem.interaction;

import com.actelion.research.chem.StereoMolecule;

public interface InteractionAtomClassifier {
	int[] classifyAtoms(StereoMolecule mol);
	String getAtomTypeName(int type);
}
