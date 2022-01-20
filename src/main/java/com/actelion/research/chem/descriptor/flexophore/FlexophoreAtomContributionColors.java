package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.chem.descriptor.flexophore.completegraphmatcher.ObjectiveBlurFlexophoreHardMatchUncovered;

public class FlexophoreAtomContributionColors {
	private int[] mARGB;
	private float[] mRadius;

	public FlexophoreAtomContributionColors(StereoMolecule mol, MolDistHist molFlexophore, MolDistHist refFlexophore) {
		if (mol != null && molFlexophore != null && refFlexophore != null) {
			mol.ensureHelperArrays(Molecule.cHelperRings);
			float[][] atomContribution = getFlexophoreAtomContributions(molFlexophore, refFlexophore, mol.getAllAtoms());
			if (atomContribution != null) {
				mARGB = new int[mol.getAtoms()];
				mRadius = new float[mol.getAtoms()];
				for (int atom=0; atom<mol.getAtoms(); atom++) {
					if (atomContribution[0][atom] > 0f) {
						final float f = 4f;
						// we stretch the range, because original values always tend to be larger than 0.75
						float nodeSimilarity = Math.max(0f, Math.min(1f, Math.round(atomContribution[0][atom] * f - f + 1f)));
//						int color = Color.HSBtoRGB(nodeSimilarity, 1f, 1f);
						int alpha = Math.round(nodeSimilarity * 255);
						mARGB[atom] = (alpha << 24) | 255;
						mRadius[atom] = atomContribution[1][atom];
						}
					}
				}
			}
		}

	public int[] getARGB() {
		return mARGB;
		}

	public float[] getRadius() {
		return mRadius;
		}

	private float[][] getFlexophoreAtomContributions(MolDistHist mdhMol, MolDistHist mdhRef, int molAtomCount) {
		int[][] molNodeAtom = mdhMol.getNodeAtoms();
		if (molNodeAtom == null)
			return null;

		DescriptorHandlerFlexophore dhFlexophore = new DescriptorHandlerFlexophore();
		ObjectiveBlurFlexophoreHardMatchUncovered objectiveBlurFlexophoreHardMatchUncovered = dhFlexophore.getObjectiveCompleteGraph();

		ModelSolutionSimilarity modelSolutionSimilarity = dhFlexophore.getBestMatch(mdhMol, mdhRef);
		if (modelSolutionSimilarity == null)
			return null;

		int heap = modelSolutionSimilarity.getSizeHeap();

		objectiveBlurFlexophoreHardMatchUncovered.setBase(mdhMol);
		objectiveBlurFlexophoreHardMatchUncovered.setQuery(mdhRef);

		float[][] atomContribution = new float[2][molAtomCount];

		for (int i=0; i<heap; i++) {
			int molNode = modelSolutionSimilarity.getIndexBaseFromHeap(i);
			int refNode = modelSolutionSimilarity.getIndexQueryFromHeap(i);

			for (int atom:molNodeAtom[molNode]) {
				atomContribution[0][atom] = modelSolutionSimilarity.getSimilarityNode(refNode);
				atomContribution[1][atom] = objectiveBlurFlexophoreHardMatchUncovered.getSimilarityHistogramsForNode(modelSolutionSimilarity, i);
				}
			}

		return atomContribution;
		}
	}