package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.chem.descriptor.flexophore.completegraphmatcher.ObjectiveBlurFlexophoreHardMatchUncovered;

import java.awt.*;

public class FlexophoreAtomContributionColors {
	private static final int[] cDiverseColor = {
			0x004900FF, 0x00FF002A, 0x0000AD66, 0x00FFFF00, 0x004FD1F8, 0x00DDA0F6, 0x00FFAB19, 0x00C0FFD0,
			0x00FF00D8, 0x00AA6B58, 0x00BCBCBC, 0x005E6B26, 0x006CFF00, 0x00FFBC97, 0x00B4C732, 0x00A60000 } ;

	private MolDistHist mMolFlexophore,mRefFlexophore;
	private int[] mRefARGB,mMolARGB;
	private float[] mRefRadius,mMolRadius;

	public FlexophoreAtomContributionColors(MolDistHist molFlexophore, MolDistHist refFlexophore, int molAtomCount, int refAtomCount) {
		// we stretch the range, because original values always tend to be larger than 0.75
		final float NODE_SIMILARITY_FACTOR = 4f;

		mMolFlexophore = molFlexophore;
		mRefFlexophore = refFlexophore;

		if (molAtomCount != 0 && refAtomCount != 0 && molFlexophore != null && refFlexophore != null) {
			int[][] molNodeAtom = molFlexophore.getNodeAtoms();
			int[][] refNodeAtom = refFlexophore.getNodeAtoms();

			if (molNodeAtom != null && refNodeAtom != null) {
				DescriptorHandlerFlexophore dhFlexophore = new DescriptorHandlerFlexophore();
				ObjectiveBlurFlexophoreHardMatchUncovered objectiveBlurFlexophoreHardMatchUncovered = dhFlexophore.getObjectiveCompleteGraph();

				ModelSolutionSimilarity modelSolutionSimilarity = dhFlexophore.getBestMatch(molFlexophore, refFlexophore);
				if (modelSolutionSimilarity != null) {
					mMolARGB = new int[molAtomCount];
					mRefARGB = new int[refAtomCount];
					mMolRadius = new float[molAtomCount];
					mRefRadius = new float[refAtomCount];

					int matchingNodeCount = modelSolutionSimilarity.getSizeHeap();

					objectiveBlurFlexophoreHardMatchUncovered.setBase(molFlexophore);
					objectiveBlurFlexophoreHardMatchUncovered.setQuery(refFlexophore);

					int[] colorList = createDiverseColorList(matchingNodeCount);

					for (int i=0; i<matchingNodeCount; i++) {
						int molNode = modelSolutionSimilarity.getIndexBaseFromHeap(i);
						int refNode = modelSolutionSimilarity.getIndexQueryFromHeap(i);

						float nodeSimilarity = modelSolutionSimilarity.getSimilarityNode(i);
						float edgeSimilarity = objectiveBlurFlexophoreHardMatchUncovered.getSimilarityHistogramsForNode(modelSolutionSimilarity, i);

						if (nodeSimilarity > 0f) {
							nodeSimilarity = Math.max(0f, Math.min(1f, Math.round(nodeSimilarity * NODE_SIMILARITY_FACTOR - NODE_SIMILARITY_FACTOR + 1f)));

							for (int atom:molNodeAtom[molNode]) {
								int alpha = Math.round(nodeSimilarity * 255);
								mMolARGB[atom] = (alpha << 24) | colorList[i];
								mMolRadius[atom] = edgeSimilarity;
							}

							for (int atom:refNodeAtom[refNode]) {
								int alpha = Math.round(nodeSimilarity * 255);
								mRefARGB[atom] = (alpha << 24) | colorList[i];
								mRefRadius[atom] = edgeSimilarity;
							}
						}
/*
						float molNodeSimilarity = modelSolutionSimilarity.getSimilarityNode(refNode);
						float refNodeSimilarity = modelSolutionSimilarity.getSimilarityNode(molNode);

						if (molNodeSimilarity > 0f) {
							molNodeSimilarity = Math.max(0f, Math.min(1f, Math.round(molNodeSimilarity * NODE_SIMILARITY_FACTOR - NODE_SIMILARITY_FACTOR + 1f)));
							for (int atom:molNodeAtom[molNode]) {
								int alpha = Math.round(molNodeSimilarity * 255);
								mMolARGB[atom] = (alpha << 24) | colorList[i];
								mMolRadius[atom] = edgeSimilarity;
							}
						}

						if (refNodeSimilarity > 0f) {
							refNodeSimilarity = Math.max(0f, Math.min(1f, Math.round(refNodeSimilarity * NODE_SIMILARITY_FACTOR - NODE_SIMILARITY_FACTOR + 1f)));
							for (int atom:refNodeAtom[refNode]) {
								int alpha = Math.round(refNodeSimilarity * 255);
								mRefARGB[atom] = (alpha << 24) | colorList[i];
								mRefRadius[atom] = edgeSimilarity;
							}
						}
 */
					}
				}
			}
		}
	}

	public MolDistHist getMolFlexophore() {
		return mMolFlexophore;
	}

	public MolDistHist getRefFlexophore() {
		return mRefFlexophore;
	}

	public int[] getMolARGB() {
		return mMolARGB;
	}

	public int[] getRefARGB() {
		return mRefARGB;
	}

	public float[] getMolRadius() {
		return mMolRadius;
	}

	public float[] getRefRadius() {
		return mRefRadius;
	}

	public static int[] createDiverseColorList(int colorCount) {
		if (colorCount <= cDiverseColor.length)
			return cDiverseColor;

		int divisor = ((colorCount & 1) == 0) ? colorCount : colorCount+1;
		final float[] s = { 1.0f, 0.4f, 0.8f, 1.0f, 1.0f, 0.6f, 0.8f, 1.0f };
		final float[] b = { 0.8f, 1.0f, 1.0f, 0.6f, 0.8f, 1.0f, 1.0f, 0.4f };
		int[] colorList = new int[colorCount];
		for (int i=0; i<colorCount; i++) {
			float hue = (float)i / (float)divisor;
			colorList[i] = Color.HSBtoRGB(hue, s[i & 7], b[i & 7]);
		}

		return colorList;
	}
}