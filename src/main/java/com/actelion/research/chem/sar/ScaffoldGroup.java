package com.actelion.research.chem.sar;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;

public class ScaffoldGroup extends ArrayList<ScaffoldData> {
	private int mRGroupCount;
	private ExitVector[] mExitVector;

	protected ScaffoldGroup(StereoMolecule query) {
		super();
		mRGroupCount = -1;
		analyzeExitVectors(query);
	}

	private void analyzeExitVectors(StereoMolecule query) {
		ArrayList<ExitVector> evList = new ArrayList<>();
		for (int atom=0; atom<query.getAtoms(); atom++) {
			if ((query.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) == 0) {
				int exitVectorCount = query.getLowestFreeValence(atom);
				for (int i=0; i<exitVectorCount; i++) {
					// If we have two exit vectors (with both one and two connAtoms in the query)
					// we assume that we can distinuish the exit vectors by topicity (any stereo criteria)
					int topicity = (exitVectorCount == 2) ? i : -1;
					evList.add(new ExitVector(atom, true, i, topicity));
				}
			}
		}
		mExitVector = evList.toArray(new ExitVector[0]);
	}

	protected int getExitVectorCount() {
		return mExitVector.length;
	}

	/**
	 * Find correct exit vector index defining the core atom and for it either the exit vectors topicity
	 * (if exit vectors are stereo-heterotop) or just an index (if exit vectors are homotop)
	 * @param queryAtom respective atom index of query
	 * @param connIndex 0-based exo-query neighbour index (not used if topicity != -1 and neighbours are stereotop
	 * @param topicity -1, if exit vector neighbours are homotop, otherise 0 or 1
	 * @return index into list of all exit vectors of scaffold
	 */
	protected int getExitVectorIndex(int queryAtom, int connIndex, int topicity) {
		for (int i=0; i<mExitVector.length; i++)
			if (mExitVector[i].getQueryAtom() == queryAtom
			 && ((topicity == -1 && mExitVector[i].getIndex() == connIndex)
			  || (topicity != -1 && mExitVector[i].getTopicity() == topicity)))
//{ System.out.println("getExitVectorIndex(queryAtom:"+queryAtom+", connIndex:"+connIndex+", topicity:"+topicity+") ev.index:"+mExitVector[i].getIndex()+" ev.topicity:"+mExitVector[i].getTopicity()+" evi:"+i);
				return i;
//}

System.out.println("getExitVectorIndex(queryAtom:"+queryAtom+", connIndex:"+connIndex+", topicity:"+topicity+") evi:"+-1);
		return -1;
	}

	protected ExitVector getExitVector(int i) {
		return mExitVector[i];
	}

	protected int getRGroupCount() {
		return mRGroupCount;
	}

	protected int assignRGroups(int firstRGroup) {
		if (mRGroupCount == -1) {
			mRGroupCount = 0;
			for (int i=0; i<mExitVector.length; i++) {
				ExitVector exitVector = mExitVector[i];
				if (exitVector.substituentVaries()) {
					exitVector.setRGroupNo(++mRGroupCount + firstRGroup - 1);
//					if (!exitVector.hasSeenBondOrder(1)) {
//						i++;
//						if (!exitVector.hasSeenBondOrder(2))
//							i++;
//					}
				}
			}
		}
		return mRGroupCount;
	}
}
