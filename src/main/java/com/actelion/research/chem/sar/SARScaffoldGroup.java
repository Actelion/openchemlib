package com.actelion.research.chem.sar;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;

/**
 * A scaffold group comprises all scaffolds that arise from a match of the same query substructure in multiple
 * processed molecules. A scaffold group may contain more than one scaffold, if the substructure contains wild card
 * elements like atom lists, or multiple allowed bond orders. A special case are bridge bonds, which cause a query
 * match to contain more atoms than the query itself. Interestingly, these atoms may also carry R-groups, which needs
 * to be handled differently from the R-groups on the core atoms. Core atoms are those atoms for which an associated
 * atom exists in the query structure. R-groups (exit vectors) on core atom are numbered consistently within all
 * scaffolds that belong to the same scaffold group, i.e. that were detected from the same query structure.
 */
public class SARScaffoldGroup {
	private int mRGroupCount;
	private ExitVector[] mExitVector;
	private ArrayList<SARScaffold> mScaffoldList;

	protected SARScaffoldGroup(StereoMolecule query) {
		super();
		mRGroupCount = -1;
		analyzeExitVectors(query);
		mScaffoldList = new ArrayList<>();
	}

	public void addScaffold(SARScaffold scaffold) {
		mScaffoldList.add(scaffold);
	}

	public ArrayList<SARScaffold> getScaffoldList() {
		return mScaffoldList;
	}

	private void analyzeExitVectors(StereoMolecule query) {
		ArrayList<ExitVector> evList = new ArrayList<>();
		for (int atom=0; atom<query.getAtoms(); atom++) {
			if ((query.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) == 0) {
				int exitVectorCount = query.getLowestFreeValence(atom);
				for (int i=0; i<exitVectorCount; i++) {
					// If we have >= two exit vectors (with both one and two connAtoms in the query)
					// we assume that we can distinuish the exit vectors by topicity (any stereo criteria)
					int topicity = (exitVectorCount >= 2) ? i : -1;
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
			 && (((topicity == -1 || mExitVector[i].getTopicity() == -1) && mExitVector[i].getIndex() == connIndex)
			  || (topicity != -1 && mExitVector[i].getTopicity() == topicity)))
//{ System.out.println("getExitVectorIndex(queryAtom:"+queryAtom+", connIndex:"+connIndex+", topicity:"+topicity+") ev.index:"+mExitVector[i].getIndex()+" ev.topicity:"+mExitVector[i].getTopicity()+" evi:"+i);
				return i;
//}

//System.out.println("getExitVectorIndex(queryAtom:"+queryAtom+", connIndex:"+connIndex+", topicity:"+topicity+") evi:"+-1);
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
