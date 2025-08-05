/*
 * Copyright 2025 Thomas Sander, openmolecules.org
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */
package com.actelion.research.chem.io.pdb.converter;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

/**
 * GraphWalker is a versatile depth-first graph walker that starts traversing
 * the molecular graph of a StereoMolecule from a given root atom exploring
 * all possible paths. It is guided by two overridable methods that define whether
 * a certain bond/atom combination qualifies to be added to the path and/or qualify
 * to finish the path. Paths can be custom scored, which allows to retrieve the best
 * scoring path once all paths have been found.
 * */
public class GraphWalker {
	private final StereoMolecule mMol;
	private final int mMaxDepth;
	private final boolean[] mIsUsedAtom;
	private final int[] mParentAtom,mConnIndex;
	private double mBestScore;
	private int[] mBestScoringPath;

	public GraphWalker(StereoMolecule mol, int maxDepth) {
		mMol = mol;
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
		mMaxDepth = maxDepth;

		mIsUsedAtom = new boolean[mMol.getAtoms()];
		mParentAtom = new int[mMol.getAtoms()];
		mConnIndex = new int[mMol.getAtoms()];
	}

	/**
	 * Starting from the given rootAtom this method locates all possible paths
	 * traversing neighbours and neighbours of neighbours whenever qualifiesAsNext()
	 * returns true. Whenever an atom was added to a path, qualifiesAsFinal() is called.
	 * If that method returns true, then calculatePathScore() is called to score the current
	 * path. Then, path detection continues until all possible paths have been scored.
	 * @param rootAtom first atom of path
	 * @param nextAtom -1 or neighbour of rootAtom, if only one direction shall be traversed
	 */
	public void find(int rootAtom, int nextAtom) {
		int level = 0;
		int parentAtom = rootAtom;
		mIsUsedAtom[rootAtom] = true;
		mParentAtom[rootAtom] = -1;
		mConnIndex[0] = 0;
		mBestScore = Double.NaN;

		while (level >= 0) {
			if (mConnIndex[level] == mMol.getConnAtoms(parentAtom)) {
				mIsUsedAtom[parentAtom] = false;
				parentAtom = mParentAtom[parentAtom];
				level--;
				continue;
			}

			int atom = mMol.getConnAtom(parentAtom, mConnIndex[level]);
			int bond = mMol.getConnBond(parentAtom, mConnIndex[level]);
			mConnIndex[level]++;

			if (!mIsUsedAtom[atom]
			 && (level != 0 || nextAtom == -1 || atom == nextAtom)
			 && ((level == 0 && qualifiesAsFirst(parentAtom, atom, bond))
			  || (level != 0 && qualifiesAsNext(parentAtom, atom, bond, level+2)))) {
				mParentAtom[atom] = parentAtom;

				if (qualifiesAsFinal(parentAtom, atom, bond, level+2)) {
					double score = calculatePathScore(atom, mParentAtom);
					if (Double.isNaN(mBestScore) || mBestScore < score) {
						mBestScore = score;
						mBestScoringPath = new int[level+2];
						int pathAtom = atom;
						for (int i=level+1; i>=0; i--) {
							mBestScoringPath[i] = pathAtom;
							pathAtom = mParentAtom[pathAtom];
						}
					}
				}

				if (level < mMaxDepth) {
					mIsUsedAtom[atom] = true;
					parentAtom = atom;
					level++;
					mConnIndex[level] = 0;
				}
			}
		}
	}

	/**
	 * Overwrite this method, if the first atom/bond to be added to the graph needs different conditions
	 * than any subsequent atom/bond. Don't override if the first atom/bond criteria are not different.
	 * @param parentAtom second front atom of current path
	 * @param atom front atom of current graph
	 * @param bond front bond of current graph connecting parentAtom and atom
	 * @return whether atom qualifies to be added via bond to the current graph ending with parentAtom
	 */
	public boolean qualifiesAsFirst(int parentAtom, int atom, int bond) {
		return qualifiesAsNext(parentAtom, atom, bond, 2);
	}

	/**
	 * Overwrite this method to control, which atoms/bonds qualify to be added to a growing graph.
	 * @param parentAtom second front atom of current path
	 * @param atom front atom of current graph
	 * @param bond front bond of current graph connecting parentAtom and atom
	 * @param size number of atoms in current graph
	 * @return whether atom qualifies to be added via bond to the current graph ending with parentAtom
	 */
	public boolean qualifiesAsNext(int parentAtom, int atom, int bond, int size) {
		return true;
	}

	/**
	 * Overwrite this method to specify additional criterion to qualifiesAsNext(), which are given for an
	 * atom/bond to qualify as the last atom/bond in a path.
	 * @param parentAtom second front atom of current path
	 * @param atom front atom of current graph
	 * @param bond front bond of current graph connecting parentAtom and atom
	 * @param size number of atoms in current graph
	 * @return whether atom qualifies as the last atom of the growing graph
	 * @return
	 */
	public boolean qualifiesAsFinal(int parentAtom, int atom, int bond, int size) {
		return true;
	}

	/**
	 * Overwrite this method if you need to score individual paths and
	 * to retrieve the best scoring path later from this PathFinder.
	 * Scores must be above 0.0 with higher scores for better paths.
	 * @param finalAtom
	 * @param graphParent contains parent atom for every path member and -1 for first atom
	 * @return
	 */
	public double calculatePathScore(int finalAtom, int[] graphParent) {
		return 0;
	}

	/**
	 * @return NaN if no scoring was done or score of best scoring path
	 */
	public double getBestScore() {
		return mBestScore;
	}

	/**
	 * If calculatePathScore() was used to score individual paths, then this method can be used$
	 * to retrieve an array of all path atoms from the root atom to the final atom.
	 * @return array of all atoms of the best scoring path
	 */
	public int[] getBestScoringPath() {
		return mBestScoringPath;
	}
}

// Usage example to find the longest path of conjugated double bonds for StereoMolecule mol:
//	private void findLongestConjugatedAtomChain(StereoMolecule mol) {
//		GraphWalker walker = new GraphWalker(mol, 12) {
//			@Override public boolean qualifiesAsNext(int parentAtom, int atom, int bond, int size) {
//				return (((size & 1) == 1) ^ (mol.getBondOrder(bond) > 1));
//			}
//
//			@Override public boolean qualifiesAsFinal(int parentAtom, int atom, int bond, int size) {
//				return (size & 1) == 0;
//			}
//
//			@Override public double calculatePathScore(int finalAtom, int[] graphParent) {
//				double score = 0.0;
//				while (finalAtom != -1) {
//					score += 1.0;
//					finalAtom = graphParent[finalAtom];
//				}
//				return score;
//			}
//		};
//
//		for (int atom=0; atom<mol.getAtoms(); atom++) {
//			if (mol.getAtomPi(atom) == 1) {
//				walker.find(atom, -1);
//				int[] bestPath = walker.getBestScoringPath();
//				System.out.print("Longest conjugated chain:");
//				for (int a : bestPath)
//					System.out.print(" "+a);
//				System.out.println();
//			}
//		}
//	}
