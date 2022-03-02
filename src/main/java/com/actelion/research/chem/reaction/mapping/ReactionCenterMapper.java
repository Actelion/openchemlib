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
 * @author Thomas Sander
 */

package com.actelion.research.chem.reaction.mapping;

import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;
import java.util.TreeMap;

/**
 * This is a helper class for the RuleEnhancedMapper, which expects most of a reaction's atoms
 * to be already mapped using an MCS or similarity guided graph-matching approach.
 * This class tries all permutations of unmapped reactant and product atoms provided their
 * atomic numbers and atom masses match. All solutions are scored considering bonds changed
 * and the best scoring solution is added to the mapping number arrays.
 * This class assumes a stoichiometrical reaction and, thus, always the same number of unmapped
 * reactant atoms and products atoms in the same atom class.
 */
public class ReactionCenterMapper {
	private static final int MAX_PERMUTATION_COUNT = 4000000;

	private ArrayList<UnmappedCenterAtoms> mAtomClasses;
	private StereoMolecule mReactant,mProduct;
	private int[] mReactantMapNo,mProductMapNo;
	private int mStartMapNo,mMapNo;

	public ReactionCenterMapper(StereoMolecule reactant, StereoMolecule product, int[] reactantMapNo, int[] productMapNo, int mapNo) {
		mReactant = reactant;
		mProduct = product;
		mReactantMapNo = reactantMapNo;
		mProductMapNo = productMapNo;
		mStartMapNo = mapNo;
		mMapNo = mapNo;

		// For every atomicNo/isotop create an UnmappedCenterAtoms class with respective reactant and product atoms
		TreeMap<Integer, UnmappedCenterAtoms> atomClassMap = new TreeMap<>();
		for (int atom=0; atom<reactant.getAtoms(); atom++) {
			if (reactantMapNo[atom] == 0) {
				int atomClass = reactant.getAtomicNo(atom) + (reactant.getAtomMass(atom) << 16);
				UnmappedCenterAtoms uca = atomClassMap.get(atomClass);
				if (uca == null) {
					uca = new UnmappedCenterAtoms();
					atomClassMap.put(atomClass, uca);
					}
				uca.addReactantAtom(atom);
				}
			}
		for (int atom=0; atom<product.getAtoms(); atom++) {
			if (productMapNo[atom] == 0) {
				int atomClass = product.getAtomicNo(atom) + (product.getAtomMass(atom) << 16);
				UnmappedCenterAtoms uca = atomClassMap.get(atomClass);
				if (uca == null) {
					uca = new UnmappedCenterAtoms();
					atomClassMap.put(atomClass, uca);
					}
				uca.addProductAtom(atom);
				}
			}

		mAtomClasses = new ArrayList<>();
		for (UnmappedCenterAtoms uca:atomClassMap.values())
			if (!uca.mapObviousAtoms())
				mAtomClasses.add(uca);
		}

	/** Tries and scores all possible mapping permutations for all of hitherto unmapped atoms.
	 * The best scoring combination is kept and its score returned.
	 * If there are no unmapped atoms, then the score of the current mapping is returned.
	 */
	public float completeAndScoreMapping() {
		// For efficient scoring we build a reactionAtomToProductAtom map,
		// which is updated with the center atom assignments for every scoring.
		MappingScorer scorer = new MappingScorer(mReactant, mProduct);
		int[] reactantToProductAtom = scorer.createReactantToProductAtomMap(mReactantMapNo, mProductMapNo);

		if (mAtomClasses.size() == 0)
			return scorer.scoreMapping(reactantToProductAtom);

		double totalPermutationCount = 1;
		for (UnmappedCenterAtoms uca:mAtomClasses)
			totalPermutationCount *= uca.getPermutationCount();
		if (totalPermutationCount > MAX_PERMUTATION_COUNT) {
			System.out.println("permutationCount exceeds maximum:"+totalPermutationCount);
			return 0;
			}

		int atomCount = 0;
		int[] cumulatedAtomCount = new int[mAtomClasses.size()];
		int[] permutationCount = new int[mAtomClasses.size()];
		for (int i=0; i<mAtomClasses.size(); i++) {
			UnmappedCenterAtoms uca = mAtomClasses.get(i);
			permutationCount[i] = uca.initializePermutations();
			cumulatedAtomCount[i] = atomCount;
			atomCount += uca.getMappableAtomCount();
			}

		float bestScore = -1e10f;
		int[] bestReactantAtom = null;
		int[] bestProductAtom = null;
		int[] permutationIndex = new int[mAtomClasses.size()];
		boolean nextPermutationAvailable = (mAtomClasses.size() != 0);
		while (nextPermutationAvailable) {
			for (int i=0; i<mAtomClasses.size(); i++)
				mAtomClasses.get(i).completeReactantToProductAtomMap(permutationIndex[i], reactantToProductAtom);

			float score = scorer.scoreMapping(reactantToProductAtom);
//System.out.print("score:"+score);
//for (int i=0; i<mAtomClasses.size(); i++)
// for (int j=0; j<mAtomClasses.get(i).mReactantAtom.length; j++)
//  System.out.print(" "+mAtomClasses.get(i).mReactantAtom[j]+":"+reactantToProductAtom[mAtomClasses.get(i).mReactantAtom[j]]);
//System.out.println();

			if (bestScore < score) {
				bestScore = score;

				bestReactantAtom = new int[atomCount];
				bestProductAtom = new int[atomCount];
				int atomOffset = 0;
				for (int i=0; i<mAtomClasses.size(); i++) {
					UnmappedCenterAtoms uca = mAtomClasses.get(i);
					uca.getReactantAtoms(permutationIndex[i], bestReactantAtom, atomOffset);
					uca.getProductAtoms(permutationIndex[i], bestProductAtom, atomOffset);
					atomOffset += uca.mMappableAtomCount;
					}
				}

			nextPermutationAvailable = false;
			for (int i=0; i<permutationIndex.length; i++) {
				permutationIndex[i]++;
				if (permutationIndex[i] < permutationCount[i]) {
					nextPermutationAvailable = true;
					break;
					}
				permutationIndex[i] = 0;
				}
			}

		if (bestScore != -1e10) {
			int mapNo = mMapNo;
			for (int i=0; i<atomCount; i++) {
				mapNo++;
				mReactantMapNo[bestReactantAtom[i]] = mapNo;
				mProductMapNo[bestProductAtom[i]] = mapNo;
				}
			}

		return bestScore;
		}

	public int getMappedAtomCount() {
		return mMapNo - mStartMapNo;
		}


	class UnmappedCenterAtoms {
		private int[] mReactantAtom = new int[0];
		private int[] mProductAtom = new int[0];
		private int mMappableAtomCount = 0;
		private ArrayList<int[]> mPermutationList;

		public void addReactantAtom(int atom) {
			mReactantAtom = addAtom(atom, mReactantAtom);
			if (mReactantAtom.length <= mProductAtom.length)
				mMappableAtomCount = mReactantAtom.length;
			}

		public void addProductAtom(int atom) {
			mProductAtom = addAtom(atom, mProductAtom);
			if (mProductAtom.length <= mReactantAtom.length)
				mMappableAtomCount = mProductAtom.length;
			}

		/**
		 * If we have only one atom of this kind, or if we have all equal un-bonded atoms
		 * on one reaction side, then we can safely map these atoms.
		 * @return true if all mappable atoms were mapped
		 */
		public boolean mapObviousAtoms() {
			if (mMappableAtomCount == 0)
				return true;

			if (mReactantAtom.length == 1 && mProductAtom.length == 1) {
				mMapNo++;
				mReactantMapNo[mReactantAtom[0]] = mMapNo;
				mProductMapNo[mProductAtom[0]] = mMapNo;
				return true;
				}

			// to qualify as equal in the context of already mapped atoms,
			// atoms must not only be symmetrical, they must also not have any neighbours.
			boolean reactantAtomsAreEqual = areEqualSingleAtoms(mReactantAtom, mReactant);
			if (reactantAtomsAreEqual && mReactantAtom.length <= mProductAtom.length) {
				for (int i=0; i<mReactantAtom.length; i++) {
					mMapNo++;
					mReactantMapNo[mReactantAtom[i]] = mMapNo;
					mProductMapNo[mProductAtom[i]] = mMapNo;
					}
				return true;
				}
			boolean productAtomsAreEqual = areEqualSingleAtoms(mProductAtom, mProduct);
			if (productAtomsAreEqual && mReactantAtom.length >= mProductAtom.length) {
				for (int i=0; i<mProductAtom.length; i++) {
					mMapNo++;
					mReactantMapNo[mReactantAtom[i]] = mMapNo;
					mProductMapNo[mProductAtom[i]] = mMapNo;
					}
				return true;
				}

			// theoretically: if we have one bonded atom on one side and multiple equivalent unbonded
			// on the other side, then we could also map the one to any of the other. This is probably a rare case.
			return false;
			}

		public double getPermutationCount() {
			int totalAtomCount = Math.max(mReactantAtom.length, mProductAtom.length);
			double permutationCount = 1;
			for (int i=totalAtomCount-mMappableAtomCount+1; i<=totalAtomCount; i++)
				permutationCount *= i;
			return permutationCount;
			}

		public int getMappableAtomCount() {
			return mMappableAtomCount;
			}

		private boolean areEqualSingleAtoms(int[] atoms, StereoMolecule mol) {
			for (int atom:atoms)
				if (mol.getConnAtoms(atom) != 0)
					return false;

			int charge = mol.getAtomCharge(atoms[0]);
			for (int i=1; i<atoms.length; i++)
				if (mol.getAtomCharge(atoms[i]) != charge)
					return false;

			return true;
			}

		public int initializePermutations() {
			mPermutationList = new ArrayList<>();    // contains pointer array from reactant to product
			int[] solution = new int[mMappableAtomCount];
			boolean[] isUsed = new boolean[Math.max(mReactantAtom.length, mProductAtom.length)];
			permute(0, isUsed, solution);
			return mPermutationList.size();
			}

		private void permute(int index, boolean[] isUsed, int[] solution) {
			for (int i=0; i<isUsed.length; i++) {
				if (!isUsed[i]) {
					isUsed[i] = true;
					solution[index] = i;
					if (index+1 == solution.length)
						mPermutationList.add(solution.clone());
					else
						permute(index+1, isUsed, solution);
					isUsed[i] = false;
					}
				}
			}

		// For the given permutation and all known reactant atoms this method
		// writes -1 or the proper product atom into the mapping array
		public void completeReactantToProductAtomMap(int permutationIndex, int[] reactantToProductAtomMap) {
			int[] permutation = mPermutationList.get(permutationIndex);
			if (mReactantAtom.length <= mProductAtom.length) {
				for (int i=0; i<mMappableAtomCount; i++)
					reactantToProductAtomMap[mReactantAtom[i]] = mProductAtom[permutation[i]];
				}
			else {
				for (int atom:mReactantAtom)
					reactantToProductAtomMap[atom] = -1;
				for (int i=0; i<mMappableAtomCount; i++)
					reactantToProductAtomMap[mReactantAtom[permutation[i]]] = mProductAtom[i];
				}
			}

		public void getReactantAtoms(int permutationIndex, int[] reactantAtom, int reactantAtomOffset) {
			if (mReactantAtom.length <= mProductAtom.length) {
				for (int i=0; i<mMappableAtomCount; i++)
					reactantAtom[reactantAtomOffset + i] = mReactantAtom[i];
				}
			else {
				int[] permutation = mPermutationList.get(permutationIndex);
				for (int i=0; i<mMappableAtomCount; i++)
					reactantAtom[reactantAtomOffset+i] = mReactantAtom[permutation[i]];
				}
			}

		public void getProductAtoms(int permutationIndex, int[] productAtom, int productAtomOffset) {
			if (mReactantAtom.length > mProductAtom.length) {
				for (int i=0; i<mMappableAtomCount; i++)
					productAtom[productAtomOffset + i] = mProductAtom[i];
				}
			else {
				int[] permutation = mPermutationList.get(permutationIndex);
				for (int i=0; i<mMappableAtomCount; i++)
					productAtom[productAtomOffset + i] = mProductAtom[permutation[i]];
				}
			}

		private int[] addAtom(int atom, int[] atoms) {
			int[] newAtoms = new int[atoms.length + 1];
			for (int i=0; i<atoms.length; i++)
				newAtoms[i] = atoms[i];
			newAtoms[atoms.length] = atom;
			return newAtoms;
			}
		}
	}
