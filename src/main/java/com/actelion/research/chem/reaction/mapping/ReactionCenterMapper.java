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
	private static final int MAX_PERMUTATION_COUNT = 20000;

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
			if (uca.mMappableAtomCount != 0)
				if (!uca.mapObviousAtoms())
					mAtomClasses.add(uca);
		}

	/** Tries and scores all possible mapping permutations for all of hitherto unmapped atoms.
	 * The best scoring combination is kept and its score returned.
	 */
	public int completeMapping() {
		int atomCount = 0;
		int mapNoAfterObvious = mMapNo;
		int[] cumulatedAtomCount = new int[mAtomClasses.size()];
		int[] permutationCount = new int[mAtomClasses.size()];
		int totalPermutationCount = 1;
		int atomClassIndex = 0;
		for (UnmappedCenterAtoms uca:mAtomClasses) {
			atomCount += uca.mappableAtomCount();
			permutationCount[atomClassIndex] = uca.initializePermutations();
			totalPermutationCount *= permutationCount[atomClassIndex];
			if (atomClassIndex != 0)
				cumulatedAtomCount[atomClassIndex] = cumulatedAtomCount[atomClassIndex-1] + uca.mappableAtomCount();
			atomClassIndex++;
			}

		if (totalPermutationCount > MAX_PERMUTATION_COUNT) {
			System.out.println("permutationCount:"+totalPermutationCount);
			return 0;
			}

		int atomOffset = 0;
		int[] reactantAtom = new int[atomCount];
		for (UnmappedCenterAtoms uca:mAtomClasses)
			atomOffset += uca.getReactantAtoms(reactantAtom, atomOffset);

		// already assign mapNos to reactant atoms
		int mapNoAfterObviousMapping = mMapNo;
		for (int i=0; i<atomCount; i++)
			mReactantMapNo[reactantAtom[i]] = ++mMapNo;

		int[] productAtom = new int[atomCount];

		int bestScore = 0;
		int[] bestProductAtom = null;
		int[] permutationIndex = new int[mAtomClasses.size()];
		boolean nextPermutationAvailable = (mAtomClasses.size() != 0);
		while (nextPermutationAvailable) {
			int index = 0;
			atomOffset = 0;
			for (UnmappedCenterAtoms uca:mAtomClasses)
				atomOffset += uca.getPermutedProductAtoms(permutationIndex[index++], productAtom, atomOffset);

			nextPermutationAvailable = false;
			for (int i=0; i<permutationIndex.length; i++) {
				permutationIndex[i]++;
				if (permutationIndex[i] < permutationCount[i]) {
					nextPermutationAvailable = true;
					break;
					}
				permutationIndex[i] = 0;
				}

			// assign product mapNos for scoring
			int mapNo = mapNoAfterObviousMapping;
			for (int i=0; i<atomCount; i++)
				mProductMapNo[productAtom[i]] = ++mapNo;

			int score = scoreCenterMapping(reactantAtom, productAtom);
			if (bestScore < score) {
				bestScore = score;
				bestProductAtom = productAtom.clone();
				}
			}

		if (bestScore != 0) {
			int mapNo = mapNoAfterObviousMapping;
			for (int i=0; i<atomCount; i++)
				mProductMapNo[bestProductAtom[i]] = ++mapNo;
			}

		return bestScore;
		}

	public int getMappedAtomCount() {
		return mMapNo - mStartMapNo;
		}

	/**
	 * @param reactantAtom list of reactant atoms to be mapped sorted by atomicNo
	 * @param productAtom list of product atoms (permuted within atomicNo sets)
	 */
	private int scoreCenterMapping(int[] reactantAtom, int[] productAtom) {
		for (int i=0; i<productAtom.length; i++) {
			int pAtom = productAtom[i];
			}

		int score = 1;
		for (int i=0; i<reactantAtom.length; i++) {
			int rAtom = reactantAtom[i];
			int pAtom = productAtom[i];
			for (int j=0; j<mReactant.getConnAtoms(rAtom); j++) {
				int rConnAtom = mReactant.getConnAtom(rAtom, j);
				for (int k=0; k<mProduct.getConnAtoms(pAtom); k++) {
					int pConnAtom = mProduct.getConnAtom(pAtom, k);
					if (mReactantMapNo[rConnAtom] != 0
					 && mReactantMapNo[rConnAtom] == mProductMapNo[pConnAtom]) {
						score++;
						break;
						}
					}
				}
			}

		return score;
		}

	class UnmappedCenterAtoms {
		private static final int MAX_ATOMS_PER_TYPE = 6;

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

		public int mappableAtomCount() {
			return mMappableAtomCount;
			}

		/**
		 * If we have only one atom of this kind, or if we have all equal un-bonded atoms
		 * on one reaction side, then we can safely map these atoms.
		 * @return true if all mappable atoms were mapped
		 */
		public boolean mapObviousAtoms() {
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
			for (int i=0; i<mMappableAtomCount; i++)
				solution[i] = i;

			if (mMappableAtomCount <= MAX_ATOMS_PER_TYPE)
				permute(solution, 0, solution.length-1);
			else
				mPermutationList.add(solution);

			return mPermutationList.size();
			}

		public int getReactantAtoms(int[] reactantAtom, int reactantAtomOffset) {
			for (int i=0; i<mMappableAtomCount; i++)
				reactantAtom[reactantAtomOffset+i] = mReactantAtom[i];

			return mMappableAtomCount;
			}

		public int getPermutedProductAtoms(int permutationIndex, int[] productAtom, int productAtomOffset) {
			int[] permutation = mPermutationList.get(permutationIndex);
			for (int i=0; i<mMappableAtomCount; i++)
				productAtom[productAtomOffset+i] = mProductAtom[permutation[i]];

			return mMappableAtomCount;
			}

		private void permute(int[] solution, int l, int r) {
			if (l == r) {
				mPermutationList.add(solution.clone());
				}
			else {
				for (int i=l; i<=r; i++) {
					int temp = solution[l];
					solution[l] = solution[i];
					solution[i] = temp;

					permute(solution, l+1, r);

					solution[i] = solution[l];
					solution[l] = temp;
				}
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
