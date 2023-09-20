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

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.CanonizerBaseValue;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.ByteArrayComparator;
import com.actelion.research.util.IntArrayComparator;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;
import java.util.TreeSet;

public class RootAtomPairSource {
	private static final int MAX_PAIR_SEQUENCES = 64;
	private static final int MAX_ENVIRONMENT_RADIUS = 7;
	private static final int MIN_ENVIRONMENT_RADIUS = 2;
	private static final int PSEUDO_MAP_NO_SKIPPED_PAIR = -99999;

	private ArrayList<RootAtomPair> mPairBuffer;
	private StereoMolecule mReactant,mProduct;
	private Canonizer mReactantCanonizer,mProductCanonizer;
	private CanonizerBaseValue[] mCanBase;
	private int mSequenceCount,mCurrentRadius,mManualMapCount,mMappableAtomCount,mCurrentEnvIndex0,mCurrentEnvIndex1,mCurrentEnvIndex2,mCurrentEnvIndex3;
	private RootAtomPairDecisionHelper mDecisionHelper;
	private boolean mIsStoichiometric;
	private int[] mReactantRank,mProductRank,mReactantFragmentNo,mProductFragmentNo,mReactantMapNo,mProductMapNo;
	private int mAtomBits,mMaxConnAtoms,mHighestReactionRank,mHighestProductRank;
	private boolean[] mReactantFragmentUsed,mProductFragmentUsed;
	private TreeMap<byte[],int[][]>[] mEnvToAtomsMap; // index: radius
	private byte[][][] mEnvKey;

	public RootAtomPairSource(StereoMolecule reactant, StereoMolecule product, int[] reactantMapNo, int[] productMapNo) {
		mReactant = reactant;
		mProduct = product;

		mReactantMapNo = reactantMapNo;
		mProductMapNo = productMapNo;

		// Reactant and product atom ranks are updated with every assigned pair in order to reflect
		// decreasing symmetry due to unsymmetrically assigned atom mapping.
		initializeAtomRanking();
		initializeRanks();

		mEnvToAtomsMap = buildEnvToAtomMaps();
		mEnvKey = getEnvKeys(mEnvToAtomsMap);

		mMappableAtomCount = 0;
		for (byte[] envKey:mEnvToAtomsMap[0].keySet()) {
			int[][] atoms = mEnvToAtomsMap[0].get(envKey);
			mMappableAtomCount += Math.min(atoms[0].length, atoms[1].length);
			}

		mIsStoichiometric = mMappableAtomCount == mReactant.getAtoms()
						 && mReactant.getAtoms() == mProduct.getAtoms();

		mReactantFragmentNo = new int[mReactant.getAtoms()];
		mProductFragmentNo = new int[mProduct.getAtoms()];
		mReactantFragmentUsed = new boolean[mReactant.getFragmentNumbers(mReactantFragmentNo, false, false)];
		mProductFragmentUsed = new boolean[mProduct.getFragmentNumbers(mProductFragmentNo, false, false)];

		mManualMapCount = assignManuallyMappedPairs();
		mCurrentRadius = MAX_ENVIRONMENT_RADIUS;
		mSequenceCount = 0;
		}

	private void initializeAtomRanking() {
		int maxAtomCount = Math.max(mReactant.getAtoms(), mProduct.getAtoms());
		mAtomBits = Canonizer.getNeededBits(maxAtomCount);

		mMaxConnAtoms = 2;
		for (int atom=0; atom<mReactant.getAtoms(); atom++)
			mMaxConnAtoms = Math.max(mMaxConnAtoms, mReactant.getConnAtoms(atom)+mReactant.getMetalBondedConnAtoms(atom));
		for (int atom=0; atom<mProduct.getAtoms(); atom++)
			mMaxConnAtoms = Math.max(mMaxConnAtoms, mProduct.getConnAtoms(atom)+mProduct.getMetalBondedConnAtoms(atom));
		int baseValueSize = Math.max(2, (62 + mAtomBits + mMaxConnAtoms * (mAtomBits+5)) / 63);

		mCanBase = new CanonizerBaseValue[maxAtomCount];
		for (int atom=0; atom<maxAtomCount; atom++)
			mCanBase[atom] = new CanonizerBaseValue(baseValueSize);

		mReactantCanonizer = new Canonizer(mReactant, Canonizer.CREATE_SYMMETRY_RANK);
		mProductCanonizer = new Canonizer(mProduct, Canonizer.CREATE_SYMMETRY_RANK);
		}

	private void initializeRanks() {
		mReactantRank = mReactantCanonizer.getSymmetryRanks().clone();
		mProductRank = mProductCanonizer.getSymmetryRanks().clone();
		for (int rank:mReactantRank)
			if (mHighestReactionRank < rank)
				mHighestReactionRank = rank;
		for (int rank:mProductRank)
			if (mHighestProductRank < rank)
				mHighestProductRank = rank;
		}

	private void reset() {
		Arrays.fill(mReactantMapNo, 0);
		Arrays.fill(mProductMapNo, 0);
		Arrays.fill(mReactantFragmentUsed, false);
		Arrays.fill(mProductFragmentUsed, false);

		initializeRanks();

		mDecisionHelper.reset();

		mManualMapCount = assignManuallyMappedPairs();

		mCurrentRadius = MAX_ENVIRONMENT_RADIUS;
		mCurrentEnvIndex0 = 0;
		mCurrentEnvIndex1 = 0;
		mCurrentEnvIndex2 = 0;
		mCurrentEnvIndex3 = 0;
		}

	public int getMappableAtomCount() {
		return mMappableAtomCount;
		}

	public int getManualMapCount() {
		return mManualMapCount;
		}

	public boolean hasNextPairSequence() {
		if (mSequenceCount++ == MAX_PAIR_SEQUENCES)
			return false;

		if (mDecisionHelper == null) {
			mDecisionHelper = new RootAtomPairDecisionHelper();
			return true;
			}

		if (mDecisionHelper.isComplete())
			return false;

		reset();

		return true;
		}

	/**
	 * RootAtomPairs are returned in similarity order. The first returned pair is that
	 * pair of atoms that is more similar than any other mutual combination of reactant
	 * and product atoms. When, however, multiple pairs are equivalent, then the choice
	 * is arbitrary. This class makes sure that symmetrically redundant choices are removed
	 * and chooses one of multiple remaining choices in a way that is different from the
	 * previous
	 * @return
	 */
	private byte[][][] classifyAtomEnvironment(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());

		byte[][][] environment = new byte[mol.getAtoms()][MAX_ENVIRONMENT_RADIUS+1][];

		int[] atomList = new int[mol.getAtoms()];
		boolean[] atomMask = new boolean[mol.getAtoms()];
		for (int rootAtom=0; rootAtom<mol.getAtoms(); rootAtom++) {
			if (rootAtom != 0)
				Arrays.fill(atomMask, false);

			int min = 0;
			int max = 0;

			// we need to mark the root atom, because otherwise close-by root atoms may end up with the same fragment
			mol.setAtomSelection(rootAtom, true);

			for (int sphere = 0; sphere<=MAX_ENVIRONMENT_RADIUS && max<mol.getAtoms(); sphere++) {
				if (max == 0) {
					atomList[0] = rootAtom;
					atomMask[rootAtom] = true;
					max = 1;
					}
				else {
					int newMax = max;
					for (int i=min; i<max; i++) {
						int atom = atomList[i];
						for (int j=0; j<mol.getConnAtoms(atom); j++) {
							int connAtom = mol.getConnAtom(atom, j);
							if (!atomMask[connAtom]) {
								atomMask[connAtom] = true;
								atomList[newMax++] = connAtom;
								}
							}
						}

					if (newMax == max)
						break;

					min = max;
					max = newMax;
					}

				if (sphere == 0) {
					environment[rootAtom][sphere] = new byte[2];
					environment[rootAtom][sphere][0] = (byte)mol.getAtomicNo(rootAtom);
					environment[rootAtom][sphere][1] = (byte)mol.getAtomMass(rootAtom);
					}
				else {
					mol.copyMoleculeByAtoms(fragment, atomMask, true, null);
					for (int atom=0; atom<fragment.getAllAtoms(); atom++) {
						fragment.setAtomCharge(atom, 0);
						fragment.setAtomRadical(atom, 0);
						}
					environment[rootAtom][sphere] = new Canonizer(fragment, Canonizer.ENCODE_ATOM_SELECTION).getIDCode().getBytes(StandardCharsets.UTF_8);
					}
				}

			mol.setAtomSelection(rootAtom, false);
			}
		return environment;
		}

	/**
	 * Builds for every radius a map environments descriptors for which we have atoms on both reaction sides.
	 * @return
	 */
	private TreeMap<byte[], int[][]>[] buildEnvToAtomMaps() {
		byte[][][] reactantAtomEnv = classifyAtomEnvironment(mReactant);
		byte[][][] productAtomEnv = classifyAtomEnvironment(mProduct);
		TreeMap<byte[], int[][]>[] envMap = new TreeMap[MAX_ENVIRONMENT_RADIUS+1];
		TreeMap<byte[], int[]>[] reactantMap = buildEnvToAtomMaps(mReactant, reactantAtomEnv);
		TreeMap<byte[], int[]>[] productMap = buildEnvToAtomMaps(mProduct, productAtomEnv);
		for (int radius=0; radius<=MAX_ENVIRONMENT_RADIUS; radius++) {
			envMap[radius] = new TreeMap<>(new ByteArrayComparator());
			for (byte[] envKey:reactantMap[radius].keySet()) {
				int[] reactantAtoms = reactantMap[radius].get(envKey);
				int[] productAtoms = productMap[radius].get(envKey);
				if (productAtoms != null) {
					int[][] atoms = new int[2][];
					atoms[0] = reactantAtoms;
					atoms[1] = productAtoms;
					envMap[radius].put(envKey, atoms);
					}
				}
			}
		return envMap;
		}

	private TreeMap<byte[], int[]>[] buildEnvToAtomMaps(StereoMolecule mol, byte[][][] atomEnv) {
		TreeMap<byte[], int[]>[] map = new TreeMap[MAX_ENVIRONMENT_RADIUS+1];
		for (int radius=0; radius<=MAX_ENVIRONMENT_RADIUS; radius++) {
			map[radius] = new TreeMap<>(new ByteArrayComparator());
			for (int atom=0; atom<mol.getAtoms(); atom++) {
				byte[] env = atomEnv[atom][radius];
				if (env != null) {
					int[] atoms = map[radius].get(env);
					atoms = (atoms == null) ? new int[1] : Arrays.copyOf(atoms, atoms.length+1);
					atoms[atoms.length-1] = atom;
					map[radius].put(env, atoms);
					}
				}
			}
		return map;
		}

	private byte[][][] getEnvKeys(TreeMap<byte[], int[][]>[] envMap) {
		byte[][][] keys = new byte[MAX_ENVIRONMENT_RADIUS+1][][];
		for (int radius=0; radius<=MAX_ENVIRONMENT_RADIUS; radius++) {
			keys[radius] = new byte[envMap[radius].size()][];
			int index = 0;
			for (byte[] key:envMap[radius].keySet())
				keys[radius][index++] = key;
			}
		return keys;
		}

	/**
	 * Assigns first mapping numbers to manually mapped atoms and add the pair list to the buffer.
	 * @return list of manually mapped atom pairs
	 */
	private int assignManuallyMappedPairs() {
		mPairBuffer = new ArrayList<>();
		int mapNo = 1;
		int maxMapNo = 0;
		for (int atom=0; atom<mProduct.getAtoms(); atom++)
			if (mProduct.getAtomMapNo(atom) != 0 && !mProduct.isAutoMappedAtom(atom))
				maxMapNo = Math.max(maxMapNo, mProduct.getAtomMapNo(atom));
		if (maxMapNo != 0) {
			int[] mapNoToProductAtom = new int[maxMapNo+1];
			for (int atom=0; atom<mProduct.getAtoms(); atom++)
				if (mProduct.getAtomMapNo(atom) != 0 && !mProduct.isAutoMappedAtom(atom))
					mapNoToProductAtom[mProduct.getAtomMapNo(atom)] = atom+1;
			for (int reactantAtom=0; reactantAtom<mReactant.getAtoms(); reactantAtom++) {
				int reactantMapNo = mReactant.getAtomMapNo(reactantAtom);
				if (reactantMapNo != 0
						&& reactantMapNo <= maxMapNo
						&& !mReactant.isAutoMappedAtom(reactantAtom)
						&& mapNoToProductAtom[reactantMapNo] != 0) {
					int productAtom = mapNoToProductAtom[reactantMapNo]-1;
					mReactantMapNo[reactantAtom] = mapNo;
					mProductMapNo[productAtom] = mapNo++;
					mPairBuffer.add(new RootAtomPair(reactantAtom, productAtom));
					}
				}
			}
		return mPairBuffer.size();
		}

	/**
	 * This creates and returns the next reactant/product atom pair to serve as starting point
	 * for building equivalent atom graphs that logically match reactant subgraphs to product subgraphs.
	 * The most plausible atom pair regarding reactant to product atom similarity is returned first.
	 * If multiple equi-plausible pairs exist, the permutation value passed to the constructor determines,
	 * which of the pairs is returned.
	 * Root atom pairs meet these conditions:<br>
	 * - they match in terms of circular fragment on reactant and products side<br>
	 * - if multiple symmetrically equivalent pairs exist, exactly one of them is marked as allowed root pair<br>
	 * - each of the reactant and product atoms have at least one unmapped neighbour<br>
	 * @return pair of currently obvious root atom pairs
	 */
	public RootAtomPair nextPair() {
		RootAtomPair pair = nextRawPair();
		while (pair != null) {
			boolean reactantAtomHasUnmappedNeighbours = false;
			for (int i=0; i<mReactant.getConnAtoms(pair.reactantAtom); i++)
				if (mReactantMapNo[mReactant.getConnAtom(pair.reactantAtom, i)] == 0)
					reactantAtomHasUnmappedNeighbours = true;

			boolean productAtomHasUnmappedNeighbours = false;
			for (int i=0; i<mProduct.getConnAtoms(pair.productAtom); i++)
				if (mProductMapNo[mProduct.getConnAtom(pair.productAtom, i)] == 0)
					productAtomHasUnmappedNeighbours = true;

			if (reactantAtomHasUnmappedNeighbours && productAtomHasUnmappedNeighbours)
				break;

			// we need to mark refused pairs as being used to avaid getting them again
			mReactantMapNo[pair.reactantAtom] = PSEUDO_MAP_NO_SKIPPED_PAIR;
			mProductMapNo[pair.productAtom] = PSEUDO_MAP_NO_SKIPPED_PAIR;

			pair = nextRawPair();
			}

		if (pair == null) { // remove pseudo map numbers once not needed anymore
			for (int i=0; i<mReactantMapNo.length; i++)
				if (mReactantMapNo[i] == PSEUDO_MAP_NO_SKIPPED_PAIR)
					mReactantMapNo[i] = 0;
			for (int i=0; i<mProductMapNo.length; i++)
				if (mProductMapNo[i] == PSEUDO_MAP_NO_SKIPPED_PAIR)
					mProductMapNo[i] = 0;
			}

		return pair;
		}

	private RootAtomPair nextRawPair() {
		while (mPairBuffer.size() != 0) {
			RootAtomPair pair = mPairBuffer.remove(0);
			if (mReactantMapNo[pair.reactantAtom] == 0 && mProductMapNo[pair.productAtom] == 0)
				return pair;
			}

		while (mCurrentRadius >= 0) {
			// We create starting pairs of reasonably similar atoms with similarity derived priority
			while (mCurrentRadius >= MIN_ENVIRONMENT_RADIUS
				&& mCurrentEnvIndex0 < mEnvKey[mCurrentRadius].length) {
				byte[] envKey = mEnvKey[mCurrentRadius][mCurrentEnvIndex0];
				int[][] atoms = mEnvToAtomsMap[mCurrentRadius].get(envKey);
				if (mReactant.getAtomicNo(atoms[0][0]) == 6) {
					RootAtomPair pair = makePairsFromSimilarAtoms(atoms[0], atoms[1]);
					if (pair != null)
						return pair;
					}
				mCurrentEnvIndex0++;    // go to next environment key once all potential pairs with current environment are depleted
				}

			while (mCurrentRadius >= MIN_ENVIRONMENT_RADIUS
			    && mCurrentEnvIndex1 < mEnvKey[mCurrentRadius].length) {
				byte[] envKey = mEnvKey[mCurrentRadius][mCurrentEnvIndex1];
				int[][] atoms = mEnvToAtomsMap[mCurrentRadius].get(envKey);
				if (mReactant.getAtomicNo(atoms[0][0]) != 6) {
					// with equal environment size, we prefer carbon root atoms TODO
					RootAtomPair pair = makePairsFromSimilarAtoms(atoms[0], atoms[1]);
					if (pair != null)
						return pair;
					}
				mCurrentEnvIndex1++;    // go to next environment key once all potential pairs with current environment are depleted
				}

			// We create a low priority starting pair if we just have one atom of a kind
			while (mIsStoichiometric && mCurrentRadius == 0
			 && mCurrentEnvIndex2 < mEnvKey[0].length) {
				byte[] envKey = mEnvKey[0][mCurrentEnvIndex2++];
				int[][] atoms = mEnvToAtomsMap[mCurrentRadius].get(envKey);
				if (atoms[0].length == 1 && atoms[1].length == 1) {
					RootAtomPair pair = tryCreatePair(atoms[0][0], atoms[1][0]);
					if (pair != null)
						return pair;
					}
				}

			// Check with decreasing sphere size down to atomicNo level, whether all reactant or product atoms
			// are in separate fragments and are equivalent. If this is the case, and if we have at least
			// as many equivalent atoms as atoms on the other side, then we just match reactant atoms
			// to product atoms in order of appearance.
			while (mCurrentEnvIndex3 < mEnvKey[mCurrentRadius].length) {
				byte[] envKey = mEnvKey[mCurrentRadius][mCurrentEnvIndex3++];
				int[][] atoms = mEnvToAtomsMap[mCurrentRadius].get(envKey);
				if (atoms[0].length == 1 && atoms[1].length == 1) {
					RootAtomPair pair = tryCreatePair(atoms[0][0], atoms[1][0]);
					if (pair != null)
						return pair;
					}
				else if ((atoms[0].length >= atoms[1].length
						&& areInDistinctEquivalentFragments(atoms[0], mReactantFragmentNo, mReactantFragmentUsed, mReactantRank))
						|| (atoms[1].length >= atoms[0].length
						&& areInDistinctEquivalentFragments(atoms[1], mProductFragmentNo, mProductFragmentUsed, mProductRank))) {
					for (int i=0; i<Math.min(atoms[0].length, atoms[1].length); i++) {
						RootAtomPair pair = tryCreatePair(atoms[0][i], atoms[1][i]);
						if (pair != null)
							return pair;
						}
					}
				}

			mCurrentRadius--;
			mCurrentEnvIndex0 = 0;
			mCurrentEnvIndex1 = 0;
			mCurrentEnvIndex2 = 0;
			mCurrentEnvIndex3 = 0;
			}

		return null;
		}

	/**
	 * Matches multiple reactant atoms to multiple product atoms, all being similar regarding
	 * their neighbourhood. Within all reactant atoms (same applies to product atoms) there
	 * may be subgroups with equal similarity ranks. All matching rank permutations are built,
	 * and one of them is chosen considering the history of choices.
	 * If an unsymmetrical matching was chosen, then ranks are consolidated such that formerly
	 * equal ranking atoms are split where atoms of these got matched to differently ranking atoms.
	 * @param reactantAtoms
	 * @param productAtoms
	 * @return
	 */
	private RootAtomPair makePairsFromSimilarAtoms(int[] reactantAtoms, int[] productAtoms) {
		int[] unusedReactantAtom = getUnmappedAtoms(reactantAtoms, mReactantMapNo);
		if (unusedReactantAtom == null)
			return null;
		int[] unusedProductAtom = getUnmappedAtoms(productAtoms, mProductMapNo);
		if (unusedProductAtom == null)
			return null;

		if (unusedReactantAtom.length == 1 && unusedProductAtom.length == 1)
			return new RootAtomPair(unusedReactantAtom[0], unusedProductAtom[0]);

		int[][] rankPairList = createDistinctRankPairs(unusedReactantAtom, unusedProductAtom);
		if (rankPairList.length == 1)
			return new RootAtomPair(unusedReactantAtom[0], unusedProductAtom[0]);

		int choice = mDecisionHelper.getNextChoice(rankPairList.length);

		int chosenReactantAtom = -1;
		for (int atom:unusedReactantAtom) {
			if (mReactantRank[atom] == rankPairList[choice][0]) {
				chosenReactantAtom = atom;
				break;
				}
			}
		int chosenProductAtom = -1;
		for (int atom:unusedProductAtom) {
			if (mProductRank[atom] == rankPairList[choice][1]) {
				chosenProductAtom = atom;
				break;
				}
			}

		mHighestReactionRank = elevateAtomRank(mReactant, mReactantRank, chosenReactantAtom, mHighestReactionRank);
		mHighestProductRank = elevateAtomRank(mProduct, mProductRank, chosenProductAtom, mHighestProductRank);

		return new RootAtomPair(chosenReactantAtom, chosenProductAtom);
		}

	private int[] getUnmappedAtoms(int[] atoms, int[] mapNo) {
		int unmappedAtomCount = 0;
		for (int atom:atoms)
			if (mapNo[atom] == 0)
				unmappedAtomCount++;
		if (unmappedAtomCount == 0)
			return null;
		int[] unmappedAtom = new int[unmappedAtomCount];
		unmappedAtomCount = 0;
		for (int atom:atoms)
			if (mapNo[atom] == 0)
				unmappedAtom[unmappedAtomCount++] = atom;
		return unmappedAtom;
		}

	private RootAtomPair tryCreatePair(int reactantAtom, int productAtom) {
		if (mReactantMapNo[reactantAtom] == 0 && mProductMapNo[productAtom] == 0)
			return createPair(reactantAtom, productAtom);

		return null;
		}

	private RootAtomPair createPair(int reactantAtom, int productAtom) {
		if (mProduct.getConnAtoms(productAtom) != 0)
			mReactantFragmentUsed[mReactantFragmentNo[reactantAtom]] = true;
		if (mReactant.getConnAtoms(reactantAtom) != 0)
			mProductFragmentUsed[mProductFragmentNo[productAtom]] = true;
		return new RootAtomPair(reactantAtom, productAtom);
		}

	private int[][] createDistinctRankPairs(int[] reactantAtom, int[] productAtom) {
		TreeSet<int[]> rankPairSet = new TreeSet<>(new IntArrayComparator());
		for (int ra:reactantAtom) {
			for (int pa:productAtom) {
				int[] pair = new int[2];
				pair[0] = mReactantRank[ra];
				pair[1] = mProductRank[pa];
				rankPairSet.add(pair);
				}
			}
		return rankPairSet.toArray(new int[0][]);
		}

	/**
	 * Increases the symmetry rank of the specified atom in regard to all other atoms with
	 * the same rank. This split of formerly equal ranks accounts for the asymmetry introduced
	 * by matching this atom to another atom on the other side of the reaction.
	 * The updated rank is then propagated recursively through all neighbours.
	 * This, of course, happens only, if the atom's initial rank is shared by at least another atom.
	 * @param mol
	 * @param atomRank reactant's or product's atom ranks to be modified in place
	 * @param atom
	 * @return
	 */
	private int elevateAtomRank(StereoMolecule mol, int[] atomRank, int atom, int maxRank) {
		int currentAtomRank = atomRank[atom];

		boolean sharedRankFound = false;
		for (int i=0; i<mol.getAtoms(); i++) {
			if (i != atom && atomRank[i] == currentAtomRank) {
				sharedRankFound = true;
				break;
				}
			}
		if (!sharedRankFound)
			return maxRank;

		for (int i=0; i<mol.getAtoms(); i++)
			if (i == atom || atomRank[i] > currentAtomRank)
				atomRank[i]++;

		int oldMaxRank;
		maxRank++;
		do {
			oldMaxRank = maxRank;
			canCalcNextBaseValues(mol, atomRank);
			maxRank = consolidateRanking(atomRank);
			} while (oldMaxRank != maxRank);

		return maxRank;
		}


	private void canCalcNextBaseValues(StereoMolecule mol, int[] atomRank) {
		int[] connRank = new int[mMaxConnAtoms];
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			// generate sorted list of ranks of neighbours
			int neighbours = mol.getConnAtoms(atom)+mol.getMetalBondedConnAtoms(atom);
			int neighbour = 0;
			for (int i=0; i<mol.getAllConnAtomsPlusMetalBonds(atom); i++) {
				if (i<mol.getConnAtoms(atom) || i>=mol.getAllConnAtoms(atom)) {
					int rank = 2 * atomRank[mol.getConnAtom(atom, i)];
					int connBond = mol.getConnBond(atom, i);
					if (mol.getBondOrder(connBond) == 2)
						if (!mol.isAromaticBond(connBond))
							rank++;        // set a flag for non-aromatic double bond
					int j;
					for (j = 0; j < neighbour; j++)
						if (rank < connRank[j])
							break;
					for (int k = neighbour; k > j; k--)
						connRank[k] = connRank[k - 1];
					connRank[j] = rank;
					neighbour++;
					}
				}

			mCanBase[atom].init(atom);
			mCanBase[atom].add(mAtomBits, atomRank[atom]);
			for (int i=neighbours; i<mMaxConnAtoms; i++)
				mCanBase[atom].add(mAtomBits + 1, 0);
			for (int i=0; i<neighbours; i++)
				mCanBase[atom].add(mAtomBits + 1, connRank[i]);
			}

		// for that reaction side with surplus base values, we fill them to be sorted to the end of the table
		for (int atom=mol.getAtoms(); atom<mCanBase.length; atom++) {
			mCanBase[atom].init(atom);
			mCanBase[atom].add(mAtomBits, mol.getAtoms()+1);
			for (int i=0; i<mMaxConnAtoms; i++)
				mCanBase[atom].add(mAtomBits + 1, 0);
			}
		}

	private int consolidateRanking(int[] atomRank) {
		int canRank = 0;
		Arrays.sort(mCanBase);
		for (int i=0; i<atomRank.length; i++) {
			if (i == 0 || mCanBase[i].compareTo(mCanBase[i-1]) != 0)
				canRank++;
			int index = mCanBase[i].getAtom();
			atomRank[index] = canRank;
			}
		return canRank;
		}

	private boolean areInDistinctEquivalentFragments(int[] atom, int[] fragmentNo, boolean[] fragmentUsed, int[] symmetryRank) {
		for (int a:atom)
			if (fragmentUsed[fragmentNo[a]])
				return false;

		for (int i=1; i<atom.length; i++)
			if (symmetryRank[atom[i]] != symmetryRank[atom[0]])
				return false;

		return true;
		}
	}

class RootAtomPairDecisionHelper {
	private RootAtomPairDecisionNode mRootNode,mCurrentNode;

	public RootAtomPairDecisionHelper() {
		mRootNode = new RootAtomPairDecisionNode(null);
		mCurrentNode = mRootNode;
		}

	public void reset() {
		mCurrentNode = mRootNode;
		}

	public boolean isComplete() {
		return mRootNode.isComplete();
		}

	public int getNextChoice(int choiceCount) {
		int choice = mCurrentNode.getNextChoice(choiceCount);
		mCurrentNode = mCurrentNode.getChild(choice);
		return choice;
		}
	}

class RootAtomPairDecisionNode {
	private int mCurrentChoice,mChoiceCount;
	private RootAtomPairDecisionNode mParent;
	private RootAtomPairDecisionNode[] mChildren;
	private boolean mIsComplete;

	public RootAtomPairDecisionNode(RootAtomPairDecisionNode parent) {
		mParent = parent;
		mChoiceCount = -1;    // no choice count known
		}

	public RootAtomPairDecisionNode getChild(int no) {
		return mChildren[no];
		}

	/**
	 * Gets the next choice of this node. If the node doesn't know the number of choices yet,
	 * it initializes to the first (highest) of possible choices.
	 * If there is no child node for the selected choice, a new uninitialized child node is created.
	 * @param choiceCount
	 * @return
	 */
	public int getNextChoice(int choiceCount) {
		if (mChoiceCount == -1) {
			mChildren = new RootAtomPairDecisionNode[choiceCount];
			mChoiceCount = choiceCount;
			mCurrentChoice = choiceCount-1;
			}
		else {
			do {
				mCurrentChoice = (mCurrentChoice == 0) ? mChoiceCount-1 : mCurrentChoice-1;
				} while (mChildren[mCurrentChoice] != null && mChildren[mCurrentChoice].isComplete());
			}
		if (mChildren[mCurrentChoice] == null)
			mChildren[mCurrentChoice] = new RootAtomPairDecisionNode(this);

		return mCurrentChoice;
		}

	public boolean isComplete() {
		if (!mIsComplete) {
			// if we never had to take a decision
			if (mChoiceCount == -1) {
				mIsComplete = true;
				}
			else {
				boolean nullChildFound = false;
				for (RootAtomPairDecisionNode child:mChildren) {
					if (child == null)
						nullChildFound = true;
					else if (!child.isComplete())
						return false;
					}
				mIsComplete = !nullChildFound;
				}
			}

		return mIsComplete;
		}
	}

class RootAtomPair implements Comparable<RootAtomPair> {
	public int reactantAtom,productAtom;

	public RootAtomPair(int reactantAtom, int productAtom) {
		this.reactantAtom = reactantAtom;
		this.productAtom = productAtom;
		}

	@Override
	public int compareTo(RootAtomPair pair) {
		return this.reactantAtom < pair.reactantAtom ? -1
			 : this.reactantAtom > pair.reactantAtom ? 1
			 : this.productAtom < pair.productAtom ? -1
			 : this.productAtom > pair.productAtom ? 1 : 0;
		}
	}
