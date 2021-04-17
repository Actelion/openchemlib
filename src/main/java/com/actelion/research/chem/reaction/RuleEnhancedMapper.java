package com.actelion.research.chem.reaction;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.ByteArrayComparator;
import com.actelion.research.util.UniqueList;

import java.util.Arrays;
import java.util.TreeMap;

public class RuleEnhancedMapper {
	private static final int MAX_EXVIRONMENT_RADIUS = 8;
	private static final int MIN_ROOT_ATOM_RADIUS = 2;
	private StereoMolecule mReactant,mProduct;
	private int mMapNo,mBestMapNo,mMappableAtomCount;
	private int[] mReactantMapNo,mProductMapNo,mBestReactantMapNo,mBestProductMapNo;
	private Canonizer mReactantCanonizer,mProductCanonizer;
	private ByteArrayComparator mSimilarityComparator;
	private byte[][][] mReactantAtomEnv,mProductAtomEnv;   // indexes: atom,radius,idcode bytes
	private TreeMap<byte[],int[]>[] mEnvToReactantAtomMap,mEnvToProductAtomMap;

	/**
	 *
	 * @param rxn reaction with proper atom coordinates
	 */
	public void map(Reaction rxn) {
		copyReactionToMolecules(rxn);
		initializeAtomEnvironments();

		RootAtomPair[] rootAtomPairs = findRootAtomPairs();

		mReactantMapNo = new int[mReactant.getAtoms()];
		mProductMapNo = new int[mProduct.getAtoms()];
		mMapNo = 0;

		mBestReactantMapNo = new int[mReactant.getAtoms()];
		mBestProductMapNo = new int[mProduct.getAtoms()];
		mBestMapNo = 0;

		for (RootAtomPair root:rootAtomPairs) {
			if (root.isFirstOfEquivalentPairs) {
				int score = 16 * mapFromRootAtoms(root);
				for (RootAtomPair pair:rootAtomPairs)
					if (pair.isNotYetMapped() && pair.isFirstOfEquivalentPairs)
						score += 16 * mapFromRootAtoms(pair);
				for (RootAtomPair pair:rootAtomPairs)
					if (pair.isNotYetMapped())
						score += 16 * mapFromRootAtoms(pair);

				if (mMapNo < mMappableAtomCount)
					score += new RuleEnhancedCenterMapper(mReactant, mProduct, mReactantMapNo, mProductMapNo, mMapNo).completeMapping();

				// for now we are happy with just one similaroty based mapping round
				System.out.println("Mapping round with root("+root.reactantAtom+","+root.productAtom+") completed:"+score);
				break;
				}
			}

		copyMapNosToReaction(rxn, mReactantMapNo, mProductMapNo);
		System.out.println("--------------------------------------------------------------");
		}

	/**
	 * This compiles a complete list of reactant/product atom pairs, which may
	 * serve as starting points, i.e. graph root atoms, for mapping a reactant
	 * subgraph to a product subgraph. Root atom pairs meet these conditions:<br>
	 * - they match in terms of circular fragment on reactant and products side<br>
	 * - if multiple symmetrically equivalent pairs exist, exactly one of them is marked as allowed root pair<br>
	 * @return
	 */
	private RootAtomPair[] findRootAtomPairs() {
		int priority = 1;
		UniqueList<RootAtomPair> list = new UniqueList<>();

		boolean[] isFirstOfReactantRank = new boolean[mReactant.getAtoms()];
		boolean[] isFirstOfProductRank = new boolean[mProduct.getAtoms()];
		boolean[] isRankUsed;

		isRankUsed = new boolean[1+mReactant.getAtoms()];
		for (int atom=0; atom<mReactant.getAtoms(); atom++) {
			int rank = mReactantCanonizer.getSymmetryRank(atom);
			if (!isRankUsed[rank]) {
				isRankUsed[rank] = true;
				isFirstOfReactantRank[atom] = true;
				}
			}
		isRankUsed = new boolean[1+mProduct.getAtoms()];
		for (int atom=0; atom<mProduct.getAtoms(); atom++) {
			int rank = mProductCanonizer.getSymmetryRank(atom);
			if (!isRankUsed[rank]) {
				isRankUsed[rank] = true;
				isFirstOfProductRank[atom] = true;
				}
			}

		// With highest priority take manually mapped seed atoms
		int maxMapNo = 0;
		for (int atom=0; atom<mProduct.getAtoms(); atom++)
			if (mProduct.getAtomMapNo(atom) != 0 && !mProduct.isAutoMappedAtom(atom))
				maxMapNo = Math.max(maxMapNo, mProduct.getAtomMapNo(atom));
		if (maxMapNo != 0) {
			int[] mapNoToProductAtom = new int[maxMapNo+1];
			for (int atom=0; atom<mProduct.getAtoms(); atom++)
				if (mProduct.getAtomMapNo(atom) != 0 && !mProduct.isAutoMappedAtom(atom))
					mapNoToProductAtom[mProduct.getAtomMapNo(atom)] = atom+1;
			for (int atom=0; atom<mReactant.getAtoms(); atom++) {
				int reactantMapNo = mReactant.getAtomMapNo(atom);
				if (reactantMapNo != 0
				 && reactantMapNo <= maxMapNo
				 && !mReactant.isAutoMappedAtom(atom)
				 && mapNoToProductAtom[reactantMapNo] != 0)
					list.add(new RootAtomPair(priority++, atom, mapNoToProductAtom[reactantMapNo]-1,true));
				}
			}

		// With second priority, we create starting pairs if we just have one atom of a kind
		for (byte[] reactantEnv:mEnvToReactantAtomMap[0].keySet()) {
			int[] productAtoms = mEnvToProductAtomMap[0].get(reactantEnv);
			if (productAtoms != null && productAtoms.length == 1) {
				int[] reactantAtoms = mEnvToReactantAtomMap[0].get(reactantEnv);
				if (reactantAtoms.length == 1)
					list.add(new RootAtomPair(priority++, reactantAtoms[0], productAtoms[0],true));
				}
			}

		for (int rootRadius=MAX_EXVIRONMENT_RADIUS-1; rootRadius>=MIN_ROOT_ATOM_RADIUS; rootRadius--) {
			for (byte[] reactantEnv:mEnvToReactantAtomMap[rootRadius].keySet()) {
				int[] productAtoms = mEnvToProductAtomMap[rootRadius].get(reactantEnv);
				if (productAtoms != null) {
					int[] reactantAtoms = mEnvToReactantAtomMap[rootRadius].get(reactantEnv);
					for (int reactantAtom:reactantAtoms)
						for (int productAtom:productAtoms)
							list.add(new RootAtomPair(priority++, reactantAtom, productAtom,
									isFirstOfReactantRank[reactantAtom] && isFirstOfProductRank[productAtom]));
					}
				}
			}

		return list.toArray(new RootAtomPair[0]);
		}

	/**
	 * Build canonical atom environment keys for all reactant and product atoms considering
	 * increasing bond count radia. Also builds maps from environment keys to respective atoms.
	 */
	private void initializeAtomEnvironments() {
		mReactantAtomEnv = classifyAtomEnvironment(mReactant);
		mProductAtomEnv = classifyAtomEnvironment(mProduct);

		mEnvToReactantAtomMap = buildEnvToAtomMaps(mReactant, mReactantAtomEnv);
		mEnvToProductAtomMap = buildEnvToAtomMaps(mProduct, mProductAtomEnv);

		mReactantMapNo = new int[mReactant.getAtoms()];
		mProductMapNo = new int[mProduct.getAtoms()];

		mMappableAtomCount = 0;
		for (byte[] reactantEnv:mEnvToReactantAtomMap[0].keySet()) {
			int[] productAtoms = mEnvToProductAtomMap[0].get(reactantEnv);
			if (productAtoms != null) {
				int[] reactantAtoms = mEnvToReactantAtomMap[0].get(reactantEnv);
				mMappableAtomCount += Math.min(productAtoms.length, reactantAtoms.length);
				}
			}

		mSimilarityComparator = new ByteArrayComparator();
		}

	/**
	 *
	 * @param reactantAtom
	 * @param productAtom
	 * @return 0: mismatch, 1:atomicNo-match, 2:directNeighbourMatch, 3:TwoShellMatch, etc.
	 */
	private int getAtomSimilarity(int reactantAtom, int productAtom) {
		for (int radius=0; radius<MAX_EXVIRONMENT_RADIUS; radius++)
			if (mSimilarityComparator.compare(mReactantAtomEnv[reactantAtom][radius], mProductAtomEnv[productAtom][radius]) != 0)
				return radius;
		return MAX_EXVIRONMENT_RADIUS;
		}

	private byte[][][] classifyAtomEnvironment(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());

		byte[][][] environment = new byte[mol.getAtoms()][MAX_EXVIRONMENT_RADIUS][];

		int[] atomList = new int[mol.getAtoms()];
		boolean[] atomMask = new boolean[mol.getAtoms()];
		for (int rootAtom=0; rootAtom<mol.getAtoms(); rootAtom++) {
			if (rootAtom != 0)
				Arrays.fill(atomMask, false);

			int min = 0;
			int max = 0;

			// we need to mark the root atom, because otherwise close-by root atoms may end up with the same fragment
			mol.setAtomSelection(rootAtom, true);

			for (int sphere=0; sphere<MAX_EXVIRONMENT_RADIUS && max<mol.getAtoms(); sphere++) {
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
					environment[rootAtom][sphere] = new Canonizer(fragment, Canonizer.ENCODE_ATOM_SELECTION).getIDCode().getBytes();
					}
				}

			mol.setAtomSelection(rootAtom, false);
			}
		return environment;
		}

	private TreeMap<byte[], int[]>[] buildEnvToAtomMaps(StereoMolecule mol, byte[][][] atomEnv) {
		TreeMap<byte[], int[]>[] map = new TreeMap[MAX_EXVIRONMENT_RADIUS];
		for (int radius=0; radius<MAX_EXVIRONMENT_RADIUS; radius++) {
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

	/**
	 * Copies all reactant/product molecules into one StereoMolecule for easier processing
	 * @param rxn
	 */
	private void copyReactionToMolecules(Reaction rxn) {
		mReactant = new StereoMolecule();
		for (int i=0; i<rxn.getReactants(); i++) {
			StereoMolecule reactant = rxn.getReactant(i);
			mReactant.addMolecule(reactant, reactant.getAtoms(), reactant.getBonds());
			}
		mReactantCanonizer = new Canonizer(mReactant, Canonizer.CREATE_SYMMETRY_RANK);

		mProduct = new StereoMolecule();
		for (int i=0; i<rxn.getProducts(); i++) {
			StereoMolecule product = rxn.getProduct(i);
			mProduct.addMolecule(product, product.getAtoms(), product.getBonds());
			}
		mProductCanonizer = new Canonizer(mProduct, Canonizer.CREATE_SYMMETRY_RANK);
		}

	/**
	 * Copies the generated mapping numbers from the two temporary molecules into the original reaction.
	 * @param rxn
	 */
	private void copyMapNosToReaction(Reaction rxn, int[] reactantMapNo, int[] productMapNo) {
		int reactantIndex = 0;
		int reactantAtom = 0;
		for (int atom=0; atom<mReactant.getAtoms(); atom++) {
			StereoMolecule reactant = rxn.getReactant(reactantIndex);
			if (reactantAtom == reactant.getAtoms()) {
				reactantAtom = 0;
				reactant = rxn.getReactant(++reactantIndex);
				}
			reactant.setAtomMapNo(reactantAtom, reactantMapNo[atom], true);
			reactantAtom++;
			}

		int productIndex = 0;
		int productAtom = 0;
		for (int atom=0; atom<mProduct.getAtoms(); atom++) {
			StereoMolecule product = rxn.getProduct(productIndex);
			if (productAtom == product.getAtoms()) {
				productAtom = 0;
				product = rxn.getProduct(++productIndex);
				}
			product.setAtomMapNo(productAtom, productMapNo[atom], true);
			productAtom++;
			}
		}

	/**
	 * Maps a many atoms a possible starting from the given root atom pair.
	 * Mapping is done by incrementally adding the most similar neighbours to
	 * the currently mapped area provided the the respective bond orders are
	 * exactly matching.
	 * @param root
	 * @return number mapped atoms
	 */
	private int mapFromRootAtoms(RootAtomPair root) {
		int[] graphAtom = new int[mReactant.getAtoms()];
		int[] productAtom = new int[mReactant.getAtoms()];

		graphAtom[0] = root.reactantAtom;
		productAtom[0] = root.productAtom;

		mMapNo++;
		mReactantMapNo[root.reactantAtom] = mMapNo;
		mProductMapNo[root.productAtom] = mMapNo;

		int bestReactantAtom = -1;
		int bestProductAtom = -1;

		int current = 0;
		int highest = 0;
		while (current <= highest) {
			// we take the best match of reaction & product neighbours as next graph members
			int bestMatch = 0;

			for (int i=0; i<mReactant.getConnAtoms(graphAtom[current]); i++) {
				int reactantCandidate = mReactant.getConnAtom(graphAtom[current], i);
				if (mReactantMapNo[reactantCandidate] != 0)
					continue;

				int bondType = mReactant.getBondTypeSimple(mReactant.getConnBond(graphAtom[current], i));

				for (int j=0; j<mProduct.getConnAtoms(productAtom[current]); j++) {
					int productCandidate = mProduct.getConnAtom(productAtom[current], j);
					if (mProductMapNo[productCandidate] == 0) {
						int candidateBond = mProduct.getConnBond(productAtom[current], j);
						if (bondType == mProduct.getBondTypeSimple(candidateBond)) {
							int match = getAtomSimilarity(reactantCandidate, productCandidate);
							if (bestMatch < match) {
								bestReactantAtom = reactantCandidate;
								bestProductAtom = productCandidate;
								bestMatch = match;
								}
							}
						}   // TODO may track multiple options and try all recursively; then end this method with best option mapped
					}
				}

			if (bestMatch != 0) {
				highest++;
				graphAtom[highest] = bestReactantAtom;
				productAtom[highest] = bestProductAtom;
				mMapNo++;
				mReactantMapNo[bestReactantAtom] = mMapNo;
				mProductMapNo[bestProductAtom] = mMapNo;
				}
			else {
				current++;
				}
			}

/*		if (mMapNo > mBestMapNo) {
			for (int atom=0; atom<mReactant.getAtoms(); atom++)
				mBestReactantMapNo[atom] = mReactantMapNo[atom];
			for (int atom=0; atom<mProduct.getAtoms(); atom++)
				mBestProductMapNo[atom] = mProductMapNo[atom];
			}

		// revert mapping numbers
		for (int i=0; i<=highest; i++) {
			mReactantMapNo[graphAtom[i]] = 0;
			mProductMapNo[productAtom[i]] = 0;
			mMapNo -= highest+1;
			}*/

		return highest+1;
		}

	class RootAtomPair implements Comparable<RootAtomPair> {
		public int priority,reactantAtom,productAtom;
		public boolean isFirstOfEquivalentPairs;

		public RootAtomPair(int priority, int reactantAtom, int productAtom, boolean isFirstOfEquivalentPairs) {
			this.priority = priority;
			this.reactantAtom = reactantAtom;
			this.productAtom = productAtom;
			this.isFirstOfEquivalentPairs = isFirstOfEquivalentPairs;
			}

		public boolean isNotYetMapped() {
			return mReactantMapNo[reactantAtom] == 0
				&& mProductMapNo[productAtom] == 0;
			}

		@Override
		public int compareTo(RootAtomPair pair) {
			return this.reactantAtom < pair.reactantAtom ? -1
				 : this.reactantAtom > pair.reactantAtom ? 1
				 : this.productAtom < pair.productAtom ? -1
				 : this.productAtom > pair.productAtom ? 1 : 0;
			}
		}
	}
