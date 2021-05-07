package com.actelion.research.chem.reaction.mapping;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.util.ByteArrayComparator;
import com.actelion.research.util.UniqueList;

import java.util.Arrays;
import java.util.TreeMap;

public class SimilarityGraphBasedReactionMapper {
	private static final int MAX_ENVIRONMENT_RADIUS = 8;
	private static final int MIN_ROOT_ATOM_RADIUS = 2;
	private StereoMolecule mReactant,mProduct;
	private int mMapNo,mGraphMapNoCount,mMappableAtomCount;
	private int[] mReactantMapNo,mProductMapNo,mReactantRingMembership,mProductRingMembership;
	private Canonizer mReactantCanonizer,mProductCanonizer;
	private ByteArrayComparator mSimilarityComparator;
	private byte[][][] mReactantAtomEnv,mProductAtomEnv;   // indexes: atom,radius,idcode bytes
	private byte[][][][] mReactantConnAtomEnv,mProductConnAtomEnv;   // indexes: atom,connIndex,radius,idcode bytes
	private TreeMap<byte[],int[]>[] mEnvToReactantAtomMap,mEnvToProductAtomMap;

	/**
	 * Entirely maps the given reaction by setting all reactant's and product's mapping numbers.
	 * If the reaction already contains manually mapped atoms, then these are untouched and used
	 * as seed atoms when building the mapping graph.
	 * @param rxn reaction with proper atom coordinates
	 */
	public void map(Reaction rxn) {
		mergeReactantsAndProducts(rxn);
		map(mReactant, mProduct, new int[mReactant.getAtoms()], new int[mProduct.getAtoms()]);
		copyMapNosToReaction(rxn, mReactantMapNo, mProductMapNo, mGraphMapNoCount);
		System.out.println("mapped:"+mMapNo+"("+mGraphMapNoCount+") of "+mMappableAtomCount+" score:"+calculateScore());
		}

	/**
	 * Entirely maps the given reaction by setting all reactant's and product's mapping numbers.
	 * @param reactant all reactants merged into one molecule; may contain seed mapping
	 * @param product all products merged into one molecule; may contain seed mapping
	 * @param reactantMapNo array to receive the created mapping numbers
	 * @param productMapNo array to receive the created mapping numbers
	 */
	public void map(StereoMolecule reactant, StereoMolecule product, int[] reactantMapNo, int[] productMapNo) {
		mReactant = reactant;
		mProduct = product;

		mReactantCanonizer = new Canonizer(mReactant, Canonizer.CREATE_SYMMETRY_RANK);
		mProductCanonizer = new Canonizer(mProduct, Canonizer.CREATE_SYMMETRY_RANK);

		initializeAtomEnvironments();
		initializeRingMembership();

		RootAtomPair[] rootAtomPairs = findRootAtomPairs();

		mReactantMapNo = reactantMapNo;
		mProductMapNo = productMapNo;
		mMapNo = 0;

		for (RootAtomPair pair:rootAtomPairs)
			if (pair.isNotYetMapped() && pair.isFirstOfEquivalentPairs)
				mapFromRootAtoms(pair);
		for (RootAtomPair pair:rootAtomPairs)
			if (pair.isNotYetMapped())
				mapFromRootAtoms(pair);

		mGraphMapNoCount = mMapNo;

		if (mMapNo < mMappableAtomCount) {
			ReactionCenterMapper centerMapper = new ReactionCenterMapper(mReactant, mProduct, mReactantMapNo, mProductMapNo, mMapNo);
			int score = centerMapper.completeMapping();
			mMapNo += centerMapper.getMappedAtomCount();
			}
		}

	/**
	 * @return number of atoms that could be mapped by graph similarity matching (no reaction center mapping)
	 */
	public int getGraphMapNoCount() {
		return mGraphMapNoCount;
		}

	public int calculateScore() {
		int score = 0;

		int[] mapNoToProduct = new int[mMapNo+1];
		Arrays.fill(mapNoToProduct, -1);
		for (int productAtom=0; productAtom<mProductMapNo.length; productAtom++)
			if (mProductMapNo[productAtom] != 0)
				mapNoToProduct[mProductMapNo[productAtom]] = productAtom;

		boolean[] productBondHandled = new boolean[mProduct.getBonds()];

		for (int reactantBond=0; reactantBond<mReactant.getBonds(); reactantBond++) {
			boolean found = false;
			int productAtom1 = mapNoToProduct[mReactantMapNo[mReactant.getBondAtom(0, reactantBond)]];
			int productAtom2 = mapNoToProduct[mReactantMapNo[mReactant.getBondAtom(1, reactantBond)]];
			if (productAtom1 != -1 && productAtom2 != -1) {
				int productBond = mProduct.getBond(productAtom1, productAtom2);
				if (productBond != -1) {
					found = true;
					productBondHandled[productBond] = true;

					// we give one point for a changed bond type, which includes delocalized as separate bond type
					if (getBondType(mReactant, reactantBond) != getBondType(mProduct, productBond))
						score--;
					}
				}

			// we give two points for removed or added bonds, no matter what type they are
			if (!found)
				score -= 2;
			}

		for (boolean pbh:productBondHandled)
			if (!pbh)
				score -= 2;

		return score;
		}

	/**
	 * This compiles a complete list of reactant/product atom pairs, which serve as starting points
	 * to build equivalent atom graphs that logically match reactant subgraphs to product subgraphs.
	 * The list is sorted with decreasing plausibility/similarity of atom pairs.
	 * Root atom pairs meet these conditions:<br>
	 * - they match in terms of circular fragment on reactant and products side<br>
	 * - if multiple symmetrically equivalent pairs exist, exactly one of them is marked as allowed root pair<br>
	 * @return pair of currently obvious root atom pairs
	 */
	private RootAtomPair[] findRootAtomPairs() {
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
					list.add(new RootAtomPair(atom, mapNoToProductAtom[reactantMapNo]-1,true));
				}
			}

		// If we don't have manual seed atoms, then with second priority we match just one pair,
		// if we find a shared atom environment with all atoms being symmetric on reactant and on product side
		for (byte[] reactantEnv:mEnvToReactantAtomMap[0].keySet()) {
			if (list.size() != 0)
				break;

			int[] productAtoms = mEnvToProductAtomMap[0].get(reactantEnv);
			if (productAtoms != null) {
				int  productAtomRank = -1;
				for (int atom:productAtoms) {
					if (productAtomRank == -1)
						productAtomRank = mProductCanonizer.getSymmetryRank(atom);
					else if (productAtomRank != mProductCanonizer.getSymmetryRank(atom)) {
						productAtomRank = -2;
						break;
						}
					}
				if (productAtomRank != -2) {
					int[] reactantAtoms = mEnvToReactantAtomMap[0].get(reactantEnv);
					int  reactantAtomRank = -1;
					for (int atom:reactantAtoms) {
						if (reactantAtomRank == -1)
							reactantAtomRank = mReactantCanonizer.getSymmetryRank(atom);
						else if (reactantAtomRank != mReactantCanonizer.getSymmetryRank(atom)) {
							reactantAtomRank = -2;
							break;
							}
						}
					if (reactantAtomRank != -2) {
						list.add(new RootAtomPair(reactantAtoms[0], productAtoms[0], true));
						}
					}
				}
			}

		// With third priority, we create starting pairs if we just have one atom of a kind or multiple
		// symmetrical ones in independent fragments
		for (byte[] reactantEnv:mEnvToReactantAtomMap[0].keySet()) {
			int[] productAtoms = mEnvToProductAtomMap[0].get(reactantEnv);
			if (productAtoms != null && productAtoms.length == 1) {
				int[] reactantAtoms = mEnvToReactantAtomMap[0].get(reactantEnv);
				if (reactantAtoms.length == 1)
					list.add(new RootAtomPair(reactantAtoms[0], productAtoms[0],true));
				}
			}

		// Check with decreasing sphere size down to atomicNo level, whether all reactant or product atoms
		// are in separate fragments and are equivalent. If this is the case, and if we have at least
		// as many equivalent atoms as atoms on the other side, then we just match reactant atoms
		// to product atoms in order of appearance.
		int[] reactantFragmentNo = new int[mReactant.getAtoms()];
		boolean[] reactantFragmentUsed = new boolean[mReactant.getFragmentNumbers(reactantFragmentNo, false, false)];
		int[] productFragmentNo = new int[mProduct.getAtoms()];
		boolean[] productFragmentUsed = new boolean[mProduct.getFragmentNumbers(productFragmentNo, false, false)];
		for (int rootRadius = MAX_ENVIRONMENT_RADIUS -1; rootRadius>=0; rootRadius--) {
			for (byte[] reactantEnv:mEnvToReactantAtomMap[rootRadius].keySet()) {
				int[] productAtoms = mEnvToProductAtomMap[rootRadius].get(reactantEnv);
				if (productAtoms != null) {
					int[] reactantAtoms = mEnvToReactantAtomMap[rootRadius].get(reactantEnv);
					if ((reactantAtoms.length >= productAtoms.length
					  && areInDistinctEquivalentFragments(reactantAtoms, reactantFragmentNo, reactantFragmentUsed, mReactantCanonizer))
					 || (productAtoms.length >= reactantAtoms.length
					  && areInDistinctEquivalentFragments(productAtoms, productFragmentNo, productFragmentUsed, mProductCanonizer))) {
						for (int i=0; i<Math.min(reactantAtoms.length, productAtoms.length); i++)
							list.add(new RootAtomPair(reactantAtoms[i], productAtoms[i], true));
						}
					}
				}
			}

		for (int rootRadius = MAX_ENVIRONMENT_RADIUS -1; rootRadius>=MIN_ROOT_ATOM_RADIUS; rootRadius--) {
			for (byte[] reactantEnv:mEnvToReactantAtomMap[rootRadius].keySet()) {
				int[] productAtoms = mEnvToProductAtomMap[rootRadius].get(reactantEnv);
				if (productAtoms != null) {
					int[] reactantAtoms = mEnvToReactantAtomMap[rootRadius].get(reactantEnv);
					for (int reactantAtom:reactantAtoms)
						for (int productAtom:productAtoms)
							list.add(new RootAtomPair(reactantAtom, productAtom,
									isFirstOfReactantRank[reactantAtom] && isFirstOfProductRank[productAtom]));
					}
				}
			}

		return list.toArray(new RootAtomPair[0]);
		}

	private boolean areInDistinctEquivalentFragments(int[] atom, int[] fragmentNo, boolean[] fragmentUsed, Canonizer canonizer) {
		int rank = canonizer.getSymmetryRank(atom[0]);
		for (int i=1; i<atom.length; i++)
			if (canonizer.getSymmetryRank(atom[i]) != rank)
				return false;

		Arrays.fill(fragmentUsed, false);
		for (int a:atom)
			if (fragmentUsed[fragmentNo[a]])
				return false;

		return true;
		}

	/**
	 * Build canonical atom environment keys for all reactant and product atoms considering
	 * increasing bond count radia. Also builds maps from environment keys to respective atoms.
	 */
	private void initializeAtomEnvironments() {
		mReactantAtomEnv = classifyAtomEnvironment(mReactant);
		mProductAtomEnv = classifyAtomEnvironment(mProduct);

		mReactantConnAtomEnv = classifyNeighbourAtomEnvironment(mReactant);
		mProductConnAtomEnv = classifyNeighbourAtomEnvironment(mProduct);

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

	private void initializeRingMembership() {
		mReactantRingMembership = initializeRingMembership(mReactant);
		mProductRingMembership = initializeRingMembership(mProduct);
		}

	private int[] initializeRingMembership(StereoMolecule mol) {
		int[] ringMembership = new int[mol.getAtoms()];
		RingCollection ringSet = mol.getRingSet();
		for (int ring=0; ring<Math.min(32, ringSet.getSize()); ring++) {
			int[] ringAtom = ringSet.getRingAtoms(ring);
			for (int atom:ringAtom)
				ringMembership[atom] |= (1 << ring);
			}
		return ringMembership;
		}

	/**
	 * @param atom1
	 * @param atom2
	 * @param atom3
	 * @return true, if atom3 is not member of any ring that is shared by atom1 and atom2
	 */
	private boolean leavesRing(int atom1, int atom2, int atom3, int[] ringMembership) {
		return (ringMembership[atom1] & ringMembership[atom2] & ~ringMembership[atom3]) != 0;
		}

	/**
	 * Compares reactantAtom to productAtom starting from and extending away from the respective root atoms.
	 * Similarity is increased, if reactantRoot and productRoot are stereo centers and reactantAtom matches
	 * productAtom in regard to the stereo configuration.
	 * @param reactantRoot
	 * @param reactantAtom
	 * @param productRoot
	 * @param productAtom
	 * @return 0: mismatch, 2:atomicNo-match, 4:directNeighbourMatch, 6:TwoShellMatch, etc.
	 */
	private int getAtomSimilarity(int reactantRoot, int reactantAtom, int productRoot, int productAtom) {
		int reactantConnIndex = -1;
		for (int i=0; i<mReactant.getConnAtoms(reactantAtom); i++) {
			if (mReactant.getConnAtom(reactantAtom, i) == reactantRoot) {
				reactantConnIndex = i;
				break;
				}
			}

		int productConnIndex = -1;
		for (int i=0; i<mProduct.getConnAtoms(productAtom); i++) {
			if (mProduct.getConnAtom(productAtom, i) == productRoot) {
				productConnIndex = i;
				break;
				}
			}

		for (int radius=0; radius<MAX_ENVIRONMENT_RADIUS; radius++)
			if (mSimilarityComparator.compare(mReactantConnAtomEnv[reactantAtom][reactantConnIndex][radius], mProductConnAtomEnv[productAtom][productConnIndex][radius]) != 0)
				return radius;

		return MAX_ENVIRONMENT_RADIUS;
		}

	private byte[][][] classifyAtomEnvironment(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());

		byte[][][] environment = new byte[mol.getAtoms()][MAX_ENVIRONMENT_RADIUS][];

		int[] atomList = new int[mol.getAtoms()];
		boolean[] atomMask = new boolean[mol.getAtoms()];
		for (int rootAtom=0; rootAtom<mol.getAtoms(); rootAtom++) {
			if (rootAtom != 0)
				Arrays.fill(atomMask, false);

			int min = 0;
			int max = 0;

			// we need to mark the root atom, because otherwise close-by root atoms may end up with the same fragment
			mol.setAtomSelection(rootAtom, true);

			for (int sphere = 0; sphere<MAX_ENVIRONMENT_RADIUS && max<mol.getAtoms(); sphere++) {
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

	private byte[][][][] classifyNeighbourAtomEnvironment(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());

		byte[][][][] environment = new byte[mol.getAtoms()][MAX_ENVIRONMENT_RADIUS][][];

		int[] atomList = new int[mol.getAtoms()];
		boolean[] atomMask = new boolean[mol.getAtoms()];
		for (int rootAtom=0; rootAtom<mol.getAtoms(); rootAtom++) {
			environment[rootAtom] = new byte[mol.getConnAtoms(rootAtom)][MAX_ENVIRONMENT_RADIUS][];
			for (int connIndex=0; connIndex<mol.getConnAtoms(rootAtom); connIndex++) {
				int exitAtom = mol.getConnAtom(rootAtom, connIndex);

				if (atomMask == null)
					atomMask = new boolean[mol.getAtoms()];
				else
					Arrays.fill(atomMask, false);

				int min = 1;
				int max = 2;

				// we need to mark the root atom, because otherwise close-by root atoms may end up with the same fragment
				mol.setAtomSelection(exitAtom, true);

				for (int sphere = 0; sphere<MAX_ENVIRONMENT_RADIUS && max<mol.getAtoms(); sphere++) {
					if (sphere == 0) {
						atomList[0] = exitAtom;
						atomMask[exitAtom] = true;
						atomList[1] = rootAtom;
						atomMask[rootAtom] = true;
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
						environment[rootAtom][connIndex][sphere] = new byte[2];
						environment[rootAtom][connIndex][sphere][0] = (byte)mol.getAtomicNo(rootAtom);
						environment[rootAtom][connIndex][sphere][1] = (byte)mol.getAtomMass(rootAtom);
						}
					else {
						mol.copyMoleculeByAtoms(fragment, atomMask, true, null);
						environment[rootAtom][connIndex][sphere] = new Canonizer(fragment, Canonizer.ENCODE_ATOM_SELECTION).getIDCode().getBytes();
						}
					}

				mol.setAtomSelection(exitAtom, false);
				}
			}
		return environment;
		}

	private TreeMap<byte[], int[]>[] buildEnvToAtomMaps(StereoMolecule mol, byte[][][] atomEnv) {
		TreeMap<byte[], int[]>[] map = new TreeMap[MAX_ENVIRONMENT_RADIUS];
		for (int radius = 0; radius<MAX_ENVIRONMENT_RADIUS; radius++) {
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
	 * Use getReactant() and getProduct() from outside to access these.
	 * @param rxn
	 */
	public void mergeReactantsAndProducts(Reaction rxn) {
		mReactant = new StereoMolecule();
		for (int i=0; i<rxn.getReactants(); i++) {
			StereoMolecule reactant = rxn.getReactant(i);
			mReactant.addMolecule(reactant, reactant.getAtoms(), reactant.getBonds());
			}

		mProduct = new StereoMolecule();
		for (int i=0; i<rxn.getProducts(); i++) {
			StereoMolecule product = rxn.getProduct(i);
			mProduct.addMolecule(product, product.getAtoms(), product.getBonds());
			}
		}

	public StereoMolecule getReactant() {
		return mReactant;
		}

	public StereoMolecule getProduct() {
		return mProduct;
		}

	/**
	 * Copies the generated mapping numbers from the merged reactant/product molecules into the original reaction.
	 * @param rxn
	 */
	public void copyMapNosToReaction(Reaction rxn, int[] reactantMapNo, int[] productMapNo, int graphMapNoCount) {
		int reactantIndex = 0;
		int reactantAtom = 0;
		for (int atom=0; atom<mReactant.getAtoms(); atom++) {
			StereoMolecule reactant = rxn.getReactant(reactantIndex);
			if (reactantAtom == reactant.getAtoms()) {
				reactantAtom = 0;
				reactant = rxn.getReactant(++reactantIndex);
				}
			reactant.setAtomMapNo(reactantAtom, reactantMapNo[atom], reactantMapNo[atom] <= graphMapNoCount);
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
			product.setAtomMapNo(productAtom, productMapNo[atom], productMapNo[atom] <= graphMapNoCount);
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
		int[] productAtom = new int[mProduct.getAtoms()];

		int[] graphParent = new int[mReactant.getAtoms()];
		int[] productParent = new int[mProduct.getAtoms()];

		graphAtom[0] = root.reactantAtom;
		graphParent[root.reactantAtom] = root.reactantAtom; // no -1 to simplify leavesRing()
		productAtom[0] = root.productAtom;
		productParent[root.productAtom] = root.productAtom; // no -1 to simplify leavesRing()

		mMapNo++;
		mReactantMapNo[root.reactantAtom] = mMapNo;
		mProductMapNo[root.productAtom] = mMapNo;

		int current = 0;
		int highest = 0;
		while (current <= highest) {
			// we take the best match of reaction & product neighbours as next graph members
			int reactantRoot = graphAtom[current];
			int productRoot = productAtom[current];
			int[][] match = new int[mReactant.getConnAtoms(reactantRoot)][mProduct.getConnAtoms(productRoot)];

			for (int i=0; i<mReactant.getConnAtoms(reactantRoot); i++) {
				int reactantCandidate = mReactant.getConnAtom(reactantRoot, i);
				if (mReactantMapNo[reactantCandidate] != 0)
					continue;

				int bondType = getBondType(mReactant, mReactant.getConnBond(reactantRoot, i));

				for (int j=0; j<mProduct.getConnAtoms(productRoot); j++) {
					int productCandidate = mProduct.getConnAtom(productRoot, j);
					if (mProductMapNo[productCandidate] == 0
					 && mReactant.getAtomicNo(reactantCandidate) == mProduct.getAtomicNo(productCandidate)) {
						int candidateBond = mProduct.getConnBond(productRoot, j);
						if (bondType == getBondType(mProduct, candidateBond)) {
							if (passesBasicRules(reactantRoot, reactantCandidate, productRoot, productCandidate)) {
								match[i][j] = 256 * getAtomSimilarity(reactantRoot, reactantCandidate, productRoot, productCandidate);
								if (match[i][j] != 0) {
									if (matchesStereo(graphParent[reactantRoot], reactantRoot, reactantCandidate, productParent[productRoot], productRoot, productCandidate))
										match[i][j]++;
									if (leavesRing(graphParent[reactantRoot], reactantRoot, reactantCandidate, mReactantRingMembership)
									 == leavesRing(productParent[productRoot], productRoot, productCandidate, mProductRingMembership))
										match[i][j] += 2;
									}
/*
if (reactantRoot == 4) {
 System.out.print("graph:"); for (int k=0; k<current; k++) System.out.print(" "+graphAtom[k]+"("+productAtom[k]+")");
 System.out.println(" match:"+match[i][j]+ " r:"+reactantCandidate+" p:"+productCandidate);
}*/
								}
							}
						// if we have a changing bond order, but the rest of the substituent has maximal similarity
						else if (getAtomSimilarity(reactantRoot, reactantCandidate, productRoot, productCandidate) == MAX_ENVIRONMENT_RADIUS) {
//							match[i][j] = 4 * MAX_ENVIRONMENT_RADIUS - 1;
							match[i][j] = 3;
							}
						}   // TODO may track multiple options and try all recursively; then end this method with best option mapped
					}
				}

			while (true) {
				int bestMatch = 0;
				int bestReactantAtom = -1;
				int bestProductAtom = -1;
				for (int i=0; i<match.length; i++) {
					int reactantCandidate = mReactant.getConnAtom(reactantRoot, i);
					if (mReactantMapNo[reactantCandidate] == 0) {
						for (int j=0; j<match[i].length; j++) {
							int productCandidate = mProduct.getConnAtom(productRoot, j);
							if (mProductMapNo[productCandidate] == 0) {
								if (bestMatch < match[i][j]) {
									bestMatch = match[i][j];
									bestReactantAtom = reactantCandidate;
									bestProductAtom = productCandidate;
									}
								}
							}
						}
					}

				if (bestMatch == 0)
					break;

				highest++;
				graphAtom[highest] = bestReactantAtom;
				graphParent[bestReactantAtom] = graphAtom[current];
				productAtom[highest] = bestProductAtom;
				productParent[bestProductAtom] = productAtom[current];
				mMapNo++;
				mReactantMapNo[bestReactantAtom] = mMapNo;
				mProductMapNo[bestProductAtom] = mMapNo;
				}

			current++;
			}

		return highest+1;
		}

	private int getBondType(StereoMolecule mol, int bond) {
		if (mol.isDelocalizedBond(bond))
			return 0;
		return mol.getBondTypeSimple(bond);
		}

	// Typically we accept the most similar neighbour pair as next atom pair in the mapping graph,
	// provided that bond order matches. In some cases, however, we apply some chemical knowledge,
	// to apply a veto, e.g. to prevent the match of the wrong oxygen atom in carboxylic acids
	// or esters.
	private boolean passesBasicRules(int reactantAtom, int reactantConn, int productAtom, int productConn) {
		// reductive ester cleavage
		if (mReactant.getAtomicNo(reactantConn) == 8
		 && mReactant.getConnAtoms(reactantConn) == 2
		 && connAtomsOfAtomicNo(mReactant, reactantAtom, 8) > connAtomsOfAtomicNo(mProduct, productAtom, 8))
			return false;
		if (mProduct.getAtomicNo(productConn) == 8
		 && mProduct.getConnAtoms(productConn) == 2
		 && connAtomsOfAtomicNo(mReactant, reactantAtom, 8) < connAtomsOfAtomicNo(mProduct, productAtom, 8))
			return false;

		// esterification / reesterification
		if (mReactant.getAtomicNo(reactantConn) == 8
		 && mReactant.getAtomPi(reactantConn) == 0
		 && hasOxo(mReactant, reactantAtom)
		 && hasOxo(mProduct, productAtom)
		 && getAtomSimilarity(reactantAtom, reactantConn, productAtom, productConn) != MAX_ENVIRONMENT_RADIUS)
			return false;

		// potential replacement of -O-R (R:H,C,any) at non-carbon atoms
		if (mReactant.getAtomicNo(reactantConn) == 8
		 && mReactant.getAtomicNo(reactantAtom) != 6
		 && getAtomSimilarity(reactantAtom, reactantConn, productAtom, productConn) != MAX_ENVIRONMENT_RADIUS)
			return false;

		// imine formation
		if (mReactant.getAtomicNo(reactantConn) == 7
		 && connAtomsOfAtomicNo(mReactant, reactantAtom, 7) < connAtomsOfAtomicNo(mProduct, productAtom, 7))
			return false;
		if (mProduct.getAtomicNo(productConn) == 7
		 && connAtomsOfAtomicNo(mReactant, reactantAtom, 7) > connAtomsOfAtomicNo(mProduct, productAtom, 7))
			return false;

		return true;
		}

	/**
	 * Returns true, if we have a stereo center (or EZ-bond) and if the
	 * @param reactantParent
	 * @param reactantAtom
	 * @param reactantConn
	 * @param productParent
	 * @param productAtom
	 * @param productConn
	 * @return
	 */
	private boolean matchesStereo(int reactantParent, int reactantAtom, int reactantConn, int productParent, int productAtom, int productConn) {
		if (mReactant.isAtomStereoCenter(reactantAtom) && mReactant.getConnAtoms(reactantAtom) == 3
		 && mProduct.isAtomStereoCenter(productAtom) && mProduct.getConnAtoms(productAtom) == 3) {
			boolean reactantInversion = (reactantParent > reactantConn);
			for (int i=0; i<mReactant.getConnAtoms(reactantAtom); i++) {
				int lastNeighbour = mReactant.getConnAtom(reactantAtom, i);
				if (lastNeighbour != reactantParent && lastNeighbour != reactantConn) {
					if (i == 1)
						reactantInversion = !reactantInversion;
					}
				}

			boolean productInversion = (productParent > productConn);
			for (int i=0; i<mProduct.getConnAtoms(productAtom); i++) {
				int lastNeighbour = mProduct.getConnAtom(productAtom, i);
				if (lastNeighbour != productParent && lastNeighbour != productConn) {
					if (i == 1)
						productInversion = !productInversion;
					}
				}

			return (reactantInversion == productInversion)
				== (mReactant.getAtomParity(reactantAtom) == mProduct.getAtomParity(productAtom));
			}

		return true;
		}

	private boolean hasOxo(StereoMolecule mol, int atom) {
		for (int i=0; i<mol.getConnAtoms(atom); i++)
			if (mol.getConnBondOrder(atom, i) == 2
			 && mol.getAtomicNo(mol.getConnAtom(atom, i)) == 8)
				return true;
		return false;
		}

	private boolean isOrganicAcid(StereoMolecule mol, int atom) {
		boolean oxoLikeFound = false;
		boolean hydroxyFound = false;
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int connAtom = mol.getConnAtom(atom, i);
			if (mol.getAtomicNo(connAtom) == 8
			 && mol.getConnAtoms(connAtom) == 1
			 && mol.getConnBondOrder(atom, i) == 1) {
				hydroxyFound = true;
				}
			if (mol.isElectronegative(connAtom)
			 && mol.getConnBondOrder(atom, i) == 2) {
				oxoLikeFound = true;
				}
			}
		return oxoLikeFound && hydroxyFound;
		}

	private int connAtomsOfAtomicNo(StereoMolecule mol, int atom, int atomicNo) {
		int count = 0;
		for (int i=0; i<mol.getConnAtoms(atom); i++)
			if (mol.getAtomicNo(mol.getConnAtom(atom, i)) == atomicNo)
				count++;

		return count;
		}

	class RootAtomPair implements Comparable<RootAtomPair> {
		public int reactantAtom,productAtom;
		public boolean isFirstOfEquivalentPairs;

		public RootAtomPair(int reactantAtom, int productAtom, boolean isFirstOfEquivalentPairs) {
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
