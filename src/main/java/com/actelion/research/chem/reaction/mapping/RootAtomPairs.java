package com.actelion.research.chem.reaction.mapping;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.ByteArrayComparator;
import com.actelion.research.util.UniqueList;

import java.util.Arrays;
import java.util.TreeMap;

public class RootAtomPairs {
	private static final int MAX_ENVIRONMENT_RADIUS = 8;
	private static final int MIN_ROOT_ATOM_RADIUS = 2;

	private StereoMolecule mReactant,mProduct;
	private Canonizer mReactantCanonizer,mProductCanonizer;
	private byte[][][] mReactantAtomEnv,mProductAtomEnv;   // indexes: atom,radius,idcode bytes
	private UniqueList<RootAtomPair> mList;
	private RootAtomPair[] mRootAtomPairs;
	private int mMappableAtomCount;
	private boolean mIsStoichiometric;
	private TreeMap<byte[],int[]>[] mEnvToReactantAtomMap,mEnvToProductAtomMap; // index: radius

	public RootAtomPairs(StereoMolecule reactant, StereoMolecule product) {
		mReactant = reactant;
		mProduct = product;
		mReactantCanonizer = new Canonizer(reactant, Canonizer.CREATE_SYMMETRY_RANK);
		mProductCanonizer = new Canonizer(product, Canonizer.CREATE_SYMMETRY_RANK);

		mReactantAtomEnv = classifyAtomEnvironment(mReactant);
		mProductAtomEnv = classifyAtomEnvironment(mProduct);
		mEnvToReactantAtomMap = buildEnvToAtomMaps(mReactant, mReactantAtomEnv);
		mEnvToProductAtomMap = buildEnvToAtomMaps(mProduct, mProductAtomEnv);

		mMappableAtomCount = 0;
		for (byte[] reactantEnv:mEnvToReactantAtomMap[0].keySet()) {
			int[] productAtoms = mEnvToProductAtomMap[0].get(reactantEnv);
			if (productAtoms != null) {
				int[] reactantAtoms = mEnvToReactantAtomMap[0].get(reactantEnv);
				mMappableAtomCount += Math.min(productAtoms.length, reactantAtoms.length);
				}
			}

		mIsStoichiometric = mMappableAtomCount == mReactant.getAtoms()
						 && mReactant.getAtoms() == mProduct.getAtoms();

		findRootAtomPairs();
		}

	public RootAtomPair[] get() {
		return mRootAtomPairs;
		}

	public int getMappableAtomCount() {
		return mMappableAtomCount;
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
					for (int atom=0; atom<fragment.getAllAtoms(); atom++) {
						fragment.setAtomCharge(atom, 0);
						fragment.setAtomRadical(atom, 0);
						}
					environment[rootAtom][sphere] = new Canonizer(fragment, Canonizer.ENCODE_ATOM_SELECTION).getIDCode().getBytes();
					}
				}

			mol.setAtomSelection(rootAtom, false);
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
	 * This compiles a complete list of reactant/product atom pairs, which serve as starting points
	 * to build equivalent atom graphs that logically match reactant subgraphs to product subgraphs.
	 * The list is sorted with decreasing plausibility/similarity of atom pairs.
	 * Root atom pairs meet these conditions:<br>
	 * - they match in terms of circular fragment on reactant and products side<br>
	 * - if multiple symmetrically equivalent pairs exist, exactly one of them is marked as allowed root pair<br>
	 * @return pair of currently obvious root atom pairs
	 */
	private void findRootAtomPairs() {
		mList = new UniqueList<>();

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

		int[] reactantFragmentNo = new int[mReactant.getAtoms()];
		int[] productFragmentNo = new int[mProduct.getAtoms()];
		boolean[] reactantFragmentUsed = new boolean[mReactant.getFragmentNumbers(reactantFragmentNo, false, false)];
		boolean[] productFragmentUsed = new boolean[mProduct.getFragmentNumbers(productFragmentNo, false, false)];

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
				 && mapNoToProductAtom[reactantMapNo] != 0) {
					int productAtom = mapNoToProductAtom[reactantMapNo]-1;
					add(atom, productAtom, (MAX_ENVIRONMENT_RADIUS+1) << 16);
					if (mProduct.getConnAtoms(productAtom) != 0)
						reactantFragmentUsed[reactantFragmentNo[atom]] = true;
					if (mReactant.getConnAtoms(atom) != 0)
						productFragmentUsed[productFragmentNo[productAtom]] = true;
					}
				}
			}

		for (int rootRadius=MAX_ENVIRONMENT_RADIUS-1; rootRadius>=MIN_ROOT_ATOM_RADIUS; rootRadius--) {
			// with equal environment size, we prefer carbon root atoms
			for (byte[] reactantEnv:mEnvToReactantAtomMap[rootRadius].keySet()) {
				int[] reactantAtoms = mEnvToReactantAtomMap[rootRadius].get(reactantEnv);
				int[] productAtoms = mEnvToProductAtomMap[rootRadius].get(reactantEnv);
				if (productAtoms != null) {
					float priority = (rootRadius << 16) + (mReactant.getAtomicNo(reactantAtoms[0]) == 6 ? 0.5f : 0f);
					for (int reactantAtom:reactantAtoms)
						for (int productAtom:productAtoms)
							add(reactantAtom, productAtom, priority);
					}
				}
			}

		// With third priority, we create starting pairs if we just have one atom of a kind
		if (mIsStoichiometric) {
			for (byte[] reactantEnv:mEnvToReactantAtomMap[0].keySet()) {
				int[] productAtoms = mEnvToProductAtomMap[0].get(reactantEnv);
				if (productAtoms != null && productAtoms.length == 1) {
					int[] reactantAtoms = mEnvToReactantAtomMap[0].get(reactantEnv);
					if (reactantAtoms.length == 1) {
						add(reactantAtoms[0], productAtoms[0], 0);
						if (mProduct.getConnAtoms(productAtoms[0]) != 0)
							reactantFragmentUsed[reactantFragmentNo[reactantAtoms[0]]] = true;
						if (mReactant.getConnAtoms(reactantAtoms[0]) != 0)
							productFragmentUsed[productFragmentNo[productAtoms[0]]] = true;
						}
					}
				}
			}

		// Check with decreasing sphere size down to atomicNo level, whether all reactant or product atoms
		// are in separate fragments and are equivalent. If this is the case, and if we have at least
		// as many equivalent atoms as atoms on the other side, then we just match reactant atoms
		// to product atoms in order of appearance.
		for (int rootRadius=MAX_ENVIRONMENT_RADIUS-1; rootRadius>=0; rootRadius--) {
			for (byte[] reactantEnv:mEnvToReactantAtomMap[rootRadius].keySet()) {
				int[] productAtoms = mEnvToProductAtomMap[rootRadius].get(reactantEnv);
				if (productAtoms != null) {
					int[] reactantAtoms = mEnvToReactantAtomMap[rootRadius].get(reactantEnv);
					if (reactantAtoms.length == 1 && productAtoms.length == 1) {
						add(reactantAtoms[0], productAtoms[0], (rootRadius << 16) - 1);
						if (mProduct.getConnAtoms(productAtoms[0]) != 0)
							reactantFragmentUsed[reactantFragmentNo[reactantAtoms[0]]] = true;
						if (mReactant.getConnAtoms(reactantAtoms[0]) != 0)
							productFragmentUsed[productFragmentNo[productAtoms[0]]] = true;
							}
					else if ((reactantAtoms.length >= productAtoms.length
							&& areInDistinctEquivalentFragments(reactantAtoms, reactantFragmentNo, reactantFragmentUsed, mReactantCanonizer))
							|| (productAtoms.length >= reactantAtoms.length
							&& areInDistinctEquivalentFragments(productAtoms, productFragmentNo, productFragmentUsed, mProductCanonizer))) {
						for (int i=0; i<Math.min(reactantAtoms.length, productAtoms.length); i++) {
							add(reactantAtoms[i], productAtoms[i], (rootRadius << 16) - 1);
							if (mProduct.getConnAtoms(productAtoms[i]) != 0)
								reactantFragmentUsed[reactantFragmentNo[reactantAtoms[i]]] = true;
							if (mReactant.getConnAtoms(reactantAtoms[i]) != 0)
								productFragmentUsed[productFragmentNo[productAtoms[i]]] = true;
							}
						}
					}
				}
			}

		mRootAtomPairs = mList.toArray(new RootAtomPair[0]);
		Arrays.sort(mRootAtomPairs, (o1, o2) -> o1.priority < o2.priority ? 1 : o1.priority > o2.priority ? -1 : 0);
		}

	private void add(int reactantAtom, int productAtom, float priority) {
		RootAtomPair pair = new RootAtomPair(reactantAtom, productAtom, priority);

		int index = mList.getIndex(pair);
		if (index == -1)
			mList.add(pair);
		else if (mList.get(index).priority < priority)
			mList.get(index).priority = priority;
		}

	private boolean areInDistinctEquivalentFragments(int[] atom, int[] fragmentNo, boolean[] fragmentUsed, Canonizer canonizer) {
		for (int a:atom)
			if (fragmentUsed[fragmentNo[a]])
				return false;

		int rank = canonizer.getSymmetryRank(atom[0]);
		for (int i=1; i<atom.length; i++)
			if (canonizer.getSymmetryRank(atom[i]) != rank)
				return false;

		return true;
		}
	}

class RootAtomPair implements Comparable<RootAtomPair> {
	public int reactantAtom,productAtom;
	public float priority;

	public RootAtomPair(int reactantAtom, int productAtom, float priority) {
		this.reactantAtom = reactantAtom;
		this.productAtom = productAtom;
		this.priority = priority;
		}

	@Override
	public int compareTo(RootAtomPair pair) {
		return this.reactantAtom < pair.reactantAtom ? -1
				: this.reactantAtom > pair.reactantAtom ? 1
				: this.productAtom < pair.productAtom ? -1
				: this.productAtom > pair.productAtom ? 1 : 0;
		}
	}
