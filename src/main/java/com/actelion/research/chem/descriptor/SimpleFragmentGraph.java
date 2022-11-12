package com.actelion.research.chem.descriptor;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.BurtleHasher;

import java.util.Arrays;

/**
 * This class contains a fragment buffer used to describe and canonicalize a molecular fragment
 * defined by certain bond and atom features and their connectivity. Its sole purpose is to
 * calculate a hash value to populate a fragment specific bit in a fingerprint.
 * By just using a small subset of the atom and bond features and neglecting any stereochemistry
 * of a normal molecule, the creation of a canonical representation is much faster than the
 * alternative of using Molecule.copyMoleculeByBonds() and a Canonizer().
 */
public class SimpleFragmentGraph {
	private static final int HASH_INIT = 13;
	private static final int MAX_CONN_ATOMS = 8;
	private static final int ATOM_INDEX_BITS = 4;    // bits needed to store either atom index or atom rank
	private static final int BOND_DESC_BITS = 2;     // bits needed to store bond descriptor
	private static final long ATOM_INDEX_MASK = 0xF;
	// MAX_CONN_ATOMS * (ATOM_INDEX_BITS + BOND_DESC_BITS) + ATOM_INDEX_BITS must not exceed 64 bits (Long.SIZE)

	private int[] mConnAtoms,mFragmentAtomFromOrig,mCanRank,mGraphAtom,mGraphIndex;
	private int[][] mConnAtom,mConnBond,mConnRank;
	private int mAtoms,mBonds;
	private byte[] mAtomDescriptor,mBondDescriptor,mBuffer;

	public SimpleFragmentGraph(int maxFragmentBonds) {
		init(maxFragmentBonds, 256);
	}

	public SimpleFragmentGraph(StereoMolecule mol, int[] bondMember, int bondCount) {
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		init(bondCount, mol.getAtoms());
		setMolecule(mol, bondMember, bondCount);
	}

	private void init(int maxFragmentBonds, int atomCount) {
		int maxFragmentAtoms = maxFragmentBonds+1;
		mAtomDescriptor = new byte[maxFragmentAtoms];
		mBondDescriptor = new byte[maxFragmentBonds];
		mConnAtoms = new int[maxFragmentAtoms];
		mCanRank = new int[maxFragmentAtoms];
		mGraphAtom = new int[maxFragmentAtoms];
		mGraphIndex = new int[maxFragmentAtoms];
		mFragmentAtomFromOrig = new int[atomCount];
		mConnAtom = new int[maxFragmentAtoms][MAX_CONN_ATOMS];
		mConnBond = new int[maxFragmentAtoms][MAX_CONN_ATOMS];
		mConnRank = new int[1+MAX_CONN_ATOMS][];
		for (int i=1; i<=MAX_CONN_ATOMS; i++)    // we keep a buffer for every connAtoms value
			mConnRank[i] = new int[i];
		mBuffer = new byte[3*maxFragmentBonds];
	}

	public void set(StereoMolecule mol, int[] bondMember, int bondCount) {
		mAtoms = 0;
		mBonds = 0;

		int maxAtomCount = mFragmentAtomFromOrig.length;
		if (maxAtomCount < mol.getAtoms()) {
			do {
				maxAtomCount *= 2;
			} while (maxAtomCount < mol.getAtoms());
			mFragmentAtomFromOrig = new int[maxAtomCount];
		}

		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		setMolecule(mol, bondMember, bondCount);
		}

	private void setMolecule(StereoMolecule mol, int[] bondMember, int bondCount) {
		boolean[] includeAtom = new boolean[mol.getAtoms()];
		Arrays.fill(mConnAtoms, 0);
		for (int i=0; i<bondCount; i++) {
			int bond = bondMember[i];
			for (int j=0; j<2; j++) {
				int atom = mol.getBondAtom(j, bond);
				if (!includeAtom[atom]) {
					includeAtom[atom] = true;
					mFragmentAtomFromOrig[atom] = mAtoms;
					mAtomDescriptor[mAtoms++] = createAtomDescriptor(mol, atom);
				}
			}

			for (int j=0; j<2; j++) {
				int atom = mol.getBondAtom(j, bond);
				int fragmentAtom = mFragmentAtomFromOrig[atom];
				mConnAtom[fragmentAtom][mConnAtoms[fragmentAtom]] = mFragmentAtomFromOrig[mol.getBondAtom(1-j, bond)];
				mConnBond[fragmentAtom][mConnAtoms[fragmentAtom]] = mBonds;
				mConnAtoms[fragmentAtom]++;
			}

			mBondDescriptor[mBonds++] = createBondDescriptor(mol, bond);
		}
	}

	public int createHashValue(int hashBits) {
		long[] baseValue = new long[mAtoms];
		for (int i=0; i<mAtoms; i++)
			baseValue[i] = (mAtomDescriptor[i] << ATOM_INDEX_BITS) + i;

		// initial ranking based on atom descriptors
		Arrays.sort(baseValue);
		int startRanks = 1;
		for (int atom=0; atom<mAtoms; atom++) {
			if (atom != 0 && ((baseValue[atom] ^ baseValue[atom-1]) & ~ATOM_INDEX_MASK) != 0)
				startRanks++;
			mCanRank[(int)(baseValue[atom] & ATOM_INDEX_MASK)] = startRanks;
		}

		while (true) {
			// calculate next base values with piggybacked atom index
			for (int atom=0; atom<mAtoms; atom++) {
				int[] connRank = mConnRank[mConnAtoms[atom]];
				for (int i = 0; i<mConnAtoms[atom]; i++)
					connRank[i] = (mBondDescriptor[mConnBond[atom][i]] << ATOM_INDEX_BITS) | mCanRank[mConnAtom[atom][i]];
				Arrays.sort(connRank);

				baseValue[atom] = mCanRank[atom];
				for (int i = 0; i<mConnAtoms[atom]; i++) {
					baseValue[atom] = (baseValue[atom] << (ATOM_INDEX_BITS + BOND_DESC_BITS)) | connRank[i];
				}

				baseValue[atom] = (baseValue[atom] << ATOM_INDEX_BITS) | atom;
			}

			// sort atom indexes by sorted neighbour rank and connection bond types
			Arrays.sort(baseValue);

			// assign new rank based on sort order
			int ranks = 1;
			for (int atom=0; atom<mAtoms; atom++) {
				if (atom != 0 && ((baseValue[atom] ^ baseValue[atom-1]) & ~ATOM_INDEX_MASK) != 0)
					ranks++;
				mCanRank[(int)(baseValue[atom] & ATOM_INDEX_MASK)] = ranks;
			}

			if (ranks == startRanks)
				break;

			startRanks = ranks;
		}

		int maxRank = -1;
		int maxAtom = -1;
		for (int atom=0; atom<mAtoms; atom++) {
			if (maxRank < mCanRank[atom]) {
				maxRank = mCanRank[atom];
				maxAtom = atom;
			}
		}

		boolean[] isUsedAtom = new boolean[mAtoms];
		boolean[] isUsedBond = new boolean[mBonds];
		mGraphAtom[0] = maxAtom;
		mGraphIndex[maxAtom] = 0;
		isUsedAtom[maxAtom] = true;
		int closureBonds = mBonds - mAtoms + 1;
		int bufferAtomIndex = 1;
		int bufferBondIndex = mAtoms;
		int bufferClosureIndex = mAtoms + mBonds - closureBonds;
		mBuffer[0] = (byte)((mAtomDescriptor[maxAtom] << ATOM_INDEX_BITS) | closureBonds);
		int maxBond = -1;
		int current = 0;
		int highest = 0;
		while (current <= highest) {
			int parent = mGraphAtom[current];
			maxRank = -1;
			maxAtom = -1;
			for (int i=0; i<mConnAtoms[parent]; i++) {
				int bond = mConnBond[parent][i];
				if (!isUsedBond[bond]) {
					int atom = mConnAtom[parent][i];
					if (maxRank < mCanRank[atom]) {
						maxRank = mCanRank[atom];
						maxAtom = atom;
						maxBond = bond;
					}
				}
			}
			if (maxRank != -1) {
				if (isUsedAtom[maxAtom]) {
					mBuffer[bufferClosureIndex++] = (byte)((mGraphIndex[parent] << ATOM_INDEX_BITS) | mGraphIndex[maxAtom]);
					mBuffer[bufferClosureIndex++] = mBondDescriptor[maxBond];
				}
				else {
					mGraphIndex[maxAtom] = ++highest;
					mGraphAtom[highest] = maxAtom;
					isUsedAtom[maxAtom] = true;
					mBuffer[bufferAtomIndex++] = mAtomDescriptor[maxAtom];
					mBuffer[bufferBondIndex++] = (byte)((mGraphIndex[parent] << BOND_DESC_BITS) | mBondDescriptor[maxBond]);
				}
				isUsedBond[maxBond] = true;
				continue;
			}
			current++;
		}

		int hash = BurtleHasher.hashlittle(mBuffer, HASH_INIT, mAtoms + mBonds + (mBonds - mAtoms + 1));
		return (hash & BurtleHasher.hashmask(hashBits));
	}

	private byte createAtomDescriptor(StereoMolecule mol, int atom) {
		return (byte)mol.getAtomicNo(atom);
	}

	private byte createBondDescriptor(StereoMolecule mol, int bond) {
		// to stay in 2 bits we pool dative with triple bonds
		if (mol.isDelocalizedBond(bond))
			return 0;
		int order = Math.min(3, mol.getBondOrder(bond));
		return (byte)(order == 0 ? 3 : order);
	}
}
