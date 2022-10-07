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
 */

package com.actelion.research.chem.descriptor;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.SortedList;

public class DescriptorHandlerAllFragmentsFP extends AbstractDescriptorHandlerLongFP<StereoMolecule> {
	private static final double CORRECTION_FACTOR = 0.68;

	public static boolean skipIDCodes;

	private static DescriptorHandlerAllFragmentsFP sDefaultInstance;

	private static final int MAX_BOND_COUNT = 6;
	private static final int HASH_BITS = 11;
	private static final int DESCRIPTOR_SIZE = (1 << HASH_BITS);
	private static final int DESCRIPTOR_LEN = DESCRIPTOR_SIZE / Long.SIZE;

	private StereoMolecule mMol;
	private boolean[] mIsAtomMember,mIsBondMember;
	private int[] mMemberBond,mMemberAtom;
	private long[] mDescriptor;
	private SortedList<BondSet> mBondSets;
	private SimpleFragmentGraph mFragmentGraph;

	public static DescriptorHandlerAllFragmentsFP getDefaultInstance() {
		synchronized(DescriptorHandlerAllFragmentsFP.class) {
			if (sDefaultInstance == null) {
				sDefaultInstance = new DescriptorHandlerAllFragmentsFP();
			}
		}
		return sDefaultInstance;
	}

	@Override
	public DescriptorInfo getInfo() {
		return DescriptorConstants.DESCRIPTOR_ALLFRAG;
	}

	@Override
	public String getVersion() {
		return DescriptorConstants.DESCRIPTOR_ALLFRAG.version;
	}

	private StereoMolecule preprocessStructure(StereoMolecule mol) {
		if (mol.isFragment()) {
			mol.ensureHelperArrays(Molecule.cHelperNeighbours);

			boolean[] isBlockedAtom = new boolean[mol.getAtoms()];
			for (int atom = 0; atom<mol.getAtoms(); atom++)
				isBlockedAtom[atom] = (mol.getAtomQueryFeatures(atom) & ~Molecule.cAtomQFNarrowing) != 0
									|| mol.getAtomList(atom) != null;

			boolean hasBlockedBond = false;
			boolean[] includeBond = new boolean[mol.getBonds()];
			for (int bond = 0; bond<mol.getBonds(); bond++) {
				includeBond[bond] = !isBlockedAtom[mol.getBondAtom(0, bond)]
						&& !isBlockedAtom[mol.getBondAtom(1, bond)]
						&& (mol.getBondQueryFeatures(bond) & ~Molecule.cBondQFNarrowing) == 0;
				hasBlockedBond |= !includeBond[bond];
			}

			if (hasBlockedBond) {
				StereoMolecule query = new StereoMolecule(mol.getAllBonds(), mol.getAllBonds());
				mol.copyMoleculeByBonds(query, includeBond, true, null);
				mol = query;
			}
		}

		return mol;
	}

	/**
	 * This descriptor requires proper up/down bonds, because it encodes stereo parities.
	 * If a passed molecule is generated from idcode parsing, make sure that coordinates
	 * and up/down/bonds are available, i.e. that the IDCodeParser was instantiated with
	 * the respective option.
	 */
	@Override
	public long[] createDescriptor(StereoMolecule mol) {
		if (mol ==null)
			return null;

//System.out.println("bit\tbonds\tidcode");

		mMol = preprocessStructure(mol);
		mMol.ensureHelperArrays(Molecule.cHelperRings);
		mFragmentGraph = new SimpleFragmentGraph(MAX_BOND_COUNT);

		mBondSets = new SortedList<>();
		mIsAtomMember = new boolean[mMol.getAtoms()];
		mIsBondMember = new boolean[mMol.getBonds()];
		mDescriptor = new long[DESCRIPTOR_LEN];
		mMemberAtom = new int[mMol.getAtoms()];
		mMemberBond = new int[mMol.getBonds()];

		for (int rootBond=0; rootBond<mMol.getBonds(); rootBond++) {
			mMemberAtom[0] = mMol.getBondAtom(0, rootBond);
			mMemberAtom[1] = mMol.getBondAtom(1, rootBond);
			mMemberBond[0] = rootBond;

			mIsAtomMember[mMemberAtom[0]] = true;
			mIsAtomMember[mMemberAtom[1]] = true;
			mIsBondMember[rootBond] = true;

			setHashBitIfNew(1);
			processOneMoreBond(2, 1);

			mIsAtomMember[mMemberAtom[0]] = false;
			mIsAtomMember[mMemberAtom[1]] = false;
			mIsBondMember[rootBond] = false;
			}

		// add bits for isolated atoms
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getConnAtoms(atom) == 0) {
				int atomicNo = mMol.getAtomicNo(atom);
				mDescriptor[atomicNo < 64 ? 0 : 1] |= (1L << (atomicNo % 64));
			}
		}

		return mDescriptor;
	}

	public void processOneMoreBond(int atomCount, int bondCount) {
		for (int i=0; i<atomCount; i++) {
			for (int j=0; j<mMol.getConnAtoms(mMemberAtom[i]); j++) {
				int candidateBond = mMol.getConnBond(mMemberAtom[i], j);
				if (!mIsBondMember[candidateBond]) {
					int candidateAtom = mMol.getConnAtom(mMemberAtom[i], j);

					mMemberBond[bondCount] = candidateBond;
					mIsBondMember[candidateBond] = true;
					boolean isAtomMember = mIsAtomMember[candidateAtom];
					if (!isAtomMember) {
						mIsAtomMember[candidateAtom] = true;
						mMemberAtom[atomCount] = candidateAtom;
						atomCount++;
						}
					bondCount++;

					setHashBitIfNew(bondCount);

					// Recursively add a connected bond to the current bond set.
					if (bondCount < MAX_BOND_COUNT)
						processOneMoreBond(atomCount, bondCount);

					bondCount--;
					mIsBondMember[candidateBond] = false;
					if (!isAtomMember) {
						atomCount--;
						mIsAtomMember[candidateAtom] = false;
					}
				}
			}
		}
	}

	/**
	 * If we have a new bond set, then create a hash code from the canonical structure and set the respective descriptor bit.
	 */
	private void setHashBitIfNew(int bondCount) {
		BondSet bondSet = new BondSet(mMemberBond, bondCount, mMol.getBonds());
		if (mBondSets.addIfNew(bondSet)) {
			mFragmentGraph.set(mMol, mMemberBond, bondCount);
			int hash = mFragmentGraph.createHashValue(HASH_BITS);
/*
StereoMolecule frag = new StereoMolecule();
boolean[] isBondMember = new boolean[mMol.getBonds()];
for (int i=0; i<bondCount; i++)
 isBondMember[mMemberBond[i]] = true;
mMol.copyMoleculeByBonds(frag, isBondMember, true, null);
System.out.println(hash+"\t"+bondCount+"\t"+new Canonizer(frag).getIDCode());
*/
			int high = DESCRIPTOR_LEN - hash / Long.SIZE - 1;
			int low = hash % Long.SIZE;
			mDescriptor[high] |= (1L << low);
		}
	}

	@Override
	public float getSimilarity(long[] o1, long[] o2) {
		return o1 == null
				|| o2 == null
				|| o1.length == 0
				|| o2.length == 0 ? 0.0f
				: normalizeValue(SSSearcherWithIndex.getSimilarityTanimoto(o1, o2));
	}

	private float normalizeValue(double value) {
		return value <= 0.0f ? 0.0f
				: value >= 1.0f ? 1.0f
				: (float)(1.0-Math.pow(1-Math.pow(value, CORRECTION_FACTOR) ,1.0/CORRECTION_FACTOR));
	}

	@Override
	public DescriptorHandler<long[], StereoMolecule> getThreadSafeCopy() {
		return new DescriptorHandlerAllFragmentsFP();
	}
}

class BondSet implements Comparable<BondSet> {
	private int[] mMask;

	public BondSet(int[] bond, int usedBonds, int bondCount) {
		mMask = new int[(bondCount + 31) / 32];
		for (int i=0; i<usedBonds; i++)
			mMask[bond[i] / 32] |= (1 << (bond[i] % 32));
	}

	public boolean equals(BondSet bs) {
		for (int i=0; i<mMask.length; i++)
			if (bs.mMask[i] != mMask[i])
				return false;
		return true;
	}

	@Override
	public int compareTo(BondSet bs) {
		for (int i=0; i<mMask.length; i++)
			if (bs.mMask[i] != mMask[i])
				return bs.mMask[i] < mMask[i] ? -1 : 1;
		return 0;
	}
}
