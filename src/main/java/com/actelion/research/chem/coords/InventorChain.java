package com.actelion.research.chem.coords;

/**
 * Created by thomas on 9/23/16.
 */
public class InventorChain {
	protected int[] mAtom;
	protected int[] mBond;

	public InventorChain(int chainLength) {
		mAtom = new int[chainLength];
		mBond = new int[chainLength];
		// last mBond[] value isn't needed if chain isn't a ring
	}

	protected int getChainLength() {
		return mAtom.length;
	}

	protected int[] getRingAtoms() {
		return mAtom;
	}

	protected int[] getRingBonds() {
		return mBond;
	}
}
