package com.actelion.research.chem.sar;

public class ExitVector {
	private int mAtom,mIndex,mSeenBondOrders,mTopicity;
	private boolean mAtomIsQueryAtom,mEmptySubstituentSeen,mSubstituentVaries;
	private String mConstantSubstituent;
	private int mRGroupNo;

	/**
	 * @param atom index of canonical core structure (if is bond bridge atom) or of original query structure
	 * @param atomIsQueryAtom
	 * @param index number to distinguish multiple exit vectors at same core atom (0-based)
	 * @param topicity -1 (no topicity) or 0 or 1
	 */
	public ExitVector(int atom, boolean atomIsQueryAtom, int index, int topicity) {
		mAtom = atom;
		mAtomIsQueryAtom = atomIsQueryAtom;
		mIndex = index;
		mTopicity = topicity;
	}

	public int getQueryAtom() {
		return mAtomIsQueryAtom ? mAtom : -1;
	}

	/**
	 * @param queryToCoreAtom
	 * @return the root atom index for the core structure
	 */
	public int getCoreAtom(int[] queryToCoreAtom) {
		return mAtomIsQueryAtom ? queryToCoreAtom[mAtom] : mAtom;
	}

	public int getIndex() {
		return mIndex;
	}

	public int getTopicity() {
		return mTopicity;
	}

	public int getRGroupNo() {
		return mRGroupNo;
	}

	public void setRGroupNo(int no) {
		mRGroupNo = no;
	}

	public boolean substituentVaries() {
		return mSubstituentVaries;
	}

	public String getConstantSubstituent() {
		return mConstantSubstituent;
	}

	public boolean hasSeenBondOrder(int bondOrder) {
		return (mSeenBondOrders & (1 << bondOrder)) != 0;
	}

	/**
	 * This method is called with all substituents at this exit vector position.
	 * After being called with all substituents from all molecules within the same scaffold
	 * or scaffols group, it accumulated the knowledge, whether this exit vector is always
	 * unsubstituted, carries always the same substituent, or carries multiple different
	 * substituents.
	 * @param substituent
	 */
	protected void checkSubstituent(String substituent, int bondOrder) {
		mSeenBondOrders |= (1 << bondOrder);
		if (!mSubstituentVaries) {
			if (substituent == null) {
				mEmptySubstituentSeen = true;
				if (mConstantSubstituent != null)
					mSubstituentVaries = true;
			}
			else {
				if (mEmptySubstituentSeen)
					mSubstituentVaries = true;
				else if (mConstantSubstituent == null)
					mConstantSubstituent = substituent;
				else if (!mConstantSubstituent.equals(substituent))
					mSubstituentVaries = true;
			}
		}
	}
}
