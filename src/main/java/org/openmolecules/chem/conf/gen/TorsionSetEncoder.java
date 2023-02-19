package org.openmolecules.chem.conf.gen;

import com.actelion.research.util.DoubleFormat;

public class TorsionSetEncoder {
	private static final long[] BITS = { 0x00L, 0x01L, 0x03L, 0x07L, 0x0fL, 0x1fL, 0x3fL, 0x7fL };

	private RigidFragment[] mRigidFragment;
	private RotatableBond[] mRotatableBond;
	private int mEncodingLongCount;
	private int[] mEncodingBitCount;
	private int[] mEncodingBitShift;
	private int[] mEncodingLongIndex;

	/**
	 * Creates a new encoder for TorsionSets (torsion and(!) rigid fragment indexes).
	 * The encoder is also used to encode elimination rules, which cover torsion indexes,
	 * but no rigid fragment indexes. Encoding of TorsionSets and EliminationRules
	 * facilitates fast matching of TorsionSets and sub-matching of EliminationRules.
	 * @param rigidFragment
	 * @param rotatableBond
	 */
	public TorsionSetEncoder(RigidFragment[] rigidFragment, RotatableBond[] rotatableBond) {
		mRigidFragment = rigidFragment;
		mRotatableBond = rotatableBond;
		mEncodingBitCount = new int[rotatableBond.length+rigidFragment.length];
		mEncodingBitShift = new int[rotatableBond.length+rigidFragment.length];
		mEncodingLongIndex = new int[rotatableBond.length+rigidFragment.length];
		int bitCount = 0;
		int longIndex = 0;
		int index = 0;
		for (RotatableBond rb:rotatableBond) {
			int bits = neededBits(rb.getTorsionCount());
			if (bitCount + bits <= 64) {
				mEncodingBitShift[index] = bitCount;
				bitCount += bits;
			}
			else {
				longIndex++;
				mEncodingBitShift[index] = 0;
				bitCount = 0;
			}
			mEncodingBitCount[index] = bits;
			mEncodingLongIndex[index] = longIndex;
			index++;
		}
		for (RigidFragment rf:rigidFragment) {
			int bits = neededBits(rf.getConformerCount());
			if (bitCount + bits <= 64) {
				mEncodingBitShift[index] = bitCount;
				bitCount += bits;
			}
			else {
				longIndex++;
				mEncodingBitShift[index] = 0;
				bitCount = 0;
			}
			mEncodingBitCount[index] = bits;
			mEncodingLongIndex[index] = longIndex;
			index++;
		}
		mEncodingLongCount = longIndex+1;
	}

	public long[] encode(int[] torsionIndex, int[] conformerIndex) {
		long[] mEncodedBits = new long[1 + mEncodingLongIndex[mEncodingLongIndex.length - 1]];

		int i = 0;
		for (int index:torsionIndex) {
			mEncodedBits[mEncodingLongIndex[i]] += (index << mEncodingBitShift[i]);
			i++;
		}
		for (int index:conformerIndex) {
			mEncodedBits[mEncodingLongIndex[i]] += (index << mEncodingBitShift[i]);
			i++;
		}

		return mEncodedBits;
	}

	public void encodeRule(int[] torsionIndex, int[] rotatableBondIndex, long[] mask, long[] data) {
		for (int rb:rotatableBondIndex) {
			data[mEncodingLongIndex[rb]] += (torsionIndex[rb] << mEncodingBitShift[rb]);
			mask[mEncodingLongIndex[rb]] += (BITS[mEncodingBitCount[rb]] << mEncodingBitShift[rb]);
		}
	}

	public int getEncodingLongCount() {
		return mEncodingLongCount;
	}

	public boolean isMaskSet(TorsionSetEliminationRule rule, int torsionAndConformerIndex) {
		return (rule.getMask()[mEncodingLongIndex[torsionAndConformerIndex]]
				& (1L << mEncodingBitShift[torsionAndConformerIndex])) != 0L;
	}

	private int[] decodeRule(TorsionSetEliminationRule rule) {
		int[] data = new int[mRotatableBond.length];
		for (int rb=0; rb<mRotatableBond.length; rb++) {
			long mask = BITS[mEncodingBitCount[rb]] << mEncodingBitShift[rb];
			if ((rule.getMask()[mEncodingLongIndex[rb]] & mask) != 0)
				data[rb] = (int)((rule.getData()[mEncodingLongIndex[rb]] & mask) >> mEncodingBitShift[rb]);
			else
				data[rb] = -1;
			}
		return data;
		}

	public String createRuleString(TorsionSetEliminationRule rule, BaseConformer baseConformer) {
		StringBuilder sb = new StringBuilder();
		int[] data = decodeRule(rule);

		sb.append(" rb:");

		int sl = sb.length();
		for (int i=0; i<data.length; i++) {
			if (data[i] != -1) {
				if (sb.length() != sl)
					sb.append(',');
				sb.append(i);
			}
		}
		sb.append(" ti:");

		sl = sb.length();
		for (int i=0; i<data.length; i++) {
			if (data[i] != -1) {
				if (sb.length() != sl)
					sb.append(',');
				sb.append(data[i]);
			}
		}
		sb.append(" (");

		sl = sb.length();
		for (int i=0; i<data.length; i++) {
			if (data[i] != -1) {
				if (sb.length() != sl)
					sb.append(',');
				sb.append(baseConformer.getTorsion(i, data[i]));
			}
		}
		sb.append(") "+ DoubleFormat.toString(rule.getCollisionIntensity(), 3));
		return sb.toString();
	}

	private int neededBits(int count) {
		int bits = 0;
		int maxIndex = count-1;
		while (maxIndex > 0) {
			maxIndex >>= 1;
			bits++;
		}
		return bits;
	}
}
