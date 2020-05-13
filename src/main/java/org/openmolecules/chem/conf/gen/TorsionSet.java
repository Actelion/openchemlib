/*
 * Copyright 2013-2020 Thomas Sander, openmolecules.org
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.gen;

import java.util.Arrays;

public class TorsionSet implements Comparable<TorsionSet> {
	private long[] mEncodedBits;
	private double mLikelihood;
	private int[] mTorsionIndex, mConformerIndex;
	private double mCollisionIntensitySum;
	private double[][] mCollisionIntensityMatrix;
	private boolean mIsUsed;

//	private int[] bitshift;
//	private int[] longIndex;

	/**
	 * Creates a new conformer description from torsion and conformer indexes.
	 *
	 * @param torsionIndex   torsion angle index for all rotatable bonds
	 * @param conformerIndex conformer index for every rigid fragment
	 * @param bitShift       bit position of torsion and conformer indexes
	 * @param longIndex      index on long array for shifted torsion and conformer indexes
	 * @param likelihood     all individual index likelyhoods multiplied
	 */
	public TorsionSet(int[] torsionIndex, int[] conformerIndex, int[] bitShift, int[] longIndex, double likelihood) {
		mTorsionIndex = torsionIndex;
		mConformerIndex = conformerIndex;
		mLikelihood = likelihood;
		mEncodedBits = new long[1 + longIndex[longIndex.length - 1]];

//		this.bitshift = bitShift;
//		this.longIndex = longIndex;

		int i = 0;
		for (int index : torsionIndex) {
			mEncodedBits[longIndex[i]] += (index << bitShift[i]);
			i++;
		}
		for (int index : conformerIndex) {
			mEncodedBits[longIndex[i]] += (index << bitShift[i]);
			i++;
		}
	}

	/**
	 * Deep-Copy constructor not including collision intensities
	 * @param ref
	 */
	public TorsionSet( TorsionSet ref ) {
//		this( Arrays.copyOf(ref.getTorsionIndexes(),ref.getTorsionIndexes().length) ,
//				Arrays.copyOf( ref.getConformerIndexes() , ref.getConformerIndexes().length  ) ,
//				Arrays.copyOf( ref.getBitshift() , ref.getBitshift().length ),
//				Arrays.copyOf( ref.getLongIndex() , ref.getLongIndex().length ) ,
//				ref.getLikelihood());
		mTorsionIndex = Arrays.copyOf(ref.mTorsionIndex, ref.mTorsionIndex.length);
		mConformerIndex = Arrays.copyOf(ref.mConformerIndex, ref.mConformerIndex.length);
		mLikelihood = ref.mLikelihood;
		mEncodedBits = Arrays.copyOf(ref.mEncodedBits, ref.mEncodedBits.length);
	}

	public int[] getTorsionIndexes() {
		return mTorsionIndex;
	}

	public int[] getConformerIndexes() {
		return mConformerIndex;
	}

	public double getLikelihood() {
		return mLikelihood;
	}

//	public int[] getBitshift() {return this.bitshift;}
//	public int[] getLongIndex() {return this.longIndex;}

	public double getCollisionIntensitySum() {
		return mCollisionIntensitySum;
	}

	public double[][] getCollisionIntensityMatrix() {
		return mCollisionIntensityMatrix;
	}

	public void setCollisionIntensity(double sum, double[][] matrix) {
		mCollisionIntensitySum = sum;
		mCollisionIntensityMatrix = matrix;
	}

	public boolean isUsed() {
		return mIsUsed;
	}

	public void setUsed() {
		mIsUsed = true;
	}

	/**
	 * Checks whether all torsion indexes of a subset of all rotatable bonds
	 * match, when comparing this TorsionSet with the torsion indexes contained
	 * in data, considering only those rotatable bonds defined by mask.
	 *
	 * @param mask
	 * @param data
	 * @return
	 */
	private boolean matches(long[] mask, long[] data) {
		for (int i = 0; i < mask.length; i++)
			if ((mEncodedBits[i] & mask[i]) != (data[i] & mask[i]))
				return false;

//		System.out.println("Eliminated mask:"+Long.toHexString(mask[0])+" rule:"+Long.toHexString(data[0])+" tSet:"+Long.toHexString(mEncodedBits[0]));
		return true;
	}

	public boolean matches(TorsionSetEliminationRule er, double tolerance) {
		return (mCollisionIntensitySum > tolerance) && matches(er.getMask(), er.getData());
	}

	/**
	 * Allows to order TorsionSets in a unique way for quick uniqueness checking against a TreeSet.
	 */
	@Override
	public int compareTo(TorsionSet ts) {
		for (int i = 0; i < mEncodedBits.length; i++)
			if (mEncodedBits[i] != ts.mEncodedBits[i])
				return (mEncodedBits[i] < ts.mEncodedBits[i]) ? -1 : 1;

		return 0;
	}

	@Override
	public boolean equals(Object ts) {
		return compareTo((TorsionSet) ts) == 0;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i : mTorsionIndex) {
			sb.append(i);
			sb.append(',');
		}
		if (mTorsionIndex.length != 0)
			sb.setLength(sb.length() - 1);
		sb.append(';');
		for (int i : mConformerIndex) {
			sb.append(i);
			sb.append(',');
		}
		if (mConformerIndex.length != 0)
			sb.setLength(sb.length() - 1);
		return sb.toString();
	}

/*
	public static List<int[][]> decodeTorsionSets(String si) {
		byte[] ea = SimpleEncoder.convertS3IDCodeToByteArray(si);
		MyByteArray a = new MyByteArray(ea, true, 32000);

		List<int[][]> decoded = new ArrayList<>();

		List<List<Byte>> splits = a.splitAtSpecificByte((byte) -2);
		for (int zi = 0; zi < splits.size(); zi++) {
			List<Byte> split = splits.get(zi);
			MyByteArray b = new MyByteArray(split);
			List<List<Byte>> bsplit = b.splitAtSpecificByte((byte) -1);
			int confis[] = new int[bsplit.get(0).size()];
			for (int zc = 0; zc < bsplit.get(0).size(); zc++) {
				confis[zc] = bsplit.get(0).get(zc);
			}
			int rotis[] = new int[bsplit.get(1).size()];
			for (int zc = 0; zc < bsplit.get(1).size(); zc++) {
				rotis[zc] = bsplit.get(1).get(zc);
			}
			decoded.add(new int[][]{confis, rotis});
		}

		return decoded;
	}
*/
	/**
	 * Encodes the conformer index array and the torsion index array of each conformer.
	 * The resulting string can be decoded with the decodeTorsionSets function.
	 * <p>
	 * In the case that we just have a single conformer of the
	 *
	 * @param list_t
	 * @return
	 *
	public static String encodeTorsionSets(List<TorsionSet> list_t) {
		List<Byte> encoded = new ArrayList<>();

		for (int zt = 0; zt < list_t.size(); zt++) {
			TorsionSet ti = list_t.get(zt);

			int ci[] = null;
			int ri[] = null;

			if (ti == null) {
				// write dummy values, they will not hurt.
				ci = new int[]{0};
				ri = new int[]{0};
			} else {
				ci = ti.getConformerIndexes();
				ri = ti.getTorsionIndexes();
			}


			for (int zi = 0; zi < ci.length; zi++) {
				encoded.add((byte) (ci[zi] % 255));
			}
			encoded.add((byte) -1);
			for (int zi = 0; zi < ri.length; zi++) {
				encoded.add((byte) (ri[zi] % 255));
			}
			if (zt < list_t.size() - 1) {
				encoded.add((byte) -2);
			}
		}

		byte ea[] = new byte[encoded.size()];
		for (int zi = 0; zi < ea.length; zi++) {
			ea[zi] = encoded.get(zi);
		}

		return SimpleEncoder.convertByteArrayToS3IDCode(MyByteArray.deflate(ea));
	}
*/
}