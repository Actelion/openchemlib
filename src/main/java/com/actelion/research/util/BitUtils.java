package com.actelion.research.util;

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

public class BitUtils {
	
    public final static long MASK_FIRST_SHORT  = 0x000000000000FFFFl;
    public final static long MASK_SEC_SHORT    = 0x00000000FFFF0000l;
    public final static long MASK_THIRD_SHORT  = 0x0000FFFF00000000l;
    public final static long MASK_FOURTH_SHORT = 0xFFFF000000000000l;

	private static int[] BIT_COUNTS = new int[65536];

	static {
		for (int i = 0; i < BIT_COUNTS.length; i++) {
			BIT_COUNTS[i] = countBitsSlow(i);
		}
	}


	/** Count the number of set bits in an int;
	 *  @returns the number of bits set in x
	 */
	static int countBitsSlow(int x) {
		int temp;

		temp = 0x55555555;
		x = (x & temp) + (x >>> 1 & temp);
		temp = 0x33333333;
		x = (x & temp) + (x >>> 2 & temp);
		temp = 0x07070707;
		x = (x & temp) + (x >>> 4 & temp);
		temp = 0x000F000F;
		x = (x & temp) + (x >>> 8 & temp);

		return (x & 0x1F) + (x >>> 16);
	}

	/**
	 *
	 * The fastest Method (by Thomas Sander)
	 * @param x
	 * @return
	 * @deprecated use Integer.bitCount() instead.
	 */
	public static int bitCount(int x) {

		int t1 = (0xFFFF0000 & x) >>> 16;
		int t2 = 0x0000FFFF & x;
		return BIT_COUNTS[t1] + BIT_COUNTS[t2];


	}
	
	public static int bitCount(long x) {
		int t1 = (int)((MASK_FOURTH_SHORT & x) >>> 48);
		int t2 = (int)((MASK_THIRD_SHORT & x) >>> 32);
		int t3 = (int)((MASK_SEC_SHORT & x) >>> 16);
		int t4 = (int)((MASK_FIRST_SHORT & x));
		return BIT_COUNTS[t1] + BIT_COUNTS[t2] + BIT_COUNTS[t3] + BIT_COUNTS[t4];
	}

	public static void setBit(int [] data, int i) {
		int ind = i / Integer.SIZE;
		int indInInt = i % Integer.SIZE;
		int mask = 1;
		mask = mask << indInInt;
		data[ind] = data[ind] | mask;
	}

	public static void unsetBit(int [] data, int i) {
		int ind = i / Integer.SIZE;
		int indInInt = i % Integer.SIZE;
		int mask = 1;
		mask = mask << indInInt;
		data[ind] = data[ind] & ~mask;
	}

	public static boolean isBitSet(int [] a, int i) {
		int ind = i / Integer.SIZE;
		int indInInt = i % Integer.SIZE;
		int mask = 1;
		mask = mask << indInInt;
		if((a[ind] & mask) != 0)
			return true;
		else
			return false;
	}

	public static boolean isValidBitIndex(int [] data, int i) {

		boolean valid = false;

		int ind = (i / Integer.SIZE);

		if(ind < data.length){
			valid=true;
		}

		return valid;
	}

}