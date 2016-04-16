package com.actelion.research.util;

import java.math.BigInteger;

/**
 * <p>Title: </p>
 * <p>Description: BitUtils </p>
 * @author Tim Tyler, Modest von Korff
 * http://mandala.co.uk/java/bitutils/BitUtils.java
 * 05.04.2005 found by MvK
 */

public class BitUtils {
	
    public final static long MASK_FIRST_SHORT  = (new BigInteger("000000000000FFFF", 16)).longValue();
    public final static long MASK_SEC_SHORT    = (new BigInteger("00000000FFFF0000", 16)).longValue();
    public final static long MASK_THIRD_SHORT  = (new BigInteger("0000FFFF00000000", 16)).longValue();
    public final static long MASK_FOURTH_SHORT = (new BigInteger("FFFF000000000000", 16)).longValue();

	private static int[] BIT_COUNTS = new int[65536];

	static {
		for (int i = 0; i < BIT_COUNTS.length; i++) {
			BIT_COUNTS[i] = bitCountSlow(i);
		}
	}


	/** Count the number of set bits in an int;
	 *  @param x the int to have its bits counted
	 *  @author Tim Tyler tt@iname.com
	 *  @returns the number of bits set in x
	 */
	static int bitCountSlow(int x) {
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
	 * The fastest Method (by Thomas Sander)
	 * @param x
	 * @return
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

	
	/** Count the number of set bits in an long;
	 *  @param x the long to have its bits counted
	 *  @author Tim Tyler tt@iname.com
	 *  @returns the number of bits set in x
	 */
	static int bitCountSlow(long x) {
		long temp;

		temp = 0x5555555555555555L;
		x = (x & temp) + (x >>> 1 & temp);
		temp = 0x3333333333333333L;
		x = (x & temp) + (x >>> 2 & temp);
		temp = 0x0707070707070707L;
		x = (x & temp) + (x >>> 4 & temp);
		temp = 0x000F000F000F000FL;
		x = (x & temp) + (x >>> 8 & temp);
		temp = 0x0000001F0000001FL;
		x = (x & temp) + (x >>> 16 & temp);

		return (int) ((x & 0x3f) + (x >>> 32));
	}


	// main method (timings)
	public static void main(String args[]) {

		testCorrectness();
		System.out.println("Finished");

	}

	public static void testCorrectness() {
		int N = 50000000;

		int a = 3;
		a = a << 30;

		
		int t1 = bitCountSlow(a);
		int t2 = bitCount(a);
		
		if (t1 != t2) {
			System.err.println("t1 = " + t1 + " and t2 " + t2);
		}

		long b = 1;
		b = b << 63;
		int t3 = bitCount(b);
		int t4 = bitCountSlow(b);
		if (t3 != t4) {
			System.err.println("t3 = " + t3 + " and t4 " + t4);
		}

		
		
		for (int x = 0; x < N; x++) {
			t1 = bitCountSlow(x);
			t2 = bitCount(x);
			if (t1 != t2) {
				System.err.println("t1 = " + t1 + " and t2 " + t2);
			}
		}

		System.out.println("Finished");

	}

	public static void testPerformance() {
		
		long start_time;
		int N = 50000000;

		System.out.println("Start");

		do {

			start_time = System.currentTimeMillis();
			for (int x = 0; x < N; x++) {
				bitCountSlow(x);
			}
			System.out.println("bitCountSlow:"
					+ (System.currentTimeMillis() - start_time));

			start_time = System.currentTimeMillis();
			for (int x = 0; x < N; x++) {
				bitCount(x);
			}
			System.out.println("bitCount:"
					+ (System.currentTimeMillis() - start_time));

			
		} while (1 == 1);
		
		
	}
}