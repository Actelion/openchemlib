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

package com.actelion.research.util;

public class DoubleFormat {
	private static String[] ZEROS = {"","0","00","000","0000","00000","000000","0000000","00000000"};
	/**
	 * Converts a double value into a String representation in scientific format
	 * rounded to 5 significant digits, e.g. 1.2345e-6.
	 * Double values that are effectively integers are turned into the scientific
	 * format only, if they consist of more than 8 digits.
	 * @param theDouble
	 * @return
	 */
	public static String toString(double theDouble) {
		return toString(theDouble, 5, true);
		}

	/**
	 * Converts a double value into a String representation in scientific format
	 * rounded to a definable number of digits, e.g. 1.2345e-6.
	 * Double values that are effectively integers are turned into the scientific
	 * format only, if they consist of more than significantDigits+3 digits.
	 * @param value
	 * @param significantDigits
	 * @return
	 */
	public static String toString(double value, int significantDigits) {
		return toString(value, significantDigits, true);
		}

	/**
	 * Converts a double value into a String representation in scientific format
	 * rounded to a definable number of digits, e.g. 1.2345e-6.
	 * Double values that are effectively integers are turned into the scientific
	 * format only, if they consist of more than significantDigits+3 digits.
	 * If skipTrailingZeros is true, then integer output is not rounded: 7173 instead of 7200.
	 * @param value
	 * @param significantDigits
	 * @param skipTrailingZeros if true then trailing zeros after a decimal point are not shown and integers are not rounded
	 * @return
	 */
	public static String toString(double value, int significantDigits, boolean skipTrailingZeros) {
		if (Double.isNaN(value))
			return "NaN";
		if (Double.isInfinite(value))
			return "Infinity";

		if (value == 0.0)
			return (skipTrailingZeros || significantDigits==1) ? "0" : "0."+zeros(significantDigits-1);

		// Determine cipher count of integer fraction in excess of significant digits.
		// For those cases, where the number is not expressed in scientific notation,
		// the increase significant digits to prevent trailing zeros through rounding.
		//		int cipherExcess = Long.toString(Math.round(value)).length() - significantDigits;
		//		if (cipherExcess > 0 && cipherExcess <= 4)
		//			significantDigits += cipherExcess;
		// Returned to original behaviour with consistent rounding no matter whether we see trailing zeros
		// in integer part that pretend precision with digits (zeros). TLS 21Jun2020

		double limit1 = 1;
		for (int i=1; i<significantDigits; i++)
			limit1 *=10;
		double limit2 = limit1 * 10;

		int exponent = 0;
		while (Math.abs(value) + 0.5 < limit1) {
			value *= 10;
			exponent--;
			}
		while (Math.abs(value) + 0.5 >= limit2) {
			value /= 10.0;
			exponent++;
			}

		return toString((long)(value+(value < 0 ? -0.5 : 0.5)), exponent, significantDigits, skipTrailingZeros);
		}

	/**
	 * Converts the value with the given exponent into a short string representation
	 * using the scientific notation if it is more compact. Trailing zeros are not shown.
	 * @param value
	 * @param exponent
	 * @return
	 */
	public static String toShortString(long value, int exponent) {
		return toString(value, exponent, 16, true);
		}

	private static String toString(long value, int exponent, int significantDigits, boolean skipTrailingZeros) {
		int noOfCiphers = 1;

		if (value == 0)
			return (skipTrailingZeros || significantDigits==1) ? "0" : "0."+zeros(significantDigits-1);

		if (value != 0) {
			while (value % 10 == 0) {
				value /= 10;
				exponent++;
				}

			noOfCiphers = 0;
			for (long temp=value; temp!=0; temp/=10)
				noOfCiphers++;
			}

		assert(significantDigits >= noOfCiphers);

		StringBuilder label = new StringBuilder();
		if (value < 0)
			label.append('-');

		int trailingZeros = Math.max(0, significantDigits-noOfCiphers);

		String cipherString = Long.toString(Math.abs(value));
		if (exponent == 0) {
			label.append(cipherString);
			if (!skipTrailingZeros && trailingZeros != 0) {
				label.append(".");
				label.append(zeros(trailingZeros));
				}
			}
		else if (exponent > 0) {
			boolean useScientificNotation = skipTrailingZeros ? (exponent > 4) : (exponent-trailingZeros > 4);
			if (useScientificNotation) {	// scientific notation for numbers above/equal 0.001
				if (noOfCiphers == 1) {
					label.append(cipherString);
					if (significantDigits > 1) {
						label.append(".");
						if (skipTrailingZeros)
							label.append("0");
						else
							label.append(zeros(trailingZeros));
						}
					}
				else {
					label.append(cipherString.substring(0, 1));
					label.append(".");
					label.append(cipherString.substring(1));
					if (!skipTrailingZeros && trailingZeros != 0)
						label.append(zeros(trailingZeros));
					}
				label.append("e");
				label.append(Integer.toString(exponent+noOfCiphers-1));
				}
			else {
				label.append(cipherString);
				label.append(zeros(exponent));
				if (!skipTrailingZeros && trailingZeros > exponent) {
					label.append(".");
					label.append(zeros(trailingZeros-exponent));
					}
				}
			}
		else {  // exponent < 0
			if (-exponent < noOfCiphers) {
				label.append(cipherString.substring(0, noOfCiphers+exponent));
				label.append(".");
				label.append(cipherString.substring(noOfCiphers+exponent));
				if (!skipTrailingZeros && trailingZeros != 0)
					label.append(zeros(trailingZeros));
				}
			else if (-exponent == noOfCiphers) {
				label.append("0.");
				label.append(cipherString);
				if (!skipTrailingZeros && trailingZeros != 0)
					label.append(zeros(trailingZeros));
				}
			else {
				if (exponent+noOfCiphers < -2) {	// scientific notation for numbers below 0.001
					if (noOfCiphers == 1) {
						label.append(cipherString);
						if (significantDigits > 1) {
							label.append(".");
							if (skipTrailingZeros)
								label.append("0");
							else
								label.append(zeros(trailingZeros));
							}
						}
					else {
						label.append(cipherString.charAt(0));
						label.append(".");
						label.append(cipherString.substring(1));
						if (!skipTrailingZeros && trailingZeros != 0)
							label.append(zeros(trailingZeros));
						}
					label.append("e-");
					label.append(Integer.toString(1-exponent-noOfCiphers));
					}
				else {
					label.append("0.");
					label.append(zeros(-exponent-noOfCiphers));
					label.append(cipherString);
					if (!skipTrailingZeros && trailingZeros != 0)
						label.append(zeros(trailingZeros));
					}
				}
			}
//System.out.println("int:"+theInt+" exp:"+exponent+" label:"+label);
		return label.toString();
		}

	private static String zeros(int i) {
		if (i<ZEROS.length)
			return ZEROS[i];
		StringBuilder b = new StringBuilder();
		while (i >= ZEROS.length) {
			b.append(ZEROS[ZEROS.length-1]);
			i -= ZEROS.length-1;
			}
		b.append(ZEROS[i]);
		return b.toString();
		}
	}