/*
* Copyright (c) 1997 - 2015
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
    /**
     * Converts a double value into a String representation in scientific format
     * rounded to 5 significant digits, e.g. 1.2345e-6.
     * Double values that are effectively integers are turned into the scientific
     * format only, if they consist of more than 8 digits.
     * @param theDouble
     * @param significantDigits
     * @return
     */
    public static String toString(double theDouble) {
        return toString(theDouble, 5);
        }

    /**
     * Converts a double value into a String representation in scientific format
     * rounded to a definable number of digits, e.g. 1.2345e-6.
     * Double values that are effectively integers are turned into the scientific
     * format only, if they consist of more than significantDigits+3 digits.
     * @param theDouble
     * @param significantDigits
     * @return
     */
    public static String toString(double theDouble, int significantDigits) {
        if (Double.isInfinite(theDouble))
            return "Infinity";
        if (Double.isNaN(theDouble))
            return "NaN";
        if (theDouble == 0.0)
            return "0";

        double limit = Math.pow(10.0, significantDigits);

        if (theDouble - (int)theDouble == 0.0
         && theDouble < 1000 * limit)
        	return ""+(int)theDouble;

        int exponent = 0;
        while (Math.abs(theDouble) + 0.5 < limit) {
            theDouble *= 10;
            exponent--;
            }
        while (Math.abs(theDouble) + 0.5 >= limit) {
            theDouble /= 10.0;
            exponent++;
            }

        return toString((int)(theDouble+(theDouble < 0 ? -0.5 : 0.5)), exponent);
        }

	public static String toString(int theInt, int exponent) {
		int noOfCiphers = 0;
        if (theInt != 0) {
            while (theInt % 10 == 0) {
                theInt /= 10;
                exponent++;
                }
    		for (int temp=theInt; temp!=0; temp/=10)
    			noOfCiphers++;
            }

		String label;

        String cipherString = ""+Math.abs(theInt);
        String sign = (theInt < 0) ? "-" : "";
		if (exponent == 0 || theInt == 0)
			label = sign+cipherString;
		else if (exponent > 0) {
			if (noOfCiphers+exponent > 5) {
                if (noOfCiphers == 1)
                    label = sign+cipherString+".0e"+exponent;
                else
                    label = sign+cipherString.substring(0, 1)+"."+cipherString.substring(1)+"e"+(exponent+noOfCiphers-1);
                }
			else {
				label = sign+cipherString;
				for (int i=0; i<exponent; i++)
					label = label.concat("0");
				}
			}
		else {  // exponent < 0
			if (-exponent < noOfCiphers)
				label = sign + cipherString.substring(0, noOfCiphers+exponent)
					  + "." + cipherString.substring(noOfCiphers+exponent);
			else if (-exponent == noOfCiphers)
				label = sign + "0."+cipherString;
		    else {
				if (exponent < -5) {
					if (noOfCiphers == 1)
						label = sign + cipherString+".0e-";
				    else
						label = sign + cipherString.charAt(0)+"."+cipherString.substring(1)+"e-";
					label = label+(1-exponent-noOfCiphers);
					}
				else {
				    label = sign + "0.";
					for (int i=0; i<-exponent-noOfCiphers; i++)
						label = label.concat("0");
					label = label+cipherString;
					}
				}
			}
//System.out.println("int:"+theInt+" exp:"+exponent+" label:"+label);
		return label;
		}
	}