/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
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