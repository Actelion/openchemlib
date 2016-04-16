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

package com.actelion.research.util.convert;

import java.util.ArrayList;
import java.util.StringTokenizer;

public class String2DoubleArray {

    private static final String DELIMITER = " \t\n\r\f,;";


    public static double [] convert(String sLine) throws Exception{

        StringTokenizer st = new StringTokenizer(sLine, DELIMITER);

        double[] dArray = new double[st.countTokens()];

        int ii = 0;
        while (st.hasMoreTokens()) {
            String sNumber = st.nextToken();
            // The formatting sign "'" is not recocnized by the Double.valueOf(...)
            // function.
            sNumber = sNumber.replaceAll("'", "");
            try {
                dArray[ii] = Double.parseDouble(sNumber);
            }
            catch (NumberFormatException ex1) {
                throw new  NumberFormatException("No number: " + sNumber + ".");
            }
            ii++;
        }

        return dArray;
    }

    public static double[] convert(ArrayList sArrList) throws Exception{

        ArrayList dArrList = new ArrayList();
        int iSize = 0;
        for (int ii = 0; ii < sArrList.size(); ii++) {
            double[] arr = convert( (String) sArrList.get(ii));
            dArrList.add(arr);
            iSize += arr.length;
        }
        double[] dArrAll = new double[iSize];
        int index = 0;
        for (int ii = 0; ii < dArrList.size(); ii++) {
            double[] arr = (double[]) dArrList.get(ii);
            for (int jj = 0; jj < arr.length; jj++) {
                dArrAll[index++] = arr[jj];
            }
        }
        return dArrAll;
    }

    public static void main(String [] args) {
        String str = "1 2 3 ";
        try {
            double[] arr = convert(str);
            for (int ii = 0; ii < arr.length; ii++) {
              System.out.println(arr[ii]);
            }

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }

    }

}