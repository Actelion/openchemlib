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
 * @author Modest v. Korff
 */

package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.calc.statistics.median.ModelMedianInteger;
import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;
import com.actelion.research.util.datamodel.IntArray;

/**
 * DistHistHelper
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 13.10.20.
 */
public class DistHistHelper {

    public static int getSpread(byte [] a){
        int start=0;
        for (int i = 0; i < a.length; i++) {
            if(a[i]>0){
                start=i;
                break;
            }
        }

        int end=0;
        for (int i = a.length-1; i >= 0; i--) {
            if(a[i]>0){
                end=i;
                break;
            }
        }

        return end-start+1;
    }

    public static int getMaxIndexNotZero(byte [] a){
        int end=0;
        for (int i = a.length-1; i >= 0; i--) {
            if(a[i]>0){
                end=i;
                break;
            }
        }

        return end;
    }

    public static int getMedianBin(byte [] a){

        int medianBin=-1;

        int sum = 0;
        for (byte b : a) {
            sum+=b;
        }
        sum /= 2;

        int s2=0;
        for (int i = 0; i < a.length; i++) {
            s2+=a[i];
            if(s2>=sum){
                medianBin=i;
                break;
            }
        }

        return medianBin;
    }

    public static RangeStatistics getRangeStatistics(MolDistHist mdh){

        RangeStatistics rangeStatisticsTotal = new RangeStatistics();

        int n = mdh.getNumPPNodes();

        int nn = ((n*n)-n)/2;
        IntArray iaMaxRange = new IntArray(nn);
        IntArray iaMedianRange = new IntArray(nn);

        for (int i = 0; i < n; i++) {

            for (int j = i+1; j < n; j++) {

                byte [] arr = mdh.getDistHist(i, j);

                int indexMax = -1;
                for (int k = arr.length-1; k >= 0; k--) {
                    if(arr[k]>0){
                        indexMax = k;
                        break;
                    }
                }

                iaMaxRange.add(indexMax);

                int indexMin = -1;
                for (int k = 0; k < arr.length; k++) {
                    if(arr[k]>0){
                        indexMin = k;
                        break;
                    }
                }

                int sum=0;
                for (int k = indexMin; k < indexMax+1; k++) {
                    sum += arr[k];
                }

                int half = sum / 2;

                sum=0;
                int indexMedian1=-1;
                for (int k = indexMin; k < indexMax+1; k++) {
                    sum += arr[k];
                    if(sum>=half){
                        indexMedian1=k;
                        break;
                    }
                }

                sum=0;
                int indexMedian2=-1;
                for (int k = indexMax; k >= indexMin; k--) {
                    sum += arr[k];
                    if(sum>=half){
                        indexMedian2=k;
                        break;
                    }
                }

                int medianRange = (int)((indexMedian1+indexMedian2)/2.0 +0.5);
                iaMedianRange.add(medianRange);
            }
        }

        ModelMedianInteger mmi = ArrayUtilsCalc.getMedian(iaMedianRange.get());

        rangeStatisticsTotal.maxRange = iaMaxRange.max();
        rangeStatisticsTotal.medianRange = mmi.median;

        return rangeStatisticsTotal;
    }

    public static int count(byte [] a){
        int c = 0;
        for (int i = 0; i < a.length; i++) {
            c += a[i];
        }
        return c;
    }

    public static byte [] normalize(double [] arrHistRaw){

        int countValuesInHistogram = 0;
        for (int i = 0; i < arrHistRaw.length; i++) {
            countValuesInHistogram += arrHistRaw[i];
        }
        byte [] arrHistPercent = new byte [arrHistRaw.length];

        for (int i = 0; i < arrHistRaw.length; i++) {
            arrHistPercent[i]= (byte)  (((arrHistRaw[i] / countValuesInHistogram) * ConstantsFlexophoreGenerator.SUM_VAL_HIST) + 0.5);
        }

//            System.out.println(StringFunctions.toString(arrHistPercent));

        return arrHistPercent;

    }
    public static byte [] normalize(byte [] arrHistRaw){

        double countValuesInHistogram = 0;
        for (int i = 0; i < arrHistRaw.length; i++) {
            countValuesInHistogram += arrHistRaw[i];
        }
        byte [] arrHistPercent = new byte [arrHistRaw.length];

        for (int i = 0; i < arrHistRaw.length; i++) {

            double v = ((arrHistRaw[i] / countValuesInHistogram) * ConstantsFlexophoreGenerator.SUM_VAL_HIST) + 0.5;


            arrHistPercent[i]= (byte)  v;
        }

//            System.out.println(StringFunctions.toString(arrHistPercent));

        return arrHistPercent;

    }


    public static class RangeStatistics {

        public int maxRange;
        public int medianRange;

    }
}
