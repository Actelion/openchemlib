/*
 * Copyright (c) 2020.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
 */

package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.calc.statistics.median.ModelMedianInteger;
import com.actelion.research.util.datamodel.IntArray;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

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


    public static class RangeStatistics {

        public int maxRange;
        public int medianRange;

    }
}
