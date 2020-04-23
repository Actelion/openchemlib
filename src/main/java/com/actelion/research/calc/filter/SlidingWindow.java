package com.actelion.research.calc.filter;

import com.actelion.research.util.ArrayUtils;

/**
 * SlidingWindow
 * <p>Copyright: Actelion Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 01.03.16.
 */
public class SlidingWindow {


    private double [] arrFilter;

    private int lenFilHalf;

    public SlidingWindow(double[] arrFilter) {

        this.arrFilter = arrFilter;

        if(arrFilter.length % 2 == 0){

            throw new RuntimeException("Odd number of filter values needed.");

        }

        lenFilHalf = arrFilter.length / 2;

    }

    public double [] filter(double [] a) {

        double[] aa = new double[a.length];

        System.arraycopy(a, 0, aa, 0, aa.length);

        int end = a.length - lenFilHalf;

        for (int i = lenFilHalf; i < end; i++) {

            double v = 0;

            for (int j = 0; j < arrFilter.length; j++) {

                int ind = i - lenFilHalf + j;

                v += a[ind] * arrFilter[j];

            }

            aa[i] = v;
        }

        return aa;
    }

    public int [] filter(int [] a) {

        int[] aa = new int[a.length];

        System.arraycopy(a, 0, aa, 0, aa.length);

        int end = a.length - lenFilHalf;

        for (int i = lenFilHalf; i < end; i++) {


            int v = 0;

            for (int j = 0; j < arrFilter.length; j++) {

                int ind = i - lenFilHalf + j;

                v += (int)(a[ind] * arrFilter[j] + 0.5);

            }

            aa[i] = v;
        }

        return aa;
    }

    public byte [] filter(byte [] a) {

        byte [] aa = new byte[a.length];

        System.arraycopy(a, 0, aa, 0, aa.length);

        int end = a.length - lenFilHalf;

        for (int i = lenFilHalf; i < end; i++) {


            byte v = 0;

            for (int j = 0; j < arrFilter.length; j++) {

                int ind = i - lenFilHalf + j;

                v += (byte)(a[ind] * arrFilter[j] + 0.5);

            }

            aa[i] = v;
        }

        return aa;
    }

    public static void main(String[] args) {


        int [] a = {0,0,0,0,100,100,0,0,0,0,0};

        double [] f = {0.25,0.5,0.25};

        SlidingWindow slidingWindow = new SlidingWindow(f);

        int [] aa = slidingWindow.filter(a);

        System.out.println(ArrayUtils.toString(aa));
    }





}
