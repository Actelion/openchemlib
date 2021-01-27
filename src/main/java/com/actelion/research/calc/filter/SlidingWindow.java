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

    private int l;
    private int l_half;

    public SlidingWindow(double[] arrFilter) {
        this.arrFilter = arrFilter;
        if(arrFilter.length % 2 == 0){
            throw new RuntimeException("Odd number of filter values needed.");
        }

        l = arrFilter.length;
        l_half = arrFilter.length/2;
    }

    public double [] filter(double [] a) {

        double [] aa = new double[a.length+l*2];
        double [] aaa = new double[a.length+l*2];
        System.arraycopy(a, 0, aa, l, a.length);
        int end = aa.length- l_half;
        for (int i = l_half; i < end; i++) {
            double v = 0;
            for (int j = 0; j < arrFilter.length; j++) {
                int ind = i - l_half + j;
                v += aa[ind] * arrFilter[j] + 0.5;
            }
            aaa[i] = v;
        }

        aa = new double[a.length];
        System.arraycopy(aaa, l, aa, 0, aa.length);
        return aa;

    }

    public int [] filter(int [] a) {

        int [] aa = new int[a.length+l*2];
        int [] aaa = new int[a.length+l*2];
        System.arraycopy(a, 0, aa, l, a.length);
        int end = aa.length- l_half;
        for (int i = l_half; i < end; i++) {
            int v = 0;
            for (int j = 0; j < arrFilter.length; j++) {
                int ind = i - l_half + j;
                v += (int)(aa[ind] * arrFilter[j] + 0.5);
            }
            aaa[i] = v;
        }

        aa = new int[a.length];
        System.arraycopy(aaa, l, aa, 0, aa.length);
        return aa;
    }

    public byte [] filter(byte [] a) {

        byte [] aa = new byte[a.length+l*2];
        byte [] aaa = new byte[a.length+l*2];
        System.arraycopy(a, 0, aa, l, a.length);
        int end = aa.length- l_half;

        for (int i = l_half; i < end; i++) {
            byte v = 0;
            for (int j = 0; j < arrFilter.length; j++) {
                int ind = i - l_half + j;
                v += (byte)(aa[ind] * arrFilter[j] + 0.5);
            }
            aaa[i] = v;
        }
        aa = new byte[a.length];
        System.arraycopy(aaa, l, aa, 0, aa.length);

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
