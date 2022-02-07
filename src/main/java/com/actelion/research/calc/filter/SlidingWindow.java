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

package com.actelion.research.calc.filter;

import com.actelion.research.util.ArrayUtils;

/**
 * SlidingWindow
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
