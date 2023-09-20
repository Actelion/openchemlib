package com.actelion.research.util;

import com.actelion.research.chem.descriptor.DescriptorEncoder;
import com.actelion.research.util.datamodel.DoubleArray;
import com.actelion.research.util.datamodel.IntArray;

import java.nio.charset.StandardCharsets;

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
public class EncoderFloatingPointNumbers {

    private static final double TINY_FACTOR_EXPONENT = -6;

    private static final double TINY_FACTOR = Math.pow(10, TINY_FACTOR_EXPONENT);

    private static final int CAPACITY_DATA = 100;

    private double min;
    private double max;
    private int precisionBits;
    private double factorPrecision;
    private double range;

    private int [] dataEncoded;

    private int ccBitCounter;
    private int ccAddedValuesCounter;

    private boolean finalized;

    private double tiny;

    /**
     * More values can be added after the constructor was called. However, the must not exceed the min and the max
     * value found in arrValue.
     * @param arrValue array with the values that will be encoded. Determines also minimum and maximum encoding value.
     * @param precisionBits
     */
    public EncoderFloatingPointNumbers(double [] arrValue, int precisionBits) {

        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;

        for (double v : arrValue) {

            if(v<min){
                min=v;
            }
            if(v>max){
                max=v;
            }
        }

        initialize(min, max, precisionBits);

        for (double v : arrValue) {
            add(v);
        }
    }

    public EncoderFloatingPointNumbers(float [] arrValue, int precisionBits) {

        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;

        for (double v : arrValue) {

            if(v<min){
                min=v;
            }
            if(v>max){
                max=v;
            }
        }

        initialize(min, max, precisionBits);

        for (double v : arrValue) {
            add(v);
        }
    }

    /**
     *
     * @param min minimum value that can be encoded.
     * @param max maximum value that can be encoded.
     * @param precisionBits precision in encoding.
     */
    private EncoderFloatingPointNumbers(double min, double max, int precisionBits) {
        initialize(min, max,  precisionBits);
    }


    private void initialize(double minValue, double maxValue, int precisionBits){

        if(precisionBits > Long.SIZE) {
            throw new RuntimeException("Maximum possible precision exceeded!");
        }

        // If we do not correct for the minimum and the maximum value we run into rounding problems with the extreme
        // values.
        double rangeValue = Math.abs(maxValue-minValue);
        int exponent;
        if(rangeValue<1E-06) 
        	exponent = (int)TINY_FACTOR_EXPONENT;
        else
        	exponent = (int)(Math.log10(rangeValue) + TINY_FACTOR_EXPONENT);

        //if(exponent <= Double.MIN_EXPONENT) {
        //   throw new RuntimeException("Minimum value out of range!");
        //}

        tiny = rangeValue * TINY_FACTOR;

        min = minValue - tiny;
        max = maxValue + tiny;

        this.precisionBits = precisionBits;

        factorPrecision = Math.pow(2,precisionBits);

        range = max-min;

        dataEncoded = new int [CAPACITY_DATA];

        ccBitCounter = 0;

        // The first byte stores the precision
        ccBitCounter += Byte.SIZE;
        // This integer stores the number of added values
        ccBitCounter += Integer.SIZE;

        // min value
        ccBitCounter += Long.SIZE;

        // range
        ccBitCounter += Long.SIZE;

        ccAddedValuesCounter = 0;

        finalized = false;
    }

    /**
     * Add a value that will be encoded.
     * @param value
     */
    private void add(double value){

        if(finalized){
            throw new RuntimeException("Already finalized!");
        }

        double delta = value-min;
        if(delta < -tiny) {
            throw new RuntimeException("Value lower than minimum!");
        }

        double max = min+range;

        if((value-max) > tiny) {
            throw new RuntimeException("Value higher than maximum!");
        }

        long encoded = getEncoded(value);

        add(encoded, precisionBits);

        ccAddedValuesCounter++;

    }

    private void add(long l, int bits){

        for (int i = 0; i < bits; i++) {
            if((1 & l)==1){
                BitUtils.setBit(dataEncoded, ccBitCounter);
            } else {
                BitUtils.unsetBit(dataEncoded, ccBitCounter);
            }

            l >>= 1;

            ccBitCounter++;

            if(!BitUtils.isValidBitIndex(dataEncoded, ccBitCounter)){

                int newlen = dataEncoded.length * 2;

                dataEncoded = IntArray.resize(dataEncoded, newlen);
            }
        }
    }

    private void set(long l, int bits, int offset){

        for (int i = 0; i < bits; i++) {

            int indexBit = offset + i;

            if((l & 1)==1){
                BitUtils.setBit(dataEncoded, indexBit);
            } else {
                BitUtils.unsetBit(dataEncoded, indexBit);
            }
            l >>= 1;
        }
    }

    private long getEncoded(double vInput) {

        double vScaled2Range = (vInput - min) / range;

        long iVScaledRange = (long)(vScaled2Range * factorPrecision);

        return iVScaledRange;

    }

    private int [] finalizeAndGet() {

        finalized = true;

        long lPrecisionBits = precisionBits;

        int offset = 0;

        set(lPrecisionBits, Byte.SIZE, offset);

        offset += Byte.SIZE;

        set(ccAddedValuesCounter, Integer.SIZE, offset);

        offset += Integer.SIZE;

        long lMin = Double.doubleToLongBits(min);

        set(lMin, Long.SIZE, offset);

        offset += Long.SIZE;

        long lRange = Double.doubleToLongBits(range);

        set(lRange, Long.SIZE, offset);

        int nInteger = ccBitCounter / Integer.SIZE + 1;

        int [] data = new int[nInteger];

        System.arraycopy(dataEncoded, 0, data, 0, nInteger);

        return data;
    }

    private String encode() {

        int [] data = finalizeAndGet();

        String strData = new String(new DescriptorEncoder().encode(data), StandardCharsets.UTF_8);

        return strData;
    }

    /**
     *
     * @param strData
     * @return array with decoded values.
     */
    public static double [] decode(String strData){

        int [] data = new DescriptorEncoder().decode(strData);

        return decode(data);
    }

    private static double [] decode(int [] data) {

        int offset = 0;

        int precisionBits = (int)decode(data, offset, Byte.SIZE);

        offset += Byte.SIZE;

        int nAddedValues =  (int)decode(data, offset, Integer.SIZE);

        offset += Integer.SIZE;

        long lMin =  decode(data, offset, Long.SIZE);

        double min = Double.longBitsToDouble(lMin);

        offset += Long.SIZE;

        long lRange =  decode(data, offset, Long.SIZE);

        double range = Double.longBitsToDouble(lRange);

        offset += Long.SIZE;

//        System.out.println("precisionBits " + precisionBits);
//        System.out.println("nAddedValues " + nAddedValues);
//        System.out.println("min " + min);
//        System.out.println("range " + range);

        Decoder decoder = new Decoder(precisionBits, min, range);

        double [] arrValue = new double[nAddedValues];

        for (int i = 0; i < nAddedValues; i++) {
            long lValue = decode(data, offset, precisionBits);

            double value = decoder.getDecoded(lValue);

            arrValue[i] = value;
            offset += precisionBits;
        }

        return arrValue;
    }

    private static long decode(int [] data, int offset, int length) {

        long l = 0;

        int start = offset+length-1;

        for (int i = start; i >= offset; i--) {
            if(BitUtils.isBitSet(data, i)){
                l |= 1;
            }
            if(i>offset) {
                l <<= 1;
            }
        }

        return l;
    }

    /**
     * Convenience method to encode a array with a given precision.
     * @param arrValue
     * @param precisionBits
     * @return
     */
    public static String encode(double [] arrValue, int precisionBits) {
        EncoderFloatingPointNumbers efpnCoeff = new EncoderFloatingPointNumbers(arrValue, precisionBits);
        return efpnCoeff.encode();
    }

    public static String encode(float [] arrValue, int precisionBits) {
        EncoderFloatingPointNumbers efpnCoeff = new EncoderFloatingPointNumbers(arrValue, precisionBits);
        return efpnCoeff.encode();
    }


    public static void main(String[] args) {

        int precisionInBits = 61;

        double min = -7.742162891457342;


        double max = 8.161857389358609;


        // double v = 0.723867952243311;
        double v = 8.161857389358609;



        // double v = 0.123;


        DoubleArray da = new DoubleArray();

        da.add(8.161857389358609);
        da.add(1000.003);
        da.add(1000000.005);
        da.add(6700000001.005);

        double [] arr = da.get();


        EncoderFloatingPointNumbers encoderFloatingPointNumbers = new EncoderFloatingPointNumbers(arr, precisionInBits);

        String strData = encoderFloatingPointNumbers.encode();

        System.out.println(strData);

        double [] arrValueDecoded = EncoderFloatingPointNumbers.decode(strData);

        for (int i = 0; i < arrValueDecoded.length; i++) {

            double delta =  arrValueDecoded[i] - arr[i];

            System.out.println(delta + "\t" + arr[i] + "\t" + arrValueDecoded[i]);

        }

    }

//    public static void main(String[] args) {
//
//        min = -7.742162891457342
//        max = 8.161857389358609
//        v = 0.723867952243311
//
//        int n = 10;
//
//        double min = -5;
//
//        double max = 34;
//
//        int precisionInBits = 19;
//
//        EncoderFloatingPointNumbers encoderFloatingPointNumbers = new EncoderFloatingPointNumbers(min, max, precisionInBits);
//
//        encoderFloatingPointNumbers.add(-4);
//        encoderFloatingPointNumbers.add(0);
//        encoderFloatingPointNumbers.add(5);
//        encoderFloatingPointNumbers.add(1.23456);
//
//        String strData = encoderFloatingPointNumbers.encode();
//
//        System.out.println(strData);
//
//        double [] arrValue = EncoderFloatingPointNumbers.decode(strData);
//
//        for (double v : arrValue) {
//            System.out.println(v);
//        }
//    }

//    public static void main(String[] args) {
//
//        int n = 10;
//
//        double min = -5;
//
//        double max = 34;
//
//        int precisionInBits = 32;
//
//        EncoderFloat encoderFloat = new EncoderFloat(min, max, precisionInBits);
//
//
//
//
//        Long lMin = Double.doubleToLongBits(min);
//
//        int offset = 0;
//
//        encoderFloat.set(lMin, Long.SIZE, offset);
//
//
//
//        long lMinDecoded = decode(encoderFloat.dataEncoded, offset, Long.SIZE);
//
//        double minFromLong = Double.longBitsToDouble(lMin);
//
//        double minFromLongDecoded = Double.longBitsToDouble(lMinDecoded);
//
//
//        System.out.println(min);
//        System.out.println(minFromLong);
//        System.out.println(minFromLongDecoded);
//
//    }

//    public static void main(String[] args) {
//
//        int n = 10;
//
//        double min = -5;
//
//        double max = 34;
//
//        int precisionInBits = 32;
//
//        EncoderFloat encoderFloat = new EncoderFloat(min, max, precisionInBits);
//
//        Random rnd = new Random();
//
//        double range = Math.abs(max-min);
//
//        for (int i = 0; i < n; i++) {
//            double  vRaw = rnd.nextDouble() * range;
//
//            double vInput = vRaw + min;
//
//            long iVScaledRange = encoderFloat.getEncoded(vInput);
//
////            System.out.println(vScaled2Range + "\t" + iVScaledRange);
//
//            double vScaled2RangeDecoded = encoderFloat.getDecoded(iVScaledRange);
//
//            System.out.println(vInput + "\t" + vScaled2RangeDecoded);
//        }
//
//        long iVScaledRange = encoderFloat.getEncoded(max);
//
//        double vScaled2RangeDecoded = encoderFloat.getDecoded(iVScaledRange);
//
//        System.out.println(max + "\t" + iVScaledRange + "\t" + vScaled2RangeDecoded);
//    }


    private static class Decoder {

        private double min;
        private double factorPrecision;
        private double range;


        public Decoder(int precisionBits, double min, double range) {
            this.min = min;
            this.range = range;

            factorPrecision = Math.pow(2,precisionBits);
        }

        private double getDecoded(long encoded){
            double vScaled2RangeDecoded = encoded / factorPrecision * range + min;

            return vScaled2RangeDecoded;
        }

    }

}
