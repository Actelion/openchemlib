package com.actelion.research.util;

import com.actelion.research.calc.Logarithm;
import com.actelion.research.chem.descriptor.DescriptorEncoder;
import com.actelion.research.util.datamodel.DoubleArray;
import com.actelion.research.util.datamodel.IntArray;

public class EncoderIntegerNumbers {

    private static final int CAPACITY_DATA = 100;

    private int min;
    private int max;
    private int bits;

    private int range;

    private int [] dataEncoded;

    private int ccBitCounter;
    private int ccAddedValuesCounter;

    private boolean finalized;


    /**
     * More values can be added after the constructor was called. However, the must not exceed the min and the max
     * value found in arrValue.
     * @param arrValue array with the values that will be encoded. Determines also minimum and maximum encoding value.
     */
    public EncoderIntegerNumbers(int [] arrValue) {

        int min = Integer.MAX_VALUE;
        int max = -Integer.MAX_VALUE;

        for (int v : arrValue) {

            if(v<min){
                min=v;
            }
            if(v>max){
                max=v;
            }
        }

        initialize(min, max);

        for (int v : arrValue) {
            add(v);
        }
    }


    private void initialize(int minValue, int maxValue){

        min = minValue;
        max = maxValue;

        range = max-min;

        this.bits = Logarithm.log2(range)+1;

        dataEncoded = new int [CAPACITY_DATA];

        ccBitCounter = 0;

        // The first byte stores the number of bits
        ccBitCounter += Byte.SIZE;
        // This integer stores the number of added values
        ccBitCounter += Integer.SIZE;

        // min value
        ccBitCounter += Integer.SIZE;

        // range
        ccBitCounter += Integer.SIZE;

        ccAddedValuesCounter = 0;

        finalized = false;
    }

    /**
     * Add a value that will be encoded.
     * @param v
     */
    private void add(int v){

        if(finalized){
            throw new RuntimeException("Already finalized!");
        }

        if(v < min) {
            throw new RuntimeException("Value lower than minimum!");
        }

        if(v > max) {
            throw new RuntimeException("Value higher than maximum!");
        }

        int valueNorm = v - min;

        for (int i = 0; i < bits; i++) {
            if((1 & valueNorm)==1){
                BitUtils.setBit(dataEncoded, ccBitCounter);
            } else {
                BitUtils.unsetBit(dataEncoded, ccBitCounter);
            }

            valueNorm >>= 1;

            ccBitCounter++;

            if(!BitUtils.isValidBitIndex(dataEncoded, ccBitCounter)){

                int newlen = dataEncoded.length * 2;

                dataEncoded = IntArray.resize(dataEncoded, newlen);
            }
        }

        ccAddedValuesCounter++;
    }

    private void set(int l, int bits, int offset){

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

    private int [] finalizeAndGet() {

        finalized = true;

        int offset = 0;

        set(bits, Byte.SIZE, offset);

        offset += Byte.SIZE;

        set(ccAddedValuesCounter, Integer.SIZE, offset);

        offset += Integer.SIZE;

        set(min, Integer.SIZE, offset);

        offset += Integer.SIZE;

        set(range, Integer.SIZE, offset);

        int nInteger = ccBitCounter / Integer.SIZE + 1;

        int [] data = new int[nInteger];

        System.arraycopy(dataEncoded, 0, data, 0, nInteger);

        return data;
    }

    private String encode() {

        int [] data = finalizeAndGet();

        String strData = new String(new DescriptorEncoder().encode(data));

        return strData;
    }

    /**
     *
     * @param strData
     * @return array with decoded values.
     */
    public static int [] decode(String strData){

        int [] data = new DescriptorEncoder().decode(strData);

        return decode(data);
    }

    private static int [] decode(int [] data) {

        int offset = 0;

        int bits = decode(data, offset, Byte.SIZE);

        offset += Byte.SIZE;

        int nAddedValues =  decode(data, offset, Integer.SIZE);

        offset += Integer.SIZE;

        int min =  decode(data, offset, Integer.SIZE);

        offset += Integer.SIZE;

        int range =  decode(data, offset, Integer.SIZE);

        offset += Integer.SIZE;

        int [] arrValue = new int[nAddedValues];

        for (int i = 0; i < nAddedValues; i++) {
            int value = decode(data, offset, bits)+min;

            arrValue[i] = value;

            offset += bits;
        }

        return arrValue;
    }

    private static int decode(int [] data, int offset, int length) {

        int l = 0;

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
     * @return
     */
    public static String encode(int [] arrValue) {
        EncoderIntegerNumbers efpnCoeff = new EncoderIntegerNumbers(arrValue);
        return efpnCoeff.encode();
    }

    public static void main(String[] args) {


        double min = -7;


        double max = 8;


        IntArray da = new IntArray();

        da.add(8);
        da.add(-7);
        da.add(1);
        da.add(2);

        int [] arr = da.get();


        EncoderIntegerNumbers encoder = new EncoderIntegerNumbers(arr);

        String strData = encoder.encode();

        System.out.println(strData);

        int [] arrValueDecoded = EncoderIntegerNumbers.decode(strData);

        for (int i = 0; i < arrValueDecoded.length; i++) {

            System.out.println(arr[i] + "\t" + arrValueDecoded[i]);

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


//    private static class Decoder {
//
//        private double min;
//        private double factorPrecision;
//        private double range;
//
//
//        public Decoder(int precisionBits, double min, double range) {
//            this.min = min;
//            this.range = range;
//
//            factorPrecision = Math.pow(2,precisionBits);
//        }
//
//        private double getDecoded(long encoded){
//            double vScaled2RangeDecoded = encoded / factorPrecision * range + min;
//
//            return vScaled2RangeDecoded;
//        }
//
//    }

}
