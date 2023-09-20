package com.actelion.research.util;

import com.actelion.research.calc.Logarithm;
import com.actelion.research.chem.descriptor.DescriptorEncoder;
import com.actelion.research.util.datamodel.IntArray;

import java.nio.charset.StandardCharsets;
import java.util.Random;

/**
 * Encodes integer numbers into a string.
 */
public class EncoderIntegerNumbers {

    private static final int CAPACITY_DATA = 100;

    private int min;
    private int max;
    private int bits;

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

        if(min == max) {
            throw new RuntimeException("Range is 0!");
        }

        initialize(min, max);

        for (int v : arrValue) {
            add(v);
        }
    }


    private void initialize(int minValue, int maxValue){

        min = minValue;
        max = maxValue;

        int range = max-min;

        this.bits = Logarithm.log2(range)+1;

        dataEncoded = new int [CAPACITY_DATA];

        ccBitCounter = 0;

        // The first byte stores the number of bits
        ccBitCounter += Byte.SIZE;
        // This integer stores the number of added values
        ccBitCounter += Integer.SIZE;

        // min value
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

        IntArray da = new IntArray();

        int n = 100;

        int min = 0;
        int max = 1000;

        Random random = new Random();
        for (int i = 0; i < n; i++) {

            int v = 0;
            if(random.nextBoolean()) {
                v = random.nextInt(max);
            } else {
                if(min>0) {
                    v = -random.nextInt(Math.abs(min));
                } else {
                    v = random.nextInt(max);
                }
            }
            da.add(v);



        }


//        da.add(18);
//        da.add(-70);
//        da.add(1);
//        da.add(2);

        int [] arr = da.get();


        EncoderIntegerNumbers encoder = new EncoderIntegerNumbers(arr);

        String strData = encoder.encode();

        System.out.println(strData);

        int [] arrValueDecoded = EncoderIntegerNumbers.decode(strData);

        for (int i = 0; i < arrValueDecoded.length; i++) {

            System.out.println(arr[i] + "\t" + arrValueDecoded[i]);

        }

    }


}
