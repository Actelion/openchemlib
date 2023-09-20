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

package com.actelion.research.chem.descriptor;

import java.nio.charset.StandardCharsets;

/**
 * DescriptorEncoder encodes int[] based descriptors
 * into byte arrays that may be used to instantiate Strings
 */

public class DescriptorEncoder {
    public static final int MAX_COUNT_VALUE = 63;
    private static final int BITS = 6;
    private static final int PAIR_BITS = 4;

    // CODE Strings must contain highest ASCII character at the end; unused characters: " ' \ `
    private static final byte[] sCode = "0123456789@ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz".getBytes(StandardCharsets.UTF_8);
    private static final byte[] sCodeMultipleMin = "!#$%&()*+,-./".getBytes(StandardCharsets.UTF_8);
    private static final byte[] sCodeMultipleMax = ":;<=>?[]^{|}~".getBytes(StandardCharsets.UTF_8);
    private static volatile int[] sDecode,sDecodeMultiple;

    private byte[]  mBytes;
    private int     mByteIndex,mAvailableBits,mTempData,mByteMask;
	private int     mTempDataLong;

    public DescriptorEncoder() {
    	if (sDecode == null) {
    		synchronized(this.getClass()) {
				if (sDecode == null) {
			        int len = 1 << BITS;
			        assert len <= sCode.length : "Error in encoding, not enough characters.";

			        sDecode = new int[sCode[sCode.length-1]+1];
			        for (int i=0; i<sCode.length; i++)
			            sDecode[sCode[i]] = i;
			        sDecodeMultiple = new int[Math.max(sCodeMultipleMin[sCodeMultipleMin.length-1],
			                                           sCodeMultipleMax[sCodeMultipleMax.length-1])+1];
			        for (int i=0; i<sCodeMultipleMin.length; i++)
			            sDecodeMultiple[sCodeMultipleMin[i]] = -i-2;
			        for (int i=0; i<sCodeMultipleMax.length; i++)
			            sDecodeMultiple[sCodeMultipleMax[i]] = i+2;
					}
    			}
    		}
        }

    /**
     * Encodes a binary fingerprint stored as int[].
     * @param data binary fingerprint
     * @return byte[] of encoded character sequence
     */
    public byte[] encode(int[] data) {
        encodeStart(32 * data.length);
        for (int i=0; i<data.length; i++)
            encodeBits(data[i], 32);
        encodeBitsEnd();
        encodeDuplicateBytes();
        return mBytes;
		}

	/**
	 * Encodes any int[][] containing positive values.
	 * @param data not necessarily rectangular 2-dimensional array without null values
	 * @return byte[] of encoded character sequence
	 */
	public byte[] encodeIntArray2D(int[][] data) {
		int maxCount = data.length;
		int maxValue = 0;
		int valueCount = 0;
		for (int i=0; i<data.length; i++) {
			if (maxCount < data[i].length)
				maxCount = data[i].length;
			for (int j = 0; j<data[i].length; j++)
				if (maxValue<data[i][j])
					maxValue = data[i][j];
			valueCount += data[i].length;
			}
		int countBits = getNeededBits(maxCount);
		int valueBits = getNeededBits(maxValue);

		encodeStart(10 + (1+data.length)*countBits + valueCount*valueBits);
		encodeBits(countBits, 5);
		encodeBits(valueBits, 5);
		encodeBits(data.length, countBits);
		for (int i=0; i<data.length; i++) {
			encodeBits(data[i].length, countBits);
			for (int j=0; j<data[i].length; j++)
				encodeBits(data[i][j], valueBits);
			}
		encodeBitsEnd();
		encodeDuplicateBytes();
		return mBytes;
	}

	/**
	 * Decodes a 2-dimensional int[][] with positive values.
	 * @param bytes encoded character sequence
	 * @return int[][] decoded 2-dimensional int array
	 */
	public int[][] decodeIntArray2D(byte[] bytes) {
		if (bytes.length == 0)
			return null;

		bytes = decodeDuplicateBytes(bytes);
		decodeStart(bytes);

		int countBits = decodeBits(5);
		int valueBits = decodeBits(5);
		int outerSize = decodeBits(countBits);

		int[][] data = new int[outerSize][];
		for (int i=0; i<outerSize; i++) {
			int innerSize = decodeBits(countBits);
			data[i] = new int[innerSize];
			for (int j=0; j<innerSize; j++)
				data[i][j] = decodeBits(valueBits);
			}

		return data;
		}

	/**
	 * Encodes a binary fingerprint stored as long[].
	 * @param data binary fingerprint
	 * @return byte[] of encoded character sequence
	 */
	public byte[] encodeLong(long[] data) {
		encodeStart(64 * data.length);
		for (int i=0; i<data.length; i++)
			encodeBits(data[i], 64);

		encodeBitsEnd();
		encodeDuplicateBytes();
		return mBytes;
		}

	/**
     * Decodes a binary fingerprint from a String into an int[].
     * @param s encoded binary fingerprint
     * @return int[] binary fingerprint
     */
	public int[] decode(String s) {
		return decode(s.getBytes(StandardCharsets.UTF_8));
		}

    /**
     * Decodes a binary fingerprint from a String into an int[].
     * @param bytes encoded binary fingerprint as byte sequence
     * @return int[] binary fingerprint
     */
	public int[] decode(byte[] bytes) {
        if (bytes.length == 0)
            return null;

        bytes = decodeDuplicateBytes(bytes);
        decodeStart(bytes);
        int[] data = new int[BITS * bytes.length / 32];
        for (int i=0; i<data.length; i++)
            data[i] = decodeBits(32);
        return data;
        }

	/**
	 * Decodes a binary fingerprint from a String into a long[].
	 * @param s encoded binary fingerprint
	 * @return int[] binary fingerprint
	 */
	public long[] decodeLong(String s) {
		return decodeLong(s.getBytes(StandardCharsets.UTF_8));
		}

	/**
	 * Decodes a binary fingerprint from a String into a long[].
	 * @param bytes encoded binary fingerprint as byte sequence
	 * @return long[] binary fingerprint
	 */
	public long[] decodeLong(byte[] bytes) {
		if (bytes.length == 0)
			return null;

		bytes = decodeDuplicateBytes(bytes);
		decodeStart(bytes);
		long[] data = new long[BITS * bytes.length / 64];
		for (int i=0; i<data.length; i++)
			data[i] = decodeBitsLong(64);

		return data;
		}

	/**
     * Encodes a fragment/hash-value count list. Count values must not
     * exceed MAX_COUNT_VALUE in order to get encoded correctly.
     * @param data list of fragment counts or hash-value counts
     * @return byte[] of encoded character sequence
     */
    public byte[] encodeCounts(byte[] data) {
        encodeStart(6 * data.length);
        for (int i=0; i<data.length; i++)
            encodeBits(data[i], 6);
        encodeBitsEnd();
        encodeDuplicateBytes();
        return mBytes;
        }

    /**
     * Decodes a list of fragment/has-value counts from an encoded String.
     * @param s count list encoded as String
     * @return byte[] array with every byte representing a count value
     */
    public byte[] decodeCounts(String s) {
        return decodeCounts(s.getBytes(StandardCharsets.UTF_8));
        }

    /**
     * Decodes a list of fragment/has-value counts from an encoded
     * byte sequence into a byte[].
     * @param bytes count list encoded as byte sequence
     * @return byte[] array with every byte representing a count value
     */
    public byte[] decodeCounts(byte[] bytes) {
        if (bytes.length == 0)
            return null;

        bytes = decodeDuplicateBytes(bytes);
        decodeStart(bytes);
        byte[] data = new byte[BITS * bytes.length / 6];
        for (int i=0; i<data.length; i++)
            data[i] = (byte)decodeBits(6);
        return data;
        }

    /**
     * Encodes an int[] containing positive values without upper limit.
     * @param data
     * @return byte[] of encoded character sequence
     */
    public byte[] encodeIntArray(int[] data) {
        int max = 0;
		int sum = 0;
		for (int v:data) {
			sum += getNeededBits(v);
			if (max < v)
				max = v;
			}

		int sizeBits = Math.min(31, getNeededBits(data.length));
		int dataBits = getNeededBits(getNeededBits(max));

		encodeStart(1+5+3+sizeBits+sum+dataBits*data.length);

        encodeBits(0, 1);	// the positive value bit; this allows to later support negative values as well

        encodeBits(sizeBits, 5);
        encodeBits(data.length, sizeBits);	// size of the array

        encodeBits(dataBits, 3);	// bits needed to store the bits needed for next data value

        for (int i=0; i<data.length; i++) {
			int bits = getNeededBits(data[i]);
			encodeBits(bits, dataBits);
			encodeBits(data[i], bits);
			}
        encodeBitsEnd();
        encodeDuplicateBytes();
        return mBytes;
        }

    /**
     * Decodes an int[] without upper limit.
     * @param s String
     * @return int[] decoded int array
     */
    public int[] decodeIntArray(String s) {
    	return s == null ? null : decodeIntArray(s.getBytes(StandardCharsets.UTF_8));
    	}

    /**
     * Decodes an int[] without upper limit.
     * @param bytes encoded character sequence
     * @return int[] decoded int array
     */
    public int[] decodeIntArray(byte[] bytes) {
        if (bytes.length == 0)
            return null;

        bytes = decodeDuplicateBytes(bytes);
        decodeStart(bytes);
        decodeBits(1);	// trash the positive value bit for now

        int sizeBits = decodeBits(5);
        int arraySize = decodeBits(sizeBits);	// size of the array

        int dataBits = decodeBits(3);

        int[] data = new int[arraySize];
        for (int i=0; i<arraySize; i++) {
        	int bits = decodeBits(dataBits);
        	data[i] = decodeBits(bits);
        	}

        return data;
    	}

    /**
     * Encodes pairs of identifying integer with corresponding count values
     * in the form of an int[n][2], where n is the number of pairs, index=0
     * refers to the ID value and index=1 refers to the count value.
     * Neither ID values nor count values must be larger than 32767.
     * @param data
     * @return
     */
    public byte[] encodePairs(int[][] data) {
    	int maxID = 0;
    	int maxCount = 0;
    	for (int[] pair:data) {
    		maxID = Math.max(maxID, pair[0]);
    		maxCount = Math.max(maxCount, pair[1]);
    		}
    	int idBits = getNeededBits(maxID);
    	int countBits = getNeededBits(maxCount);
    	int requiredBits = getNeededBits(BITS)+2*PAIR_BITS+data.length*(idBits+countBits);
    	encodeStart(requiredBits);
    	encodeBits(mBytes.length * BITS - requiredBits, getNeededBits(BITS));	// store number of unused bits at the end of encoding
    	encodeBits(idBits, PAIR_BITS);
    	encodeBits(countBits, PAIR_BITS);
    	for (int[] pair:data) {
        	encodeBits(pair[0], idBits);
        	encodeBits(pair[1], countBits);
    		}
        encodeBitsEnd();
        return mBytes;
    	}

    /**
     * Decodes pairs of identifying integer with corresponding count values
     * in the form of an int[n][2], where n is the number of pairs, index=0
     * refers to the ID value and index=1 refers to the count value.
     * @param bytes
     * @return
     */
    public int[][] decodePairs(byte[] bytes) {
        if (bytes.length == 0)
            return null;

        decodeStart(bytes);
        int unusedBits = decodeBits(getNeededBits(BITS));
        int idBits = decodeBits(PAIR_BITS);
        int countBits = decodeBits(PAIR_BITS);
        if (idBits + countBits == 0)
        	return new int[0][2];
        int length = (BITS * bytes.length - getNeededBits(BITS) - 2*PAIR_BITS - unusedBits) / (idBits+countBits);
        int[][] data = new int[length][2];
        for (int i=0; i<data.length; i++) {
            data[i][0] = decodeBits(idBits);
            data[i][1] = decodeBits(countBits);
        	}
        return data;
    	}

    /**
     * Decodes pairs of identifying integer with corresponding count values
     * in the form of an int[n][2], where n is the number of pairs, index=0
     * refers to the ID value and index=1 refers to the count value.
     * @param s
     * @return
     */
    public int[][] decodePairs(String s) {
    	return decodePairs(s.getBytes(StandardCharsets.UTF_8));
    	}

	private int getNeededBits(int no) {
		int bits = 0;
		while (no > 0) {
			no >>= 1;
			bits++;
			}
		return bits;
		}

	private void encodeStart(int bitCount) {
        mBytes = new byte[(bitCount + BITS - 1) / BITS];
        mAvailableBits = BITS;
        mByteIndex = 0;
        }

    private void encodeBits(int data, int bits) {
        int mask = (bits == 0) ? 0 : 1 << (bits - 1);
        while (mask != 0) {
            if (mAvailableBits == 0) {
                mBytes[mByteIndex] = sCode[mBytes[mByteIndex]];
                mByteIndex++;
                mAvailableBits = BITS;
                }
            mBytes[mByteIndex] <<= 1;
            if ((data & mask) != 0)
                mBytes[mByteIndex] |= 1;
            mask >>>= 1;
            mAvailableBits--;
            }
        }

	private void encodeBits(long data, int bits) {
		long mask = (bits == 0) ? 0L : 1L << (bits - 1);
		while (mask != 0) {
			if (mAvailableBits == 0) {
				mBytes[mByteIndex] = sCode[mBytes[mByteIndex]];
				mByteIndex++;
				mAvailableBits = BITS;
			}
			mBytes[mByteIndex] <<= 1;
			if ((data & mask) != 0)
				mBytes[mByteIndex] |= 1;
			mask >>>= 1;
			mAvailableBits--;
			}
    	}

	private void encodeBitsEnd() {
        mBytes[mByteIndex] <<= mAvailableBits;
        mBytes[mByteIndex] = sCode[mBytes[mByteIndex]];
        }

    private void decodeStart(byte[] bytes) {
        mBytes = bytes;
        mByteIndex = 0;
        mTempData = sDecode[mBytes[0]];
		mTempDataLong = sDecode[mBytes[0]];
        mByteMask = 1 << (BITS - 1);
        }

    private int decodeBits(int bits) {
        int data = 0;
        while (bits != 0) {
            if (mByteMask == 0) {
                mByteIndex++;
                mTempData = sDecode[mBytes[mByteIndex]];
                mByteMask = 1 << (BITS - 1);
                }
            data <<= 1;
            if ((mTempData & mByteMask) != 0)
                data |= 1;
            mByteMask >>>= 1;
            bits--;
            }
        return data;
        }

	private long decodeBitsLong(int bits) {
		long data = 0L;
		while (bits != 0) {
			if (mByteMask == 0) {
				mByteIndex++;
				mTempDataLong = sDecode[mBytes[mByteIndex]];
				mByteMask = 1 << (BITS - 1);
			}
			data <<= 1;
			if ((mTempDataLong & (long)mByteMask) != 0)
				data |= 1L;
			mByteMask >>>= 1;
			bits--;
			}
		return data;
		}

	private byte[] decodeDuplicateBytes(byte[] bytes) {
        int length = bytes.length;
        for (int i=0; i<bytes.length; i++) {
            if (sDecodeMultiple[bytes[i]] != 0)
                length += Math.abs(sDecodeMultiple[bytes[i]]) - 1;
            }
        if (length == bytes.length)
            return bytes;

        byte[] newBytes = new byte[length];
        int oldIndex = 0;
        int newIndex = 0;
        while (oldIndex<bytes.length) {
            if (sDecodeMultiple[bytes[oldIndex]] != 0) {
                int count = Math.abs(sDecodeMultiple[bytes[oldIndex]]);
                byte code = (sDecodeMultiple[bytes[oldIndex]] < 0) ? sCode[0] : sCode[sCode.length-1];
                for (int i=0; i<count; i++)
                    newBytes[newIndex++] = code;
                oldIndex++;
                }
            else
                newBytes[newIndex++] = bytes[oldIndex++];
            }
        return newBytes;
        }

    private void encodeDuplicateBytes() {
        int length = encodeDuplicateBytes(sCode[0], sCodeMultipleMin, mBytes.length);
        length = encodeDuplicateBytes(sCode[sCode.length-1], sCodeMultipleMax, length);
        byte[] oldBytes = mBytes;
        mBytes = new byte[length];
        for (int i=0; i<length; i++)
            mBytes[i] = oldBytes[i];
        }

    private int encodeDuplicateBytes(byte code, byte[] replacement, int length) {
        int oldIndex = 0;
        int newIndex = 0;
        while (oldIndex<length) {
            if (mBytes[oldIndex] == code) {
                int count = 1;
                for (int i=oldIndex+1; i<length && mBytes[i] == code; i++)
                    count++;
                while (count > replacement.length+1) {
                    mBytes[newIndex++] = replacement[replacement.length-1];
                    oldIndex += replacement.length+1;
                    count -= replacement.length+1;
                    }
                if (count > 1) {
                    mBytes[newIndex++] = replacement[count-2];
                    oldIndex += count;
                    continue;
                    }
                }

            mBytes[newIndex++] = mBytes[oldIndex++];
            }
        return newIndex;
        }
    }
