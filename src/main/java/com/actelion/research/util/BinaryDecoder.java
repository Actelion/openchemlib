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

package com.actelion.research.util;

import java.io.*;

public class BinaryDecoder {
	private BufferedReader	mReader;
	private int				mAvailableBufferBits,mData,mDataBitCount,
							mBitsPerCharacter,mBaseCharacter;
	private int[]			mMask;

	/**
	 * Convenience method to directly decode a String-encoded byte array.
	 * @param encodedBytes
	 * @param dataBitCount count of bits per byte
	 * @return
	 */
	public static byte[] toBytes(String encodedBytes, int dataBitCount) {
		BinaryDecoder decoder = new BinaryDecoder(new BufferedReader(new StringReader(encodedBytes)));
		try {
			byte[] data = new byte[decoder.initialize(dataBitCount)];
			for (int i=0; i<data.length; i++)
				data[i] = (byte)decoder.read();
			return data;
			} catch (IOException ioe) {}
		return null;
		}

	/**
	 * Convenience method to directly decode a String-encoded int array.
	 * @param encodedBytes
	 * @param dataBitCount count of bits per byte
	 * @return
	 */
	public static int[] toInts(String encodedBytes, int dataBitCount) {
		BinaryDecoder decoder = new BinaryDecoder(new BufferedReader(new StringReader(encodedBytes)));
		try {
			int[] data = new int[decoder.initialize(dataBitCount)];
			for (int i=0; i<data.length; i++)
				data[i] = decoder.read();
			return data;
			} catch (IOException ioe) {}
		return null;
		}

	public BinaryDecoder(BufferedReader reader) {
		this(reader, 6, 64);
		}

	public BinaryDecoder(BufferedReader reader, int bitsPerCharacter, int baseCharacter) {
		mReader = reader;
		mBitsPerCharacter = bitsPerCharacter;
		mBaseCharacter = baseCharacter;
		mAvailableBufferBits = 0;

		mMask = new int[bitsPerCharacter+1];
		for (int i=1; i<=bitsPerCharacter; i++)
			mMask[i] = (mMask[i-1] << 1) | 1;
		}

	public int initialize(int dataBitCount) throws IOException {
		mDataBitCount = 32;
		int totalByteCount = read();
		mDataBitCount = dataBitCount;
		return totalByteCount;
		}

	public int read() throws IOException {
		int data = 0;
		int neededBits = mDataBitCount;
		while (neededBits != 0) {
			if (mAvailableBufferBits == 0) {
				while ((mData = mReader.read() - mBaseCharacter) < 0);
				mAvailableBufferBits = mBitsPerCharacter;
				}
			int bitsToRead = Math.min(mAvailableBufferBits, neededBits);
			int dataBits = mData & (mMask[bitsToRead] << (mAvailableBufferBits - bitsToRead));
			data |= (mAvailableBufferBits > neededBits) ?
						dataBits >> (mAvailableBufferBits - neededBits)
					  : dataBits << (neededBits - mAvailableBufferBits);
			neededBits -= bitsToRead;
			mAvailableBufferBits -= bitsToRead;
			}
		return data;
		}
	}
