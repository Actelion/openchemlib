/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * 3. Neither the name of the copyright holder nor the
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
 * @author Thomas Sander
 */

package com.actelion.research.util;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StringWriter;

public class BinaryEncoder {
	private BufferedWriter	mWriter;
	private int				mFreeBufferBits,mData,mCharsWritten,mDataBitCount,
							mBitsPerCharacter,mBaseCharacter,mLineLength;
	private int[]			mMask;

	/**
	 * Convenience method to directly encode an int array into a String
	 * Uses 6 bits per character using 64 characters from '@' (ASCII 64).
	 * @param data
	 * @param dataBitCount
	 * @return
	 */
	public static String toString(int[] data, int dataBitCount) {
		StringWriter writer = new StringWriter();
		BinaryEncoder encoder = new BinaryEncoder(new BufferedWriter(writer));
		try {
			encoder.initialize(dataBitCount, data.length);
			for (int d:data)
				encoder.write(d);
			encoder.finalize();
			} catch (IOException ioe) {}
		return writer.toString();
		}

	/**
	 * Convenience method to directly encode a byte array into a String
	 * Uses 6 bits per character using 64 characters from '@' (ASCII 64).
	 * @param data
	 * @param dataBitCount
	 * @return
	 */
	public static String toString(byte[] data, int dataBitCount) {
		StringWriter writer = new StringWriter();
		BinaryEncoder encoder = new BinaryEncoder(new BufferedWriter(writer));
		try {
			encoder.initialize(dataBitCount, data.length);
			for (int d:data)
				encoder.write(d);
			encoder.finalize();
			} catch (IOException ioe) {}
		return writer.toString();
		}

	/**
	 * Creates an encoder to make a BufferedWriter write byte/int arrays as Strings.
	 * This version uses 6 bits per character using 64 characters from '@' (ASCII 64).
	 * It inserts a platform specific newline/linefeed after every 80 characters.
	 * @param writer
	 */
	public BinaryEncoder(BufferedWriter writer) {
		this(writer, 6, 64, 80);
		}

	/**
	 * Instantiates an encoder for byte/int arrays to String.
	 * @param writer
	 * @param bitsPerCharacter no of bits to be used per character
	 * @param baseCharacter lowest character, which encodes value '0'
	 * @param lineLength -1 (no LF/NL) or count of character after which a NL/LF is inserted
	 */
	public BinaryEncoder(BufferedWriter writer, int bitsPerCharacter,
									  int baseCharacter, int lineLength) {
		mWriter = writer;
		mBitsPerCharacter = bitsPerCharacter;
		mBaseCharacter = baseCharacter;
		mLineLength = lineLength;
		mData = 0;
		mCharsWritten = 0;
		mFreeBufferBits = bitsPerCharacter;

		mMask = new int[bitsPerCharacter+1];
		for (int i=1; i<=bitsPerCharacter; i++)
			mMask[i] = (mMask[i-1] << 1) | 1;
		}

	/**
	 * Initializes encoder for writing one binary array
	 * @param dataBitCount count of bits to be written per data value
	 * @param totalByteCount count of data values to be encoded
	 * @throws IOException
	 */
	public void initialize(int dataBitCount, int totalByteCount) throws IOException {
		mDataBitCount = 32;
		write(totalByteCount);
		mDataBitCount = dataBitCount;
		}

	/**
	 * Encodes and writes one data value
	 * @param data
	 * @throws IOException
	 */
	public void write(int data) throws IOException {
		int remainingSourceBits = mDataBitCount;
		while (remainingSourceBits != 0) {
			int bitsToWrite = Math.min(mFreeBufferBits, remainingSourceBits);
			int dataBits = data & (mMask[bitsToWrite] << (remainingSourceBits - bitsToWrite));
			mData |= (remainingSourceBits > mFreeBufferBits) ?
						dataBits >> (remainingSourceBits - mFreeBufferBits)
					  : dataBits << (mFreeBufferBits - remainingSourceBits);
			remainingSourceBits -= bitsToWrite;
			mFreeBufferBits -= bitsToWrite;
			if (mFreeBufferBits == 0) {
				mWriter.write(mBaseCharacter+mData);
				mFreeBufferBits = mBitsPerCharacter;
				mData = 0;
				mCharsWritten++;
				if (mLineLength != -1 && (mCharsWritten % mLineLength == 0))
					mWriter.newLine();
				}
			}
		}

	/**
	 * Writes remaining bits from buffer, if there are any and possibly adds final NL/LF.
	 */
	public void finalize() throws IOException {
		if (mFreeBufferBits < mBitsPerCharacter) {
			mWriter.write(mBaseCharacter+mData);
			if (mLineLength != -1)
				mWriter.newLine();
			}
		else if (mLineLength != -1 && (mCharsWritten % mLineLength != 0)) {
			mWriter.newLine();
			}
		mWriter.flush();
		}
	}
