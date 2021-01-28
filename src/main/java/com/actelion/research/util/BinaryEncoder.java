/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
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
