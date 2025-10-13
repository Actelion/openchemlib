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


public class LittleEndianDataInputStream implements DataInput
{
	private DataInputStream dataInputStream; // to get at high level readFully methods of DataInputStream
	private InputStream inputStream; // to get at the low-level read methods of InputStream
	private byte[] workBuffer; // work array for buffering input

	public LittleEndianDataInputStream(InputStream inputStream)
	{
		this.inputStream = inputStream;
		this.dataInputStream = new DataInputStream(inputStream);
		workBuffer = new byte[8];
	}

	public final short readShort() throws IOException
	{
		dataInputStream.readFully(workBuffer, 0, 2);

		return (short)(((workBuffer[1] & 0xff) << 8) | (workBuffer[0] & 0xff));
	}

	public final int readUnsignedShort() throws IOException
	{
		dataInputStream.readFully(workBuffer, 0, 2);

		return (((workBuffer[1] & 0xff) << 8) | (workBuffer[0] & 0xff));
	}

	public final char readChar() throws IOException
	{
		dataInputStream.readFully(workBuffer, 0, 2);

		return (char)(((workBuffer[1] & 0xff) << 8) | (workBuffer[0] & 0xff));
	}

	public final int readInt() throws IOException
	{
		dataInputStream.readFully(workBuffer, 0, 4);

		return ((workBuffer[3]) << 24) | ((workBuffer[2] & 0xff) << 16) | ((workBuffer[1] & 0xff) << 8) | (workBuffer[0] & 0xff);
	}

	public final long readLong() throws IOException
	{
		dataInputStream.readFully(workBuffer, 0, 8);

		return ((long)(workBuffer[7]) << 56) /* long cast needed or shift done modulo 32 */
		| ((long)(workBuffer[6] & 0xff) << 48) | ((long)(workBuffer[5] & 0xff) << 40) | ((long)(workBuffer[4] & 0xff) << 32)
		| ((long)(workBuffer[3] & 0xff) << 24) | ((long)(workBuffer[2] & 0xff) << 16) | ((long)(workBuffer[1] & 0xff) << 8)
		| (long)(workBuffer[0] & 0xff);
	}


	public final float readFloat() throws IOException
	{
		return Float.intBitsToFloat(readInt());
	}

	public final double readDouble() throws IOException
	{
		return Double.longBitsToDouble(readLong());
	}


	public final int read(byte[] b, int off, int len) throws IOException
	{
		return inputStream.read(b, off, len);
	}

	public final void readFully(byte[] b) throws IOException
	{
		dataInputStream.readFully(b, 0, b.length);
	}

	public final void readFully(byte[] b, int off, int len)
		throws IOException
	{
        if (len < 0)
            System.err.println("Error: Negative number of bytes to read");
		else
            dataInputStream.readFully(b, off, len);
	}


	public final int skipBytes(int n) throws IOException
	{
		return dataInputStream.skipBytes(n);
	}

	public final boolean readBoolean() throws IOException
	{
		return dataInputStream.readBoolean();
	}

	public final byte readByte() throws IOException
	{
		return dataInputStream.readByte();
	}

	public final int readUnsignedByte() throws IOException
	{
		return dataInputStream.readUnsignedByte();
	}

	@Deprecated
	public final String readLine() throws IOException
	{
		return dataInputStream.readLine();
	}

	public final String readUTF() throws IOException
	{
		return dataInputStream.readUTF();
	}

	public static String readUTF(DataInput in) throws IOException
	{
		return DataInputStream.readUTF(in);
	}

	public final void close() throws IOException
	{
		dataInputStream.close();
	}


	public int available() throws IOException
	{
		return inputStream.available();
	}

	public void mark(int i)
	{
		inputStream.mark(i);
	}

	public void reset() throws IOException
	{
		inputStream.reset();
	}

	public int read(byte[] buffer) throws IOException
	{
		return read(buffer,0,buffer.length);
	}
}
