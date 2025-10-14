package com.actelion.research.util;

import java.io.DataInputStream;
import java.nio.*;

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


import java.io.*;


public class EndianInputStream implements DataInput
{
    private ByteBuffer buffer = ByteBuffer.allocate(8);

    public EndianInputStream(ByteOrder order, byte[] bytes) throws IOException
    {
        buffer = ByteBuffer.allocate(bytes.length);
        buffer.mark();
        buffer.put(bytes);
        buffer.reset();
        buffer.order(order);
    }

    public final short readShort() throws IOException
    {
        return buffer.getShort();
    }

    public final int readUnsignedShort() throws IOException
    {
        return readShort();
    }

    public final char readChar() throws IOException
    {
        return buffer.getChar();
    }

    public final int readInt() throws IOException
    {
        return buffer.getInt();
    }

    public final long readLong() throws IOException
    {
        return buffer.getLong();
    }


    public final float readFloat() throws IOException
    {
        return buffer.getFloat();
    }

    public final double readDouble() throws IOException
    {
        return buffer.getDouble();
    }

    @Override
    public String readLine() throws IOException
    {
        StringBuilder sb = new StringBuilder();
        buffer.mark();
        boolean ok = true;
        while(ok) {
            byte b = readByte();
            switch(b) {
                case -1:
                case '\n':
                    ok = false;
                    break;
                case '\r':
                    break;
                default: sb.append((char)b);

            }
        }
        return sb.length() == 0 ? null : sb.toString();
    }


    public final int read(byte[] b, int off, int len) throws IOException
    {
        try {
            ByteBuffer byteBuffer = buffer.get(b, off, len);
            return len;
        } catch (Exception e) {
            return -1;
        }
        //return in.read(b, off, len);
    }

    public final void readFully(byte[] b) throws IOException
    {
        buffer.get(b,0,b.length);
    }

    public final void readFully(byte[] b, int off, int len)
            throws IOException
    {
        if (len < 0)
            System.err.println("Error: Negative number of bytes to read");
        else
            buffer.get(b, off, len);
    }


    public final int skipBytes(int n) throws IOException
    {
        buffer.get(new byte[n],0,n);
        return n;
    }

    public final boolean readBoolean() throws IOException
    {
        return readByte() == (byte)1;
    }

    public final byte readByte() throws IOException
    {
        try {
            return buffer.get();
        } catch (BufferUnderflowException e) {
            return -1;
        }
    }

    public final int readUnsignedByte() throws IOException
    {
        return readByte();
    }

    public final String readUTF() throws IOException
    {
        return DataInputStream.readUTF(this);
    }


    public final void close() throws IOException
    {
    }

    public void mark()
    {
        buffer.mark();
    }

    public void reset()
    {
        buffer.reset();
    }
}