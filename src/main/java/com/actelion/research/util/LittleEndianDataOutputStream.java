package com.actelion.research.util;
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

import java.io.DataOutput;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;

public class LittleEndianDataOutputStream implements DataOutput
{
    private DataOutputStream delegate;
    private byte[] buffer;

    public LittleEndianDataOutputStream(OutputStream out)
    {
        this.delegate = new DataOutputStream(out);
        buffer = new byte[8];
    }

    public final void writeShort(int v) throws IOException
    {
        buffer[0] = (byte) v;
        buffer[1] = (byte) (v >> 8);
        delegate.write(buffer, 0, 2);
    }

    public final void writeChar(int v) throws IOException
    {
        // same code as writeShort
        buffer[0] = (byte) v;
        buffer[1] = (byte) (v >> 8);
        delegate.write(buffer, 0, 2);
    }

    public final void writeInt(int v) throws IOException
    {
        buffer[0] = (byte) v;
        buffer[1] = (byte) (v >> 8);
        buffer[2] = (byte) (v >> 16);
        buffer[3] = (byte) (v >> 24);
        delegate.write(buffer, 0, 4);
    }

    public final void writeLong(long v) throws IOException
    {
        buffer[0] = (byte) v;
        buffer[1] = (byte) (v >> 8);
        buffer[2] = (byte) (v >> 16);
        buffer[3] = (byte) (v >> 24);
        buffer[4] = (byte) (v >> 32);
        buffer[5] = (byte) (v >> 40);
        buffer[6] = (byte) (v >> 48);
        buffer[7] = (byte) (v >> 56);
        delegate.write(buffer, 0, 8);
    }

    public final void writeFloat(float v) throws IOException
    {
        writeInt(Float.floatToIntBits(v));
    }

    public final void writeDouble(double v) throws IOException
    {
        writeLong(Double.doubleToLongBits(v));
    }

    public final void writeChars(String s) throws IOException
    {
        int len = s.length();

        for (int i = 0; i < len; i++) {
            writeChar(s.charAt(i));
        }
    }

    public final synchronized void write(int b) throws IOException
    {
        delegate.write(b);
    }

    public final synchronized void write(byte[] b, int off, int len)
        throws IOException
    {
        delegate.write(b, off, len);
    }

    public void flush() throws IOException
    {
        delegate.flush();
    }

    public final void writeBoolean(boolean v) throws IOException
    {
        delegate.writeBoolean(v);
    }

    public final void writeByte(int v) throws IOException
    {
        delegate.writeByte(v);
    }

    public final void writeBytes(String s) throws IOException
    {
        delegate.writeBytes(s);
    }

    public final void writeUTF(String str) throws IOException
    {
        delegate.writeUTF(str);
    }

    public final int size()
    {
        return delegate.size();
    }

    public final void write(byte[] b) throws IOException
    {
        delegate.write(b, 0, b.length);
    }

    public final void close() throws IOException
    {
        delegate.close();
    }
}
