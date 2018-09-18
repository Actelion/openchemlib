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

public class ByteArray
{
    private byte data[];

    ByteArray(byte data[])
    {
        this.data = data;
    }

    byte[] byteValue()
    {
        return data;
    }

    int compareTo(ByteArray o)
    {
        int len = data.length;
        int lenb = o.data.length;
        for (int i = 0;; i++) {
            int a = 0, b = 0;
            if (i < len) {
                a = ((int) data[i]) & 0xff;
            } else if (i >= lenb) {
                return 0;
            }
            if (i < lenb) {
                b = ((int) o.data[i]) & 0xff;
            }
            if (a > b) {
                return 1;
            }
            if (b > a) {
                return -1;
            }
        }
    }

    int compareTo(ByteArray o,int maxlen)
    {
        int len = Math.min(data.length,maxlen);
        int lenb = Math.min(o.data.length,maxlen);
        if (maxlen > data.length || maxlen > o.data.length)
            throw new RuntimeException("Cannot compare bytearray!");
        
        for (int i = 0; i < maxlen; i++) {
            int a = 0, b = 0;
            a = ((int) data[i]) & 0xff;
            b = ((int) o.data[i]) & 0xff;
            if (a > b) {
                return 1;
            }
            if (b > a) {
                return -1;
            }
        }
        return 0;
    }
    
    public static int compareArrays(byte[] rhs,  byte[] lhs)
    {
        ByteArray r = new ByteArray(rhs);
        ByteArray l = new ByteArray(lhs);
        return r.compareTo(l);
    }
    public static int compareArrays(byte[] rhs,  byte[] lhs, int maxlen)
    {
        ByteArray r = new ByteArray(rhs);
        ByteArray l = new ByteArray(lhs);
        return r.compareTo(l,maxlen);
    }
}
