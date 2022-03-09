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

import java.io.Serializable;
import java.util.Comparator;

public class ByteArrayArrayComparator implements Comparator<byte[][]>,Serializable {
    static final long serialVersionUID = 0x20131015;
 
    public int compare(byte[][] b1, byte[][] b2) {
        if (b1 == null)
            return (b2 == null) ? 0 : 1;
        if (b2 == null)
            return -1;
        for (int i=0; i<b1.length; i++) {
            if (b2.length == i)
                return -1;
            if (b1[i] == null && b2[i] == null)
            	continue;
            if (b1[i] == null)
                return 1;
            if (b2[i] == null)
                return -1;
            for (int j=0; j<b1[i].length; j++) {
                if (b2[i].length == j)
                    return 1;
                if (b1[i][j] != b2[i][j])
                    return (b1[i][j] < b2[i][j]) ? -1 : 1;
            	}
            if (b2[i].length > b1[i].length)
            	return -1;
            }
        return (b2.length > b1.length) ? -1 : 0;
		}
	}
