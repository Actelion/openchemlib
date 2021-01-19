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
