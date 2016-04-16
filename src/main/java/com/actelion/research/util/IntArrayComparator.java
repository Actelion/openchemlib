package com.actelion.research.util;

import java.io.Serializable;
import java.util.Comparator;

public class IntArrayComparator  implements Comparator<int[]>,Serializable {
    static final long serialVersionUID = 0x20141022;
 
    public int compare(int[] ia1, int[] ia2) {
        if (ia1 == null)
            return (ia2 == null) ? 0 : 1;
        if (ia2 == null)
            return -1;
        for (int i=0; i<ia1.length; i++) {
            if (ia2.length == i)
                return 1;
            if (ia1[i] != ia2[i])
                return (ia1[i] < ia2[i]) ? -1 : 1;
            }
        return (ia2.length > ia1.length) ? -1 : 0;
		}
	}
