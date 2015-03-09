/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
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

package com.actelion.research.chem.descriptor;

import java.util.Arrays;

import com.actelion.research.chem.*;

abstract public class AbstractDescriptorHandlerFP<U extends Object> implements DescriptorHandler<int[], U> {
    protected static final int[] FAILED_OBJECT = new int[0];

    public String encode(int[] o) {
        return calculationFailed(o) ? FAILED_STRING
                                    : new String(new DescriptorEncoder().encode(o));
    	}

    public int[] decode(String s) {
        return s == null ?               null
             : s.equals(FAILED_STRING) ? FAILED_OBJECT
             :                           new DescriptorEncoder().decode(s);
    	}

    public int[] decode(byte[] bytes) {
        return bytes == null ?           			null
             : Arrays.equals(bytes, FAILED_BYTES) ? FAILED_OBJECT
             :                        				new DescriptorEncoder().decode(bytes);
    	}

    public boolean calculationFailed(int[] o) {
        return o==null || o.length == 0;
    	}

    public float getSimilarity(int[] o1, int[] o2) {
        return o1 == null
            || o2 == null
            || o1.length == 0
            || o2.length == 0 ? 0.0f
               : SSSearcherWithIndex.getSimilarityTanimoto((int[])o1, (int[])o2);
    	}
	}
