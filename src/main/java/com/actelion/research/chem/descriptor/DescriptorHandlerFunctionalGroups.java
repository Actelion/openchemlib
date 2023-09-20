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

package com.actelion.research.chem.descriptor;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.FunctionalGroupClassifier;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.TreeSet;

public class DescriptorHandlerFunctionalGroups implements DescriptorHandler<int[][], StereoMolecule> {
    private static final double CORRECTION_FACTOR = 0.7;

    private static final int[][] FAILED_OBJECT = { { 1, 0 } };

    private static DescriptorHandlerFunctionalGroups sDefaultInstance;

    public static DescriptorHandlerFunctionalGroups getDefaultInstance() {
    	synchronized(DescriptorHandlerFunctionalGroups.class) {
    		if (sDefaultInstance == null) {
        		sDefaultInstance = new DescriptorHandlerFunctionalGroups();
        		}
        	}
        return sDefaultInstance;
    	}

    public boolean calculationFailed(int[][] d) {
        return d==null || d.length == 1 && d[0][1] == 0;
        }

    /**
     * This descriptor contains non-hashed counts of predefined overlapping
     * small fragments of organic functional groups. Areas of a molecule that
     * consist of aliphatic carbon atoms without pi-electrons or hetero atoms
     * are neglected by the descriptor generation.
     */
    public int[][] createDescriptor(StereoMolecule mol) {
	    if (mol ==null)
		    return null;

        FunctionalGroupClassifier fgc = new FunctionalGroupClassifier(mol);
		return fgc.getOrganicFunctionalGroupCounts();
        }

    public int[][] decode(String s) {
        return s == null ?               null
             : s.equals(FAILED_STRING) ? FAILED_OBJECT
             :                           new DescriptorEncoder().decodePairs(s);
        }

    public int[][] decode(byte[] bytes) {
        return bytes == null ?               		null
             : Arrays.equals(bytes, FAILED_BYTES) ? FAILED_OBJECT
             :                           			new DescriptorEncoder().decodePairs(bytes);
        }

    public String encode(int[][] d) {
        return calculationFailed(d) ? FAILED_STRING
             : new String(new DescriptorEncoder().encodePairs(d), StandardCharsets.UTF_8);
        }

    public DescriptorInfo getInfo() {
        return DescriptorConstants.DESCRIPTOR_OrganicFunctionalGroups;
        }

    public String getVersion() {
        return DescriptorConstants.DESCRIPTOR_OrganicFunctionalGroups.version;
        }

    public float getSimilarity(int[][] dl1, int[][] dl2) {
        if (dl1 == null || dl2 == null)
            return Float.NaN;

        if (dl1.length == 0 && dl2.length == 0)
        	return 1f;
        if (dl1.length == 0 || dl2.length == 0)
        	return 0f;

        TreeSet<Match> matchList = new TreeSet<Match>();
        int i1 = 0;
        int i2 = -1;
        for (int[] d1:dl1) {
            i2 = 0;
            for (int[] d2:dl2) {
            	int matchLevel = FunctionalGroupClassifier.getFunctionalGroupEquivalenceLevel(d1[0], d2[0]);
            	if (matchLevel != -1) {
                	for (int i=0; i<d1[1]; i++) {
                    	for (int j=0; j<d2[1]; j++) {
                    		matchList.add(new Match(i1+i, i2+j, matchLevel));
                    		}
                		}
            		}
                i2 += d2[1];
            	}
            i1 += d1[1];
            }

        float total = i1 + i2;
        float matching = 0;
        boolean[] used1 = new boolean[i1];
        boolean[] used2 = new boolean[i2];
        for (Match match:matchList) {
        	if (!used1[match.fg1] && !used2[match.fg2]) {
        		float m = 1.0f - 0.1f * match.level;
        		matching += m;
        		total -= m;
        		used1[match.fg1] = true;
        		used2[match.fg2] = true;
        		}
        	}

        return normalizeValue(matching/total);
        }

	private float normalizeValue(double value) {
		return value <= 0.0f ? 0.0f
			 : value >= 1.0f ? 1.0f
			 : (float)(1.0-Math.pow(1-Math.pow(value, CORRECTION_FACTOR) ,1.0/CORRECTION_FACTOR));
		}

    public DescriptorHandler<int[][], StereoMolecule> getThreadSafeCopy() {
		return new DescriptorHandlerFunctionalGroups();
    	}

    private class Match implements Comparable<Match> {
    	int fg1,fg2,level;

    	Match(int fg1, int fg2, int level) {
    		this.fg1 = fg1;
    		this.fg2 = fg2;
    		this.level = level;
    		}

		@Override
		public int compareTo(Match o) {
			if (level != o.level)
				return (level < o.level) ? -1 : 1;

			int v = (fg1 < fg2) ? (fg1 << 10) + fg2 : (fg2 << 10) + fg1;
			int ov = (o.fg1 < o.fg2) ? (o.fg1 << 10) + o.fg2 : (o.fg2 << 10) + o.fg1;
			return (v < ov) ? -1 :  (v > ov) ? 1 : 0;
			}
    	}
	}
