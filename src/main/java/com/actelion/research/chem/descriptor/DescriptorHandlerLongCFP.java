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

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.BurtleHasher;

import java.util.Arrays;

public class DescriptorHandlerLongCFP extends AbstractDescriptorHandlerLongFP<StereoMolecule> {
    private static final double CORRECTION_FACTOR = 0.6;

    private static DescriptorHandlerLongCFP sDefaultInstance;
    
    private static final int SPHERE_COUNT = 5;
    private static final int HASH_BITS = 10;
    private static final int HASH_INIT = 13;
    private static final int DESCRIPTOR_SIZE = (1 << HASH_BITS);

    public static DescriptorHandlerLongCFP getDefaultInstance() {
    	synchronized(DescriptorHandlerLongCFP.class) {
    		if (sDefaultInstance == null) {
                sDefaultInstance = new DescriptorHandlerLongCFP();
        	}
        }
        return sDefaultInstance;
    }

    @Override
    public DescriptorInfo getInfo() {
        return DescriptorConstants.DESCRIPTOR_HashedCFp;
    }

    @Override
    public String getVersion() {
        return DescriptorConstants.DESCRIPTOR_HashedCFp.version;
    }

    /**
     * This descriptor requires proper up/down bonds, because it encodes stereo parities.
     * If a passed molecule is generated from idcode parsing, make sure that coordinates
     * and up/down/bonds are available, i.e. that the IDCodeParser was instantiated with
     * the respective option.
     */
    @Override
    public long[] createDescriptor(StereoMolecule mol) {
        if (mol ==null)
            return null;

        mol.ensureHelperArrays(Molecule.cHelperRings);
        StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());

//System.out.println("atom\tsphere\tbit\tidcode");
        final int len = DESCRIPTOR_SIZE / Long.SIZE;
        long[] data = new long[len];

		int[] atomList = new int[mol.getAtoms()];
        boolean[] atomMask = new boolean[mol.getAtoms()];
        for (int rootAtom=0; rootAtom<mol.getAtoms(); rootAtom++) {
            if (rootAtom != 0)
                Arrays.fill(atomMask, false);

            int min = 0;
            int max = 0;
            for (int sphere=0; sphere<SPHERE_COUNT && max<mol.getAtoms(); sphere++) {
                if (max == 0) {
                    atomList[0] = rootAtom;
                    atomMask[rootAtom] = true;
                    max = 1;
                    }
                else {
                    int newMax = max;
                    for (int i=min; i<max; i++) {
                        int atom = atomList[i];
                        for (int j=0; j<mol.getConnAtoms(atom); j++) {
                            int connAtom = mol.getConnAtom(atom, j);
                            if (!atomMask[connAtom]) {
                                atomMask[connAtom] = true;
                                atomList[newMax++] = connAtom;
                                }
                            }
                        }
                    min = max;
                    max = newMax;
                    }

                mol.copyMoleculeByAtoms(fragment, atomMask, true, null);

                // take fragment as it is
                String idcode = new Canonizer(fragment).getIDCode();
                int h = BurtleHasher.hashlittle(idcode, HASH_INIT);
                h = (h & BurtleHasher.hashmask(HASH_BITS));

	            int index = len - h / Long.SIZE - 1;
	            int bitNo = h % 32;     // we need this strange 32 bit block handling to be compatible with the older 32-bit version
	            if (h % 64 >= 32)
	            	bitNo += 32;
	            data[index] |= (1L << bitNo);
//System.out.println(rootAtom+"\t"+sphere+"\t"+h+"\t"+idcode);
                }
            }

        return data;
        }

    @Override
    public float getSimilarity(long[] o1, long[] o2) {
        return o1 == null
            || o2 == null
            || o1.length == 0
            || o2.length == 0 ? 0.0f
        : normalizeValue(SSSearcherWithIndex.getSimilarityTanimoto(o1, o2));
    }
    
	private float normalizeValue(double value) {
		return value <= 0.0f ? 0.0f
			 : value >= 1.0f ? 1.0f
			 : (float)(1.0-Math.pow(1-Math.pow(value, CORRECTION_FACTOR) ,1.0/CORRECTION_FACTOR));
	}

    @Override
    public DescriptorHandler<long[], StereoMolecule> getThreadSafeCopy() {
    	return this;
    }

}
