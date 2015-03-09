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

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.BurtleHasher;
import com.actelion.research.util.datamodel.IntVec;

public class DescriptorHandlerHashedCFp extends AbstractDescriptorHandlerFP<StereoMolecule> {
    private static final double CORRECTION_FACTOR = 0.6;

    private static DescriptorHandlerHashedCFp sDefaultInstance;
    
    private static final int SPHERE_COUNT = 5;
    private static final int HASH_BITS = 10;
    private static final int HASH_INIT = 13;
    private static final int DESCRIPTOR_SIZE = (1 << HASH_BITS);

    public static DescriptorHandlerHashedCFp getDefaultInstance() {
    	synchronized(DescriptorHandlerHashedCFp.class) {
    		if (sDefaultInstance == null) {
                sDefaultInstance = new DescriptorHandlerHashedCFp();
        	}
        }
        return sDefaultInstance;
    }

    public DescriptorInfo getInfo() {
        return DescriptorConstants.DESCRIPTOR_HashedCFp;
    }

    public String getVersion() {
        return "2.1";
    }

    /**
     * This descriptor requires proper up/down bonds, because it encodes stereo parities.
     * If a passed molecule is generated from idcode parsing, make sure that coordinates
     * and up/down/bonds are available, i.e. that the IDCodeParser was instantiated with
     * the respective option.
     */
    public int[] createDescriptor(StereoMolecule mol) {
        mol.ensureHelperArrays(Molecule.cHelperRings);
        StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());

        // byte[] descriptor = new byte[DESCRIPTOR_SIZE];
        int len = DESCRIPTOR_SIZE / Integer.SIZE;
        IntVec iv = new IntVec(len); 

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
                iv.setBit(h);

//System.out.println("atom:"+rootAtom+"\tsphere:"+sphere+"\thash:"+h+"\t"+idcode);
                }
            }

        return iv.get();
        }

    public float getSimilarity(int[] o1, int[] o2) {
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

    public DescriptorHandler<int[], StereoMolecule> getDeepCopy() {
    	return new DescriptorHandlerHashedCFp();
    }

}
