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
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.BurtleHasher;

public class DescriptorHandlerSkeletonSpheres implements DescriptorHandler<byte[], StereoMolecule> {
    private static final double CORRECTION_FACTOR = 0.7;

    private static final byte[] FAILED_OBJECT = new byte[0];
    private static final int MAX_SPHERE_COUNT = 5;
    private static final int EXACT_SPHERE_COUNT = 4;
    private static final int SKELETON_SPHERE_COUNT = 5;
    private static final int HASH_BITS = 10;
    private static final int HASH_INIT = 13;
    private static final int DESCRIPTOR_SIZE = (1 << HASH_BITS);

    private static DescriptorHandlerSkeletonSpheres sDefaultInstance;

    public static DescriptorHandlerSkeletonSpheres getDefaultInstance() {
    	synchronized(DescriptorHandlerSkeletonSpheres.class) {
    		if (sDefaultInstance == null) {
        		sDefaultInstance = new DescriptorHandlerSkeletonSpheres();
        		}
        	}
        return sDefaultInstance;
    	}

    public boolean calculationFailed(byte[] o) {
        return o==null || o.length == 0;
        }

    /**
     * This descriptor requires proper up/down bonds, because it encodes stereo parities.
     * If a passed molecule is generated from idcode parsing, make sure that coordinates
     * and up/down/bonds are available, i.e. that the IDCodeParser was instantiated with
     * the respective option.
     */
    public byte[] createDescriptor(StereoMolecule mol) {
        mol.ensureHelperArrays(Molecule.cHelperRings);
        StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());

		byte[] descriptor = new byte[DESCRIPTOR_SIZE];

//System.out.println("descriptor skeleton spheres:");
        int[] atomList = new int[mol.getAtoms()];
        boolean[] atomMask = new boolean[mol.getAtoms()];
        for (int rootAtom=0; rootAtom<mol.getAtoms(); rootAtom++) {
            if (rootAtom != 0)
                Arrays.fill(atomMask, false);

            int min = 0;
            int max = 0;

            for (int sphere=0; sphere<MAX_SPHERE_COUNT && max<mol.getAtoms(); sphere++) {
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
                if (sphere < EXACT_SPHERE_COUNT) {
                    String idcode = new Canonizer(fragment).getIDCode();
                    int h = BurtleHasher.hashlittle(idcode, HASH_INIT);
                    h = (h & BurtleHasher.hashmask(HASH_BITS));
                    if (descriptor[h] < DescriptorEncoder.MAX_COUNT_VALUE)
                    	descriptor[h]++;
//System.out.println("atom:"+rootAtom+"\tfragment\tradius:"+sphere+"\thash:"+h+"\t"+idcode);
                    }

                // take atomic no reduced fragment skeleton also
                if (sphere < SKELETON_SPHERE_COUNT) {
                    for (int atom=0; atom<fragment.getAllAtoms(); atom++)
                        fragment.setAtomicNo(atom, 6);
                    String idcode = new Canonizer(fragment).getIDCode();
                    int h = BurtleHasher.hashlittle(idcode, HASH_INIT);
                    h = (h & BurtleHasher.hashmask(HASH_BITS));
                    if (descriptor[h] < DescriptorEncoder.MAX_COUNT_VALUE)
                    	descriptor[h]++;
//System.out.println("atom:"+rootAtom+"\tskeleton\tradius:"+sphere+"\thash:"+h+"\t"+idcode);
                    }
                }
            }

        return descriptor;
        }

    public byte[] decode(String s) {
        return s == null ?               null
             : s.equals(FAILED_STRING) ? FAILED_OBJECT
             :                           new DescriptorEncoder().decodeCounts(s);
        }

    public byte[] decode(byte[] bytes) {
        return bytes == null ?               		null
             : Arrays.equals(bytes, FAILED_BYTES) ? FAILED_OBJECT
             :                           			new DescriptorEncoder().decodeCounts(bytes);
        }

    public String encode(byte[] o) {
        return calculationFailed(o) ? FAILED_STRING
             : new String(new DescriptorEncoder().encodeCounts(o));
        }

    public DescriptorInfo getInfo() {
        return DescriptorConstants.DESCRIPTOR_SkeletonSpheres;
        }

    public String getVersion() {
        return "1.1";
        }

    public float getSimilarity(byte[] d1, byte[] d2) {
        if (d1 == null || d2 == null)
            return Float.NaN;

        int total = 0;
        int matching = 0;
        for (int i=0; i<d1.length; i++) {
            total += Math.max(d1[i], d2[i]);
            matching += Math.min(d1[i], d2[i]);
            }
/*
if (((double)matching/(double)total) > 0.8) {
    System.out.print("i:");
    for (int i=0; i<d1.length; i++)
        if (d1[i] != d2[i])
            System.out.print(" "+i);
    System.out.println();
    System.out.print("d1:");
    for (int i=0; i<d1.length; i++)
        if (d1[i] != d2[i])
            System.out.print(" "+d1[i]);
    System.out.println();
    System.out.print("d2:");
    for (int i=0; i<d2.length; i++)
        if (d1[i] != d2[i])
            System.out.print(" "+d2[i]);
    System.out.println();
    }
*/        
        return normalizeValue((double)matching/(double)total);
        }

	private float normalizeValue(double value) {
		return value <= 0.0f ? 0.0f
			 : value >= 1.0f ? 1.0f
			 : (float)(1.0-Math.pow(1-Math.pow(value, CORRECTION_FACTOR) ,1.0/CORRECTION_FACTOR));
		}

    public DescriptorHandler<byte[], StereoMolecule> getDeepCopy() {
		return new DescriptorHandlerSkeletonSpheres();
    	}
	}
