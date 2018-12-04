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
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.util.BurtleHasher;

import java.util.Arrays;

public class DescriptorHandlerReactionFP extends AbstractDescriptorHandlerLongFP<Reaction> {
	public static final String cVersion = "1.0.0";

	public static final int REACTION_CENTER_LONG_COUNT = 8;
	public static final float REACTION_CENTER_WEIGHT = 0.8f;	// the percentage of reaction center similarity to contribute to the total similarity
	public static final float PERIPHERY_WEIGHT = 1f - REACTION_CENTER_WEIGHT;

	private static final int SPHERE_COUNT = 5;
	private static final int HASH_BITS = 10;
	private static final int HASH_INIT = 13;
	private static final int DESCRIPTOR_SIZE = (1 << HASH_BITS);

    private static DescriptorHandlerReactionFP sDefaultInstance;

    public static DescriptorHandlerReactionFP getDefaultInstance() {
        if (sDefaultInstance == null) {
        	synchronized(DescriptorHandlerReactionFP.class) {
        		sDefaultInstance = new DescriptorHandlerReactionFP();
        	}
        }
        return sDefaultInstance;
    }

    public DescriptorInfo getInfo() {
        return DescriptorConstants.DESCRIPTOR_ReactionFP;
    }

    public String getVersion() {
        return DescriptorConstants.DESCRIPTOR_ReactionFP.version;
    }

    public long[] createDescriptor(Reaction rxn) {
	    if (rxn == null)
		    return null;

	    boolean[] isReactionCenterMapNo = rxn.getReactionCenterMapNos();
		if (isReactionCenterMapNo == null)
			return FAILED_OBJECT;

		final int len = DESCRIPTOR_SIZE / Long.SIZE;
		long[] data = new long[len];

//System.out.println("mol\tisProduct\tatom\tsphere\tisRC\thash\tindex\tbit\tidcode");
		for (int m=0; m<rxn.getMolecules(); m++) {
			StereoMolecule mol = rxn.getMolecule(m);
			mol.ensureHelperArrays(Molecule.cHelperRings);
			boolean[] isReactionCenterAtom = new boolean[mol.getAllAtoms()];
			int reactionCenterAtomCount = rxn.getReactionCenterAtoms(m, isReactionCenterMapNo, isReactionCenterAtom, null);

			StereoMolecule fragment = new StereoMolecule(mol.getAllAtoms(), mol.getAllBonds());

			int[] atomList = new int[mol.getAllAtoms()];
			boolean[] atomMask = new boolean[mol.getAllAtoms()];
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
								if (!atomMask[connAtom]
								 && (isReactionCenterAtom[rootAtom] || !isReactionCenterAtom[connAtom])) {
									atomMask[connAtom] = true;
									atomList[newMax++] = connAtom;
									}
								}
							}

						if (max == newMax)
							break;

						min = max;
						max = newMax;
						}

					mol.copyMoleculeByAtoms(fragment, atomMask, true, null);

					// take fragment as it is
					String idcode = new Canonizer(fragment).getIDCode();

					// we want product fragments to create a different set of hash codes
					if (m >= rxn.getReactants())
						idcode = idcode.concat("P");

					// For reaction center atoms smaller fragments have higher weights, which is achieved by a higher number of representing bits
					int weight = isReactionCenterAtom[rootAtom] ? SPHERE_COUNT - sphere : 1;
					for (int i=0; i<weight; i++) {
						int h = BurtleHasher.hashlittle(idcode.concat(Integer.toString(i)), HASH_INIT);
						h = (h & BurtleHasher.hashmask(HASH_BITS-1));		// '-1' divides the full bit span into two equal parts for reaction center and non-reaction center

						// we put the non-reaction-center bits at the end; high h values translate into low indexes
						if (isReactionCenterAtom[rootAtom])
							h += (1 << (HASH_BITS-1));

						int index = len - h / Long.SIZE - 1;
						data[index] |= (1L << (h % Long.SIZE));
//System.out.println(""+m+"\t"+Boolean.toString(m >= rxn.getReactants())+"\t"+rootAtom+"\t"+sphere+"\t"+isReactionCenterAtom[rootAtom]+"\t"+h+"\t"+index+"\t"+(h % Long.SIZE)+"\t"+idcode);
						}
					}
				}
			}

		return data;
//		return (descriptor == null) ? FAILED_OBJECT : descriptor;
        }

	public float getSimilarity(long[] o1, long[] o2) {
		if (o1 == null
		 || o2 == null
		 || o1.length == 0
		 || o2.length == 0)
			return 0.0f;

		float reactionCenterSimilarity = getSimilarityTanimoto(o1, o2, 0, REACTION_CENTER_LONG_COUNT);
		float nonReactionCenterSimilarity = getSimilarityTanimoto(o1, o2, REACTION_CENTER_LONG_COUNT, o1.length);
		return REACTION_CENTER_WEIGHT * reactionCenterSimilarity + PERIPHERY_WEIGHT * nonReactionCenterSimilarity;
	}

	public float getReactionCenterSimilarity(long[] o1, long[] o2) {
		if (o1 == null
		 || o2 == null
		 || o1.length == 0
		 || o2.length == 0)
			return 0.0f;

		return getSimilarityTanimoto(o1, o2, 0, REACTION_CENTER_LONG_COUNT);
		}

	public float getPeripherySimilarity(long[] o1, long[] o2) {
		if (o1 == null
		 || o2 == null
		 || o1.length == 0
		 || o2.length == 0)
			return 0.0f;

		return getSimilarityTanimoto(o1, o2, REACTION_CENTER_LONG_COUNT, o1.length);
		}

	public static float getSimilarityTanimoto(final long[] index1, final long[] index2, int i1, int i2) {
		int sharedKeys = 0;
		int allKeys = 0;
		for (int i=i1; i<i2; i++) {
			sharedKeys += Long.bitCount(index1[i] & index2[i]);
			allKeys += Long.bitCount(index1[i] | index2[i]);
		}
		return allKeys == 0 ? 1f : (float)sharedKeys/(float)allKeys;
	}

	public DescriptorHandler<long[], Reaction> getThreadSafeCopy() {
		return this;
	}
}
