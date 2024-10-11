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

package com.actelion.research.chem.forcefield.mmff;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.RingCollection;

import java.util.Arrays;

/**
 * MMFF molecule is a wrapper class for the ExtendedMolecule. It holds some
 * additional data such as a cache of the atom types, whether the molecule is
 * valid for MMFF and the ring mmff aromaticity property.
 */
public final class MMFFMolecule extends StereoMolecule  {
    private final RingBoolean[] mRingArom;
    private final int[] mAtomTypes;
    private final int[] mHydrogenMap;

    public MMFFMolecule(StereoMolecule mol) throws BadAtomTypeException,BadRingAromException {
    	super(mol);
        mHydrogenMap = getHandleHydrogenMap();
        RingCollection rings = getRingSet();
        mRingArom = new RingBoolean[rings.getSize()];
		Arrays.fill(mRingArom, RingBoolean.NOT_SET);

        boolean allset = false, changed = true;
        while (!allset && changed) {
            allset = true;
            changed = false;
            for (int r=0; r<rings.getSize(); r++) {
                if (mRingArom[r] == RingBoolean.NOT_SET) {
                    mRingArom[r] = com.actelion.research.chem.forcefield.mmff.type.Atom.ringIsMMFFAromatic(this, r);

                    if (mRingArom[r] != RingBoolean.NOT_SET)
                        changed = true;
                }

                if (mRingArom[r] == RingBoolean.NOT_SET)
                    allset = false;
            }
        }

        if (!allset)
            throw new BadRingAromException();

        // Assign the atom types to the atom type cache.
        mAtomTypes = new int[getAllAtoms()];
        for (int i=0; i<mAtomTypes.length; i++) {

        	mAtomTypes[i] = -1;
        	mAtomTypes[i] = com.actelion.research.chem.forcefield.mmff.type.Atom.getType(this, i);

            if (mAtomTypes[i] == 0)
                throw new BadAtomTypeException("Couldn't assign an atom type "
                        +"to atom "+i+" ("+getAtomLabel(i)+")");
        }
    }

    /**
     * Get the MMFF atom type of an atom. This returns the cached value.
     *  @param a The atom index in the molecule.
     *  @return The MMFF atom type.
     */
    public int getAtomType(int a) {
        return mAtomTypes[a];
    }
    
    public int[] getHydrogenMap() {
    	return mHydrogenMap;
    }

    /**
     * Determine if a ring is aromatic according to MMFF criteria. Only
     * designed to work with rings of size 5 and 6. Returns the cached value.
     *  @param r The ring index in the molecule.
     *  @return True if the ring is aromatic, false otherwise.
     */
    public boolean ringIsMMFFAromatic(int r) {
        return mRingArom[r] == RingBoolean.TRUE;
    }

    /**
     * Returns true if the given ring has had its MMFF aromaticity flag set.
     *  @param r The ring index in the molecule.
     *  @return True if the ring has had its flag set, false otherwise.
     */
    public boolean isSetRingMMFFAromaticity(int r) {
        return mRingArom[r] != RingBoolean.NOT_SET;
    }


}
