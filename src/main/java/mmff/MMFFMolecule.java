/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of ActelionMMFF94.
 * 
 * ActelionMMFF94 is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * ActelionMMFF94 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with ActelionMMFF94.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Paolo Tosco,Daniel Bergmann
 */

package mmff;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.RingCollection;

/**
 * MMFF molecule is a wrapper class for the ExtendedMolecule. It holds some
 * additional data such as a cache of the atom types, whether the molecule is
 * valid for MMFF and the ring mmff aromaticity property.
 */
public final class MMFFMolecule extends StereoMolecule  {
    private RingBoolean[] mRingArom;
    private int[] mAtomTypes;
    private int[] mHydrogenMap;

    public MMFFMolecule(StereoMolecule mol) throws BadAtomTypeException,
                                                 BadRingAromException
    {
        super(mol);
        mHydrogenMap = mol.getHandleHydrogenMap();
        mol.ensureHelperArrays(StereoMolecule.cHelperRings);

        RingCollection rings = mol.getRingSet();
        mRingArom = new RingBoolean[rings.getSize()];
        for (int i=0; i<mRingArom.length; i++)
            mRingArom[i] = RingBoolean.NOT_SET;

        boolean allset = false, changed = true;
        while (!allset && changed) {
            allset = true;
            changed = false;
            for (int r=0; r<rings.getSize(); r++) {
                if (mRingArom[r] == RingBoolean.NOT_SET) {
                    mRingArom[r] = mmff.type.Atom.ringIsMMFFAromatic(this, r);

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
        mAtomTypes = new int[mol.getAllAtoms()];
        for (int i=0; i<mAtomTypes.length; i++) {
        	mAtomTypes[i] = -1;
        	mAtomTypes[i] = mmff.type.Atom.getType(this, i);

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
        return mRingArom[r] == RingBoolean.TRUE ? true : false;
    }

    /**
     * Returns true if the given ring has had its MMFF aromaticity flag set.
     *  @param r The ring index in the molecule.
     *  @return True if the ring has had its flag set, false otherwise.
     */
    public boolean isSetRingMMFFAromaticity(int r) {
        return mRingArom[r] == RingBoolean.NOT_SET ? false : true;
    }


}
