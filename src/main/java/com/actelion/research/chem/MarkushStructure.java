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

package com.actelion.research.chem;

import com.actelion.research.chem.coords.CoordinateInventor;

import java.util.ArrayList;

public class MarkushStructure implements java.io.Serializable {
	static final long serialVersionUID = 0x20080515;

	private static final int MAX_RGROUP_DEPTH = 5;

	private	ArrayList<StereoMolecule>   mCoreList;
	private	ArrayList<StereoMolecule>   mRGroupList;
	private StereoMolecule[][]          mSubstituent;
	private int                         mCoreIndex;
	private int[]                       mSubstituentCount,mSubstituentIndex;
	private int[][]                     mCoreRGroup;
	private int[][][]                   mSubstituentRGroup;

	public MarkushStructure() {
	    mCoreList = new ArrayList<StereoMolecule>();
	    mRGroupList = new ArrayList<StereoMolecule>();
		}

	public MarkushStructure(StereoMolecule[] fragment, int coreCount) {
		this();
		if (fragment != null) {
			for (int i=0; i<coreCount; i++)
			    mCoreList.add(fragment[i]);
			for (int i=coreCount; i<fragment.length; i++)
			    mRGroupList.add(fragment[i]);
			}
		}

	public StereoMolecule getCoreStructure(int no) {
		return mCoreList.get(no);
		}


	public int getCoreCount() {
		return mCoreList.size();
		}


	public StereoMolecule getRGroup(int no) {
		return mRGroupList.get(no);
		}


	public int getRGroupCount() {
		return mRGroupList.size();
		}


	public void addCore(StereoMolecule core) {
		mCoreList.add(core);
		}


	/**
	 * This adds a substituent list representing one R-group as multiple
	 * fragments within one StereoMolecule. Every fragment is supposed to
	 * contain one atom with atomicNo=0 that is considered the attachment point.
	 * Fragments may contain atoms representing other R-groups (R1-R16).
	 * The first call of this method adds R1, the second R2, etc.
	 * @param substituents multiple fragments representing R-group substituents
	 */
	public void addRGroup(StereoMolecule substituents) {
	    mRGroupList.add(substituents);
		}

	/**
	 * Check the validity of a defined Markush structure.
	 * This method must be called before the enumeration.
	 * @return a warning if defined R-groups are unreferenced
     * @throws Exception if an error prevents proper enumeration
	 */
    public String validate() throws Exception {
        if (mCoreList.size() == 0)
            throw new Exception("You didn't define a core structure.");

        int rGroupFlags = 0;
        int maxRGroup = -1;
        int rGroupCount = 0;
        mCoreRGroup = new int[mCoreList.size()][];
		for (int i=0; i<mCoreList.size()+mRGroupList.size(); i++) {
		    StereoMolecule mol = (i<mCoreList.size()) ? mCoreList.get(i) : mRGroupList.get(i-mCoreList.size());
		    mol.ensureHelperArrays(Molecule.cHelperNeighbours);
            ArrayList<Integer> rGroupList = new ArrayList<Integer>();
			for (int atom=0; atom<mol.getAllAtoms(); atom++) {
				int rGroupIndex = getRGroupIndex(mol.getAtomicNo(atom));
				if (rGroupIndex != -1) {
				    if (mol.getConnAtoms(atom) != 1)
                        throw new Exception("R"+(rGroupIndex+1)+" has multiple neighbours.");
				    maxRGroup = Math.max(maxRGroup, rGroupIndex);
                    int flag = (1 << rGroupIndex);
                    if ((rGroupFlags & flag) == 0) {
                        rGroupFlags |= flag;
                        rGroupCount++;
                        }
                    if (i<mCoreList.size())
                        rGroupList.add(new Integer(rGroupIndex));
				    }
				}
            if (i<mCoreList.size()) {
    	        mCoreRGroup[i] = new int[rGroupList.size()];
    	        for (int j=0; j<rGroupList.size(); j++)
    	            mCoreRGroup[i][j] = rGroupList.get(j).intValue();
                }
			}
        if (maxRGroup >= mRGroupList.size())
            throw new Exception("Not all used R-groups are defined.");

        String warning = null;
        if (rGroupCount < mRGroupList.size())
            warning = "One or more defined R-group(s) are unreferenced.";

        mSubstituent = new StereoMolecule[mRGroupList.size()][];
        mSubstituentRGroup = new int[mRGroupList.size()][][];
        for (int i=0; i<mRGroupList.size(); i++) {
            mSubstituent[i] = mRGroupList.get(i).getFragments();
            mSubstituentRGroup[i] = new int[mSubstituent[i].length][];
            for (int j=0; j<mSubstituent[i].length; j++) {
                mSubstituent[i][j].ensureHelperArrays(Molecule.cHelperNeighbours);
                boolean attachmentPointFound = false;
                ArrayList<Integer> rGroupList = new ArrayList<Integer>();
                for (int atom=0; atom<mSubstituent[i][j].getAllAtoms(); atom++) {
                    if (mSubstituent[i][j].getAtomicNo(atom) == 0) {
                        if (attachmentPointFound)
                            throw new Exception("Variation "+(j+1)+" of R"+(i+1)+" has multiple attachment points.");
                        if (mSubstituent[i][j].getConnAtoms(atom) != 1)
                            throw new Exception("The attachment point of variation "+(j+1)+" of R"+(i+1)+" has multiple neighbours.");
                        attachmentPointFound = true;
                        }
                    else {
                        int rGroupIndex = getRGroupIndex(mSubstituent[i][j].getAtomicNo(atom));
                        if (rGroupIndex != -1)
                            rGroupList.add(new Integer(rGroupIndex));
                        }
                    }
                if (!attachmentPointFound)
                    throw new Exception("Variation "+(j+1)+" of R"+(i+1)+" has no attachment point ('?'-atom).");
                mSubstituentRGroup[i][j] = new int[rGroupList.size()];
                for (int k=0; k<rGroupList.size(); k++)
                    mSubstituentRGroup[i][j][k] = rGroupList.get(k).intValue();
                }
            }

        int rGroupReferenced = 0;
        for (int i=0; i<mCoreRGroup.length; i++)
            for (int j=0; j<mCoreRGroup[i].length; j++)
                rGroupReferenced |= (1 << mCoreRGroup[i][j]);

        int rGroupReferencedPreviousLevel = rGroupReferenced;
        for (int level=1; level<=MAX_RGROUP_DEPTH; level++) {
            int rGroupReferencedCurrentLevel = 0;
            for (int rGroup=0; rGroup<=maxRGroup; rGroup++) {
                if ((rGroupReferencedPreviousLevel & (1 << rGroup)) != 0) {
                    for (int i=0; i<mSubstituentRGroup[rGroup].length; i++) {
                        for (int j=0; j<mSubstituentRGroup[rGroup][i].length; j++) {
                            rGroupReferencedCurrentLevel |= (1 << mSubstituentRGroup[rGroup][i][j]);
                            rGroupReferenced |= (1 << mSubstituentRGroup[rGroup][i][j]);
                            }
                        }
                    }
                }
            if (rGroupReferencedCurrentLevel == 0)
                break;
            if (level == MAX_RGROUP_DEPTH)
                throw new Exception("R-Groups are recursively definined or maximum depth of indirection ("+MAX_RGROUP_DEPTH+") exceeded");

            rGroupReferencedPreviousLevel = rGroupReferencedCurrentLevel;
            }

        mSubstituentCount = new int[mRGroupList.size()];
        for (int rGroup=0; rGroup<mRGroupList.size(); rGroup++)
            mSubstituentCount[rGroup] = ((rGroupReferenced & (1 << rGroup)) == 0) ?
                    0 : mSubstituent[rGroup].length;
        mSubstituentIndex = new int[mRGroupList.size()];

        mCoreIndex = 0;
		return warning;
		}

    /**
     * After calling validate() this method may be called until null is returned
     * to construct one by one a new representation of the Markush structure.
     */
    public StereoMolecule getNextEnumeration() {
        if (mCoreIndex == mCoreList.size())
            return null;

        StereoMolecule mol = createCurrentEnumeration();

        boolean incrementCoreIndex = true;
        for (int rGroup=0; rGroup<mRGroupList.size(); rGroup++) {
            if (mSubstituentCount[rGroup] > 1) {
                if (mSubstituentIndex[rGroup]+1 < mSubstituentCount[rGroup]) {
                    mSubstituentIndex[rGroup]++;
                    incrementCoreIndex = false;
                    break;
                    }
                else {
                    mSubstituentIndex[rGroup] = 0;
                    }
                }
            }
        if (incrementCoreIndex)
            mCoreIndex++;

        return mol;
        }

    private int getRGroupIndex(int atomicNo) {
        if (atomicNo < 129 || atomicNo > 144)
            return -1;
        return (atomicNo >= 142) ? atomicNo - 142 : atomicNo - 126;
        }

    private StereoMolecule createCurrentEnumeration() {
        StereoMolecule mol = new StereoMolecule();

        mol.addMolecule(mCoreList.get(mCoreIndex));
        for (int atom=0; atom<mol.getAllAtoms(); atom++) {
            int rGroup = getRGroupIndex(mol.getAtomicNo(atom));
            if (rGroup != -1)
                mol.addSubstituent(mSubstituent[rGroup][mSubstituentIndex[rGroup]],
                                   getAttachmentAtom(mol, atom));
            }
        for (int atom=0; atom<mol.getAllAtoms(); atom++)
            mol.setAtomSelection(atom, getRGroupIndex(mol.getAtomicNo(atom)) != -1);
        mol.deleteSelectedAtoms();

        int coreAtoms = mCoreList.get(mCoreIndex).getAllAtoms() - mCoreRGroup[mCoreIndex].length;
        for (int atom=0; atom<coreAtoms; atom++)
            mol.setAtomMarker(atom, true);
        new CoordinateInventor(CoordinateInventor.MODE_REMOVE_HYDROGEN | CoordinateInventor.MODE_PREFER_MARKED_ATOM_COORDS).invent(mol);
        
        return mol;
        }

    private int getAttachmentAtom(Molecule mol, int atom) {
        for (int bond=0; bond<mol.getAllBonds(); bond++) {
            if (mol.getBondAtom(0, bond) == atom)
                return mol.getBondAtom(1, bond);
            if (mol.getBondAtom(1, bond) == atom)
                return mol.getBondAtom(0, bond);
            }
        return -1;
        }
    }
