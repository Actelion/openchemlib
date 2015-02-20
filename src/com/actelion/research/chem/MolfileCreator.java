/*
* Copyright (c) 1997 - 2015
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

import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

public class MolfileCreator {
    private static final float TARGET_AVBL = 1.5f;
	
    private StringBuilder mMolfile;
    
    public MolfileCreator(ExtendedMolecule mol) {
        this(mol, true);
        }
    
    /**
     * This constructor is needed for the reactionfile writer, since we cannot allow scaling on individual molecule within a reaction
     * @param mol
     * @param scale
     */
    public MolfileCreator(ExtendedMolecule mol, boolean scale) {
        mol.ensureHelperArrays(Molecule.cHelperParities);

        boolean isRacemic = true;
        for (int atom=0; atom<mol.getAtoms(); atom++) {
            if (mol.getAtomParity(atom) != Molecule.cAtomParityNone
             && mol.getAtomParity(atom) != Molecule.cAtomParityUnknown
             && mol.getAtomESRType(atom) != Molecule.cESRTypeAnd) {
                isRacemic = false;
                break;
                }
            }
        int maxESRGroup = -1;
        if (isRacemic) {
            int[] esrGroupCount = new int[Molecule.cESRMaxGroups];
            int maxGroupCount = 0;
            for (int atom=0; atom<mol.getAtoms(); atom++) {
                if (mol.getAtomParity(atom) != Molecule.cAtomParityNone
                 && mol.getAtomParity(atom) != Molecule.cAtomParityUnknown
                 && mol.getAtomESRType(atom) == Molecule.cESRTypeAnd) {
                    int group = mol.getAtomESRGroup(atom);
                    esrGroupCount[group]++;
                    if (maxGroupCount < esrGroupCount[group]) {
                        maxGroupCount = esrGroupCount[group];
                        maxESRGroup = group;
                        }
                    break;
                    }
                }
            }
        
        mMolfile = new StringBuilder(32768);
        String name = (mol.getName() != null) ? mol.getName() : "";
        mMolfile.append(name+"\n");
        mMolfile.append("Actelion Java MolfileCreator 1.0\n\n");

        appendThreeDigitInt(mol.getAllAtoms());
        appendThreeDigitInt(mol.getAllBonds());
        mMolfile.append("  0  0");
        appendThreeDigitInt((!isRacemic) ? 1 : 0);
        mMolfile.append("  0  0  0  0  0999 V2000\n");

        boolean hasCoordinates = (mol.getAllAtoms() == 1);
        for(int atom=1; atom<mol.getAllAtoms(); atom++) {
            if (mol.getAtomX(atom) != mol.getAtomX(0)
             || mol.getAtomY(atom) != mol.getAtomY(0)
             || mol.getAtomZ(atom) != mol.getAtomZ(0)) {
                hasCoordinates = true;
                break;
            }
        }

        float grafac = 1.0f;

        if (hasCoordinates && scale) {
            float avbl = mol.getAverageBondLength();
            if (avbl != 0.0f) {
            	if (avbl < 1.0f || avbl > 3.0f)
            		grafac = TARGET_AVBL / avbl;
            	}
            else { // make the minimum distance between any two atoms twice as long as TARGET_AVBL
                float minDistance = Float.MAX_VALUE;
                for (int atom1=1; atom1<mol.getAllAtoms(); atom1++) {
                    for (int atom2=0; atom2<atom1; atom2++) {
                    	float dx = mol.getAtomX(atom2) - mol.getAtomX(atom1);
                    	float dy = mol.getAtomY(atom2) - mol.getAtomY(atom1);
                    	float dz = mol.getAtomZ(atom2) - mol.getAtomZ(atom1);
                    	float distance = dx*dx + dy*dy + dz*dz;
                        if (minDistance > distance)
                            minDistance = distance;
                        }
                    }
                grafac = 2.0f * TARGET_AVBL / minDistance;
	            }
	        }

        for (int atom=0; atom<mol.getAllAtoms(); atom++) {
            if (hasCoordinates) {
                appendTenDigitDouble(grafac * mol.getAtomX(atom));
                appendTenDigitDouble(grafac * -mol.getAtomY(atom));
                appendTenDigitDouble(grafac * -mol.getAtomZ(atom));
                }
            else {
                mMolfile.append("    0.0000    0.0000    0.0000");
                }
            if (mol.getAtomList(atom) != null)
                mMolfile.append(" L  ");
            else if ((mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0)
                mMolfile.append(" A  ");
            else {
                String atomLabel = mol.getAtomLabel(atom);
                mMolfile.append(" "+atomLabel);
                if (atomLabel.length() == 1)
                    mMolfile.append("  ");
                else if (atomLabel.length() == 2)
                    mMolfile.append(" ");
                }

            mMolfile.append(" 0  0  0");	// massDif, charge, parity

            int hydrogenFlags = Molecule.cAtomQFHydrogen & mol.getAtomQueryFeatures(atom);
            if (hydrogenFlags == 0)
                mMolfile.append("  0");
            else if (hydrogenFlags == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen))
                mMolfile.append("  3");	// at least 2 hydrogens
            else if (hydrogenFlags == Molecule.cAtomQFNot0Hydrogen)
                mMolfile.append("  2");	// at least 1 hydrogens
            else if (hydrogenFlags == (Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen))
                mMolfile.append("  1");	// no hydrogens
            else if (hydrogenFlags == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen))
                mMolfile.append("  2");	// use at least 1 hydrogens as closest match for exactly one

            mMolfile.append(((mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFMatchStereo) != 0) ? "  1" : "  0");

            int valence = mol.getAtomAbnormalValence(atom);
            if (valence == -1)
                mMolfile.append("  0");
            else if (valence == 0)
                mMolfile.append(" 15");
            else
                appendThreeDigitInt(valence);

            mMolfile.append("  0  0  0");
            appendThreeDigitInt(mol.getAtomMapNo(atom));
            mMolfile.append("  0  0\n");
            }

        for (int bond=0; bond<mol.getAllBonds(); bond++) {
            int order,stereo;
            switch (mol.getBondType(bond)) {
            case Molecule.cBondTypeSingle:		order = 1; stereo = 0; break;
            case Molecule.cBondTypeDouble:		order = 2; stereo = 0; break;
            case Molecule.cBondTypeTriple:		order = 3; stereo = 0; break;
            case Molecule.cBondTypeDown:		order = 1; stereo = 6; break;
            case Molecule.cBondTypeUp:			order = 1; stereo = 1; break;
            case Molecule.cBondTypeCross:		order = 2; stereo = 3; break;
            case Molecule.cBondTypeDelocalized:	order = 4; stereo = 0; break;
            default:							order = 1; stereo = 0; break; }

            if (isRacemic && (stereo == 1 || stereo == 6)) {
                if (mol.getAtomESRGroup(mol.getBondAtom(0, bond)) != maxESRGroup)
                    stereo = 0;
                }
                    // if query features cannot be expressed exactly stay on the loosely defined side
            int bondType = mol.getBondQueryFeatures(bond) & Molecule.cBondQFBondTypes;
            if (bondType != 0) {
                if (bondType == Molecule.cBondQFDelocalized)
                    order = 4;	// aromatic
                else if (bondType == (Molecule.cBondQFSingle | Molecule.cBondQFDouble))
                    order = 5;	// single or double
                else if (bondType == (Molecule.cBondQFSingle | Molecule.cBondQFDelocalized))
                    order = 6;	// single or aromatic
                else if (bondType == (Molecule.cBondQFDouble | Molecule.cBondQFDelocalized))
                    order = 7;	// single or double
                else
                    order = 8;	// any
                }

            int ringState = mol.getBondQueryFeatures(bond) & Molecule.cBondQFRingState;
            int topology = (ringState == 0) ? 0 : (ringState == Molecule.cBondQFRing) ? 1 : 2;

            appendThreeDigitInt(1 + mol.getBondAtom(0, bond));
            appendThreeDigitInt(1 + mol.getBondAtom(1, bond));
            appendThreeDigitInt(order);
            appendThreeDigitInt(stereo);
            mMolfile.append("  0");
            appendThreeDigitInt(topology);
            mMolfile.append("  0\n");
            }

        int no = 0;
        for (int atom=0; atom<mol.getAllAtoms(); atom++)
            if (mol.getAtomCharge(atom) != 0)
                no++;

        if (no != 0) {
            mMolfile.append("M  CHG");
            appendThreeDigitInt(no);
            for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                if (mol.getAtomCharge(atom) != 0) {
                    mMolfile.append(" ");
                    appendThreeDigitInt(atom + 1);
                    int charge = mol.getAtomCharge(atom);
                    if (charge < 0) {
                        mMolfile.append("  -");
                        charge = -charge;
                        }
                    else
                        mMolfile.append("   ");
                    mMolfile.append((char)('0' + charge));
                    }
                }
            mMolfile.append("\n");
            }

        no = 0;
        for (int atom=0; atom<mol.getAllAtoms(); atom++)
            if (!mol.isNaturalAbundance(atom))
                no++;
        if (no != 0) {
            mMolfile.append("M  ISO");
            appendThreeDigitInt(no);

            for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                if (!mol.isNaturalAbundance(atom)) {
                    mMolfile.append(" ");
                    appendThreeDigitInt(atom + 1);
                    mMolfile.append(" ");
                    appendThreeDigitInt(mol.getAtomMass(atom));
                    }
                }

            mMolfile.append("\n");
            }

        no = 0;
        for (int atom=0; atom<mol.getAllAtoms(); atom++)
            if (mol.getAtomRadical(atom) != 0)
                no++;
        if (no != 0) {
            mMolfile.append("M  RAD");
            appendThreeDigitInt(no);

            for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                if (mol.getAtomRadical(atom) != 0) {
                    mMolfile.append(" ");
                    appendThreeDigitInt(atom + 1);
                    switch (mol.getAtomRadical(atom)) {
                    case Molecule.cAtomRadicalStateS:
                        mMolfile.append("   1");
                        break;
                    case Molecule.cAtomRadicalStateD:
                        mMolfile.append("   2");
                        break;
                    case Molecule.cAtomRadicalStateT:
                        mMolfile.append("   3");
                        break;
                        }
                    }
                }

            mMolfile.append("\n");
            }

        if (mol.isFragment()) {
            no = 0;
            for (int atom=0; atom<mol.getAllAtoms(); atom++)
                if ((mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFRingState) != 0)
                    no++;
            if (no != 0) {
                mMolfile.append("M  RBD");
                appendThreeDigitInt(no);

                for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                    int ringFeatures = mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFRingState;
                    if (ringFeatures != 0) {
                        mMolfile.append(" ");
                        appendThreeDigitInt(atom + 1);
                        switch (ringFeatures) {
                            case Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds:
                                mMolfile.append("  -1");
                                break;
                            case Molecule.cAtomQFNotChain:
                                mMolfile.append("   1");	// any ring atom; there is no MDL equivalent
                                break;
                            case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds:
                                mMolfile.append("   2");
                                break;
                            case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot4RingBonds:
                                mMolfile.append("   3");
                                break;
                            case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds:
                                mMolfile.append("   4");
                                break;
                            }
                        }
                    }
                mMolfile.append("\n");
                }

            for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                int[] atomList = mol.getAtomList(atom);
                if (atomList != null) {
                    mMolfile.append("M  ALS ");
                    appendThreeDigitInt(atom+1);
                    appendThreeDigitInt(atomList.length);
                    mMolfile.append(((mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0) ? " T " : " F ");
                    for (int i=0; i<atomList.length; i++) {
                        String label = Molecule.cAtomLabel[atomList[i]];
                        switch (label.length()) {
                        case 1:
                            mMolfile.append(label+"   ");
                            break;
                        case 2:
                            mMolfile.append(label+"  ");
                            break;
                        case 3:
                            mMolfile.append(label+" ");
                            break;
                        default:
                            mMolfile.append("   ?");
                            break;
                            }
                        }
                    mMolfile.append("\n");
                    }
                }

            no = 0;
            for (int atom=0; atom<mol.getAllAtoms(); atom++)
                if ((mol.getAtomQueryFeatures(atom) & (Molecule.cAtomQFMoreNeighbours | Molecule.cAtomQFNoMoreNeighbours)) != 0)
                    no++;
            if (no != 0) {
                mMolfile.append("M  SUB");
                appendThreeDigitInt(no);

                for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                    int substitution = mol.getAtomQueryFeatures(atom) & (Molecule.cAtomQFMoreNeighbours | Molecule.cAtomQFNoMoreNeighbours);
                    if (substitution != 0) {
                        mMolfile.append(" ");
                        appendThreeDigitInt(atom + 1);
                        if ((substitution & Molecule.cAtomQFMoreNeighbours) != 0)
                            mMolfile.append("   "+(mol.getAllConnAtoms(atom)+1));
                        else
                            mMolfile.append("  -2");
                        }
                    }
                mMolfile.append("\n");
                }
            }

        mMolfile.append("M  END\n");
        }


    public String getMolfile() {
        return mMolfile.toString();
        }


    public void writeMolfile(Writer theWriter) throws IOException {
        theWriter.write(mMolfile.toString());
        }


    private void appendThreeDigitInt(int data) {
        if (data < 0 || data > 999) {
            mMolfile.append("  ?");
            return;
            }

        boolean digitFound = false;
        for (int i=0; i<3; i++) {
            int theChar = data / 100;
            if (theChar==0) {
                if (i==2 || digitFound)
                    mMolfile.append((char)'0');
                else
                    mMolfile.append((char)' ');
                }
            else {
                mMolfile.append((char)('0' + theChar));
                digitFound = true;
                }
            data = 10 * (data % 100);
            }
        }

    DecimalFormat df = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.ENGLISH)); //English local ('.' for the dot)
    private void appendTenDigitDouble(double theDouble) {
    	String val = df.format(theDouble);  
    	for(int i=val.length(); i<10; i++) mMolfile.append(' ');   
    	mMolfile.append(val);  
        }
	}

