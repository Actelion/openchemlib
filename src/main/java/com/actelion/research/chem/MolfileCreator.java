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
 * @author Thomas Sander
 */

package com.actelion.research.chem;

import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

public class MolfileCreator {
    private static final float TARGET_AVBL = 1.5f;

    private StringBuilder mBuilder;
    private DecimalFormat mDoubleFormat;

    /**
     * This creates a new molfile version 2 from the given molecule.
     * If the average bond length is smaller than 1.0 or larger than 3.0,
     * then all coordinates are scaled to achieve an average bond length of 1.5.
     * @param mol
     */
    public MolfileCreator(ExtendedMolecule mol) {
        this(mol, true);
        }

    /**
     * This creates a new molfile version 2 from the given molecule.
     * If allowScaling==true and the average bond length is smaller than 1.0 or larger than 3.0,
     * then all coordinates are scaled to achieve an average bond length of 1.5.
     * @param mol
     * @param allowScaling
     */
    public MolfileCreator(ExtendedMolecule mol, boolean allowScaling) {
        this(mol, allowScaling, new StringBuilder(32768));
        }

    /**
     * This creates a new molfile version 2 from the given molecule.
     * If allowScaling==true and the average bond length is smaller than 1.0 or larger than 3.0,
     * then all coordinates are scaled to achieve an average bond length of 1.5.
     * If a StringBuilder is given, then the molfile will be appended to that.
     * @param mol
     * @param allowScaling
     * @param builder null or StringBuilder to append to
     */
    public MolfileCreator(ExtendedMolecule mol, boolean allowScaling, StringBuilder builder) {
        this(mol, allowScaling, 0.0, builder);
    }

    /**
     * This creates a new molfile version 2 from the given molecule.
     * If allowScaling==true and the average bond length is smaller than 1.0 or larger than 3.0,
     * then all coordinates are scaled to achieve an average bond length of 1.5.
     * If scalingFactor is given, then the molecule is scaled accordingly independent of the average bond length.
     * If a StringBuilder is given, then the molfile will be appended to that.
     * @param mol
     * @param allowScaling if false, then no scaling is performed
     * @param scalingFactor if not 0.0 then the molecule is scaled by this factor
     * @param builder null or StringBuilder to append to
     */
    public MolfileCreator(ExtendedMolecule mol, boolean allowScaling, double scalingFactor, StringBuilder builder) {
		mDoubleFormat = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.ENGLISH)); //English local ('.' for the dot)
        final String nl = System.lineSeparator();

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

        mBuilder = (builder == null) ? new StringBuilder() : builder;

        String name = (mol.getName() != null) ? mol.getName() : "";
        mBuilder.append(name+nl);
        mBuilder.append("Actelion Java MolfileCreator 1.0"+nl+nl);

        appendThreeDigitInt(mol.getAllAtoms());
        appendThreeDigitInt(mol.getAllBonds());
        mBuilder.append("  0  0");
        appendThreeDigitInt((!isRacemic) ? 1 : 0);
        mBuilder.append("  0  0  0  0  0999 V2000"+nl);

        boolean hasCoordinates = (mol.getAllAtoms() == 1);
        for(int atom=1; atom<mol.getAllAtoms(); atom++) {
            if (mol.getAtomX(atom) != mol.getAtomX(0)
             || mol.getAtomY(atom) != mol.getAtomY(0)
             || mol.getAtomZ(atom) != mol.getAtomZ(0)) {
                hasCoordinates = true;
                break;
            }
        }

        double grafac = 1.0;

        if (hasCoordinates) {
        	if (scalingFactor != 0) {
                grafac = scalingFactor;
                }
            else if (allowScaling) {
                double avbl = mol.getAverageBondLength();
                if (avbl != 0.0f) {
                    if (avbl < 1.0f || avbl > 3.0f)
                        grafac = TARGET_AVBL / avbl;
                    }
                else { // make the minimum distance between any two atoms twice as long as TARGET_AVBL
                    double minDistance = Double.MAX_VALUE;
                    for (int atom1=1; atom1<mol.getAllAtoms(); atom1++) {
                        for (int atom2=0; atom2<atom1; atom2++) {
                            double dx = mol.getAtomX(atom2) - mol.getAtomX(atom1);
                            double dy = mol.getAtomY(atom2) - mol.getAtomY(atom1);
                            double dz = mol.getAtomZ(atom2) - mol.getAtomZ(atom1);
                            double distance = dx*dx + dy*dy + dz*dz;
                            if (minDistance > distance)
                                minDistance = distance;
                            }
                        }
                    grafac = 2.0f * TARGET_AVBL / minDistance;
                    }
                }
            }

        for (int atom=0; atom<mol.getAllAtoms(); atom++) {
            if (hasCoordinates) {
                appendTenDigitDouble(grafac * mol.getAtomX(atom));
                appendTenDigitDouble(grafac * -mol.getAtomY(atom));
                appendTenDigitDouble(grafac * -mol.getAtomZ(atom));
                }
            else {
                mBuilder.append("    0.0000    0.0000    0.0000");
                }
            if (mol.getAtomList(atom) != null)
                mBuilder.append(" L  ");
            else if ((mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0)
                mBuilder.append(" A  ");
            else if ((mol.getAtomicNo(atom) >= 129 && mol.getAtomicNo(atom) <= 144) || mol.getAtomicNo(atom) == 154)
                mBuilder.append(" R# ");
            else {
                String atomLabel = mol.getAtomLabel(atom);
                mBuilder.append(" "+atomLabel);
                if (atomLabel.length() == 1)
                    mBuilder.append("  ");
                else if (atomLabel.length() == 2)
                    mBuilder.append(" ");
                }

            mBuilder.append(" 0  0  0");	// massDif, charge, parity

            long hydrogenFlags = Molecule.cAtomQFHydrogen & mol.getAtomQueryFeatures(atom);
            if (hydrogenFlags == 0)
                mBuilder.append("  0");
            else if (hydrogenFlags == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen))
                mBuilder.append("  3");	// at least 2 hydrogens
            else if (hydrogenFlags == Molecule.cAtomQFNot0Hydrogen)
                mBuilder.append("  2");	// at least 1 hydrogens
            else if (hydrogenFlags == (Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen))
                mBuilder.append("  1");	// no hydrogens
            else if (hydrogenFlags == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen))
                mBuilder.append("  2");	// use at least 1 hydrogens as closest match for exactly one

            mBuilder.append(((mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFMatchStereo) != 0) ? "  1" : "  0");

            int valence = mol.getAtomAbnormalValence(atom);
            if (valence == -1)
                mBuilder.append("  0");
            else if (valence == 0)
                mBuilder.append(" 15");
            else
                appendThreeDigitInt(valence);

            mBuilder.append("  0  0  0");
            appendThreeDigitInt(mol.getAtomMapNo(atom));
            mBuilder.append("  0  0"+nl);
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
            case Molecule.cBondTypeMetalLigand:	order = 8; stereo = 0; break;
            default:							order = 1; stereo = 0; break; }

            if (isRacemic && (stereo == 1 || stereo == 6)) {
                int atom = mol.getBondAtom(0, bond);
                if (mol.getAtomESRType(atom) == Molecule.cESRTypeOr)
                    stereo = 0; // we interpret 'either' bonds as racemic and don't use it for OR atoms
                else if (mol.getAtomESRType(atom) == Molecule.cESRTypeAnd
                      && mol.getAtomESRGroup(atom) != maxESRGroup)
                    stereo = 4; // we use 'either' bonds for all racemic atoms that are not in the largest AND group
                }

            // if query features cannot be expressed exactly stay on the loosely defined side
            int bondType = mol.getBondQueryFeatures(bond) & Molecule.cBondQFBondTypes;
            if (bondType != 0) {
                if (bondType == Molecule.cBondTypeDelocalized)
                    order = 4;	// aromatic
                else if (bondType == (Molecule.cBondTypeSingle | Molecule.cBondTypeDouble))
                    order = 5;	// single or double
                else if (bondType == (Molecule.cBondTypeSingle | Molecule.cBondTypeDelocalized))
                    order = 6;	// single or aromatic
                else if (bondType == (Molecule.cBondTypeDouble | Molecule.cBondTypeDelocalized))
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
            mBuilder.append("  0");
            appendThreeDigitInt(topology);
            mBuilder.append("  0"+nl);
            }

        int no = 0;
        for (int atom=0; atom<mol.getAllAtoms(); atom++)
            if (mol.getAtomCharge(atom) != 0)
                no++;

        if (no != 0) {
            int count = 0;
            for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                if (mol.getAtomCharge(atom) != 0) {
                    if (count == 0) {
                        mBuilder.append("M  CHG");
                        appendThreeDigitInt(Math.min(8, no));
                        }
                    mBuilder.append(" ");
                    appendThreeDigitInt(atom + 1);
                    int charge = mol.getAtomCharge(atom);
                    if (charge < 0) {
                        mBuilder.append("  -");
                        charge = -charge;
                        }
                    else
                        mBuilder.append("   ");
                    mBuilder.append((char)('0' + charge));
                    no--;
                    if (++count == 8 || no == 0) {
                        count = 0;
                        mBuilder.append(nl);
                        }
                    }
                }
            }

        no = 0;
        for (int atom=0; atom<mol.getAllAtoms(); atom++)
            if (!mol.isNaturalAbundance(atom))
                no++;

        if (no != 0) {
            int count = 0;
            for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                if (!mol.isNaturalAbundance(atom)) {
                    if (count == 0) {
                        mBuilder.append("M  ISO");
                        appendThreeDigitInt(Math.min(8, no));
                        }
                    mBuilder.append(" ");
                    appendThreeDigitInt(atom + 1);
                    mBuilder.append(" ");
                    appendThreeDigitInt(mol.getAtomMass(atom));
                    no--;
                    if (++count == 8 || no == 0) {
                        count = 0;
                        mBuilder.append(nl);
                        }
                    }
                }
            }

        no = 0;
        for (int atom=0; atom<mol.getAllAtoms(); atom++)
            if (mol.getAtomRadical(atom) != 0)
                no++;

        if (no != 0) {
            int count = 0;
            for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                if (mol.getAtomRadical(atom) != 0) {
                    if (count == 0) {
                        mBuilder.append("M  RAD");
                        appendThreeDigitInt(Math.min(8, no));
                        }
                    mBuilder.append(" ");
                    appendThreeDigitInt(atom + 1);
                    switch (mol.getAtomRadical(atom)) {
                    case Molecule.cAtomRadicalStateS:
                        mBuilder.append("   1");
                        break;
                    case Molecule.cAtomRadicalStateD:
                        mBuilder.append("   2");
                        break;
                    case Molecule.cAtomRadicalStateT:
                        mBuilder.append("   3");
                        break;
                        }
                    no--;
                    if (++count == 8 || no == 0) {
                        count = 0;
                        mBuilder.append(nl);
                        }
                    }
                }
            }

        no = 0;
        for (int atom=0; atom<mol.getAllAtoms(); atom++)
            if ((mol.getAtomicNo(atom) >= 129 && mol.getAtomicNo(atom) <= 144) || mol.getAtomicNo(atom) == 154)
                no++;

        if (no != 0) {
            int count = 0;
            for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                int atomicNo = mol.getAtomicNo(atom);
                if ((atomicNo >= 129 && atomicNo <= 144) || atomicNo == 154) {
                    if (count == 0) {
                        mBuilder.append("M  RGP");
                        appendThreeDigitInt(Math.min(8, no));
                        }
                    mBuilder.append(" ");
                    appendThreeDigitInt(atom + 1);
                    mBuilder.append(" ");
                    appendThreeDigitInt(atomicNo == 154 ? 0 : atomicNo >= 142 ? atomicNo - 141 : atomicNo - 125);
                    no--;
                    if (++count == 8 || no == 0) {
                        count = 0;
                        mBuilder.append(nl);
                        }
                    }
                }
            }

        if (mol.isFragment()) {
            no = 0;
            for (int atom=0; atom<mol.getAllAtoms(); atom++)
                if ((mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFRingState) != 0)
                    no++;

            if (no != 0) {
                int count = 0;
                for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                    long ringFeatures = mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFRingState;
                    if (ringFeatures != 0) {
                        if (count == 0) {
                            mBuilder.append("M  RBC");
                            appendThreeDigitInt(Math.min(8, no));
                            }
                        mBuilder.append(" ");
                        appendThreeDigitInt(atom + 1);
                        if (ringFeatures == (Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds))
                                mBuilder.append("  -1");
                        else if (ringFeatures == Molecule.cAtomQFNotChain)
                                mBuilder.append("   1");	// any ring atom; there is no MDL equivalent
                        else if (ringFeatures == (Molecule.cAtomQFNotChain | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds))
                                mBuilder.append("   2");
                        else if (ringFeatures == (Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot4RingBonds))
                                mBuilder.append("   3");
                        else if (ringFeatures == (Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds))
                                mBuilder.append("   4");
                        no--;
                        if (++count == 8 || no == 0) {
                            count = 0;
                            mBuilder.append(nl);
                            }
                        }
                    }
                }

            for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                int[] atomList = mol.getAtomList(atom);
                if (atomList != null) {
                    mBuilder.append("M  ALS ");
                    appendThreeDigitInt(atom+1);
                    appendThreeDigitInt(atomList.length);
                    mBuilder.append(((mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0) ? " T " : " F ");
                    for (int i=0; i<atomList.length; i++) {
                        String label = Molecule.cAtomLabel[atomList[i]];
                        switch (label.length()) {
                        case 1:
                            mBuilder.append(label+"   ");
                            break;
                        case 2:
                            mBuilder.append(label+"  ");
                            break;
                        case 3:
                            mBuilder.append(label+" ");
                            break;
                        default:
                            mBuilder.append("   ?");
                            break;
                            }
                        }
                    mBuilder.append(nl);
                    }
                }

            no = 0;
            for (int atom=0; atom<mol.getAllAtoms(); atom++)
                if ((mol.getAtomQueryFeatures(atom) & (Molecule.cAtomQFMoreNeighbours | Molecule.cAtomQFNoMoreNeighbours)) != 0)
                    no++;

            if (no != 0) {
                int count = 0;
                for (int atom=0; atom<mol.getAllAtoms(); atom++) {
                    long substitution = mol.getAtomQueryFeatures(atom) & (Molecule.cAtomQFMoreNeighbours | Molecule.cAtomQFNoMoreNeighbours);
                    if (substitution != 0) {
                        if (count == 0) {
                            mBuilder.append("M  SUB");
                            appendThreeDigitInt(Math.min(8, no));
                            }
                        mBuilder.append(" ");
                        appendThreeDigitInt(atom + 1);
                        if ((substitution & Molecule.cAtomQFMoreNeighbours) != 0)
                            mBuilder.append("   "+(mol.getAllConnAtoms(atom)+1));
                        else
                            mBuilder.append("  -2");
                        no--;
                        if (++count == 8 || no == 0) {
                            count = 0;
                            mBuilder.append(nl);
                            }
                        }
                    }
                }
            }

        mBuilder.append("M  END"+nl);
        }


    /**
     * If a pre-filled StringBuilder was passed to the constructor, then this returns
     * the original content with the appended molfile.
     * @return
     */
    public String getMolfile() {
        return mBuilder.toString();
        }


    public void writeMolfile(Writer theWriter) throws IOException {
        theWriter.write(mBuilder.toString());
        }


    private void appendThreeDigitInt(int data) {
        if (data < 0 || data > 999) {
            mBuilder.append("  ?");
            return;
            }

        boolean digitFound = false;
        for (int i=0; i<3; i++) {
            int theChar = data / 100;
            if (theChar==0) {
                if (i==2 || digitFound)
                    mBuilder.append((char)'0');
                else
                    mBuilder.append((char)' ');
                }
            else {
                mBuilder.append((char)('0' + theChar));
                digitFound = true;
                }
            data = 10 * (data % 100);
            }
        }

    private void appendTenDigitDouble(double theDouble) {
    	String val = mDoubleFormat.format(theDouble);
    	for(int i=val.length(); i<10; i++) mBuilder.append(' ');
    	mBuilder.append(val);
        }
	}

