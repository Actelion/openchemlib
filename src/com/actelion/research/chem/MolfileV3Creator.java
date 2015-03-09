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


package com.actelion.research.chem;

import java.io.IOException;
import java.io.Writer;

/**
 * This class generates an MDL molfile version 3.0 from a StereoMolecule
 * as described by MDL in 'CTFile Formats June 2005'.
 * Since the MDL enhanced stereo recognition concept doesn't include
 * support for axial chirality as bond property, we added object type 'BONDS'
 * to the internal collection types STEABS,STERAC and STEREL, in order to
 * properly encode ESR assignments of axial stereo bonds, e.g. BINAP kind of
 * stereo bonds.
 * @author sandert
 *
 */
public class MolfileV3Creator
{
    private StringBuilder mMolfile;
    private static final float TARGET_AVBL = 1.5f;
    private static final float PRECISION_FACTOR = 10000f;

    private float mScalingFactor = 1.0f;

    public MolfileV3Creator(ExtendedMolecule mol) {
        this(mol, true);
    	}

    // this constructor is needed for the reactionfile writer, since we cannot allow scaling on individual molecule within a reaction
    public MolfileV3Creator(ExtendedMolecule mol, boolean scale) {
        mol.ensureHelperArrays(Molecule.cHelperParities);

        mMolfile = new StringBuilder(32768);
        String name = (mol.getName() != null) ? mol.getName() : "";
        mMolfile.append(name + "\n");
        mMolfile.append("Actelion Java MolfileCreator 2.0\n\n");
        mMolfile.append("  0  0  0  0  0  0              0 V3000\n");

        boolean hasCoordinates = (mol.getAllAtoms() == 1);
        for(int atom=1; atom<mol.getAllAtoms(); atom++) {
            if (mol.getAtomX(atom) != mol.getAtomX(0)
             || mol.getAtomY(atom) != mol.getAtomY(0)
             || mol.getAtomZ(atom) != mol.getAtomZ(0)) {
                hasCoordinates = true;
                break;
            }
        }

        mScalingFactor = 1.0f;

        if (hasCoordinates && scale) {
        	// Calculate a reasonable molecule size for ISIS-Draw default settings.
            float avbl = mol.getAverageBondLength();
            if (avbl != 0.0f) {
            	// 0.84 seems to be the average bond distance in ISIS Draw 2.5 with the default setting of 0.7 cm standard bond length.
                // grafac = 0.84f / mol.getAverageBondLength();

            	if (avbl < 1.0f || avbl > 3.0f)
            		mScalingFactor = TARGET_AVBL / avbl;
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
                mScalingFactor = 2.0f * TARGET_AVBL / minDistance;
	            }
	        }

        writeBody(mol, hasCoordinates);
        mMolfile.append("M  END\n");
    	}

    public static String writeCTAB(ExtendedMolecule mol, boolean hasCoordinates) {
        MolfileV3Creator mf = new MolfileV3Creator();
        mol.ensureHelperArrays(Molecule.cHelperParities);
        mf.writeBody(mol, hasCoordinates);
        return mf.getMolfile();
    	}

    private MolfileV3Creator() {
        mMolfile = new StringBuilder(32768);
    	}

    private void writeBody(ExtendedMolecule mol, boolean hasCoordinates) {
        mMolfile.append("M  V30 BEGIN CTAB\n");
        mMolfile.append("M  V30 COUNTS " + mol.getAllAtoms() + " " + mol.getAllBonds() + " 0 0 0\n");
        mMolfile.append("M  V30 BEGIN ATOM\n");

        for (int atom=0; atom<mol.getAllAtoms(); atom++) {
            mMolfile.append("M  V30 " + (atom + 1));

            if (mol.getAtomList(atom) != null) {
                // mMolfile.append(" L");
                int[] atomList = mol.getAtomList(atom);
                boolean notlist = (mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0;
                mMolfile.append(notlist ? " NOT[" : " [");
                for (int i = 0;i < atomList.length;i++) {
                    if(i > 0) {
                        mMolfile.append(",");
                    	}
                    String label = Molecule.cAtomLabel[atomList[i]];
                    switch (label.length()) {
                        case 1:
                            mMolfile.append(label);
                            break;
                        case 2:
                            mMolfile.append(label);
                            break;
                        case 3:
                            mMolfile.append(label);
                            break;
                        default:
                            mMolfile.append("?");
                            break;
                    	}
                	}
                mMolfile.append("]");
            	}
            else if((mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0) {
                mMolfile.append(" A");
            	}
            else {
                mMolfile.append(" " + mol.getAtomLabel(atom));
            	}

            if (hasCoordinates) {
                mMolfile.append(" " + ((double)((int)(PRECISION_FACTOR * mScalingFactor * mol.getAtomX(atom))) / PRECISION_FACTOR));
                mMolfile.append(" " + ((double)((int)(PRECISION_FACTOR * mScalingFactor * -mol.getAtomY(atom))) / PRECISION_FACTOR));
                mMolfile.append(" " + ((double)((int)(PRECISION_FACTOR * mScalingFactor * -mol.getAtomZ(atom))) / PRECISION_FACTOR));
            	}
            else {
                mMolfile.append(" 0 0 0");
            	}

            mMolfile.append(" " + mol.getAtomMapNo(atom));

            if (mol.getAtomCharge(atom) != 0) {
                mMolfile.append(" CHG=" + mol.getAtomCharge(atom));
            	}

            if (mol.getAtomRadical(atom) != 0) {
                mMolfile.append(" RAD=");
                switch (mol.getAtomRadical(atom)) {
                    case Molecule.cAtomRadicalStateS:
                        mMolfile.append("1");
                        break;
                    case Molecule.cAtomRadicalStateD:
                        mMolfile.append("2");
                        break;
                    case Molecule.cAtomRadicalStateT:
                        mMolfile.append("3");
                        break;
                	}
            	}

            if (mol.getAtomParity(atom) == Molecule.cAtomParity1
             || mol.getAtomParity(atom) == Molecule.cAtomParity2) {
                mMolfile.append(" CFG=");
                if (mol.getAtomParity(atom) == Molecule.cAtomParity1) {
                    mMolfile.append("1");
                	}
                else {
                    mMolfile.append("2");
                	}
            	}

            if (mol.getAtomMass(atom) != 0) {
                mMolfile.append(" MASS=" + mol.getAtomMass(atom));
            	}

            int valence = mol.getAtomAbnormalValence(atom);
            if (valence != -1) {
                mMolfile.append(" VAL=" + ((valence == 0) ? "-1" : valence));
            	}

            int hydrogenFlags = Molecule.cAtomQFHydrogen & mol.getAtomQueryFeatures(atom);
            if (hydrogenFlags == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen)) {
                mMolfile.append(" HCOUNT=2"); // at least 2 hydrogens
            	}
            else if(hydrogenFlags == Molecule.cAtomQFNot0Hydrogen) {
                mMolfile.append(" HCOUNT=1"); // at least 1 hydrogens
            	}
            else if(hydrogenFlags == (Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen)) {
                mMolfile.append(" HCOUNT=-1"); // no hydrogens
            	}
            else if(hydrogenFlags == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen)) {
                mMolfile.append(" HCOUNT=1"); // use at least 1 hydrogens as closest match for exactly one
            	}

            int substitution = mol.getAtomQueryFeatures(atom) & (Molecule.cAtomQFMoreNeighbours | Molecule.cAtomQFNoMoreNeighbours);
            if (substitution != 0) {
                if ((substitution & Molecule.cAtomQFMoreNeighbours) != 0) {
                    mMolfile.append(" SUBST=" + (mol.getAllConnAtoms(atom) + 1));
                	}
                else {
                    mMolfile.append(" SUBST=-1");
                	}
            	}

            int ringFeatures = mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFRingState;
            if (ringFeatures != 0) {
                switch(ringFeatures) {
                    case Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds:
                        mMolfile.append(" RBCNT=-1");
                        break;
                    case Molecule.cAtomQFNotChain:
                        mMolfile.append(" RBCNT=2"); // any ring atom; there is no MDL equivalent
                        break;
                    case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds:
                        mMolfile.append(" RBCNT=2");
                        break;
                    case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot4RingBonds:
                        mMolfile.append(" RBCNT=3");
                        break;
                    case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds:
                        mMolfile.append(" RBCNT=4");
                        break;
                	}
            	}

            mMolfile.append("\n");
        	}

        mMolfile.append("M  V30 END ATOM\n");
        mMolfile.append("M  V30 BEGIN BOND\n");

        for (int bond=0; bond<mol.getAllBonds(); bond++) {
            mMolfile.append("M  V30 " + (bond + 1));

            int order,stereo;
            switch (mol.getBondType(bond)) {
                case Molecule.cBondTypeSingle:
                    order = 1;
                    stereo = 0;
                    break;
                case Molecule.cBondTypeDouble:
                    order = 2;
                    stereo = 0;
                    break;
                case Molecule.cBondTypeTriple:
                    order = 3;
                    stereo = 0;
                    break;
                case Molecule.cBondTypeDown:
                    order = 1;
                    stereo = 3;
                    break;
                case Molecule.cBondTypeUp:
                    order = 1;
                    stereo = 1;
                    break;
                case Molecule.cBondTypeCross:
                    order = 2;
                    stereo = 2;
                    break;
                case Molecule.cBondTypeDelocalized:
                    order = 4;
                    stereo = 0;
                    break;
                default:
                    order = 1;
                    stereo = 0;
                    break;
            	}

            // if query features cannot be expressed exactly stay on the loosely defined side
            int bondType = mol.getBondQueryFeatures(bond) & Molecule.cBondQFBondTypes;
            if (bondType != 0) {
                if (bondType == Molecule.cBondQFDelocalized) {
                    order = 4; // aromatic
                	}
                else if (bondType == (Molecule.cBondQFSingle | Molecule.cBondQFDouble)) {
                    order = 5; // single or double
                	}
                else if (bondType == (Molecule.cBondQFSingle | Molecule.cBondQFDelocalized)) {
                    order = 6; // single or aromatic
                	}
                else if (bondType == (Molecule.cBondQFDouble | Molecule.cBondQFDelocalized)) {
                    order = 7; // single or double
                	}
                else {
                    order = 8; // any
                	}
            	}

            mMolfile.append(" " + order
                            + " " + (mol.getBondAtom(0,bond) + 1)
                            + " " + (mol.getBondAtom(1,bond) + 1));

            if (stereo != 0) {
                mMolfile.append(" CFG=" + stereo);
            	}

            int ringState = mol.getBondQueryFeatures(bond) & Molecule.cBondQFRingState;
            int topology = (ringState == 0) ? 0 : (ringState == Molecule.cBondQFRing) ? 1 : 2;

            if (topology != 0) {
                mMolfile.append(" TOPO=" + topology);
            	}

            mMolfile.append("\n");
        	}

        mMolfile.append("M  V30 END BOND\n");

        boolean paritiesFound = false;
        int absAtomsCount = 0;
        int[] orAtomsCount = new int[Molecule.cESRMaxGroups];
        int[] andAtomsCount = new int[Molecule.cESRMaxGroups];
        for(int atom = 0;atom < mol.getAtoms();atom++) {
            if (mol.getAtomParity(atom) == Molecule.cAtomParity1
             || mol.getAtomParity(atom) == Molecule.cAtomParity2) {
                paritiesFound = true;
                int type = mol.getAtomESRType(atom);
                if (type == Molecule.cESRTypeAnd) {
                    andAtomsCount[mol.getAtomESRGroup(atom)]++;
                	}
                else if (type == Molecule.cESRTypeOr) {
                    orAtomsCount[mol.getAtomESRGroup(atom)]++;
                	}
                else {
                    absAtomsCount++;
                	}
            	}
        	}

        int absBondsCount = 0;
        int[] orBondsCount = new int[Molecule.cESRMaxGroups];
        int[] andBondsCount = new int[Molecule.cESRMaxGroups];
        for (int bond=0; bond<mol.getBonds(); bond++) {
            if (mol.getBondOrder(bond) != 2
             && (mol.getBondParity(bond) == Molecule.cBondParityEor1
              || mol.getBondParity(bond) == Molecule.cBondParityZor2)) {
                paritiesFound = true;
                int type = mol.getBondESRType(bond);
                if(type == Molecule.cESRTypeAnd) {
                    andBondsCount[mol.getBondESRGroup(bond)]++;
                	}
                else if (type == Molecule.cESRTypeOr) {
                    orBondsCount[mol.getBondESRGroup(bond)]++;
                	}
                else {
                    absBondsCount++;
                	}
            	}
        	}

        if(paritiesFound) {
            mMolfile.append("M  V30 BEGIN COLLECTION\n");
            if(absAtomsCount != 0) {
                mMolfile.append("M  V30 MDLV30/STEABS ATOMS=(" + absAtomsCount);
                for(int atom = 0;atom < mol.getAtoms();atom++) {
                    if((mol.getAtomParity(atom) == Molecule.cAtomParity1
                        || mol.getAtomParity(atom) == Molecule.cAtomParity2)
                       && mol.getAtomESRType(atom) == Molecule.cESRTypeAbs) {
                        mMolfile.append(" " + (atom + 1));
                    	}
                	}
                mMolfile.append(")\n");
            	}
            if(absBondsCount != 0) {
                mMolfile.append("M  V30 MDLV30/STEABS BONDS=(" + absBondsCount);
                for(int bond = 0;bond < mol.getBonds();bond++) {
                    if(mol.getBondOrder(bond) != 2
                     && (mol.getBondParity(bond) == Molecule.cBondParityEor1
                      || mol.getBondParity(bond) == Molecule.cBondParityZor2)
                     && mol.getBondESRType(bond) == Molecule.cESRTypeAbs) {
                        mMolfile.append(" " + (bond + 1));
                    	}
                	}
                mMolfile.append(")\n");
            	}
            for(int group = 0;group < Molecule.cESRMaxGroups;group++) {
                if(orAtomsCount[group] != 0) {
                    mMolfile.append("M  V30 MDLV30/STEREL" + (group + 1) + " ATOMS=(" + orAtomsCount[group]);
                    for(int atom = 0;atom < mol.getAtoms();atom++) {
                        if((mol.getAtomParity(atom) == Molecule.cAtomParity1
                            || mol.getAtomParity(atom) == Molecule.cAtomParity2)
                           && mol.getAtomESRType(atom) == Molecule.cESRTypeOr
                           && mol.getAtomESRGroup(atom) == group) {
                            mMolfile.append(" " + (atom + 1));
                        	}
                    	}
                    mMolfile.append(")\n");
                	}
                if(andAtomsCount[group] != 0) {
                    mMolfile.append("M  V30 MDLV30/STERAC" + (group + 1) + " ATOMS=(" + andAtomsCount[group]);
                    for(int atom = 0;atom < mol.getAtoms();atom++) {
                        if((mol.getAtomParity(atom) == Molecule.cAtomParity1
                            || mol.getAtomParity(atom) == Molecule.cAtomParity2)
                           && mol.getAtomESRType(atom) == Molecule.cESRTypeAnd
                           && mol.getAtomESRGroup(atom) == group) {
                            mMolfile.append(" " + (atom + 1));
                        	}
                    	}
                    mMolfile.append(")\n");
                	}
                if(orBondsCount[group] != 0) {
                    mMolfile.append("M  V30 MDLV30/STEREL" + (group + 1) + " BONDS=(" + orBondsCount[group]);
                    for(int bond = 0;bond < mol.getBonds();bond++) {
                        if(mol.getBondOrder(bond) != 2
                         && (mol.getBondParity(bond) == Molecule.cBondParityEor1
                          || mol.getBondParity(bond) == Molecule.cBondParityZor2)
                         && mol.getBondESRType(bond) == Molecule.cESRTypeOr
                         && mol.getBondESRGroup(bond) == group) {
                            mMolfile.append(" " + (bond + 1));
                        	}
                    	}
                    mMolfile.append(")\n");
                	}
                if(andBondsCount[group] != 0) {
                    mMolfile.append("M  V30 MDLV30/STERAC" + (group + 1) + " BONDS=(" + andBondsCount[group]);
                    for(int bond = 0;bond < mol.getBonds();bond++) {
                        if(mol.getBondOrder(bond) != 2
                         && (mol.getBondParity(bond) == Molecule.cBondParityEor1
                          || mol.getBondParity(bond) == Molecule.cBondParityZor2)
                         && mol.getBondESRType(bond) == Molecule.cESRTypeAnd
                         && mol.getBondESRGroup(bond) == group) {
                            mMolfile.append(" " + (bond + 1));
                        	}
                    	}
                    mMolfile.append(")\n");
                	}
            	}
            mMolfile.append("M  V30 END COLLECTION\n");
        	}

        mMolfile.append("M  V30 END CTAB\n");
    	}

    public String getMolfile() {
        return mMolfile.toString();
    	}

    public float getScalingFactor() {
        return mScalingFactor;
    	}

    public void writeMolfile(Writer theWriter) throws IOException {
        theWriter.write(mMolfile.toString());
	    }
	}
