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
 * @author Modest v. Korff, Thomas Sander
 */

package com.actelion.research.chem.io.pdb.converter;

import com.actelion.research.chem.*;
import com.actelion.research.chem.io.pdb.parser.AtomRecord;

import java.util.*;

public class AminoAcids {
	private static final String[][] AA_TEMPLATES = {
            { "gGX`BDdwMULPGzILwXM[jD", "Alanine", "ala" },
            { "dctd@BE]ADf{UYjjihp`GzBfMvCS]Plw^OMtbK]hrwUj}tRfwnbXp", "Arginine", "arg" },
            { "diEL@BDDyInvZjZL`OtiL[lFfzaYn|^{iFLO]Hi`", "Asparagine", "asn" },
            { "diFB@BANEInvZjZLHA~eIc]`twTKMwcw]Hqa{iEL", "Aspartic Acid", "asp" },
            { "gNxhMV@aI[jihj@?SHF{ac]PinpP", "Cysteine", "Cys" },
            { "defB@BAAeInufjihr@?QdqnpZ[jEf{qyndQ{mFMO]hi`","Glutamic Acid", "glu" },
            { "deeL@BdDEInufjihp`GzLfMvCS]Plw^OMtbO]hqi{mEL", "Glutamine", "gln" },
            { "gJX`BDdvu@OtbYnpP", "Glycine", "gly" },
            { "dmwD@ByPQInvVUZjejL`OtyL[lFfzaYn|^{iFLO]Hii{mFLo]hi`", "Histidine", "his" },
            { "diFD@BADf{ejjdrU@_iRXwXMMuBw]xqn{oELO]Hq`", "Isoleucine", "ile" },
            { "diFD@BADf{Yjjhr@?RdqnpZ[jEf{q{ndTp}tcF@", "Leucine", "leu" },
            { "deeD@BdDR[mUjjjL`OtYL[lFfzaYn|^[iDV{QenkP", "Lysine", "lys" },
            { "diFD`JxPBDivzjihI@?RdAndX[oEF{QqnhR[lD", "Methionine", "met" },
            { "dcND@BADf{YU]Zj@@cHC}ASF{AinhV[oGnzQSCwRLZ^{QSKwZL[Vzm@", "Phenylalanine", "phe" },
            { "daFD@BADfyVyjjhr@?PdqnpZ[jEfzQyn|P", "Proline", "pro" },
            { "gNy`BDtf{ZjfHC}Lf[lFmuBv{q@", "Serine", "ser" },
            { "dazL@BAFR[nZjdrT`_hRXwXMMuBw]xqn{oEL", "Threonine", "thr" },
            { "foAP`@BZ@aInvYWejsfjiB@bFB@OttfF{AhwTKF{qywRJXW]Hqi]vbfUwZN[W]hqc]uZfmwUnYw]Di`", "Tryptophane", "trp" },
            { "dknL@BACR[me]]Zj@BHr@?RTqnpZ[jEf{q{ndTp}tcFgntTr}vcFunkS[hd", "Tyrosine", "tyr" },
            { "dazD@BADf{fjjL`OtIL[lFfza[n|Tw]wcF@", "Valine", "val" }
    };

    private static TreeMap<String,StereoMolecule> sShortLabelMap;

    public static StereoMolecule getStructure(String label) {
        return ensureAAMap().get(label.toLowerCase());
    }

    public static Molecule3D createResidue(String label, List<AtomRecord> atomRecordList) {
        StereoMolecule mol = getStructure(label);
        if (mol != null) {
            Molecule3D aminoAcid = new Molecule3D(mol);
            Map<String,AtomRecord> recordMap = new HashMap<>();
            for(AtomRecord record : atomRecordList)
                recordMap.put(record.getAtomName(), record);

            for (int atom=0; atom<aminoAcid.getAllAtoms(); atom++) {
                String atomName = mol.getAtomLabel(atom);
                String customLabel = mol.getAtomCustomLabel(atom);
                if (customLabel != null && customLabel.startsWith("]"))
                    atomName = atomName.concat(customLabel.substring(1));

                AtomRecord record = recordMap.get(atomName);
                if (record != null) {
                    Coordinates coords3d = new Coordinates(record.getX(), record.getY(), record.getZ());
                    aminoAcid.setAtomName(atom, record.getAtomName());
                    aminoAcid.setAtomAmino(atom, record.getResName());
                    aminoAcid.setAtomSequence(atom, record.getSerialId());
                    aminoAcid.setResSequence(atom, record.getResNum());
                    aminoAcid.setAtomAmino(atom, record.getResName());
                    aminoAcid.setAtomChainId(atom, record.getChainID());
                    aminoAcid.setAtomX(atom, coords3d.x);
                    aminoAcid.setAtomY(atom, coords3d.y);
                    aminoAcid.setAtomZ(atom, coords3d.z);
                }
            };

            for (int atom=0; atom<aminoAcid.getAllAtoms(); atom++)
                if (aminoAcid.getAtomName(atom) == null)
                    aminoAcid.markAtomForDeletion(atom);
            aminoAcid.deleteMarkedAtomsAndBonds();

            return aminoAcid;
        }
        return null;
    }

    private static TreeMap<String,StereoMolecule> ensureAAMap() {
        if (sShortLabelMap == null) {
            sShortLabelMap = new TreeMap<>();
            for (String[] template : AA_TEMPLATES) {
                StereoMolecule mol = new IDCodeParserWithoutCoordinateInvention().getCompactMolecule(template[0]);
                mol.setName(template[1]);
                sShortLabelMap.put(template[2], mol);
            }
        }
        return sShortLabelMap;
    }
}
