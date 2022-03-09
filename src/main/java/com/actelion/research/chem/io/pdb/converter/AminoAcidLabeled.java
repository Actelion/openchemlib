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
 * @author Modest v. Korff
 */

package com.actelion.research.chem.io.pdb.converter;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.io.pdb.parser.AtomRecord;

import java.lang.StringBuilder;
import java.util.Map;
import java.util.Optional;
import java.util.stream.IntStream;

/**
 * AminoAcidLabeled
 * Created by korffmo1 on 17.04.18.
 */
public class AminoAcidLabeled {

    private StereoMolecule mol;

    private String name;

    // Three letter, upper case
    private String abbreviation;

    public AminoAcidLabeled(StereoMolecule mol, String name, String abbreviation) {
        this.mol = mol;
        this.name = name;
        this.abbreviation = abbreviation;
    }
    
    //TODO: set Atom Names, res Names, res numbers!
    public Molecule3D createResidue(Map<String,AtomRecord> recordMap) {
    	Molecule3D residue;
    	Molecule3D aminoAcid = new Molecule3D(mol);
    	try {
    	if(mol.getAtoms()-recordMap.size()!=1) {
    		throw new RuntimeException();
    	}
    	IntStream.range(0,mol.getAtoms()).forEach(atom -> {
    		Optional<String> customLabel = Optional.ofNullable(mol.getAtomCustomLabel(atom));
    		StringBuilder sb = new StringBuilder(customLabel.orElse(" "));
			sb.setCharAt(0, mol.getAtomLabel(atom).charAt(0));
			AtomRecord record = recordMap.get(sb.toString());
			Coordinates coords3d = new Coordinates(record.getX(),record.getY(),record.getZ());
			aminoAcid.setAtomName(atom, record.getAtomName());
			aminoAcid.setAtomAmino(atom, record.getResName());
			aminoAcid.setAtomSequence(atom,record.getSerialId());
			aminoAcid.setResSequence(atom, record.getResNum());
			aminoAcid.setAtomAmino(atom, record.getResName());
			aminoAcid.setAtomChainId(atom, record.getChainID());
			aminoAcid.setAtomX(atom, coords3d.x);
			aminoAcid.setAtomY(atom, coords3d.y);
			aminoAcid.setAtomZ(atom, coords3d.z);
    	});
    	residue = aminoAcid;
    	}
    	catch(Exception e) {
    		residue = null;
    	}

    	return residue;
    }
    


    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("AminoAcidLabeled{");
        sb.append("name='").append(name).append('\'');
        sb.append(", abbreviation='").append(abbreviation).append('\'');
        sb.append('}');
        return sb.toString();
    }
}
