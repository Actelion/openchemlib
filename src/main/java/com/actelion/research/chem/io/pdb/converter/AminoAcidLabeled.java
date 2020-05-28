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
 * <p>Copyright: Idorsia Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Idorsia Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
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
