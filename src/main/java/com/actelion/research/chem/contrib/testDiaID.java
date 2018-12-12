package com.actelion.research.chem.contrib;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.contrib.DiastereotopicAtomID;
import com.actelion.research.chem.contrib.HydrogenHandler;

public class testDiaID {

    public static void main(String[] args) {
        StereoMolecule mol = new StereoMolecule();
        SmilesParser parser = new SmilesParser();
        try {
            parser.parse(mol, "C1CCCCC1C");
        } catch (Exception e) {
            System.out.println(e.toString());
        }

        String[] ids = DiastereotopicAtomID.getAtomIds(mol);
        System.out.println(ids[0]);

        HydrogenHandler.addImplicitHydrogens(mol);
        ids = DiastereotopicAtomID.getAtomIds(mol);
        System.out.println(ids[0]);

        System.out.println("");

        for (String id : ids) {
            System.out.println(id);
        }


    }

}

