package com.actelion.research.chem.contrib;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.contrib.DiastereotopicAtomID;
import com.actelion.research.chem.contrib.HydrogenHandler;

public class testChiralDia {

    public static void main(String[] args) {
        StereoMolecule mol = new StereoMolecule();
        SmilesParser parser = new SmilesParser();
        try {
            parser.parse(mol, "C1CCCCC1C");
        } catch (Exception e) {
            System.out.println(e.toString());
        }

        HydrogenHandler.addImplicitHydrogens(mol);

        String[] ids = DiastereotopicAtomID.getAtomIds(mol);
        for (String id : ids) {
            System.out.println(id);
        }
        /*
        System.out.println("With hydrogens");
        HydrogenHandler.addImplicitHydrogens(mol);
        ids = DiastereotopicAtomID.getAtomIds(mol);


        for (String id : ids) {
            System.out.println(id);
        }
*/

    }

}

