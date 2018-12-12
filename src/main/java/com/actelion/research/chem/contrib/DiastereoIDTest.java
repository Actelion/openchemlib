package com.actelion.research.chem.contrib;

import java.io.IOException;
import java.util.TreeMap;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.contrib.HydrogenHandler;
import com.actelion.research.chem.contrib.DiastereotopicAtomID;
import com.actelion.research.chem.io.SDFileParser;
import java.io.File;


public class DiastereoIDTest {


    public static void main(String[] args) throws IOException {

        // we need to testDiaID if the diastereotopic protons generation works
        // starting from a 2D molecule without hydrogen and from a 3D molecule generate by Moloc (with added hydrogen)

        String sdf2d="./src/main/java/com/actelion/research/chem/contrib/test/test-diastereo.sdf";
        String sdf3d="./src/main/java/com/actelion/research/chem/contrib/test/test-diastereo-3d.sdf";

        TreeMap<String,StereoMolecule> mol2d=new TreeMap<String,StereoMolecule>();
        TreeMap<String,StereoMolecule> mol2dh=new TreeMap<String,StereoMolecule>();
        TreeMap<String,StereoMolecule> mol3dh=new TreeMap<String,StereoMolecule>();


        // String file=new File(sdf2d).getCanonicalPath();
        // loading the 2D molecules

        SDFileParser parser=new SDFileParser(sdf2d);
        while (parser.next()) {
            StereoMolecule molecule=parser.getMolecule();
            String idCode=new Canonizer(molecule).getIDCode();
            mol2d.put(idCode, molecule);
            // we add the hydrogens ...
            StereoMolecule moleculeAllH=molecule.getCompactCopy();
            HydrogenHandler.addImplicitHydrogens(moleculeAllH);
            String idCodeAllH=new Canonizer(moleculeAllH).getIDCode();
            mol2dh.put(idCodeAllH, moleculeAllH);

            if (! idCode.equals(idCodeAllH)) {
                System.out.println("ERROR : 2D: "+idCode+" - "+"2DH: "+idCodeAllH);
            }
        }

        // loading the 3D molecules
        parser=new SDFileParser(sdf3d);
        while (parser.next()) {
            StereoMolecule molecule=parser.getMolecule();
            String idCode=new Canonizer(molecule).getIDCode();
            mol3dh.put(idCode, molecule);
        }

        System.out.println("SIZE (should be the same): "+mol2d.size()+" - "+mol2dh.size()+" - "+mol3dh.size());

        // FIRST REPORT: we check if the ids are the same in the 3 sets !
        for (String id : mol2d.keySet()) {
            System.out.println(id+" - "+mol2dh.containsKey(id)+" - "+mol3dh.containsKey(id));
        }

        System.out.println("ID from 3D molecule not present in 2D: ");
        for (String id : mol3dh.keySet()) {
            if (! mol2d.containsKey(id)) {
                System.out.println(id);
            }
        }

        // we check if the molecules have the same number of different atoms !
 /*
        for (String id : mol2d.keySet()) {
            System.out.println("---> "+id);

            String[] dia2=DiastereotopicAtomID.getAtomIds(mol2d);
            String[] dia2h=DiastereotopicAtomID.getAtomIds(mol2dh);
            String[] dia3h=DiastereotopicAtomID.getAtomIds(mol3dh);


            if ((dia2h.length==dia3h.length) && (dia2.length<dia2h.length)) {
                System.out.println("Size of diastereotopic protons: OK");
            } else {
                System.out.println("Size of diastereotopic protons: ERROR");
                System.out.println(dia2.length+" - "+dia2h.length+" - "+dia3h.length);
            }

            for (String idDia2 : dia2) {
                //	System.out.println("Dia2: "+idDia2);
            }

            for (String idDia2h : dia2h) {
                //	System.out.println("Dia2 - H: "+idDia2h);
            }

            /*
            // we check is all the atoms ID of the mol2d are present in mol2dh and mol3dh
            for (String atomID : dia2h) {
                if (mol2dh.getNumberAtoms(atomID)!=mol3dh.getNumberAtoms(atomID)) {
                    System.out.println("ERROR in count dia atoms in mol2D-H and mol3D-H: "+atomID+" - "+eMol2dh.getNumberAtoms(atomID)+" - "+eMol3dh.getNumberAtoms(atomID));
                }
                if (eMol2d.getNumberAtoms(atomID)!=eMol2dh.getNumberAtoms(atomID)) {
                    if (! eMol2dh.getAtomLabel(atomID).equals("H")) {
                        System.out.println("ERROR in count dia atoms in mol2D-H and mol2D: "+atomID+" - "+eMol2dh.getNumberAtoms(atomID)+" - "+eMol2d.getNumberAtoms(atomID));
                    }
                }
            }



        }
*/

    }
}
