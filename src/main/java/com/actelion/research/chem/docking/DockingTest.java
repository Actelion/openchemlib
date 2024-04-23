package com.actelion.research.chem.docking;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.docking.DockingEngine.DockingResult;
import com.actelion.research.chem.docking.DockingEngine.ScoringFunction;
import com.actelion.research.chem.io.Mol2FileParser;

public class DockingTest {

	public static final String DUDE_PLUS_DIR = "C:\\Users\\wahljo1\\benchmark_phesa\\DUD_E-PLUS\\DUDE-plus-rev1\\ace";
	
	 
	
	public static void main(String[] args) throws Exception {
		String receptor = DUDE_PLUS_DIR +"/" + "protein.mol2";
		String ligand = DUDE_PLUS_DIR +"/"+ "ligand.mol2";
		Mol2FileParser m2fp = new Mol2FileParser();

		Molecule3D lig = new Molecule3D(m2fp.load(new File(ligand)));
		Canonizer can = new Canonizer(lig);
		DockingUtils.repairMolecule3D(lig);
		StereoMolecule canMol = can.getCanMolecule(true);
		lig = new Molecule3D(canMol);
		lig.ensureHelperArrays(Molecule.cHelperCIP);

		int[] ligAtomTypes = new int[lig.getAtoms()];


		//mols.values().forEach(e -> e.forEach(i -> rec.addMolecule(i)));
		Molecule3D rec = m2fp.load(receptor);
		rec.normalizeAmbiguousBonds();
		can = new Canonizer(rec);
		rec = new Molecule3D(rec);
		rec.ensureHelperArrays(StereoMolecule.cHelperCIP);
		DockingUtils.repairMolecule3D(rec);
		
		DockingEngine engine = new DockingEngine(rec, lig,50,10,8.0,ScoringFunction.CHEMPLP);
		StereoMolecule toDock = new StereoMolecule(lig);
		toDock.ensureHelperArrays(Molecule.cHelperParities);
		DockingResult result = engine.dockMolecule(toDock);

		System.out.println(result.getScore());
		System.out.println(result.getContributions());
		
		String delim = System.getProperty("line.separator");
		File outfile = new File("docked_random.sdf");
		   FileWriter fileWriter;
		try {
			
		    fileWriter = new FileWriter(outfile);
		    MolfileCreator mfc = new MolfileCreator(result.getPose());
		    mfc.writeMolfile(fileWriter);
		    fileWriter.write(delim);
		    fileWriter.write("$$$$");
		    fileWriter.write(delim);
		    fileWriter.flush();
			fileWriter.close();
			System.out.println("written");
		} catch (IOException e) {
			// TODO Auto-generated catch block
		e.printStackTrace();
		}
			

		
	}
}
