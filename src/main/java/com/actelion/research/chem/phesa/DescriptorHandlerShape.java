package com.actelion.research.chem.phesa;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParserWithoutCoordinateInvention;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.descriptor.DescriptorConstants;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorInfo;
import org.openmolecules.chem.conf.gen.ConformerGenerator;
import java.util.ArrayList;



/** 
 * @version: 1.0, February 2018
 * Author: JW
 * a molecular shape descriptor is a shape-ensemble (molecular shapes of generated conformers of a molecule)
 * it also contains information about the atom-connectivity,elements,... so for any molecular shape, the corresponding
 * molecule can be reconstructed 
 * there are different modi for calculating the shape similarity:
 * 0: only takes into account molecule shape for alignment and overlap calculation
 * 1: alignment solely based on shape, but overlap calculation incorporates pharmacophore overlaps
 * 2: both alignment and overlap calculation take a combined score of shape and pharmacophore overlap
 * 3: only the pharmacophore overlap is used for both aligment and overlap calculation -> this is the fastest method!
 * 19 April 2018: performance enhancement by using a cutoff for the calculation of atom-atom overlaps and preculated exp-values with linear interpolation
 * TODO: add Tversky index
 * July 2019: various improvements in the Code, moved to DD_core
*/


public class DescriptorHandlerShape implements DescriptorHandler<PheSAMolecule,StereoMolecule> {

		
	private static final long SEED = 123456789;
	
	private static final int CONFORMATIONS = 50;


	private static DescriptorHandlerShape INSTANCE;
	
	public static final PheSAMolecule FAILED_OBJECT = new PheSAMolecule();


	private  boolean singleBaseConformation; // take conformation of base molecule as is and don't generate conformers
	
	private double[][] transforms;// = ShapeAlignment.initialTransform(2);
	
	private StereoMolecule[] previousAlignment;// = new StereoMolecule[2];

	// Maximum number of tries to generate conformers with the torsion rule based conformer generator from Thomas Sander
	
	

	
	private ConformerSetGenerator conformerGenerator;
	
	public DescriptorHandlerShape() {
		this(false);
	}
	
	public DescriptorHandlerShape(boolean useSingleBaseConformation) {
		singleBaseConformation = useSingleBaseConformation;
		init();

	}
	
	public PheSAMolecule createDescriptor(ConformerSet fullSet) {
		ConformerSet confSet = fullSet.getSubset(CONFORMATIONS);
		init();
		
		ArrayList<MolecularVolume> molecularVolumes = new ArrayList<MolecularVolume>(); 
		
		StereoMolecule mol = confSet.first().toMolecule();
		MolecularVolume refMolVol = new MolecularVolume(mol);
		MolecularVolume molVol;
		
		for(Conformer conformer : confSet) {

            if(conformer==null) {

                break;

            } else {
				molVol = new MolecularVolume(refMolVol,conformer);
				PheSAAlignment.preProcess(conformer, molVol);
				molecularVolumes.add(molVol);
            }
        }
 
        return new PheSAMolecule(mol,molecularVolumes);
	}

	
	public void init() {
		transforms = PheSAAlignment.initialTransform(2);
		previousAlignment = new StereoMolecule[2];
		conformerGenerator = new ConformerSetGenerator(CONFORMATIONS,ConformerGenerator.STRATEGY_LIKELY_RANDOM,false,SEED);
		
	}
	

	
	/**
	 * the ShapeDescriptor consists of a whole ensemble of MolecularVolumes (MolecularGaussians),
	 * obtained from a conformational search algorithm
	*/
	
	public PheSAMolecule createDescriptor(StereoMolecule mol) {
		StereoMolecule shapeMolecule = new StereoMolecule(mol);
		boolean has3Dcoordinates = false;
		for (int atom=1; atom<mol.getAllAtoms(); atom++) {
			if (Math.abs(mol.getAtomZ(atom) - mol.getAtomZ(0)) > 0.1) {
				has3Dcoordinates = true;
				break;
				}
			}
		shapeMolecule.stripSmallFragments();
		new Canonizer(shapeMolecule);
		shapeMolecule.ensureHelperArrays(StereoMolecule.cHelperCIP);

		ConformerSet confSet = new ConformerSet();
		if (!singleBaseConformation) {
			confSet = conformerGenerator.generateConformerSet(shapeMolecule);
		}
		
		else if(!has3Dcoordinates) {
			return FAILED_OBJECT;
		}
			
		else {
			if(shapeMolecule.getAllAtoms()-shapeMolecule.getAtoms()>0) {
				ConformerGenerator.addHydrogenAtoms(shapeMolecule);
			}
			confSet.add(new Conformer(shapeMolecule)); //take input conformation
			
			}

		return this.createDescriptor(confSet);

	}
	/**
	 * calculates the Shape- and/or Pharmacophore similarity of a query molecule with a base molecule
	 * 
	 */
	
	public float getSimilarity(PheSAMolecule query, PheSAMolecule base) {
		StereoMolecule[] bestPair = {query.getMolecule(),base.getMolecule()};
		double similarity = PheSAAlignmentOptimizer.align(query,base,bestPair);
		this.setPreviousAlignment(bestPair);
		return (float)similarity;
	}


	
	public StereoMolecule[] getPreviousAlignment() {
		return this.previousAlignment;
	}
	
	public void setPreviousAlignment(StereoMolecule[] previousAlignment) {
		this.previousAlignment = previousAlignment;
	}
	
	

	

	

	
	public String getVersion() {
		return DescriptorConstants.DESCRIPTOR_ShapeAlign.version;
	}
	
	public DescriptorInfo getInfo() {
		return DescriptorConstants.DESCRIPTOR_ShapeAlign;
	}
	

	

	
	public String encode(PheSAMolecule o) {

		if(calculationFailed(o)){
			return FAILED_STRING;
		}

		ArrayList<MolecularVolume> molVols = null;
		PheSAMolecule shapeMol;

		if(o instanceof PheSAMolecule){
			
			shapeMol = (PheSAMolecule)o;
			molVols = shapeMol.getVolumes();


		} else {
			return FAILED_STRING;
		}
		
		StringBuilder shapeString = new StringBuilder();
		int nrOfMolVols = molVols.size();
		shapeString.append(Integer.toString(nrOfMolVols));
		shapeString.append("   ");
		shapeString.append(molVols.get(0).encodeFull());
		shapeString.append("   ");
		
		for(int i=1;i<nrOfMolVols;i++) {
			shapeString.append(molVols.get(i).encodeCoordsOnly());
			shapeString.append("   ");
		}
			
		shapeString.append("   ");
		StereoMolecule mol = shapeMol.getConformer(shapeMol.getVolumes().get(0));
		Canonizer can = new Canonizer(mol, Canonizer.COORDS_ARE_3D);
		mol = can.getCanMolecule(true);
		String idcoords = can.getEncodedCoordinates(true);
		String idcode = can.getIDCode();
		shapeString.append(idcode);
		shapeString.append("   ");
		shapeString.append(idcoords);
		return shapeString.toString();
	}
	
	public PheSAMolecule decode(String s) {
		try {
			return s == null || s.length() == 0 ? null
					: s.equals(FAILED_STRING) ? FAILED_OBJECT
					:                           getDecodedObject(s);
		} catch (RuntimeException e1) {
			return FAILED_OBJECT;
		}
	}
	
	private PheSAMolecule getDecodedObject(String string64)  {
		String[] splitted = string64.split("   ");
		String idcode = splitted[splitted.length-2];
		String idcoords = splitted[splitted.length-1];
		StereoMolecule mol = new StereoMolecule();
		IDCodeParserWithoutCoordinateInvention parser = new IDCodeParserWithoutCoordinateInvention();
		parser.parse(mol, idcode, idcoords);
		mol.ensureHelperArrays(Molecule.cHelperCIP);
		ArrayList<MolecularVolume> molVols = new ArrayList<MolecularVolume>();
		int nrOfMolVols = Integer.decode(splitted[0]);
		MolecularVolume refMolVol = MolecularVolume.decodeFull(splitted[1], mol);
		molVols.add(refMolVol);
		for(int i=2;i<nrOfMolVols+1;i++) {
			molVols.add(MolecularVolume.decodeCoordsOnly(splitted[i], refMolVol));
			
		}

		PheSAMolecule shapeMol = new PheSAMolecule(mol,molVols);
		return shapeMol;
		
	}
	
	public PheSAMolecule decode(byte[] arr) {

		return decode(new String(arr));
		
	}
	
	public boolean calculationFailed(PheSAMolecule o) {


			return o.getVolumes().size()==0;
	
	}
	
	public DescriptorHandlerShape getThreadSafeCopy() {

		DescriptorHandlerShape dhs = new DescriptorHandlerShape();

		return dhs;
	}

	public static DescriptorHandlerShape getDefaultInstance(){

		if(INSTANCE==null){
			INSTANCE = new DescriptorHandlerShape();
		}

		return INSTANCE;
	}

}