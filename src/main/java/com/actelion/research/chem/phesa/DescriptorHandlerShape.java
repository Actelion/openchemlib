package com.actelion.research.chem.phesa;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.*;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer.PheSASetting;
import com.actelion.research.chem.alignment3d.transformation.Rotation;
import com.actelion.research.chem.alignment3d.transformation.Transformation;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.alignment3d.transformation.Translation;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.descriptor.DescriptorConstants;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorInfo;
import com.actelion.research.chem.phesaflex.FlexibleShapeAlignment;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;



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

		

	private static final int CONFORMATIONS = 200;
	public static final int SIZE_CUTOFF = 200;


	private static DescriptorHandlerShape INSTANCE;
	
	public static final PheSAMolecule FAILED_OBJECT = new PheSAMolecule();


	private  boolean singleBaseConformation; // take conformation of base molecule as is and don't generate conformers
	
	private List<Transformation> preProcessTransformations;
	
	private StereoMolecule[] previousAlignment;// = new StereoMolecule[2];
	
	private double[] previousPhesaResult;
	
	private PheSASetting phesaSetting;

	protected int maxConfs;
	
	protected double ppWeight;
	
	protected boolean flexible;
	
	protected ThreadMaster threadMaster;
	
	// Maximum number of tries to generate conformers with the torsion rule based conformer generator from Thomas Sander
	
	

	
	private ConformerSetGenerator conformerGenerator;
	
	public DescriptorHandlerShape() {
		this(false,CONFORMATIONS,0.5);
	}
	
	public DescriptorHandlerShape(boolean useSingleBaseConformation) {
		this(useSingleBaseConformation,CONFORMATIONS,0.5);
	}
	
	public DescriptorHandlerShape(boolean useSingleBaseConformation, double ppWeight) {
		this(useSingleBaseConformation,CONFORMATIONS,ppWeight);
	}
	
	public DescriptorHandlerShape(int maxConfs,double ppWeight) {
		this(false,maxConfs,ppWeight);
	}

	/**
	 *
	 * @param ppWeight similarity weight for the pharmacophore in total similarity.
	 */
	public DescriptorHandlerShape(double ppWeight) {
		this(false,CONFORMATIONS,ppWeight);
	}
	
	public DescriptorHandlerShape(boolean useSingleBaseConformation,int maxConfs, double ppWeight) {
		flexible = false;
		singleBaseConformation = useSingleBaseConformation;
		this.maxConfs = maxConfs;
		this.ppWeight = ppWeight;
		init();
		conformerGenerator = new ConformerSetGenerator(maxConfs);
		conformerGenerator.setThreadMaster(threadMaster);
		preProcessTransformations = new ArrayList<>();
		phesaSetting = new PheSASetting();

	}
	
	public PheSASetting getPhesaSetting() {
		return phesaSetting;
	}

	public void setPhesaSetting(PheSASetting phesaSetting) {
		this.phesaSetting = phesaSetting;
	}

	public void setThreadMaster(ThreadMaster tm) {
		this.threadMaster = tm;
	}
		
	
	public List<Transformation> getPreProcessTransformations() {
		return preProcessTransformations;
	}

	public PheSAMolecule createDescriptor(ConformerSet confSet) {
		preProcessTransformations = new ArrayList<>();
		try {

			init();
			
			ArrayList<MolecularVolume> molecularVolumes = new ArrayList<MolecularVolume>(); 
			
			StereoMolecule mol = confSet.first().toMolecule();
			if(mol.getAtoms()>SIZE_CUTOFF)
				return FAILED_OBJECT;
			MolecularVolume refMolVol = new MolecularVolume(mol);
			MolecularVolume molVol;
			
			for(Conformer conformer : confSet) {
	            if(conformer==null) {
	
	                break;
	
	            } else {
					molVol = new MolecularVolume(refMolVol,conformer);
					Coordinates com = molVol.getCOM();
					Rotation rotation = molVol.preProcess(conformer);
					TransformationSequence transformation = new TransformationSequence();
					transformation.addTransformation(rotation.getInvert());
					transformation.addTransformation(new Translation(new double[] {com.x,com.y,com.z}));
					preProcessTransformations.add(transformation);
					molecularVolumes.add(molVol);
	            }
	        }
			 mol = confSet.first().toMolecule();

        return new PheSAMolecule(mol,molecularVolumes);
		}
		catch(Exception e) {
			return FAILED_OBJECT;
		}
	}

	
	public void init() {
		previousAlignment = new StereoMolecule[2];
		
	}
	

	
	/**
	 * the ShapeDescriptor consists of a whole ensemble of MolecularVolumes (MolecularGaussians),
	 * obtained from a conformational search algorithm
	*/
	
	public PheSAMolecule createDescriptor(StereoMolecule mol) {
		StereoMolecule shapeMolecule = new StereoMolecule(mol);
		shapeMolecule.ensureHelperArrays(StereoMolecule.cHelperCIP);
		boolean has3Dcoordinates = false;
		for (int atom=1; atom<mol.getAllAtoms(); atom++) {
			if (Math.abs(mol.getAtomZ(atom) - mol.getAtomZ(0)) > 0.1) {
				has3Dcoordinates = true;
				break;
				}
			}
		shapeMolecule.stripSmallFragments();
		Canonizer can = new Canonizer(shapeMolecule);
		shapeMolecule = can.getCanMolecule(true);

		ConformerSet confSet = new ConformerSet();
		if (!singleBaseConformation) {
			confSet = conformerGenerator.generateConformerSet(shapeMolecule);
		}
		
		else if(!has3Dcoordinates) {
			return FAILED_OBJECT;
		}
			
		else {
			if(shapeMolecule.getAllAtoms()-shapeMolecule.getAtoms()==0) {
				System.err.println("missing hydrogens in 3D structure");
				return FAILED_OBJECT;
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
		double[] result = PheSAAlignmentOptimizer.align(query,base,bestPair,phesaSetting);
		this.setPreviousAlignment(bestPair);
		this.setPreviousPheSAResult(result);
		if(flexible) {
			FlexibleShapeAlignment fsa = new FlexibleShapeAlignment(bestPair[0],bestPair[1]);
			fsa.setSettings(phesaSetting);
			result = fsa.align();
			this.setPreviousPheSAResult(result);
		}
		
		return (float)result[0];
	}



	public StereoMolecule[] getPreviousAlignment() {
		return this.previousAlignment;
	}

	/***
	 * additional output:
	 * element 0: total similarity (identical to getSimilarity(...))
	 * element 1: pharmacophore similarity
	 * element 2: shape similarity
	 * element 3: contribution to similarity that originates from additional volumes (incl/excl)
	 * @return
	 */
	public double[] getPreviousPheSAResult() {
		return this.previousPhesaResult;
	}
	
	public void setPreviousAlignment(StereoMolecule[] previousAlignment) {
		this.previousAlignment = previousAlignment;
	}
	
	public void setPreviousPheSAResult(double[] previousPhesaResult) {
		this.previousPhesaResult = previousPhesaResult;
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
		//StereoMolecule mol = shapeMol.getConformer(shapeMol.getVolumes().get(0));
		StereoMolecule mol = new StereoMolecule(shapeMol.getMolecule());
		Canonizer can = new Canonizer(mol, Canonizer.COORDS_ARE_3D);
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

		return decode(new String(arr, StandardCharsets.UTF_8));
		
	}
	
	public boolean calculationFailed(PheSAMolecule o) {


			return o.getVolumes().size()==0;
	
	}
	
	public DescriptorHandlerShape getThreadSafeCopy() {

		DescriptorHandlerShape dhs = new DescriptorHandlerShape();
		dhs.ppWeight = ppWeight;
		dhs.flexible = flexible;
		dhs.maxConfs = maxConfs;

		return dhs;
	}

	public static DescriptorHandlerShape getDefaultInstance(){

		if(INSTANCE==null){
			INSTANCE = new DescriptorHandlerShape();
		}

		return INSTANCE;
	}
	
	public void setMaxConfs(int maxConfs) {
		this.maxConfs = maxConfs;
		init();
	}
	
	public void setFlexible(boolean flexible) {
		this.flexible = flexible;
	}

	public boolean isFlexible() {
		return flexible;
	}

	public double getPpWeight() {
		return ppWeight;
	}
}