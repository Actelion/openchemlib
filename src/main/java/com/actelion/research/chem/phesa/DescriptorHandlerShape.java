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


public class DescriptorHandlerShape implements DescriptorHandler<ShapeMolecule,StereoMolecule> {

		
	private static final long SEED = 123456789;
	
	private static final int CONFORMATIONS = 20;


	private static DescriptorHandlerShape INSTANCE;
	
	public static final ShapeMolecule FAILED_OBJECT = new ShapeMolecule();


	private  boolean singleBaseConformation; // take conformation of base molecule as is and don't generate conformers
	
	private double[][] transforms;// = ShapeAlignment.initialTransform(2);
	
	private StereoMolecule[] previousAlignment;// = new StereoMolecule[2];

	// Maximum number of tries to generate conformers with the torsion rule based conformer generator from Thomas Sander
	
	private static final int MAX_NUM_TRIES = 10000;
	

	
	private ConformerSetGenerator conformerGenerator;
	
	public DescriptorHandlerShape() {
		this(false);
	}
	
	public DescriptorHandlerShape(boolean useSingleBaseConformation) {
		singleBaseConformation = useSingleBaseConformation;
		init();

	}
	
	public ShapeMolecule createDescriptor(ConformerSet fullSet) {
		ConformerSet confSet = fullSet.getSubset(CONFORMATIONS);
		init();
		
		ArrayList<MolecularVolume> molecularVolumes = new ArrayList<MolecularVolume>(); 
		
		StereoMolecule conf = confSet.first().toMolecule();
		MolecularVolume molVol;
		
		for(Conformer conformer : confSet) {

            if(conformer==null) {

                break;

            } else {
				conf = conformer.toMolecule(null);
				new Canonizer(conf);
				conf.ensureHelperArrays(StereoMolecule.cHelperCIP);
				molVol = new MolecularVolume(conf);
				ShapeAlignment.preProcess(conf, molVol);
				molecularVolumes.add(molVol);
            }
        }
 
        return new ShapeMolecule(conf,molecularVolumes);
	}

	
	public void init() {
		transforms = ShapeAlignment.initialTransform(2);
		previousAlignment = new StereoMolecule[2];
		conformerGenerator = new ConformerSetGenerator(20,ConformerGenerator.STRATEGY_LIKELY_RANDOM,false,SEED);
		
	}
	

	
	/**
	 * the ShapeDescriptor consists of a whole ensemble of MolecularVolumes (MolecularGaussians),
	 * obtained from a conformational search algorithm
	*/
	
	public ShapeMolecule createDescriptor(StereoMolecule mol) {
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
			confSet = conformerGenerator.generateConformerSet(mol);
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
	
	public float getSimilarity(ShapeMolecule query, ShapeMolecule base) {
		StereoMolecule[] bestPair = {base.getMolecule(),query.getMolecule()};
		ArrayList<MolecularVolume>  molVolQuery = query.getVolumes();
		ArrayList<MolecularVolume>  molVolBase = base.getVolumes();
		float maxSimilarity=0.0f;
		double Oaa = 0.0;
		double Obb = 0.0;
		double ppOaa = 0.0;
		double ppObb = 0.0;
		double [] alignment = {1.0,0.0,0.0,0.0,0.0,0.0,0.0};
		
		for(MolecularVolume queryVol: molVolQuery) {
			for(MolecularVolume baseVol: molVolBase) {
				MolecularVolume baseVolume = new MolecularVolume(baseVol);
				MolecularVolume queryVolume = new MolecularVolume(queryVol);
				ShapeAlignment shapeAlign = new ShapeAlignment(queryVolume, baseVolume);
				Oaa = shapeAlign.getSelfAtomOverlapRef();
				Obb = shapeAlign.getSelfAtomOverlapFit();
				ppOaa = shapeAlign.getSelfPPOverlapRef();
				ppObb = shapeAlign.getSelfPPOverlapFit();
				EvaluableOverlap eval = new EvaluableOverlap(shapeAlign, new double[7]);
				ShapeOptimizerLBFGS opt = new ShapeOptimizerLBFGS(200,0.001);
				double ppScaling = 1.0;

				//iterate over all initial alignments (necessary since optimizer just finds next local minimum, so we need different initial guesses
				
				for(double [] transform:transforms) { //iterate over all initial alignments (necessary since optimizer just finds next local minimum, so we need different initial guesses
					eval.setState(transform);
					double[] bestTransform = opt.optimize(eval);
					double atomOverlap = 0.0;
					double ppOverlap = 0.0;
					float similarity = 0.0f;
					atomOverlap = shapeAlign.getTotalAtomOverlap(bestTransform);
					ppOverlap = shapeAlign.getTotalPPOverlap(bestTransform);
					float atomSimilarity = (float)(atomOverlap/(Oaa+Obb-atomOverlap));
					float ppSimilarity = 0.0f;
					if(shapeAlign.getRefMolGauss().getPPGaussians().size()==0 && shapeAlign.getMolGauss().getPPGaussians().size()==0 )
						ppSimilarity = 1.0f;
					else ppSimilarity=(float)(ppOverlap/(ppOaa+ppObb-ppOverlap));
					if(ppSimilarity>1.0) //can happen because of weights
						ppSimilarity = 1.0f;
					similarity = (1.0f/(1+(float)ppScaling))* (atomSimilarity + (float)ppScaling*ppSimilarity) ;
					if (similarity>maxSimilarity) {
						maxSimilarity = similarity;
						alignment = bestTransform;
						bestPair[1] = base.getConformer(baseVol);
						bestPair[0] = query.getConformer(queryVol);
				}
				}
				
			}
		}
		ShapeAlignment.rotateMol(bestPair[1], alignment);
		

		this.setPreviousAlignment(bestPair);
		return maxSimilarity;	
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
	

	

	
	public String encode(ShapeMolecule o) {

		if(calculationFailed(o)){
			return FAILED_STRING;
		}

		ArrayList<MolecularVolume> molVols = null;
		ShapeMolecule shapeMol;

		if(o instanceof ShapeMolecule){
			
			shapeMol = (ShapeMolecule)o;
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
		StereoMolecule mol = shapeMol.getMolecule();
		Canonizer can = new Canonizer(mol, Canonizer.COORDS_ARE_3D);
		mol = can.getCanMolecule(true);
		String idcoords = can.getEncodedCoordinates(true);
		String idcode = can.getIDCode();
		shapeString.append(idcode);
		shapeString.append("   ");
		shapeString.append(idcoords);
		return shapeString.toString();
	}
	
	public ShapeMolecule decode(String s) {
		try {
			return s == null || s.length() == 0 ? null
					: s.equals(FAILED_STRING) ? FAILED_OBJECT
					:                           getDecodedObject(s);
		} catch (RuntimeException e1) {
			return FAILED_OBJECT;
		}
	}
	
	private ShapeMolecule getDecodedObject(String string64)  {
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

		ShapeMolecule shapeMol = new ShapeMolecule(mol,molVols);
		return shapeMol;
		
	}
	
	public ShapeMolecule decode(byte[] arr) {

		return decode(new String(arr));
		
	}
	
	public boolean calculationFailed(ShapeMolecule o) {


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