package com.actelion.research.chem.phesa;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.IDCodeParserWithoutCoordinateInvention;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer.SimilarityMode;
import com.actelion.research.chem.alignment3d.transformation.ExponentialMap;
import com.actelion.research.chem.alignment3d.transformation.Quaternion;
import com.actelion.research.chem.alignment3d.transformation.Transformation;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.alignment3d.transformation.Translation;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.optimization.OptimizerLBFGS;
import java.util.Arrays;
import java.util.Base64;
import java.util.Base64.Decoder;
import java.util.Base64.Encoder;
import java.util.stream.IntStream;

/** 
 * @version: 1.0, February 2018
 * Author: J. Wahl
 * this class provides functionalities to calculate the overlap between two molecules
 * 
*/


public class PheSAAlignment {

	private ShapeVolume refMolGauss;
	private ShapeVolume molGauss;
	private double ppWeight;
	public enum axis {X,Y,Z};
	public static final double TVERSKY_COEFFICIENT = 0.95;


	
	
	
	public PheSAAlignment(StereoMolecule refMol, StereoMolecule mol,double ppWeight) {
		this.ppWeight = ppWeight;
		this.refMolGauss = new MolecularVolume(refMol);
		this.molGauss = new MolecularVolume(mol);
	}
	
	public PheSAAlignment(StereoMolecule refMol, StereoMolecule mol) {
		this(refMol,mol,0.5);
	}
	
	public PheSAAlignment(MolecularVolume refMolGauss, MolecularVolume molGauss) {
		this(refMolGauss,molGauss,0.5);
	}
	
	public PheSAAlignment(ShapeVolume refMolGauss, ShapeVolume molGauss,double ppWeight) {
		this.ppWeight = ppWeight;
		this.refMolGauss= refMolGauss;
		this.molGauss = molGauss;
	}
	


	public ShapeVolume getRefMolGauss() {
		return refMolGauss;
	}

	public ShapeVolume getMolGauss() {
		return molGauss;
	}

	
	public static void rotateMolAroundAxis180(Conformer conf,axis a) {
			IntStream.range(0,conf.getSize()).forEach(i -> {
				Coordinates coords = conf.getCoordinates(i);
				rotateCoordsAroundAxis180(coords,a);
			});
	}
	
	public static void rotateCoordsAroundAxis180(Coordinates coords,axis a) {
		if (a == axis.X) {
			coords.y = -coords.y;
			coords.z = -coords.z;

		}
		else if (a == axis.Y) {
			coords.x = -coords.x;
			coords.z = -coords.z;
		}
		
		else  {
			coords.x = -coords.x;
			coords.y = -coords.y;

		}

	}
	

	
	/**
	 * .
	 * generate initial orientations of the molecule: 
     * mode1: 4 orientations: initial orientation and 180 degree rotation about each axis
	 * mode2: mode1 and 90 degree rotations about each axis
	 * a transformation vector consists of 7 elements: the first 4 elements form a Quaternion and describe the rotation
	 * the last three elements are the translation vector
	 * @param mode
	 * @return
	 */
	public static double[][] initialTransform(int mode) {
		double[][] transforms1 = {
				{0.00, 0.00, 0.00,0.0,0.0,0.0},
				{3.14, 0.00, 0.00,0.0,0.0,0.0},
				{0.00, 3.14, 0.00,0.0,0.0,0.0},
				{0.00, 0.00, 3.14,0.0,0.0,0.0}
		};
		double[][] transforms2 = {
				{0.00, 0.00, 0.00,0.0,0.0,0.0},
				{3.14, 0.00, 0.00,0.0,0.0,0.0},
				{0.00, 3.14, 0.00,0.0,0.0,0.0},
				{0.00, 0.00, 3.14,0.0,0.0,0.0},
				{1.57, 0.00, 0.00,0.0,0.0,0.0},
				{0.00, 1.57, 0.00,0.0,0.0,0.0},
				{0.00, 0.00, 1.57,0.0,0.0,0.0},
				{-1.21, -1.21, 1.21,0.0,0.0,0.0},
				{-1.21, 1.21, -1.21,0.0,0.0,0.0},
				{1.21, 1.21, -1.21,0.0,0.0,0.0},
				{-1.21, -1.21, -1.21,0.0,0.0,0.0},
				{1.21, -1.21, -1.21,0.0,0.0,0.0}
		};
		
	
		switch(mode){
		case 1:
			return transforms1;
		case 2:
			return transforms2;
		
			
		default:
		
			double [][] transform = {{0.00, 0.00, 0.00,0.0,0.0,0.0}};
			return transform;
		}
	
			
	}
		
	
	public double getSelfAtomOverlapRef(){
		
		return refMolGauss.getSelfAtomOverlap();
	}
	
	public double getSelfAtomOverlapFit(){
		
		return molGauss.getSelfAtomOverlap();
	}
	
	public double getSelfPPOverlapRef(){
		
		return refMolGauss.getSelfPPOverlap();
	}
	
	public double getSelfPPOverlapFit(){
		
		return molGauss.getSelfPPOverlap();
	}
		

	
	
	public static void rotateMol(Conformer conf,Quaternion rotor, double[] transl) {

		double normFactor = 1/rotor.normSquared();

		
		int nrOfAtoms = conf.getSize();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords = conf.getCoordinates(i);
			double[][] m = rotor.getRotMatrix().getArray();
			coords.rotate(m);
			
			coords.scale(normFactor);
			coords.add(transl[0],transl[1],transl[2]);

		}
		
	}
	
	public static void rotateMol(StereoMolecule mol,Quaternion rotor, double[] transl) {

		double normFactor = 1/rotor.normSquared();

		
		int nrOfAtoms = mol.getAllAtoms();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords = mol.getCoordinates(i);
			double[][] m = rotor.getRotMatrix().getArray();
			coords.rotate(m);	
			coords.scale(normFactor);
			coords.add(transl[0],transl[1],transl[2]);

		}
		
	}
	
	public static void rotateMol(StereoMolecule mol,double[][] m) {
		int nrOfAtoms = mol.getAllAtoms();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords = mol.getCoordinates(i);
			coords.rotate(m);
		}
		
	}
	
	public static void rotateMol(Conformer conf,double[][] m) {
		int nrOfAtoms = conf.getMolecule().getAllAtoms();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords = conf.getCoordinates(i);
			coords.rotate(m);
		}
		
	}
	

	
	public static void translateMol(StereoMolecule mol,double[] translate) {
		int nrOfAtoms = mol.getAllAtoms();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords = mol.getCoordinates(i);
			coords.x += translate[0];
			coords.y += translate[1];
			coords.z += translate[2];
		}
		
	}
	
	public static void multiplyMatrix(double[][] r, double[][] s, double[][] rs) {
		rs[0][0] = r[0][0]*s[0][0] + r[0][1]*s[1][0] + r[0][2]*s[2][0];
		rs[0][1] = r[0][0]*s[0][1] + r[0][1]*s[1][1] + r[0][2]*s[2][1];
		rs[0][2] = r[0][0]*s[0][2] + r[0][1]*s[1][2] + r[0][2]*s[2][2];
		
		rs[1][0] = r[1][0]*s[0][0] + r[1][1]*s[1][0] + r[1][2]*s[2][0];
		rs[1][1] = r[1][0]*s[0][1] + r[1][1]*s[1][1] + r[1][2]*s[2][1];
		rs[1][2] = r[1][0]*s[0][2] + r[1][1]*s[1][2] + r[1][2]*s[2][2];
		
		rs[2][0] = r[2][0]*s[0][0] + r[2][1]*s[1][0] + r[2][2]*s[2][0];
		rs[2][1] = r[2][0]*s[0][1] + r[2][1]*s[1][1] + r[2][2]*s[2][1];
		rs[2][2] = r[2][0]*s[0][2] + r[2][1]*s[1][2] + r[2][2]*s[2][2];
	}
	
	public static void multiplyInverseMatrix(double[][] r, double[][] s, double[][] rs) {
		rs[0][0] = r[0][0]*s[0][0] + r[0][1]*s[0][1] + r[0][2]*s[0][2];
		rs[0][1] = r[0][0]*s[1][0] + r[0][1]*s[1][1] + r[0][2]*s[1][2];
		rs[0][2] = r[0][0]*s[2][0] + r[0][1]*s[2][1] + r[0][2]*s[2][2];
		
		rs[1][0] = r[1][0]*s[0][0] + r[1][1]*s[0][1] + r[1][2]*s[0][2];
		rs[1][1] = r[1][0]*s[1][0] + r[1][1]*s[1][1] + r[1][2]*s[1][2];
		rs[1][2] = r[1][0]*s[2][0] + r[1][1]*s[2][1] + r[1][2]*s[2][2];
		
		rs[2][0] = r[2][0]*s[0][0] + r[2][1]*s[0][1] + r[2][2]*s[0][2];
		rs[2][1] = r[2][0]*s[1][0] + r[2][1]*s[1][1] + r[2][2]*s[1][2];
		rs[2][2] = r[2][0]*s[2][0] + r[2][1]*s[2][1] + r[2][2]*s[2][2];
	}
	
	public static void getRotationMatrix(double theta, Coordinates axis, double[][] r) {
		double x = axis.x;
		double y = axis.y;
		double z = axis.z;
		double c = Math.cos(theta);
		double s = Math.sin(theta);
		double t = 1-c;
		r[0][0] = c+x*x*t;
		r[1][0] = x*y*t-z*s;
		r[2][0] = x*z*t+y*s;
		r[0][1] = x*y*t+z*s;
		r[1][1] = c+y*y*t;
		r[2][1] = y*z*t-x*s;
		r[0][2] = z*x*t-y*s;
		r[1][2] = z*y*t+x*s;
		r[2][2] = c+z*z*t;

	}
	
	
	public double[] findAlignment(double[][] initialTransforms, TransformationSequence optimizedTransform) {
		return findAlignment(initialTransforms,optimizedTransform,true);
	}
	
	public double[] findAlignment(double[][] initialTransforms, TransformationSequence optimizedTransform, boolean optimize) {
		return findAlignment(initialTransforms,optimizedTransform,optimize,SimilarityMode.TANIMOTO);
	}
	
	public double[] findAlignment(double[][] initialTransforms, TransformationSequence optimizedTransform,boolean optimize, SimilarityMode simMode) {
		boolean tversky = true;
		if(simMode==SimilarityMode.TANIMOTO)
			tversky=false;
		double tverskyCoeff = simMode==SimilarityMode.TVERSKY ? TVERSKY_COEFFICIENT : 1.0-TVERSKY_COEFFICIENT;
		double Oaa = getSelfAtomOverlapRef();
		double Obb = getSelfAtomOverlapFit();
		double ppOaa = getSelfPPOverlapRef();
		double ppObb = getSelfPPOverlapFit();
		EvaluableOverlap eval = new EvaluableOverlap(this, new double[6],ppWeight);
		OptimizerLBFGS opt = new OptimizerLBFGS(200,0.001);
		double maxSimilarity = 0.0;
		double maxPPSimilarity = 0.0;
		double maxVolSimilarity = 0.0;
		double maxShapeSimilarity = 0.0;
		double[] bestTransform = new double[6];
		for(double [] transform:initialTransforms) { 
			double ppSimilarity = 0.0;//iterate over all initial alignments (necessary since optimizer just finds next local minimum, so we need different initial guesses
			double atomSimilarity = 0.0;
			double volSimilarity = 0.0;
			double[] currentTransform;
			eval.setState(transform);
			if(optimize) {
				currentTransform = opt.optimize(eval);
			}
			else
				currentTransform = transform;
			double atomOverlap = 0.0;
			double ppOverlap = 0.0;
			double similarity = 0.0;
			ppOverlap = refMolGauss.getTotalPPOverlap(currentTransform,molGauss);
			if(getRefMolGauss().getPPGaussians().size()==0 && getMolGauss().getPPGaussians().size()==0 )
				ppSimilarity = 1.0;
			else {
				if(tversky)
					ppSimilarity = ppOverlap/(tverskyCoeff*ppObb+(1.0-tverskyCoeff)*ppOaa);
				else
					ppSimilarity=(ppOverlap/(ppOaa+ppObb-ppOverlap));
			}
			double correctionFactor = refMolGauss.getPPGaussians().size()/refMolGauss.getPPGaussians().stream().mapToDouble(g -> g.getWeight()).sum();
			ppSimilarity*=correctionFactor;
			if(ppSimilarity>1.0) //can happen because of weights
				ppSimilarity = 1.0f;
			double[] result = refMolGauss.getTotalAtomOverlap(currentTransform,molGauss);
			atomOverlap = result[0];
			double additionalVolOverlap = result[1];
			if(tversky)
				atomSimilarity = atomOverlap/(tverskyCoeff*Obb+(1.0-tverskyCoeff)*Oaa);
			else
				atomSimilarity = atomOverlap/(Oaa+Obb-atomOverlap);
			if(!tversky && atomSimilarity>1.0) //can happen because of weights
				atomSimilarity = 1.0f;
			volSimilarity = (additionalVolOverlap/atomOverlap);
			similarity = (1.0-ppWeight)*atomSimilarity + ppWeight*ppSimilarity;
			if (similarity>maxSimilarity) {
				maxSimilarity = similarity;
				maxVolSimilarity = volSimilarity;
				maxShapeSimilarity = atomSimilarity;
				maxPPSimilarity = ppSimilarity;
				bestTransform = currentTransform;
			}
		}
		ExponentialMap eMap = new ExponentialMap(bestTransform[0],bestTransform[1],bestTransform[2]);
		Quaternion rotor = eMap.toQuaternion();
		Translation translate = new Translation(new double[] {bestTransform[3],bestTransform[4],bestTransform[5]});
		TransformationSequence transformation = new TransformationSequence(rotor);
		transformation.addTransformation(translate);
		for(Transformation trans : transformation.getTransformations())
			optimizedTransform.addTransformation(trans);
		if(!tversky && maxSimilarity>1.0) // can happen because of manually placed inclusion spheres
			maxSimilarity = 1.0;
		return Arrays.stream(new double[] {maxSimilarity,maxPPSimilarity,maxShapeSimilarity,maxVolSimilarity}).toArray();
		}
		
		
	/*
	public static void rotateMol(Conformer conf, Matrix rotMat) {
		int nrOfAtoms = conf.getSize();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords1 = conf.getCoordinates(i);
			coords1.rotate(rotMat.getArray());
		}

	}
	

	
	
	
	public static void rotateMol(StereoMolecule mol, double[] transform) {
		Quaternion rotor = new Quaternion(transform[0],transform[1], transform [2], transform[3]);
		double[] translate =  {transform[4], transform[5], transform[6]};
		PheSAAlignment.rotateMol(mol, rotor, translate);

	}
	
	public static void rotateMol(StereoMolecule mol, Matrix rotMat) {
		int nrOfAtoms = mol.getAllAtoms();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords1 = mol.getCoordinates(i);
			coords1.rotate(rotMat.getArray());
		}

	}
	
	
	
	public static void rotateMol(Conformer conf, double[] transform) {
		Quaternion rotor = new Quaternion(transform[0],transform[1], transform [2], transform[3]);
		double[] translate =  {transform[4], transform[5], transform[6]};
		PheSAAlignment.rotateMol(conf, rotor, translate);

	}
	*/
	
	public static class PheSAResult implements Comparable <PheSAResult>{
		private StereoMolecule refMol;
		private StereoMolecule fitMol;
		private StereoMolecule fitInput;
		
		private double sim;
		private double[] contributions;
		private static final String DELIMITER = ";";
		
		public PheSAResult(StereoMolecule refMol, StereoMolecule fitInput, StereoMolecule fitMol, double sim) {
			this.refMol = refMol;
			this.fitMol = fitMol;
			this.sim = sim;
			this.contributions = new double[4];
			this.fitInput = fitInput;
		}
		
		public void setFitInput(StereoMolecule fitInput) {
			this.fitInput = fitInput;
		}


		public StereoMolecule getRefMol() {
			return refMol;
		}
		
		public StereoMolecule getFitMol() {
			return fitMol;
		}
		
		public double getSim() {
			return sim;
		}
		
		public void setContributions(double[] contributions) {
			this.contributions = contributions;
		}
		
		public double[] getContributions() {
			return contributions;
		}
		
		public String encode() {
			Encoder encoder = Base64.getEncoder();
			StringBuilder sb = new StringBuilder();
			Canonizer can = new Canonizer(refMol, Canonizer.COORDS_ARE_3D);
			String idcoords = can.getEncodedCoordinates(true);
			String idcode = can.getIDCode();
			sb.append(idcode);
			sb.append(DELIMITER);
			sb.append(idcoords);
			sb.append(DELIMITER);
			Canonizer can2 = new Canonizer(fitMol, Canonizer.COORDS_ARE_3D);
			String idcoords2 = can2.getEncodedCoordinates(true);
			String idcode2 = can2.getIDCode();
			sb.append(idcode2);
			sb.append(DELIMITER);
			sb.append(idcoords2);
			sb.append(DELIMITER);
			sb.append(encoder.encodeToString(EncodeFunctions.doubleToByteArray(sim)));
			sb.append(DELIMITER);
			sb.append(encoder.encodeToString(EncodeFunctions.doubleArrayToByteArray(contributions)));
			sb.append(DELIMITER);
			sb.append(fitInput.getIDCode());
			
			return sb.toString();
		}
		
		public static PheSAResult decode(String resultString) {
			Decoder decoder = Base64.getDecoder();
			String[] s = resultString.split(DELIMITER);
			String idcode = s[0];
			String idcoords = s[1];
			StereoMolecule refMol = new StereoMolecule();
			IDCodeParserWithoutCoordinateInvention parser = new IDCodeParserWithoutCoordinateInvention();
			parser.parse(refMol, idcode, idcoords);
			refMol.ensureHelperArrays(Molecule.cHelperCIP);
			idcode = s[2];
			idcoords = s[3];
			StereoMolecule fitMol = new StereoMolecule();
			parser = new IDCodeParserWithoutCoordinateInvention();
			parser.parse(fitMol, idcode, idcoords);
			fitMol.ensureHelperArrays(Molecule.cHelperCIP);
			double sim = EncodeFunctions.byteArrayToDouble(decoder.decode(s[4].getBytes()));
			double[] contributions = EncodeFunctions.byteArrayToDoubleArray(decoder.decode(s[5].getBytes()));
			StereoMolecule fitInput = new StereoMolecule();
			new IDCodeParser().parse(fitInput,s[6]);
			PheSAResult pheSAResult = new PheSAResult(refMol,fitInput,fitMol,sim);
			pheSAResult.setContributions(contributions);
			return pheSAResult;
		}

		@Override
		public int compareTo(PheSAResult o) {
				return Double.compare(sim, o.sim);
		}
	}
	
	
	

	
}

	
	
	

	

	
	
	
	


