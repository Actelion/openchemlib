package com.actelion.research.chem.alignment3d;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.transformation.Rotation;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.PheSAAlignment;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.actelion.research.chem.phesa.ShapeVolume;
import com.actelion.research.chem.phesa.pharmacophore.PPTriangle;
import com.actelion.research.chem.phesa.pharmacophore.PPTriangleCreator;
import com.actelion.research.chem.phesa.pharmacophore.PPTriangleMatcher;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;
import com.actelion.research.chem.phesa.pharmacophore.pp.PPGaussian;

import java.util.*;
import java.util.stream.Collectors;

public class PheSAAlignmentOptimizer {
	
	public static final int TRIANGLE_OPTIMIZATIONS = 50;
	public static final int PMI_OPTIMIZATIONS = 10;
	private static final double EXIT_VECTOR_WEIGHT = 10.0;
	private static final int BEST_RESULT_SIZE = 20;
	
	public enum SimilarityMode {REFTVERSKY,TVERSKY, TANIMOTO
		}
	
	private PheSAAlignmentOptimizer() {}
	
	
	public static double alignTwoMolsInPlace(StereoMolecule refMol, StereoMolecule fitMol, ThreadMaster tm) {
		return alignTwoMolsInPlace(refMol, fitMol, 0.5, tm);
	}
	
	public static double alignTwoMolsInPlace(StereoMolecule refMol, StereoMolecule fitMol, double ppWeight, ThreadMaster tm) {
		PheSASetting setting = new PheSASetting();
		setting.setPpWeight(ppWeight);
		double similarity = 0.0;
		MolecularVolume refVol = new MolecularVolume(refMol);
		MolecularVolume fitVol = new MolecularVolume(fitMol);
		Coordinates origCOM = new Coordinates(refVol.getCOM());
		Conformer refConf = new Conformer(refMol);
		Conformer fitConf = new Conformer(fitMol);
		Rotation rotation = refVol.preProcess(refConf);
		rotation = rotation.getInvert();
		fitVol.preProcess(fitConf);
		AlignmentResult bestSolution = createAlignmentSolutions(Collections.singletonList(refVol), Collections.singletonList(fitVol),setting, tm).get(0);
		similarity = bestSolution.getSimilarity();
		
		for(int a=0;a<fitMol.getAllAtoms();a++) {
			fitMol.setAtomX(a, fitConf.getX(a));
			fitMol.setAtomY(a, fitConf.getY(a));
			fitMol.setAtomZ(a, fitConf.getZ(a));
		}
		bestSolution.getTransform().apply(fitMol);
		rotation.apply(fitMol);
		fitMol.translate(origCOM.x, origCOM.y, origCOM.z);

		return similarity;
		
		
		
	}
	
	public static List<AlignmentResult> alignToNegRecImg(ShapeVolume ref,
														 List<? extends ShapeVolume> fitVols,
														 PheSASetting setting, ThreadMaster tm) {
		for(ShapeVolume shapeVol : fitVols) {
			shapeVol.removeRings();
		}
		List<AlignmentResult> alignmentSolutions = createAlignmentSolutions(Collections.singletonList(ref), fitVols,setting, tm);
		List<AlignmentResult> results = new ArrayList<>();
		int counter = 0;
		for(AlignmentResult solution : alignmentSolutions) {
			if(counter++>=BEST_RESULT_SIZE) {
				break;
			}
			results.add(solution);
		
		}
		return results;
	}

	
	public static List<AlignmentResult> createAlignmentSolutions(List<? extends ShapeVolume> refVols,
																 List<? extends ShapeVolume> fitVols,
																 PheSASetting setting, ThreadMaster tm) {
		int pmiOpti = setting.nrOptimizationsPMI;
		for(ShapeVolume molVol : refVols) {
			for(PPGaussian ppg : molVol.getPPGaussians()) {
				if(ppg.getPharmacophorePoint().getFunctionalityIndex()==PharmacophoreCalculator.EXIT_VECTOR_ID) {
					ppg.setWeight(EXIT_VECTOR_WEIGHT);
				}
			}
		}
	
		List<AlignmentResult> triangleSolutions = new ArrayList<AlignmentResult>();
		if(setting.useTriangle)
			triangleSolutions = getBestTriangleAlignments(refVols,fitVols,setting);
		List<AlignmentResult> alignmentSolutions = new ArrayList<>(triangleSolutions);
		List<AlignmentResult> pmiSolutions = new ArrayList<AlignmentResult>();
		for(int i=0;i<refVols.size() && (tm==null || !tm.threadMustDie());i++) {
			ShapeVolume refVol = refVols.get(i);
			for(int j=0;j<fitVols.size() && (tm==null || !tm.threadMustDie());j++) {
				ShapeVolume fitVol = new ShapeVolume(fitVols.get(j));
				PheSAAlignment shapeAlignment = new PheSAAlignment(refVol,fitVol, setting.ppWeight);
				TransformationSequence pmiTransformation = new TransformationSequence();
				double[][] transforms = PheSAAlignment.initialTransform(2);
				double[] r = shapeAlignment.findAlignment(transforms,pmiTransformation,false,setting.getSimMode());
				AlignmentResult pmiAlignment = new AlignmentResult(r[0],pmiTransformation,i,j);
				pmiAlignment.setSimilarityContributions(r);
				pmiSolutions.add(pmiAlignment);
			}
		}
		//optimize best PMI alignments
		List<AlignmentResult> sortedPMISolutions = pmiSolutions.stream()
				.sorted(Comparator.reverseOrder())
				.collect(Collectors.toList());
		int counter = 0;
		for(AlignmentResult pmiAlignment : sortedPMISolutions) {
			if(counter++>pmiOpti || (tm!=null && tm.threadMustDie()))
				break;
			ShapeVolume refVol = refVols.get(pmiAlignment.getRefConformerIndex());
			ShapeVolume fitVol = new ShapeVolume(fitVols.get(pmiAlignment.getConformerIndex()));
			TransformationSequence bestTransform = pmiAlignment.getTransform();
			fitVol.transform(bestTransform);
			PheSAAlignment shapeAlignment = new PheSAAlignment(refVol,fitVol,setting.ppWeight);
			TransformationSequence optimizedTransform = new TransformationSequence();
			double[] r = shapeAlignment.findAlignment(PheSAAlignment.initialTransform(0),optimizedTransform,true,setting.simMode);
			pmiAlignment.getTransform().addTransformation(optimizedTransform);
			AlignmentResult optimizedPMIAlignment = new AlignmentResult(r[0],pmiAlignment.getTransform(),pmiAlignment.getRefConformerIndex(),pmiAlignment.getConformerIndex());
			optimizedPMIAlignment.setSimilarityContributions(r);
			alignmentSolutions.add(optimizedPMIAlignment);
			
		}
		return alignmentSolutions.stream()
				.sorted(Comparator.reverseOrder())
				.collect(Collectors.toList());
		
	}
	
	public static double[] align(PheSAMolecule refShape, PheSAMolecule fitShape, StereoMolecule[] bestAlignment, PheSASetting setting, ThreadMaster tm) {
		double[] result = new double[4]; //overall sim, ppSimilarity and additional volume similarity contribution
		
		List<AlignmentResult> alignmentSolutions = createAlignmentSolutions(refShape.getVolumes(),fitShape.getVolumes(),setting, tm);
		AlignmentResult bestResult = alignmentSolutions.get(0);
		StereoMolecule refMol = refShape.getConformer(bestResult.getRefConformerIndex());
		StereoMolecule fitMol = fitShape.getConformer(bestResult.getConformerIndex());
		bestResult.getTransform().apply(fitMol);
		result = bestResult.getSimilarityContributions();
		bestAlignment[0] = refMol;
		bestAlignment[1] = fitMol;
		int n1 = refShape.getVolumes().get(0).getExitVectorGaussians().size();
		int n2 = fitShape.getVolumes().get(0).getExitVectorGaussians().size();
		if(!refShape.getVolumes().get(0).getExitVectorGaussians().isEmpty() ||
				!refShape.getVolumes().get(0).getExitVectorGaussians().isEmpty()) {
			// there are exit vectors 
			
			if(n1!=n2) { //different number of exit vectors --> no match
				result[0] = 0.0;
				result[1] = 0.0;
				result[2] = 0.0;
				result[3] = 0.0;
			}
			
		}
	
		return result;
	}
	

	
	private static List<AlignmentResult> getBestTriangleAlignments(List<? extends ShapeVolume> refVols, List<? extends ShapeVolume> fitVols, PheSASetting setting) {
		List<AlignmentResult> triangleResults = new ArrayList<AlignmentResult>();
		for(int i=0;i<refVols.size();i++) {
			ShapeVolume refVol = refVols.get(i);
			Map<Integer,ArrayList<PPTriangle>> refTriangles = PPTriangleCreator.create(refVol.getPPGaussians(), refVol.getCOM());
			for(int j=0;j<fitVols.size();j++) {
				ShapeVolume fitVol = fitVols.get(j);
				Map<Integer,ArrayList<PPTriangle>> fitTriangles = PPTriangleCreator.create(fitVol.getPPGaussians(),fitVol.getCOM());
				triangleResults.addAll(PPTriangleMatcher.getMatchingTransforms(refTriangles, fitTriangles,i,j,setting.useDirectionality));
			}
		}

		List<AlignmentResult> sortedTriangleResults = triangleResults.stream()
				.sorted(Comparator.reverseOrder())
				.collect(Collectors.toList());
		List<AlignmentResult> optimizedResults = new ArrayList<AlignmentResult>();
		if(!triangleResults.isEmpty()) { // found triangle alignments

			double[][] alignments = PheSAAlignment.initialTransform(0);
			int counter = 0;
			for(AlignmentResult result: sortedTriangleResults) {
				if(counter++>setting.nrOptimizationsTriangle)
					break;
				ShapeVolume refVol = refVols.get(result.getRefConformerIndex());
				ShapeVolume fitVol = new ShapeVolume(fitVols.get(result.getConformerIndex()));
				TransformationSequence bestTransform = result.getTransform();
				fitVol.transform(bestTransform);
				PheSAAlignment shapeAlignment = new PheSAAlignment(refVol,fitVol,setting.ppWeight);
				TransformationSequence optTransform = new TransformationSequence();
				double[] r = shapeAlignment.findAlignment(alignments,optTransform,true, setting.simMode);
				bestTransform.addTransformation(optTransform);
				AlignmentResult optimizedResult = new AlignmentResult(r[0], bestTransform, result.getRefConformerIndex(), result.getConformerIndex());
				optimizedResult.setSimilarityContributions(r);
				optimizedResults.add(optimizedResult);
			}
		}

		return optimizedResults.stream()
				.sorted(Comparator.reverseOrder())
				.collect(Collectors.toList());
	}
	
	private static double[][] createSubAlignments(ShapeVolume refVol, double[][] baseTransforms) {
		final long seed = 12345L;
		final int maxPoints = 10;
		final int points = Math.min(refVol.getAtomicGaussians().size(), maxPoints);
		Random rnd = new Random(seed);
		List<double[]> transforms = new ArrayList<>(); 
		for(int i=0;i<points;i++) {
			int index = rnd.nextInt(refVol.getAtomicGaussians().size());
			Coordinates c = refVol.getAtomicGaussians().get(index).getCenter();
			for(double[] t : baseTransforms) {
				transforms.add(new double[] {t[0],t[1],t[2],c.x,c.y,c.z});
			}
		}
		transforms.addAll(Arrays.asList(baseTransforms));
		return transforms.toArray(new double[0][]);
	}
	

	
	public static class AlignmentResult implements Comparable<AlignmentResult>{
		private final double similarity;
		private final TransformationSequence transformation;
		private final int conformerIndex;
		private final int refConformerIndex;
		private double[] similarityContributions;
		
		public AlignmentResult(double similarity, TransformationSequence transformation, int refConformerIndex, int conformerIndex) {
			this.similarity = similarity;
			this.transformation = transformation;
			this.refConformerIndex = refConformerIndex;
			this.conformerIndex = conformerIndex;
		}
		
		
		public double[] getSimilarityContributions() {
			return similarityContributions;
		}


		public void setSimilarityContributions(double[] similarityContributions) {
			this.similarityContributions = similarityContributions;
		}


		public TransformationSequence getTransform() {
			return transformation;
		}

		
		public double getSimilarity() {
			return similarity;
		}
		
		public int getConformerIndex() {
			return conformerIndex;
		}
		
		public int getRefConformerIndex() {
			return refConformerIndex;
		}


		@Override
		public int compareTo(AlignmentResult o) {
			return Double.compare(similarity, o.similarity);
		}
		
		
		
	}

	/**
	 * @author wahljo1
	 *
	 */
	public static class PheSASetting {
		public double getPpWeight() {
			return ppWeight;
		}
		public void setPpWeight(double ppWeight) {
			this.ppWeight = ppWeight;
		}
		public SimilarityMode getSimMode() {
			return simMode;
		}
		public void setSimMode(SimilarityMode simMode) {
			this.simMode = simMode;
		}
		public boolean isUseDirectionality() {
			return useDirectionality;
		}
		public void setUseDirectionality(boolean useDirectionality) {
			this.useDirectionality = useDirectionality;
		}

		public int getNrOptimizationsPMI() {
			return nrOptimizationsPMI;
		}
		public void setNrOptimizationsPMI(int nrOptimizationsPMI) {
			this.nrOptimizationsPMI = nrOptimizationsPMI;
		}
		public int getNrOptimizationsTriangle() {
			return nrOptimizationsTriangle;
		}
		public void setNrOptimizationsTriangle(int nrOptimizationsTriangle) {
			this.nrOptimizationsTriangle = nrOptimizationsTriangle;
		}

		private double ppWeight;
		private SimilarityMode simMode; 
		private boolean useDirectionality;
		private boolean useTriangle;
		public boolean isUseTriangle() {
			return useTriangle;
		}
		public void setUseTriangle(boolean useTriangle) {
			this.useTriangle = useTriangle;
		}

		private int nrOptimizationsPMI;
		private int nrOptimizationsTriangle;
		private static final String DELIMITER1 ="_";
		
		public PheSASetting() {
			this.ppWeight = 0.5;
			this.simMode = SimilarityMode.TANIMOTO;
			this.useDirectionality = true;
			this.useTriangle = true;
			this.nrOptimizationsPMI = PheSAAlignmentOptimizer.PMI_OPTIMIZATIONS;
			this.nrOptimizationsTriangle = PheSAAlignmentOptimizer.TRIANGLE_OPTIMIZATIONS;
		}
		
		public String toString() {
			return Double.toString(ppWeight) +
					DELIMITER1 +
					simMode.toString() +
					DELIMITER1 +
					Boolean.toString(useDirectionality) +
					DELIMITER1 +
					Boolean.toString(useTriangle) +
					DELIMITER1 +
					Integer.toString(nrOptimizationsPMI) +
					DELIMITER1 +
					Integer.toString(nrOptimizationsTriangle);
		}
		
		public static PheSASetting fromString(String encoded) {
			PheSASetting setting = new PheSASetting();
			String[] split = encoded.split(DELIMITER1);
			setting.setPpWeight(Double.parseDouble(split[0]));
			setting.setSimMode(SimilarityMode.valueOf(split[1]));
			setting.setUseDirectionality(Boolean.parseBoolean(split[2]));
			setting.setUseTriangle(Boolean.parseBoolean(split[3]));
			setting.setNrOptimizationsPMI(Integer.parseInt(split[4]));
			setting.setNrOptimizationsTriangle(Integer.parseInt(split[5]));
			return setting;
		}
	}
}
