package com.actelion.research.chem.phesa;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.phesa.pharmacophore.PPTriangle;
import com.actelion.research.chem.phesa.pharmacophore.PPTriangleCreator;
import com.actelion.research.chem.phesa.pharmacophore.PPTriangleMatcher;
import com.actelion.research.chem.phesa.pharmacophore.PPTriangleMatcher.AlignmentResult;

public class PheSAAlignmentOptimizer {
	
	private static int OPTIMIZATIONS = 20;
	
	private PheSAAlignmentOptimizer() {}
	
	
	public static double alignTwoMolsInPlace(StereoMolecule refMol, StereoMolecule fitMol) {
		double similarity = 0.0;
		MolecularVolume refVol = new MolecularVolume(refMol);
		MolecularVolume fitVol = new MolecularVolume(fitMol);
		Coordinates origCOM = new Coordinates(refVol.getCOM());
		Conformer refConf = new Conformer(refMol);
		Conformer fitConf = new Conformer(fitMol);
		Matrix rotation = PheSAAlignment.preProcess(refConf, refVol);
		rotation = rotation.getTranspose();
		PheSAAlignment.preProcess(fitConf, fitVol);
		
		Map<Integer,ArrayList<PPTriangle>> refTriangles = PPTriangleCreator.create(refVol.getPPGaussians(), refVol.getCOM());
		Map<Integer,ArrayList<PPTriangle>> fitTriangles = PPTriangleCreator.create(fitVol.getPPGaussians(), fitVol.getCOM());
		List<AlignmentResult> results = PPTriangleMatcher.getMatchingTransforms(refTriangles, fitTriangles,0,0);
		double[] bestTransformTriangle = new double[7];
		double[] bestTransformPMI = new double[7];
		
		StereoMolecule bestMatch = fitMol;
		
		double bestScoreTriangle = 0.0;
		if(results.size()!=0) { // found triangle alignments
			results.sort((c1,c2) -> {
				return Double.compare(c1.getSimilarity(),c2.getSimilarity());
			});
		
			AlignmentResult[] bestAlignments = new AlignmentResult[Math.min(results.size(),OPTIMIZATIONS)];
			for(int i=1;i<bestAlignments.length+1;i++) 
				bestAlignments[i-1] = results.get(results.size()-i);
			double[][] rotate = new double[3][3];
			double[][] alignments = PheSAAlignment.initialTransform(0);
			MolecularVolume refVol2 = new MolecularVolume(refVol);
			for(AlignmentResult result: bestAlignments) {
				MolecularVolume fitVol2 = new MolecularVolume(fitVol);
				double[][] bestTransform = result.getTransform();
				rotate[0] = bestTransform[0];
				rotate[1] = bestTransform[1];
				rotate[2] = bestTransform[2];
				StereoMolecule fitMol2 = fitConf.toMolecule(null);
				fitMol2.ensureHelperArrays(Molecule.cHelperNeighbours);
				PheSAAlignment.translateMol(fitMol2, bestTransform[3]);
				PheSAAlignment.rotateMol(fitMol2, rotate);
				PheSAAlignment.translateMol(fitMol2, bestTransform[4]);
				fitVol2.update(fitMol2);
				PheSAAlignment shapeAlignment = new PheSAAlignment(refVol2,fitVol2);
				double[] r = shapeAlignment.findAlignment(alignments);
				if(r[0]>bestScoreTriangle) {
					bestScoreTriangle = r[0];
					bestTransformTriangle[0] = r[1];
					bestTransformTriangle[1] = r[2];
					bestTransformTriangle[2] = r[3];
					bestTransformTriangle[3] = r[4];
					bestTransformTriangle[4] = r[5];
					bestTransformTriangle[5] = r[6];
					bestTransformTriangle[6] = r[7];
					bestMatch = fitMol2;
				}
			}
		}
		PheSAAlignment shapeAlignment = new PheSAAlignment(refVol,fitVol);
		double[] r = shapeAlignment.findAlignment(PheSAAlignment.initialTransform(2),true);
		if(r[0]>bestScoreTriangle) { //alignment found by PMI initial alignment is better
			similarity = r[0];
			bestTransformPMI = new double[] {r[1],r[2], r[3], r[4], r[5], r[6], r[7]};
			bestMatch = fitConf.toMolecule();
			PheSAAlignment.rotateMol(bestMatch, bestTransformPMI);
			for(int a=0;a<fitMol.getAllAtoms();a++) {
				fitMol.setAtomX(a, bestMatch.getAtomX(a));
				fitMol.setAtomY(a, bestMatch.getAtomY(a));
				fitMol.setAtomZ(a, bestMatch.getAtomZ(a));
			}
		}
		
		else { 
			similarity = bestScoreTriangle;
			PheSAAlignment.rotateMol(bestMatch, bestTransformTriangle);
			for(int a=0;a<fitMol.getAllAtoms();a++) {
				fitMol.setAtomX(a, bestMatch.getAtomX(a));
				fitMol.setAtomY(a, bestMatch.getAtomY(a));
				fitMol.setAtomZ(a, bestMatch.getAtomZ(a));
			}
		
		}
	PheSAAlignment.rotateMol(fitMol, rotation);
	fitMol.translate(origCOM.x, origCOM.y, origCOM.z);

	return similarity;
		
		
		
	}
	
	
	public static double align(PheSAMolecule refShape, PheSAMolecule fitShape, StereoMolecule[] bestAlignment) {
		StereoMolecule[] bestPairTriangle = new StereoMolecule[2];
		double[] bestTransformTriangle = new double[7];
		double bestScoreTriangle = getBestTriangleAlignment(refShape,fitShape,bestPairTriangle,bestTransformTriangle);
		double bestScorePMI = 0.0; 
		double[] bestTransformPMI = new double[7];
		MolecularVolume[] bestPairPMI = new MolecularVolume[2];
		for(int i=0;i<refShape.getVolumes().size();i++) {
			MolecularVolume refVol = new MolecularVolume(refShape.getVolumes().get(i));
			for(int j=0;j<fitShape.getVolumes().size();j++) {
				MolecularVolume fitVol = new MolecularVolume(fitShape.getVolumes().get(j));
				PheSAAlignment shapeAlignment = new PheSAAlignment(refVol,fitVol);
				double[] r = shapeAlignment.findAlignment(PheSAAlignment.initialTransform(1),false);
				if(r[0]>bestScorePMI) {
					bestScorePMI = r[0];
					bestTransformPMI = new double[] {r[1],r[2], r[3], r[4], r[5], r[6], r[7]};
					bestPairPMI[0] = refVol;
					bestPairPMI[1] = fitVol;
				}
				
			}
		}
		//optimize best PMI alignment
		double similarity = 0.0;
		PheSAAlignment shapeAlignment = new PheSAAlignment(bestPairPMI[0],bestPairPMI[1]);
		double[] r = shapeAlignment.findAlignment(new double[][] {bestTransformPMI},true);
		if(r[0]>bestScoreTriangle) { //alignment found by PMI initial alignment is better
			similarity = r[0];
			bestAlignment[0] = refShape.getConformer(bestPairPMI[0]);
			bestAlignment[1] = fitShape.getConformer(bestPairPMI[1]);
			PheSAAlignment.rotateMol(bestAlignment[1], bestTransformPMI);
		}
		
		else { 
			similarity = bestScoreTriangle;
			bestAlignment[0] = bestPairTriangle[0];
			bestAlignment[1] = bestPairTriangle[1];
		}
			
	return similarity;
	}
	
	public static double align(PheSAMolecule fitShape, MolecularVolume refVol, MolecularVolume fitVol, StereoMolecule aligned) {
		Coordinates[] origCoords = new Coordinates[aligned.getAllAtoms()];
		IntStream.range(0,aligned.getAllAtoms()).forEach(i -> origCoords[i] = new Coordinates(aligned.getCoordinates(i)));
		double[] bestTransformTriangle = new double[7];
		double bestScoreTriangle = getTriangleAlignment(fitShape,new MolecularVolume(refVol),new MolecularVolume(fitVol),bestTransformTriangle,aligned);
		double bestScorePMI = 0.0; 
		double[] bestTransformPMI = new double[7];
		PheSAAlignment shapeAlignment = new PheSAAlignment(new MolecularVolume(refVol),new MolecularVolume(fitVol));
		double[] r = shapeAlignment.findAlignment(PheSAAlignment.initialTransform(2),true);
		bestScorePMI = r[0];
		bestTransformPMI = new double[] {r[1],r[2], r[3], r[4], r[5], r[6], r[7]};

				
		double similarity = 0.0;

		if(bestScorePMI>bestScoreTriangle) { //alignment found by PMI initial alignment is better
			IntStream.range(0,aligned.getAllAtoms()).forEach(i -> {	
			aligned.setAtomX(i, origCoords[i].x);
			aligned.setAtomY(i, origCoords[i].y);
			aligned.setAtomZ(i, origCoords[i].z);
			});
			similarity = r[0];
			PheSAAlignment.rotateMol(aligned, bestTransformPMI);
		}
		
		else { 
			similarity = bestScoreTriangle;
		}

		return similarity;
	}
	
	private static double getBestTriangleAlignment(PheSAMolecule refShape, PheSAMolecule fitShape, StereoMolecule[] bestPairTriangle, double[] bestTransformTriangle) {
		List<AlignmentResult> results = new ArrayList<AlignmentResult>();
		for(int i=0;i<refShape.getVolumes().size();i++) {
			MolecularVolume refVol = refShape.getVolumes().get(i);
			Map<Integer,ArrayList<PPTriangle>> refTriangles = PPTriangleCreator.create(refVol.getPPGaussians(), refVol.getCOM());
			for(int j=0;j<fitShape.getVolumes().size();j++) {
				MolecularVolume fitVol = fitShape.getVolumes().get(j);
				Map<Integer,ArrayList<PPTriangle>> fitTriangles = PPTriangleCreator.create(fitVol.getPPGaussians(),fitVol.getCOM());
				results.addAll(PPTriangleMatcher.getMatchingTransforms(refTriangles, fitTriangles,i,j));
			}
		}
		double bestScoreTriangle = 0.0;
		if(results.size()!=0) { // found triangle alignments
			results.sort((c1,c2) -> {
				return Double.compare(c1.getSimilarity(),c2.getSimilarity());
			});
		
			AlignmentResult[] bestAlignments = new AlignmentResult[Math.min(results.size(),OPTIMIZATIONS)];
			for(int i=1;i<bestAlignments.length+1;i++) 
				bestAlignments[i-1] = results.get(results.size()-i);
			double[][] rotate = new double[3][3];
			double[][] alignments = PheSAAlignment.initialTransform(0);
			for(AlignmentResult result: bestAlignments) {
				MolecularVolume refVol = new MolecularVolume(refShape.getVolumes().get(result.getRefConformerIndex()));
				MolecularVolume fitVol = new MolecularVolume(fitShape.getVolumes().get(result.getConformerIndex()));
				double[][] bestTransform = result.getTransform();
				rotate[0] = bestTransform[0];
				rotate[1] = bestTransform[1];
				rotate[2] = bestTransform[2];
				StereoMolecule fitMol = fitShape.getConformer(fitVol);
				PheSAAlignment.translateMol(fitMol, bestTransform[3]);
				PheSAAlignment.rotateMol(fitMol, rotate);
				PheSAAlignment.translateMol(fitMol, bestTransform[4]);
				fitVol.update(fitMol);
				PheSAAlignment shapeAlignment = new PheSAAlignment(refVol,fitVol);
				double[] r = shapeAlignment.findAlignment(alignments);
				if(r[0]>bestScoreTriangle) {
					bestScoreTriangle = r[0];
					bestTransformTriangle[0] = r[1];
					bestTransformTriangle[1] = r[2];
					bestTransformTriangle[2] = r[3];
					bestTransformTriangle[3] = r[4];
					bestTransformTriangle[4] = r[5];
					bestTransformTriangle[5] = r[6];
					bestTransformTriangle[6] = r[7];
					bestPairTriangle[1] = fitShape.getConformer(fitVol);
					bestPairTriangle[0] = refShape.getConformer(refVol);
			}
			}
		}
		if(bestScoreTriangle>0.0)
			PheSAAlignment.rotateMol(bestPairTriangle[1], bestTransformTriangle);
	
		return bestScoreTriangle;
	}
	
	private static double getTriangleAlignment(PheSAMolecule fitShape,MolecularVolume refVol_, MolecularVolume fitVol_, double[] bestTransformTriangle, StereoMolecule aligned) {
		MolecularVolume refVol = new MolecularVolume(refVol_);
		MolecularVolume fitVol = new MolecularVolume(fitVol_);
		List<AlignmentResult> results = new ArrayList<AlignmentResult>();
		Map<Integer,ArrayList<PPTriangle>> refTriangles = PPTriangleCreator.create(refVol.getPPGaussians(), refVol.getCOM());
		Map<Integer,ArrayList<PPTriangle>> fitTriangles = PPTriangleCreator.create(fitVol.getPPGaussians(),fitVol.getCOM());
		results.addAll(PPTriangleMatcher.getMatchingTransforms(refTriangles, fitTriangles,0,0));
		double bestScoreTriangle = 0.0;
		if(results.size()!=0) { // found triangle alignments
			results.sort((c1,c2) -> {
				return Double.compare(c1.getSimilarity(),c2.getSimilarity());
			});
		
			AlignmentResult[] bestAlignments = new AlignmentResult[Math.min(results.size(),OPTIMIZATIONS)];
			for(int i=1;i<bestAlignments.length+1;i++) 
				bestAlignments[i-1] = results.get(results.size()-i);
			double[][] rotate = new double[3][3];
			double[][] alignments = PheSAAlignment.initialTransform(0);
			double[][] bestTransform = new double[4][4];
			for(AlignmentResult result: bestAlignments) {
				double[][] transform = result.getTransform();
				rotate[0] = transform[0];
				rotate[1] = transform[1];
				rotate[2] = transform[2];
				StereoMolecule fitMol = fitShape.getConformer(fitVol);
				PheSAAlignment.translateMol(fitMol, transform[3]);
				PheSAAlignment.rotateMol(fitMol, rotate);
				PheSAAlignment.translateMol(fitMol, transform[4]);
				fitVol.update(fitMol);
				PheSAAlignment shapeAlignment = new PheSAAlignment(refVol,fitVol);
				double[] r = shapeAlignment.findAlignment(alignments);
				if(r[0]>bestScoreTriangle) {
					bestTransform = transform;
					bestScoreTriangle = r[0];
					bestTransformTriangle[0] = r[1];
					bestTransformTriangle[1] = r[2];
					bestTransformTriangle[2] = r[3];
					bestTransformTriangle[3] = r[4];
					bestTransformTriangle[4] = r[5];
					bestTransformTriangle[5] = r[6];
					bestTransformTriangle[6] = r[7];

			}
			}
			if(bestScoreTriangle>0.0) {
				rotate[0] = bestTransform[0];
				rotate[1] = bestTransform[1];
				rotate[2] = bestTransform[2];
				PheSAAlignment.translateMol(aligned, bestTransform[3]);
				PheSAAlignment.rotateMol(aligned, rotate);
				PheSAAlignment.translateMol(aligned, bestTransform[4]);
				PheSAAlignment.rotateMol(aligned, bestTransformTriangle);
			}
		}
		return bestScoreTriangle;
	}

}
