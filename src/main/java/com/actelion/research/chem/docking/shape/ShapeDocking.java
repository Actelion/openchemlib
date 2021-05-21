package com.actelion.research.chem.docking.shape;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer.AlignmentResult;
import com.actelion.research.chem.alignment3d.transformation.Transformation;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.DescriptorHandlerShapeOneConf;
import com.actelion.research.chem.phesa.PheSAAlignment;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.actelion.research.chem.phesa.ShapeVolume;

public class ShapeDocking {
	private static final double DEFAULT_PP_WEIGHT = 0.5;
	private Transformation transformation;
	private ShapeVolume negRecImage;
	private DescriptorHandlerShape dhs;
	
	/**
	 * 
	 * @param negRecImage
	 * @param transformation: transformation used during creation of ShapeVolume of NegRecImage (translate COM to (0,0,0)
	 */
	
	public ShapeDocking(ShapeVolume negRecImage, Transformation transformation) {
		this.negRecImage = negRecImage;
		this.transformation = transformation;
		dhs = new DescriptorHandlerShape(500,0.5);
	}
	
	public List<StereoMolecule> dock(StereoMolecule candidate) {
		List<StereoMolecule> aligned = new ArrayList<>();
		PheSAMolecule candidateShape = dhs.createDescriptor(candidate);
		List<AlignmentResult> results = PheSAAlignmentOptimizer.alignToNegRecImg(negRecImage, candidateShape.getVolumes(), DEFAULT_PP_WEIGHT, true);
		for(AlignmentResult result: results) {
			if(result.getSimilarity()==0.0) {
				return aligned;
			}
			StereoMolecule alignedMol = candidateShape.getConformer(candidateShape.getVolumes().get(result.getConformerIndex()));
			result.getTransform().apply(alignedMol);
			this.transformation.apply(alignedMol);
			aligned.add(alignedMol);
		}
		return aligned;

		/*
		StereoMolecule aligned = candidateShape.getConformer(candidateShape.getVolumes().get(bestIndeces[1]));
		PheSAAlignment.rotateMol(aligned, bestTransform);
		for(int a=0;a<aligned.getAllAtoms();a++) {
			transformation.applyRotation(aligned.getCoordinates(a));
			transformation.applyTranslation(aligned.getCoordinates(a));
		}
		System.out.println(Arrays.toString(bestTransform));
		System.out.println(result[0]);
		return aligned;
		*/
	}
}
