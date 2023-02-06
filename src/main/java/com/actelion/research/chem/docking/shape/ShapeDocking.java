package com.actelion.research.chem.docking.shape;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer.AlignmentResult;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer.PheSASetting;
import com.actelion.research.chem.alignment3d.transformation.Transformation;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.PheSAMolecule;
import com.actelion.research.chem.phesa.ShapeVolume;

import java.util.ArrayList;
import java.util.List;

public class ShapeDocking {
	private static final double DEFAULT_PP_WEIGHT = 0.5;
	private Transformation transformation;
	private ShapeVolume negRecImage;
	private DescriptorHandlerShape dhs;
	private ThreadMaster threadMaster;
	private PheSASetting phesaSetting;
	
	/**
	 * 
	 * @param negRecImage
	 * @param transformation: transformation used during creation of ShapeVolume of NegRecImage (translate COM to (0,0,0)
	 */
	
	public ShapeDocking(ShapeVolume negRecImage, Transformation transformation) {
		this.negRecImage = negRecImage;
		this.transformation = transformation;
		dhs = new DescriptorHandlerShape(500,0.5);
		phesaSetting = new PheSASetting();
		phesaSetting.setNrOptimizationsPMI(2*PheSAAlignmentOptimizer.PMI_OPTIMIZATIONS);
		phesaSetting.setNrOptimizationsTriangle(2*PheSAAlignmentOptimizer.TRIANGLE_OPTIMIZATIONS);
		phesaSetting.setUseDirectionality(false);
		phesaSetting.setSimMode(PheSAAlignmentOptimizer.SimilarityMode.TVERSKY);
	}
	
	public PheSASetting getPhesaSetting() {
		return phesaSetting;
	}

	public void setThreadMaster(ThreadMaster tm) {
		// TODO use the ThreadMaster!!!
		this.threadMaster = tm;
	}
	
	public List<StereoMolecule> dock(StereoMolecule candidate) {
		List<StereoMolecule> aligned = new ArrayList<>();
		PheSAMolecule candidateShape = dhs.createDescriptor(candidate);
		List<AlignmentResult> results = PheSAAlignmentOptimizer.alignToNegRecImg(negRecImage, candidateShape.getVolumes(), phesaSetting);
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


	}
}
