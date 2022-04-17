package com.actelion.research.chem.phesa.pharmacophore;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer.AlignmentResult;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;

public class PPTriangleMatcher {
	
	private static final double CUTOFF = 2.0; //if lengths of edges of two triangles differ by more than 2.5A, it is
											  // not considered as a match
	
	private static final double SCORE_CUTOFF = 0.3;
	private static final double SCORE_CUTOFF_DIREC = 0.6;
	
	private PPTriangleMatcher() {}
	
	public static List<AlignmentResult> getMatchingTransforms(Map<Integer,ArrayList<PPTriangle>> triangleSetRef, 
			Map<Integer,ArrayList<PPTriangle>> triangleSetFit, int refConformerId, int conformerId) {
		return getMatchingTransforms(triangleSetRef, 
				triangleSetFit, refConformerId, conformerId, true);
	}
	
	public static List<AlignmentResult> getMatchingTransforms(Map<Integer,ArrayList<PPTriangle>> triangleSetRef, 
			Map<Integer,ArrayList<PPTriangle>> triangleSetFit, int refConformerId, int conformerId, boolean useDirectionality) {
		double cutoff = SCORE_CUTOFF;
		if(!useDirectionality)
			cutoff=SCORE_CUTOFF_DIREC;
		List<AlignmentResult> results = new ArrayList<AlignmentResult>();
		for(int tHash : triangleSetRef.keySet())  {
			if(!triangleSetFit.containsKey(tHash))
				continue;
			List<PPTriangle> fitTriangles = triangleSetFit.get(tHash);
			List<PPTriangle> refTriangles = triangleSetRef.get(tHash);
			for(PPTriangle refTriangle : refTriangles) {
				for(PPTriangle fitTriangle : fitTriangles) {
					if(doEdgeLengthsMatch(refTriangle,fitTriangle)) {
						TransformationSequence transform = new TransformationSequence();
						double score = refTriangle.getMatchingTransform(fitTriangle, transform,useDirectionality);
						if(score>cutoff)
							results.add(new AlignmentResult(score, transform,refConformerId,conformerId));
					}
				}
			}
		}
	
		return results;
		
	}
	
	
	private static boolean doEdgeLengthsMatch(PPTriangle triangle1, PPTriangle triangle2) {
		double[] d1 = triangle1.getEdgeLengths();	
		double[] d2 = triangle2.getEdgeLengths();
		if(Math.abs(d1[0]-d2[0])<CUTOFF && Math.abs(d1[1]-d2[1])<CUTOFF && Math.abs(d1[2]-d2[2])<CUTOFF) {
			return true;
		}
		else
			return false;
	}
	


}
