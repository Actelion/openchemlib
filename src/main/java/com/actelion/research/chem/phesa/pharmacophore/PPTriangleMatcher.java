package com.actelion.research.chem.phesa.pharmacophore;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class PPTriangleMatcher {
	
	private static final double CUTOFF = 2.5; //if lengths of edges of two triangles differ by more than 2.5A, it is
											  // not considered as a match
	
	private PPTriangleMatcher() {}
	
	public static List<AlignmentResult> getMatchingTransforms(Map<Integer,ArrayList<PPTriangle>> triangleSetRef, 
			Map<Integer,ArrayList<PPTriangle>> triangleSetFit, int refConformerId, int conformerId) {
		List<AlignmentResult> results = new ArrayList<AlignmentResult>();
		for(int tHash : triangleSetRef.keySet())  {
			List<PPTriangle> fitTriangles = triangleSetFit.get(tHash);
			if(fitTriangles!=null) {
				List<PPTriangle> refTriangles = triangleSetRef.get(tHash);
				for(PPTriangle refTriangle : refTriangles) {
					for(PPTriangle fitTriangle : fitTriangles) {
						if(doEdgeLengthsMatch(refTriangle,fitTriangle)) {
							double[][] transform = new double[5][3];
							double score = refTriangle.getMatchingTransform(fitTriangle, transform);
							results.add(new AlignmentResult(score, transform,refConformerId,conformerId));
						}
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
	
	public static class AlignmentResult{
		private double similarity;
		private double[][] transform;
		private int conformerIndex;
		private int refConformerIndex;
		
		public AlignmentResult(double similarity, double[][] transform, int refConformerIndex, int conformerIndex) {
			this.similarity = similarity;
			this.transform = transform;
			this.refConformerIndex = refConformerIndex;
			this.conformerIndex = conformerIndex;
		}
		
		
		public double[][] getTransform() {
			return transform;
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
		
		
		
	}

}
