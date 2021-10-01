package com.actelion.research.chem.phesa.pharmacophore;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.phesa.pharmacophore.pp.IPharmacophorePoint;
import com.actelion.research.chem.phesa.pharmacophore.pp.PPGaussian;

public class PPTriangleCreator {
	
	private static final double SIDE_LENGTH_CUTOFF = 2.5;
	private static final int MAX_TRIANGLES = 1000;

	
	private PPTriangleCreator() {
	}
	
	public static Map<Integer,ArrayList<PPTriangle>> create(List<PPGaussian> ppGaussians, Coordinates com) {
		List<IPharmacophorePoint> pharmacophorePoints = new ArrayList<IPharmacophorePoint>();
		for(PPGaussian ppGaussian : ppGaussians)
			pharmacophorePoints.add(ppGaussian.getPharmacophorePoint());
		pharmacophorePoints.sort((p1,p2) -> {
			int f1 = p1.getFunctionalityIndex();
			int f2 = p2.getFunctionalityIndex();
			return Integer.compare(f1, f2);	
		});
		List<Integer> toSkip = new ArrayList<Integer>();
		int n = pharmacophorePoints.size();
		double[][] distMat = new double[n][n];
		for(int i=0;i<n;i++) {
			IPharmacophorePoint pp1 = pharmacophorePoints.get(i);
			for(int j=i+1;j<n;j++) {
				IPharmacophorePoint pp2 = pharmacophorePoints.get(j);
				distMat[i][j] = pp1.getCenter().distance(pp2.getCenter());
				if(distMat[i][j]<0.01) 
					toSkip.add(j);
			}
		}
		PPTriangle triangle;
		Map<Integer,ArrayList<PPTriangle>> triangles = new HashMap<Integer,ArrayList<PPTriangle>>();
		int counter = 0;
		for(int i=0;i<n && counter<MAX_TRIANGLES;i++) {
			if(toSkip.contains(i))
				continue;
			IPharmacophorePoint pp1 = pharmacophorePoints.get(i);
			for(int j=i+1;j<n;j++) {
				if(toSkip.contains(j))
					continue;
				IPharmacophorePoint pp2 = pharmacophorePoints.get(j);
				for(int k=j+1;k<n;k++) {
					if(toSkip.contains(k))
						continue;
					IPharmacophorePoint pp3 = pharmacophorePoints.get(k);
					double d12 =  distMat[i][j];
					double d13 = distMat[i][k];
					double d23 = distMat[j][k];
					if(d12<SIDE_LENGTH_CUTOFF  || d13<SIDE_LENGTH_CUTOFF  || d23<SIDE_LENGTH_CUTOFF )
						continue; //skip triangles with very short sides
					triangle = new PPTriangle(pp1,pp2,pp3,d12,d13,d23,com);
					counter++;
					int key = triangle.getHash();
					if(triangles.containsKey(key))
						triangles.get(key).add(triangle);
					else  {
						ArrayList<PPTriangle> triangleList = new ArrayList<PPTriangle>(); 
						triangleList.add(triangle);
						triangles.put(key, triangleList);
					}
					
					
				}
			}
			
		}	

		return triangles;
	}
	

}
