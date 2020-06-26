package com.actelion.research.chem.docking;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.SingularValueDecomposition;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.phesa.AtomicGaussian;
import com.actelion.research.chem.phesa.PheSAAlignment;

public class DockingUtils {
	
	private DockingUtils() {};
	
	public static Coordinates getCOM(Conformer conf) {
		int counter = 0;
		Coordinates com = new Coordinates();
		for(Coordinates coords:conf.getCoordinates()) {
			com.add(coords);
			counter++;
		}
		com.scale(1.0/counter);
		return com;
	}
	
	public static Matrix createInitialOrientation(Conformer conf) {
		Matrix m = calculateMassCovarianceMatrix(conf);
		SingularValueDecomposition svd = new SingularValueDecomposition(m.getArray(),null,null);
		Matrix u = new Matrix(svd.getU());
		double det = u.det();
		if(det<0) {
			u.set(0,1,-u.get(0, 1));
			u.set(1,1,-u.get(1, 1));
			u.set(2,1,-u.get(2, 1));
		}
		PheSAAlignment.rotateMol(conf,u);
		return u;
		
	}
	
	public static Matrix calculateMassCovarianceMatrix(Conformer conf) {
		Matrix massMatrix = new Matrix(3,3); 
		int counter = 0;
		for (Coordinates coords : conf.getCoordinates()){
			counter++;
			double value = coords.x*coords.x;
			massMatrix.addToElement(0,0,value);
			value = coords.x*coords.y;
			massMatrix.addToElement(0,1,value);
			value = coords.x*coords.z;
			massMatrix.addToElement(0,2,value);
			value = coords.y*coords.y;
			massMatrix.addToElement(1,1,value);
			value = coords.y*coords.z;
			massMatrix.addToElement(1,2,value);
			value = coords.z*coords.z;
			massMatrix.addToElement(2,2,value);	
		}
		massMatrix.set(0,0,massMatrix.get(0,0)/counter);
		massMatrix.set(0,1,massMatrix.get(0,1)/counter);
		massMatrix.set(0,2,massMatrix.get(0,2)/counter);
		massMatrix.set(1,1,massMatrix.get(1,1)/counter);
		massMatrix.set(1,2,massMatrix.get(1,2)/counter);
		massMatrix.set(2,2,massMatrix.get(2,2)/counter);
		massMatrix.set(1,0,massMatrix.get(0,1));
		massMatrix.set(2,0,massMatrix.get(0,2));
		massMatrix.set(2,1,massMatrix.get(1,2));
		
		return massMatrix;
	}
	

}
