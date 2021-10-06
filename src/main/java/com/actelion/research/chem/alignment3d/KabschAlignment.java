package com.actelion.research.chem.alignment3d;

import java.util.Arrays;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.SingularValueDecomposition;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.conf.Conformer;


import java.util.stream.IntStream;

/**
 * described in: DOI 10.1002/jcc.20110  "Using Quaternions to Calculate RMSD"
 * @author wahljo1
 * the first coordinates are rigid, the second coordinate set is rotated
 */

public class KabschAlignment {
	
	private Coordinates[] coords1;
	private Coordinates[] coords2;
	private Coordinates com1, com2;
	int [][] mapping;
	
	public KabschAlignment(Coordinates[] coords1, Coordinates[] coords2, int[][] mapping) {
		this.coords1 = coords1;
		this.coords2 = coords2;
		this.mapping = mapping;
	}
	/**
	 * first conformer is the reference, second conformer is rotated and translated
	 * to minimize the RMSD
	 * @param conf1
	 * @param conf2
	 */
	public KabschAlignment(Conformer conf1, Conformer conf2) {
		this(conf1.getCoordinates(), conf2.getCoordinates(),
				IntStream.range(0, conf1.getMolecule().getAtoms()).mapToObj(e -> new int[]{e,e}).toArray(int[][]::new));
	}
	
	
	private Coordinates getCOM(Coordinates[] coords) {
		int counter = 0;
		Coordinates com = new Coordinates();
		for(Coordinates c:coords) {
			com.add(c);
			counter++;
		}
		com.scale(1.0/counter);
		return com;
	
	}
	
	public void align() {
		align(new Coordinates(), new Matrix(3,3), new Coordinates());
	}
	
	public void align(Coordinates trans1, Matrix rot, Coordinates trans2) {
		
		com1 = getCOM(coords1);
		com2 = getCOM(coords2);
		for(Coordinates c1 : coords1)
			c1.sub(com1);
		
		for(Coordinates c2 : coords2)
			c2.sub(com2);
		
		Coordinates t1 = com2.scaleC(-1.0);
		trans1.x = t1.x;
		trans1.y = t1.y;
		trans1.z = t1.z;
		
		trans2.x = com1.x;
		trans2.y = com1.y;
		trans2.z = com1.z;
		
		Matrix m = new Matrix(3,3);
		double [][] c1 = Arrays.stream(coords1).map(e -> new double[] {e.x,e.y,e.z}).toArray(double[][]::new);
		double [][] c2 = Arrays.stream(coords2).map(e -> new double[] {e.x,e.y,e.z}).toArray(double[][]::new);
		for(int i=0;i<3;i++) {
			for(int j=0;j<3;j++) {
				double rij = 0.0;
				for(int[] map : mapping) {
					rij+= c2[map[1]][i]* c1[map[0]][j];
				}
				m.set(i, j, rij);
			}
		}
		
		SingularValueDecomposition svd = new SingularValueDecomposition(m.getArray(),null,null);

		Matrix u = new Matrix(svd.getU());
		Matrix v = new Matrix(svd.getV());
		Matrix ut = u.getTranspose();
		Matrix vut = v.multiply(ut);
		double det = vut.det();

		Matrix ma = new Matrix(3,3);
		ma.set(0,0,1.0);
		ma.set(1,1,1.0);
		ma.set(2,2,det);
		
		
		Matrix r = ma.multiply(ut);
		r = v.multiply(r);
		assert(r.det()>0.0);
		r = r.getTranspose();
		
		
		
	    for(Coordinates c : coords2) {
	    	c.rotate(r.getArray());
	    	c.add(com1);
	    }
	    
	    for(Coordinates c : coords1) {
	    	c.add(com1);
	    }
	    
	    rot.set(r.getArray());
	    

		
	}
	
	

}

