package com.actelion.research.chem.phesa.pharmacophore;

import java.util.Arrays;

import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.phesa.PheSAAlignment;

public class PPTriangle {
	

	private double[] d = new double[3]; //d12, d13, d23
	private Coordinates[] c = new Coordinates[3]; //d12, d13, d23
	private int[] f = new int[3]; //f1, f2, f3
	private double[][] initialRot; 
	private double[] initialTranslate;
	private double[][] s; //final rotation for alignment
	private double[][] u;
	private double[][] us; //product of final rotation and inverse of rotation performed to align the reference Triangle to the xy-plane 
	private Coordinates molCom; //center of mass of the associated MolecularVolume
	private Coordinates com; // center of mass of the triangle
	private Coordinates[] dirs = new Coordinates[3]; // directionalities
	
	public PPTriangle (IPharmacophorePoint pp1, IPharmacophorePoint pp2, IPharmacophorePoint pp3, double d12, double d13, double d23, Coordinates molCom) {
		this.molCom = molCom;
		f[0] = pp1.getFunctionalityIndex();
		f[1] = pp2.getFunctionalityIndex();	
		f[2] = pp3.getFunctionalityIndex();	
		d[0] = d12;
		d[1] = d13;
		d[2] = d23;
		c[0] = new Coordinates(pp1.getCenter());
		c[1] = new Coordinates(pp2.getCenter());
		c[2] = new Coordinates(pp3.getCenter());
		dirs[0] = pp1.getDirectionality();
		dirs[1] = pp2.getDirectionality();
		dirs[2] = pp3.getDirectionality();
		canonizeOrder();

	}
	
	private void canonizeOrder() {
		if(f[0]!=f[1] && f[0]!=f[2] && f[1]!=f[2]) {
			return;
		}
		else if(f[0]==f[1] && f[1]==f[2]) { //three equal features
			if(d[0]<=d[1]) { //d12<=d13
				if(d[1]<=d[2])  //d13<=d23
					return;
				else { // d12<=d13, d23<d13
					if(d[0]<=d[2])  //d12<=d23<d13  -->swap f1 and f2
						swap(0,1);
					
					else { // d23<d12<d13
						swap(0,1);
						swap(1,2);
					}
				}
			}
			else { // d13<d12
				if(d[2]<=d[1])  //d23<=d13<d12
					swap(0,2);
				else { //d13<d12, d13<d23
					if(d[0]<=d[2])  //d13<d12<=d23
						swap(1,2);
					else { //d13<d23<d12
						swap(0,1);
						swap(0,2);
					}
				}
			}	
		}
	}
	
	private void swap(int i, int j) {
		int fiold = f[i];
		double diold = d[i];
		Coordinates dirold = new Coordinates(dirs[i]);
		dirs[i] = new Coordinates(dirs[j]);
		f[i] = f[j];
		d[i] = d[j];
		f[j] = fiold;
		d[j] = diold;
		dirs[j] = dirold;
		if(i==0 && j==1) {
			Coordinates ciold = new Coordinates(c[1]);
			c[1] = new Coordinates(c[2]); 
			c[2] = new Coordinates(ciold); 
		}
		
		else if(i==0 && j==2) {
			Coordinates ciold = new Coordinates(c[0]);
			c[0] = new Coordinates(c[2]); 
			c[2] = new Coordinates(ciold); 
		}
		
		else if(i==1 && j==2) {
			Coordinates ciold = new Coordinates(c[0]);
			c[0] = new Coordinates(c[1]); 
			c[1] = new Coordinates(ciold); 
		}
			

	}
	
	public int getHash() { // triangles with equal hash match in their functionality
		return f[0] + 10*f[1] + 100*f[2];
	}
	
	public double[] getEdgeLengths() {
		return d;
	}
	
	private void preAlign() {  //move com to origin, rotate triangle into xy-plane (c1 should lie on positive x-axis)
		s = new double[3][3];
		us = new double[3][3];
		com = (c[0].addC(c[1]).addC(c[2])).scaleC(0.33333);
		c[0].sub(com);
		c[1].sub(com);
		c[2].sub(com);
		Coordinates xInt = c[0].unitC();
		Coordinates v = c[1].unitC();
		Coordinates zInt = xInt.cross(v).unitC();
		Coordinates yInt = zInt.cross(xInt).unitC();
		// align three internal axis to lab-frame coord system -> direction cosine
		double[][] m = new double[3][3];
		m[0][0] = xInt.x;
		m[0][1] = xInt.y;
		m[0][2] = xInt.z;
		m[1][0] = yInt.x;
		m[1][1] = yInt.y;
		m[1][2] = yInt.z;
		m[2][0] = zInt.x;
		m[2][1] = zInt.y;
		m[2][2] = zInt.z;
		PheSAAlignment.rotateCoords(c[0], m);
		PheSAAlignment.rotateCoords(c[1], m);
		PheSAAlignment.rotateCoords(c[2], m);
		initialRot = m;
		initialTranslate = new double[]{-com.x,-com.y,-com.z};

		
	}
	
	public double[][] getInitialRot() {
		return initialRot;
	}
	
	public double[] getInitialTranslate() {
		return initialTranslate;
	}
	
	
	public double getMatchingTransform(PPTriangle fitTriangle, double[][] ur) {
		if(initialRot==null && initialTranslate==null)
			preAlign();
		if(fitTriangle.initialRot==null && fitTriangle.initialTranslate==null)
			fitTriangle.preAlign();
		Coordinates a1 = c[0];
		Coordinates a2 = c[1];
		Coordinates a3 = c[2];
		Coordinates b1 = fitTriangle.c[0];
		Coordinates b2 = fitTriangle.c[1];
		Coordinates b3 = fitTriangle.c[2];

		
		double ppFit = ((4-2*Math.min(a1.subC(b1).dist(),2)) + 
				(4-2*Math.min(a2.subC(b2).dist(),2)) + 
					(4-2*Math.min(a3.subC(b3).dist(),2)))/12.0;
		
		
		Matrix mu = new Matrix(initialRot);
		
		u = mu.getTranspose().getArray();
		

		
		PheSAAlignment.multiplyMatrix(u, fitTriangle.initialRot, ur);
		

		ur[3][0] = fitTriangle.initialTranslate[0];  
		ur[3][1] = fitTriangle.initialTranslate[1];  
		ur[3][2] = fitTriangle.initialTranslate[2];  
		
		ur[4][0] = -initialTranslate[0];  
		ur[4][1] = -initialTranslate[1];  
		ur[4][2] = -initialTranslate[2];  
		
		
		//first move com to origin
		
		Coordinates fitComNew = new Coordinates();
		fitComNew.x = fitTriangle.molCom.x-fitTriangle.com.x;
		fitComNew.y = fitTriangle.molCom.y-fitTriangle.com.y;
		fitComNew.z = fitTriangle.molCom.z-fitTriangle.com.z;
		
		Coordinates fitComNewRot = new Coordinates();
		//rotate
		fitComNewRot.x = ur[0][0]*(fitComNew.x)+
				ur[0][1]*(fitComNew.y)+
				ur[0][2]*(fitComNew.z);
		
		fitComNewRot.y = ur[1][0]*(fitComNew.x)+
				ur[1][1]*(fitComNew.y)+
				ur[1][2]*(fitComNew.z);
		
		fitComNewRot.z = ur[2][0]*(fitComNew.x)+
				ur[2][1]*(fitComNew.y)+
				ur[2][2]*(fitComNew.z);
		
		//move to refTriangle
		

		//compare length and directionalities of vectors to triangle origin a and b
		Coordinates a = com.subC(molCom); 
		Coordinates b = fitComNewRot.scale(-1.0); 

		double ra = a.dist();
		double rb = b.dist();
		double Sdirec = a.unitC().dot(b.unitC());
		if(Sdirec<0)
				Sdirec = 0.0;
		double a_b = a.subC(b).dist();

		double comScore = (1-((1-Math.exp(-0.25*Math.sqrt(ra*rb)))*(1-Sdirec)))*Math.exp(-0.125*(a_b*a_b));
		Coordinates fitDir1 = new Coordinates(fitTriangle.dirs[0]);
		Coordinates fitDir2 = new Coordinates(fitTriangle.dirs[1]);
		Coordinates fitDir3 = new Coordinates(fitTriangle.dirs[2]);
		PheSAAlignment.rotateCoords(fitDir1, ur);
		PheSAAlignment.rotateCoords(fitDir2, ur);
		PheSAAlignment.rotateCoords(fitDir3, ur);
		double dirScore = 0.33333*(Math.max(0,fitDir1.dot(dirs[0])) + Math.max(0,fitDir2.dot(dirs[1])) + Math.max(0,fitDir3.dot(dirs[2])));
		return ppFit*dirScore*comScore;
	
	}
	

	
	
	
	
	

}
