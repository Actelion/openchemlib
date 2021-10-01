package com.actelion.research.chem.alignment3d.transformation;
import java.util.Random;

import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;


/**
 * @author J.Wahl, February 2018
 * Describes rotation using quaternion formulation
 * convention: q0 -> scalar
 * q1,q2,q3 -> vector part
 * http://www.cs.cmu.edu/~spiff/exp-map/
 * Grassi,98
 * --> we take the transpose because of different convention!
 */

public class Quaternion {
	
	private double q0;
	private double q1;
	private double q2;
	private double q3;
	
	public Quaternion(double q0, double q1, double q2, double q3){
		
		this.q0 = q0;
		this.q1 = q1;
		this.q2 = q2;
		this.q3 = q3;
	}
	
	public Quaternion(Coordinates axis, double angle) {
		if(angle>Math.PI) {
			int n = (int)((angle)/Math.PI);
			angle-=n*Math.PI;
		}
		if(angle<-Math.PI) {
			int n = (int)((Math.abs(angle))/Math.PI);
			angle+=n*Math.PI;
		}
		double c = Math.cos(angle*0.5);
		double s = Math.sin(angle*0.5);
		this.q0 = c;
		this.q1 = s*axis.x;
		this.q2 = s*axis.y;
		this.q3 = s*axis.z;
		normalize();	
			
	}
	
	public void normalize() {
		double s = normSquared();
		double a = Math.sqrt(s);
		q0 = q0 * 1.0/a;
		q1 = q1 * 1.0/a;
		q2 = q2 * 1.0/a;
		q3 = q3 * 1.0/a;
		
	}
	
	public void setQ0(double q0){
		this.q0 = q0;
	}
	public void setQ1(double q1){
		this.q1 = q1;
	}
	public void setQ2(double q2){
		this.q2 = q2;
	}
	public void setQ3(double q3){
		this.q3 = q3;
	}
	
	public double getQ0(){
		return this.q0;
	}
	
	public double getQ1(){
		return this.q1;
	}
	
	public double getQ2(){
		return this.q2;
	}
	
	public double getQ3(){
		return this.q3;
	}
	
	public double normSquared(){
		return (q0*q0+q1*q1+q2*q2+q3*q3);
	}
	
	/**
	 * convert Quaternion to a rotation matrix
	 * @return
	 */
	
	public void multiply(Quaternion r) {
		q0 = r.q0*q0 - r.q1*q1 - r.q2*q2 - r.q3*q3;
		q1 = r.q0*q1 + r.q1*q0 - r.q2*q3 + r.q3*q2;
		q2 = r.q0*q2 + r.q1*q3 + r.q2*q0 - r.q3*q1;
		q3 = r.q0*q3 - r.q1*q2 + r.q2*q1 + r.q3*q0;
	}
	/**
	 * from Grassia: q0 corresponds to qw (scalar part)
	 * @return
	 */
	public Matrix getRotMatrix(){ 
		Matrix rotMat = new Matrix(3,3);
		double[][] data = rotMat.getArray();
		double q1q1 = q1*q1;
		double q2q2 = q2*q2;
		double q0q1 = q0*q1;
		double q1q2 = q1*q2;
		double q3q3 = q3*q3;
		double q0q3 = q0*q3;
		double q1q3 = q1*q3;
		double q0q2 = q0*q2;
		double q2q3 = q2*q3;
		data[0][0] = 1.0-2*(q2q2+q3q3);
		data[1][0] = 2*(q1q2-q0q3);
		data[2][0] = 2*(q1q3+q0q2);
		data[0][1] = 2*(q1q2+q0q3);
		data[1][1] = 1.0-2*(q1q1+q3q3);
		data[2][1] = 2*(q2q3-q0q1);
		data[0][2] = 2*(q1q3-q0q2);
		data[1][2] = 2*(q2q3+q0q1);
		data[2][2] = 1.0-2*(q1q1+q2q2);

		
		
		return rotMat;
	}
	
	/**
	 * get a random quaternion, from: https://stackoverflow.com/questions/31600717/how-to-generate-a-random-quaternion-quickly
	 * @return
	 */
	public static Quaternion getRandomRotation() {
		Random r = new Random();
		double q0 = 0.0;
		double q1 = 0.0;
		double q2 = 0.0;
		double q3 = 0.0;
		double s1 = 0.0;
		double s2 = 0.0;
		boolean notFulfilled = true;
		while(notFulfilled) {
			q0 = -1 + 2 * r.nextDouble();
			q1 = -1 + 2 * r.nextDouble();
			s1 = q0*q0 + q1*q1;
			if(s1<1.0)
				notFulfilled=false;
		}
		notFulfilled = true;
		while(notFulfilled) {
			q2 = -1 + 2 * r.nextDouble();
			q3 = -1 + 2 * r.nextDouble();
			s2 = q2*q2 + q3*q3;
			if(s2<1.0)
				notFulfilled=false;
		}
		double s = Math.sqrt((1-s1)/s2);
		q2 *= s;
		q3 *= s;
		
		return new Quaternion(q0,q1,q2,q3);
	}
	
	@Override 
	public String toString() {
		String s = q0 + " " + q1 + " " + q2 + " " + q3;
		return s;
	}
	
	

}
