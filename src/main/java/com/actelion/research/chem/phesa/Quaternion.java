package com.actelion.research.chem.phesa;
import com.actelion.research.calc.Matrix;


/**
 * @author J.Wahl, February 2018
 * Describes rotation using quaternion formulation
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
	
	public Matrix getRotMatrix(){ 
		Matrix rotMat = new Matrix(3,3);
		double[][] data = rotMat.getArray();
		double e2 = q1*q1;
		double n2 = q2*q2;
		double s2 = q3*q3;
		double x2 = q0*q0;
		double xe = q0*q1;
		double ne = q1*q2;
		double xs = q0*q3;
		double es = q1*q3;
		double nx = q0*q2;
		double ns = q2*q3;
		data[0][0] = e2-n2-s2+x2;
		data[0][1] = 2*(ne+xs);
		data[0][2] = 2*(es-nx);
		data[1][0] = 2*(ne-xs);
		data[1][1] = -e2+n2-s2+x2;
		data[1][2] = 2*(ns+xe);
		data[2][0] = 2*(es+nx);
		data[2][1] = 2*(ns-xe);
		data[2][2] = -e2-n2+s2+x2;
		
		
		return rotMat;
	}
	
	

}
