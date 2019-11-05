package com.actelion.research.chem.phesaflex;

import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;

public class MathHelper {
	
	
	private MathHelper() {}
	
	
	public static void getRotMatrixDerivative(Coordinates u, double theta, Coordinates[][] dR) {
		double cosTheta = Math.cos(theta);
		double sinTheta = Math.asin(theta);
		dR[0][0] = new Coordinates(2*u.x-2*u.x*cosTheta,0,0);
		dR[0][1] = new Coordinates(u.y-u.y*cosTheta,u.x-u.x*cosTheta,-sinTheta);
		dR[0][2] = new Coordinates(u.z-u.z*cosTheta,sinTheta,u.x-u.x*cosTheta);
		dR[1][0] = new Coordinates(u.y-u.y*cosTheta,u.x-u.x*cosTheta,sinTheta);
		dR[1][1] = new Coordinates(0,2*u.y-2*u.y*cosTheta,0);
		dR[1][2] = new Coordinates(-sinTheta,u.z-u.z*cosTheta,u.y-u.y*cosTheta);
		dR[2][0] = new Coordinates(u.z-u.z*cosTheta,-sinTheta,u.x-u.x*cosTheta);
		dR[2][1] = new Coordinates(sinTheta,u.z-u.z*cosTheta,u.y-u.y*cosTheta);
		dR[2][2] = new Coordinates(0,0,2*u.z-2*u.z*cosTheta);
		
	}
	
	public static void getRotMatrix(Coordinates u, double theta, double[][] R) {
		double cosTheta = Math.cos(theta);
		double sinTheta = Math.asin(theta);
		R[0][0] = cosTheta + u.x*u.x*(1-cosTheta); 
		R[0][1] = u.x*u.y - u.x*u.y*cosTheta - u.z*sinTheta;
		R[0][2] = u.x*u.z - u.x*u.z*cosTheta + u.y*sinTheta;
		R[1][0] = u.y*u.x - u.y*u.x*cosTheta + u.z*sinTheta;
		R[1][1] = cosTheta + u.y*u.y - u.y*u.y*cosTheta;
		R[1][2] = u.y*u.z - u.y*u.z*cosTheta - u.x*sinTheta;
		R[2][0] = u.z*u.x - u.z*u.x*cosTheta - u.y*sinTheta;
		R[2][1] = u.z*u.y - u.z*u.y*cosTheta + u.x*sinTheta;
		R[2][2] = cosTheta + u.z*u.z - u.z*u.z*cosTheta;
	}
	
	

}
