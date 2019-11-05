package com.actelion.research.chem.phesa;

import com.actelion.research.util.EncoderFloatingPointNumbers;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;

import com.actelion.research.chem.PeriodicTable;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.Coordinates;

/** 
 * @version: 1.0, February 2018
 * Author: J. Wahl
 * basic class to describe Gaussian functions used for the calculation of Molecular Volumes
 * Gaussian functions have a center (3D coordinates), a width and a height 
 * this class provides functionalities for calculating higher order overlaps of Gaussians

*/

public abstract class Gaussian3D {
	public static final double DIST_CUTOFF = 10.0;
	protected int atomId;
	protected int atomicNo;
	protected Coordinates center;
	protected double coeff;
	protected double alpha;
	protected double volume;
	protected double weight;

	
	public Gaussian3D(int atomId, int atomicNo, Coordinates center, double weight){
		this.weight = weight;
		this.atomId = atomId;
		this.atomicNo = atomicNo;
		this.center  = center;
		this.coeff = calculateHeight();
		this.alpha = calculateWidth();
		this.volume = calculateVolume();
	}
	
	public Gaussian3D(Gaussian3D original){
		this.atomId = original.atomId;
		this.atomicNo = original.atomicNo;
		this.coeff = original.coeff;
		this.center = new Coordinates(original.center);
		this.alpha = original.alpha;
		this.volume = original.volume;
		this.weight = original.weight;
	}
	
	public Gaussian3D() {}
	

	
	public abstract double calculateHeight();
	
	public abstract double calculateWidth();
	
	public double calculateVolume() {
		double vdwR = PeriodicTable.getElement(atomicNo).getVDWRadius();		
		double volume = (4.0 * Math.PI/3.0) * vdwR*vdwR*vdwR;
		return volume;
	}
	
	public double getHeight() {
		return this.coeff;
	}
	
	public void setHeight(double height) {
		this.coeff = height;
	}
	
	public double getWidth() {
		return this.alpha;
	}
	
	public void setWidth(double width) {
		this.alpha = width;
	}
	
	public double getVolume() {
		return this.volume;
	}
	
	public void setVolume(double volume) {
		this.volume = volume;
	}
	
	public Coordinates getCenter() {
		return this.center;
	}
	
	public void setCenter(Coordinates center) {
		this.center = center;
	}
	
	public int getAtomicNo() {
		return atomicNo;
	}
	
	public void setAtomicNo(int atomicNo) {
		this.atomicNo = atomicNo;
	}
	
	public int getAtomId() {
		return atomId;
	}
	
	public void setAtomId(int atomId) {
		this.atomId = atomId;
	}
	
	public double getWeight() {
		return weight;
	}

	public void setWeight(double weight) {
		this.weight = weight;
	}
	
	public Coordinates getRotatedCenter(double[][] rotMatrix, double scaleFactor, double[] translation) {
		Coordinates centerCoords = this.getCenter();
		Coordinates centerModCoords = new Coordinates();
		centerModCoords.x = centerCoords.x*rotMatrix[0][0] + centerCoords.y*rotMatrix[0][1] + centerCoords.z*rotMatrix[0][2];
		centerModCoords.y = centerCoords.x*rotMatrix[1][0] + centerCoords.y*rotMatrix[1][1] + centerCoords.z*rotMatrix[1][2];
		centerModCoords.z = centerCoords.x*rotMatrix[2][0] + centerCoords.y*rotMatrix[2][1] + centerCoords.z*rotMatrix[2][2];
		//centerModCoords = this.getCenter().rotateC(rotMatrix); //we operate on the transformed coordinates of the molecule to be fitted
		centerModCoords.scale(scaleFactor); // scale by the inverse squared norm of the quaternion, necessary if quaternion is not a unit quaternion
		centerModCoords.add(translation[0], translation[1], translation[2]);
		return centerModCoords;
	}
	

		

	
	public double getVolumeOverlap(Gaussian3D g2,Coordinates c2, double distCutoff) {
		double alphaSum = getWidth() + g2.getWidth();
		double Vij = 0.0;
		double Kij=0.0;

		double dx = getCenter().x-c2.x;
		double dy = getCenter().y-c2.y;
		double dz = getCenter().z-c2.z;
		double Rij2 = dx*dx+dy*dy+dz*dz;
		if(Rij2<distCutoff) {
			double c = -( getWidth() * g2.getWidth()* Rij2)/alphaSum;
			Kij = getHeight()*g2.getHeight()*QuickMathCalculator.getInstance().quickExp(c); 
			double factor = QuickMathCalculator.getInstance().getPrefactor(getAtomicNo(),g2.getAtomicNo());
			Vij = weight*factor*Kij;
			
		}
		return Vij;
	}
		
	public double getVolumeOverlap(Gaussian3D g2) {
		return getVolumeOverlap(g2,DIST_CUTOFF);
	}
	
	public double getVolumeOverlap(Gaussian3D g2, double distCutoff) {
		return getVolumeOverlap(g2,g2.getCenter(),distCutoff);
	}
	
	public void updateCoordinates(StereoMolecule mol) {
		center = new Coordinates(mol.getCoordinates(atomId));

	}
	
	public void updateCoordinates(Conformer conf) {
		center = new Coordinates(conf.getCoordinates(atomId));

	}
	
	public void updateAtomIndeces(int[] map) {
		atomId = map[atomId];
	}
	
	abstract public String encode();
	

	
}

