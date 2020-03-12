package com.actelion.research.chem.phesa;

import com.actelion.research.util.EncoderFloatingPointNumbers;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;

import com.actelion.research.chem.PeriodicTable;
import com.actelion.research.chem.Coordinates;

/** 
 * @version: 1.0, February 2018
 * Author: J. Wahl
 * basic class to describe Gaussian functions used for the calculation of Molecular Volumes
 * Gaussian functions have a center (3D coordinates), a width and a height 
 * this class provides functionalities for calculating higher order overlaps of Gaussians

*/

public class AtomicGaussian extends Gaussian3D {


	
	public AtomicGaussian(int atomId,int atomicNo,Coordinates center){
		super(atomId,atomicNo,center,1.0);
	}
	
	public AtomicGaussian(AtomicGaussian original){
		super(original);

	}
	
	private AtomicGaussian(String encodedGaussian) {
		decode(encodedGaussian);
	}
	
	public static AtomicGaussian fromString(String encodedGaussian) {
		return new AtomicGaussian(encodedGaussian);
	}


	
	
	@Override
	public String encode() { //encodes all information of a atomicGaussian using the Base64 encoder
		double[] coords = new double[] {center.x,center.y,center.z};
		StringBuilder molVolString = new StringBuilder();
		molVolString.append(Integer.toString(atomicNo));
		molVolString.append(" ");
		molVolString.append(Integer.toString(atomId));
		molVolString.append(" ");
		molVolString.append(EncoderFloatingPointNumbers.encode(new double[]{weight},13));
		molVolString.append(" ");
		molVolString.append(EncoderFloatingPointNumbers.encode(coords,13));
		return molVolString.toString();
	}
	
	public void decode(String string64)  {
		String[] strings = string64.split(" ");
		atomicNo = Integer.decode(strings[0]);
		atomId = Integer.decode(strings[1]);
		double [] w = EncoderFloatingPointNumbers.decode(strings[2]);
		double [] coords = EncoderFloatingPointNumbers.decode(strings[3]);
		alpha = calculateWidth(); //the width of the Gaussian depends on the atomic radius of the atom
		volume = calculateVolume();
		coeff = calculateHeight();
		center = new Coordinates(coords[0],coords[1],coords[2]);
		weight = w[0];
		
	}

	@Override
	public double calculateHeight() {
		return MolecularVolume.p;
	}

	@Override
	public double calculateWidth() {
		double vdwR = PeriodicTable.getElement(atomicNo).getVDWRadius();
		return MolecularVolume.alpha_pref/(vdwR*vdwR);
	}


}

