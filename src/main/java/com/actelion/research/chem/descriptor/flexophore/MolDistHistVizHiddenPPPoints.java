/*
 * Copyright (c) 2020.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
 */

package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.chem.Molecule3D;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashSet;

public class MolDistHistVizHiddenPPPoints extends MolDistHistViz implements Serializable {

	private static final long serialVersionUID = 22102012;
	
	public static final String TAG_FLEXOPHORE_ALL = "FlexophoreWithHidden";
	
	public static final String TAG_MOLDISTHISTVIZ_HIDDENPPPOINTS = "MolDistHistVizHiddenPPPoints";
	
	/**
	 * The native object MolDistHistVizHiddenPPPoints.
	 * For the map in FFMolecule etc.
	 */
	public static final String TAG_FLEXOPHORE_OBJECT = "MolDistHistVizHiddenPPPointsObject";
	
	public static final String TAG_FLEXOPHORE_OBJECTLIST = "ListMolDistHistVizHiddenPPPointsObject";
	
	
	private byte [] arrMapHiddenPPPoints;
	
	private HashSet<Byte> hsHiddenIndex;


	
	public MolDistHistVizHiddenPPPoints(MolDistHistVizHiddenPPPoints mdhvhpp) {
		super(mdhvhpp);
		
		init();
		
		hsHiddenIndex.addAll(mdhvhpp.hsHiddenIndex);
		
		hidden2Map();
		
	}
	
	public MolDistHistVizHiddenPPPoints(MolDistHistViz mdhv) {
		super(mdhv);
		init();
	}
	
	public MolDistHistVizHiddenPPPoints(MolDistHist mdh) {
		super(mdh);
		init();
	}
	
	private void init(){
		hsHiddenIndex = new HashSet<Byte>();
		
		int size = super.getNumPPNodes();
		
		arrMapHiddenPPPoints = new byte[size];
		
		hidden2Map();
	}
	
	private void hidden2Map(){
		
		int size = super.getNumPPNodes();
		
		int cc=0;
		for (byte i = 0; i < size; i++) {
			if(!hsHiddenIndex.contains(i)){
				arrMapHiddenPPPoints[cc++]=i;
			}
		}

		
	}
	
	public void addInevitablePharmacophorePoint(int indexNodeAbsolute){
		
		if(isHiddenPharmacophorePointAbsolute((byte)indexNodeAbsolute)) {
			removeHiddenPharmacophorePoint((byte)indexNodeAbsolute);
		}
		
		super.addInevitablePharmacophorePoint(indexNodeAbsolute);
		
	}

	public void removeInevitablePharmacophorePoint(int indexNodeAbsolute){
		
		// System.out.println("MolDistHistVizHiddenPPPoints removeInevitablePharmacophorePoint " + indexNodeAbsolute);
		
		super.removeInevitablePharmacophorePoint(indexNodeAbsolute);
	}

	public void addHiddenPharmacophorePoint(byte indexPPNodeAbsolute){
		
		// System.out.println("MolDistHistVizHiddenPPPoints addHiddenPharmacophorePoint " + indexPPNodeAbsolute);
		
		hsHiddenIndex.add(indexPPNodeAbsolute);
		
		hidden2Map();
	}
	
	public void removeHiddenPharmacophorePoint(byte indexPPNodeAbsolute){
		
		// System.out.println("MolDistHistVizHiddenPPPoints removeHiddenPharmacophorePoint " + indexPPNodeAbsolute);

		hsHiddenIndex.remove(indexPPNodeAbsolute);
		
		hidden2Map();

	}
	
	public void resetHiddenPharmacophorePoints(){
		
		hsHiddenIndex.clear();
		
		hidden2Map();
	}

	public double getRelMaxDistInHistSkipHidden(int indexAt1, int indexAt2) {
		return super.getRelMaxDistInHist(arrMapHiddenPPPoints[indexAt1], arrMapHiddenPPPoints[indexAt2]);
	}
	
	public PPNodeViz getNodeSkipHidden(int indexNode) {
		return super.getNode(arrMapHiddenPPPoints[indexNode]);
	}
	
	public byte [] getDistHistSkipHidden(int indexAt1, int indexAt2, byte [] arr) {
		return super.getDistHist(arrMapHiddenPPPoints[indexAt1], arrMapHiddenPPPoints[indexAt2], arr);
	}
	
	public byte [] getDistHistSkipHidden(int indexAt1, int indexAt2) {
		return super.getDistHist(arrMapHiddenPPPoints[indexAt1], arrMapHiddenPPPoints[indexAt2]);
	}
	
	public boolean isInevitablePharmacophorePointSkipHidden(int indexNode) {
		return super.isInevitablePharmacophorePoint(arrMapHiddenPPPoints[indexNode]);
	}
	
	public int getNumInevitablePharmacophorePointsSkipHidden() {
		
		int inevitable=0;
		
		int size = getSizeSkipHidden();
		for (int i = 0; i < size; i++) {
			if(isInevitablePharmacophorePointSkipHidden(i)){
				inevitable++;
			}
		}
		
		return inevitable;
		
	}
	
	public int getSizeSkipHidden() {
		return super.getNumPPNodes()-hsHiddenIndex.size();
	}
	
//	public int getSizeAbsolute() {
//		return super.getSize();
//	}
	
//	public PPNodeViz getNodeTotal(int indexNodeAbsolute) {
//		return super.getNode(indexNodeAbsolute);
//	}
	
	/**
	 * @return the hsHiddenIndex
	 */
	protected HashSet<Byte> getHashSetHiddenIndex() {
		return hsHiddenIndex;
	}

	public boolean isHiddenPharmacophorePointAbsolute(byte indexNodeAbsolute) {
		
		// System.out.println("MolDistHistVizHiddenPPPoints isHiddenPharmacophorePointAbsolute " + indexNodeAbsolute + " " + hsHiddenIndex.contains(indexNodeAbsolute));

		return hsHiddenIndex.contains(indexNodeAbsolute);
	}

	public boolean isInevitablePharmacophorePointAbsolute(int indexNodeAbsolute) {
		return super.isInevitablePharmacophorePoint(indexNodeAbsolute);
	}
	

	/**
	 * 
	 * @return deep copy.
	 */
	public MolDistHistViz getMolDistHistVizSkipHidden(){
		
//		System.out.println("MolDistHistVizHiddenPPPoints getMolDistHistVizSkipHidden()");
		
		Molecule3D ffCopy = finalizeMolecule(molecule3D);
		
		int size = getSizeSkipHidden();
		
		MolDistHistViz mdhv = new MolDistHistViz(size, ffCopy);
		
		
		int sizeAbsolute = getNumPPNodes();
		
		for (int i = 0; i < sizeAbsolute; i++) {
			
			if(!isHiddenPharmacophorePointAbsolute((byte)i)){
			
				PPNodeViz node = new PPNodeViz(getNode(i));
				
				int indexNNew = mdhv.addNode(node);
				
				if(isInevitablePharmacophorePointAbsolute(i)){
					mdhv.addInevitablePharmacophorePoint(indexNNew);
				}
				
//				System.out.println(i + " " + node.toString());
			
			}
		}
		
		for (int i = 0; i < size; i++) {
			
			for (int j = i+1; j < size; j++) {
				
				byte [] arrDistHist = getDistHistSkipHidden(i, j);
				
				mdhv.setDistHist(i, j, arrDistHist);
			}
		}
		
		mdhv.realize();
				
		return mdhv;
		
	}
	
	public String toStringHiddenAndInevitable() {
		
		StringBuilder sb = new StringBuilder();
		
		sb.append(toStringInevitable());
		sb.append("\n");
		sb.append("Num hidden " + hsHiddenIndex.size());
		sb.append("\n");
		sb.append("Index hidden ");
		
		for (byte indexHidden : hsHiddenIndex) {
			sb.append(indexHidden + " ");
		}
		sb.append("\n");
		sb.append("Map hidden " + Arrays.toString(arrMapHiddenPPPoints));

		return sb.toString();
	}
	

}
