/*
 * Copyright (c) 1997 - 2016
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the the copyright holder nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Modest v. Korff
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
	
	public void addMandatoryPharmacophorePoint(int indexNodeAbsolute){
		
		if(isHiddenPharmacophorePointAbsolute((byte)indexNodeAbsolute)) {
			removeHiddenPharmacophorePoint((byte)indexNodeAbsolute);
		}
		
		super.addMandatoryPharmacophorePoint(indexNodeAbsolute);
		
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
		return super.isMandatoryPharmacophorePoint(arrMapHiddenPPPoints[indexNode]);
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
	public HashSet<Byte> getHashSetHiddenIndex() {
		return hsHiddenIndex;
	}

	public boolean isHiddenPharmacophorePointAbsolute(byte indexNodeAbsolute) {
		
		// System.out.println("MolDistHistVizHiddenPPPoints isHiddenPharmacophorePointAbsolute " + indexNodeAbsolute + " " + hsHiddenIndex.contains(indexNodeAbsolute));

		return hsHiddenIndex.contains(indexNodeAbsolute);
	}

	public boolean isInevitablePharmacophorePointAbsolute(int indexNodeAbsolute) {
		return super.isMandatoryPharmacophorePoint(indexNodeAbsolute);
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
					mdhv.addMandatoryPharmacophorePoint(indexNNew);
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
