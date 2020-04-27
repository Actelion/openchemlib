package com.actelion.research.chem.descriptor.flexophore.entity;

import com.actelion.research.util.datamodel.ByteVec;
import com.actelion.research.util.datamodel.IntArray;

/**
 * 
 * Linker
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Jan 15, 2013 MvK Start implementation
 */
public class Linker {
	
	private static final int CAPACITY_ATOMINDEX = 12;
	
	private int id;
	
	private byte [] arrDistanceHistogram;
	
	private String idcode;
	
	private IntArray iaOriginalAtomIndex;
	
	private int idFlexophorePoint1;
	
	private int idFlexophorePoint2;
	
	private int hash;
	
	public Linker(int id) {
		this.id = id;
		
		hash = -1;
		iaOriginalAtomIndex = new IntArray(CAPACITY_ATOMINDEX);
	}
	
	
	public void calculateHashCode() {
		
		if(arrDistanceHistogram == null || idcode==null){
			hash=-1;
			return;
		}
		
		int hashHist = ByteVec.getHashCode(arrDistanceHistogram);
		
		int hashIdCode = idcode.hashCode();
		
		hash = hashHist ^  hashIdCode;
	}
	
	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}

	public int hashCode() {
		
		return hash;
	}

	public void addOriginalAtomIndex(int atomIndex){
		iaOriginalAtomIndex.add(atomIndex);
	}
	
	public void addOriginalAtomIndex(int [] arrAtomIndex){
		iaOriginalAtomIndex.add(arrAtomIndex);
	}
	
	
	public int [] getOriginalAtomIndex(){
		return iaOriginalAtomIndex.get();
	}

	/**
	 * @return the arrDistanceHistogram
	 */
	public byte[] getDistanceHistogram() {
		return arrDistanceHistogram;
	}


	/**
	 * @param arrDistanceHistogram the arrDistanceHistogram to set
	 */
	public void setDistanceHistogram(byte[] arrDistanceHistogram) {
		this.arrDistanceHistogram = arrDistanceHistogram;
		calculateHashCode();
	}


	/**
	 * @return the idcode
	 */
	public String getIdCode() {
		return idcode;
	}


	/**
	 * @param idcode the idcode to set
	 */
	public void setIdCode(String idcode) {
		this.idcode = idcode;
		calculateHashCode();
	}
	
	
	public boolean equals(Object obj) {
		boolean equal = true;
		
		if(!(obj instanceof Linker)){
			return false;
		}
		
		Linker l = (Linker)obj;
		
		if((arrDistanceHistogram == null) && (l.arrDistanceHistogram != null)){
			return false;
		} else if((arrDistanceHistogram != null) && (l.arrDistanceHistogram == null)){
			return false;
		} else if((idcode == null) && (l.idcode != null)){
			return false;
		} else if((idcode != null) && (l.idcode == null)){
			return false;
		} 
		
		if(arrDistanceHistogram != null){
			
			if(l.arrDistanceHistogram != null){
				
				if(arrDistanceHistogram.length == l.arrDistanceHistogram.length) {
					
					for (int i = 0; i < arrDistanceHistogram.length; i++) {
						if(arrDistanceHistogram[i] != l.arrDistanceHistogram[i]){
							equal = false;
							break;
						}
					}
				}
			}
		}
		
		if(idcode != null){
			if(l.idcode != null){
				if(!idcode.equals(l.idcode)){
					equal = false;
				}
			}
		}
		
		return equal;
	}
	
	/**
	 * @return the idFlexophorePoint1
	 */
	public int getIdFlexophorePoint1() {
		return idFlexophorePoint1;
	}


	/**
	 * @param idFlexophorePoint1 the idFlexophorePoint1 to set
	 */
	public void setIdFlexophorePoint1(int idFlexophorePoint1) {
		this.idFlexophorePoint1 = idFlexophorePoint1;
	}


	/**
	 * @return the idFlexophorePoint2
	 */
	public int getIdFlexophorePoint2() {
		return idFlexophorePoint2;
	}


	/**
	 * @param idFlexophorePoint2 the idFlexophorePoint2 to set
	 */
	public void setIdFlexophorePoint2(int idFlexophorePoint2) {
		this.idFlexophorePoint2 = idFlexophorePoint2;
	}

	

}