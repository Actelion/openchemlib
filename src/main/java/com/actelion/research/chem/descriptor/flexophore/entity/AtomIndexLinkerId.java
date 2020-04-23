package com.actelion.research.chem.descriptor.flexophore.entity;

/**
 * AtomIndexLinkerId
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Jan 17, 2013 MvK Start implementation
 */
public class AtomIndexLinkerId {

	private int atomIndexConnectionFlexophorePoint;
	
	private int atomIndexFirstLinkerAtom;
	
	private int linkerId;
	
	/**
	 * 
	 */
	public AtomIndexLinkerId() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param atomIndexConnectionFlexophorePoint
	 * @param atomIndexFirstLinkerAtom
	 * @param linkerId
	 */
	public AtomIndexLinkerId(int atomIndexConnectionFlexophorePoint,
			int atomIndexFirstLinkerAtom, int linkerId) {
		super();
		this.atomIndexConnectionFlexophorePoint = atomIndexConnectionFlexophorePoint;
		this.atomIndexFirstLinkerAtom = atomIndexFirstLinkerAtom;
		this.linkerId = linkerId;
	}

	/**
	 * @return the atomIndexConnectionFlexophorePoint
	 */
	public int getAtomIndexConnectionFlexophorePoint() {
		return atomIndexConnectionFlexophorePoint;
	}

	/**
	 * @param atomIndexConnectionFlexophorePoint the atomIndexConnectionFlexophorePoint to set
	 */
	public void setAtomIndexConnectionFlexophorePoint(
			int atomIndexConnectionFlexophorePoint) {
		this.atomIndexConnectionFlexophorePoint = atomIndexConnectionFlexophorePoint;
	}

	/**
	 * @return the atomIndexFirstLinkerAtom
	 */
	public int getAtomIndexFirstLinkerAtom() {
		return atomIndexFirstLinkerAtom;
	}

	/**
	 * @param atomIndexFirstLinkerAtom the atomIndexFirstLinkerAtom to set
	 */
	public void setAtomIndexFirstLinkerAtom(int atomIndexFirstLinkerAtom) {
		this.atomIndexFirstLinkerAtom = atomIndexFirstLinkerAtom;
	}

	/**
	 * @return the linkerId
	 */
	public int getLinkerId() {
		return linkerId;
	}

	/**
	 * @param linkerId the linkerId to set
	 */
	public void setLinkerId(int linkerId) {
		this.linkerId = linkerId;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("AtomIndexLinkerId [atomIndexConnectionFlexophorePoint=");
		sb.append(atomIndexConnectionFlexophorePoint);
		sb.append(", atomIndexFirstLinkerAtom=");
		sb.append(atomIndexFirstLinkerAtom);
		sb.append(", linkerId=");
		sb.append(linkerId);
		sb.append("]");
		return sb.toString();
	}
	
	
	
	

}
