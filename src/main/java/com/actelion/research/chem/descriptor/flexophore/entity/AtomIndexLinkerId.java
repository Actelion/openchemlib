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
 */

package com.actelion.research.chem.descriptor.flexophore.entity;

/**
 * AtomIndexLinkerId
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
