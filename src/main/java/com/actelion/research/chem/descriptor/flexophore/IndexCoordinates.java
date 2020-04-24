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

import com.actelion.research.chem.Coordinates;

public class IndexCoordinates {

	private int index;
	
	private int indexOriginalAtom;
	
	private Coordinates coord;
	
	/**
	 * 
	 */
	public IndexCoordinates() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param index
	 * @param coord
	 */
	public IndexCoordinates(int index, int indexOriginalAtom, Coordinates coord) {
		super();
		this.index = index;
		this.indexOriginalAtom = indexOriginalAtom;
		this.coord = coord;
	}

	/**
	 * @return the index
	 */
	public int getIndex() {
		return index;
	}

	/**
	 * @param index the index to set
	 */
	public void setIndex(int index) {
		this.index = index;
	}
	
	/**
	 * @return the indexOriginalAtom
	 */
	public int getIndexOriginalAtom() {
		return indexOriginalAtom;
	}

	/**
	 * @param indexOriginalAtom the indexOriginalAtom to set
	 */
	public void setIndexOriginalAtom(int indexOriginalAtom) {
		this.indexOriginalAtom = indexOriginalAtom;
	}


	/**
	 * @return the coord
	 */
	public Coordinates getCoord() {
		return coord;
	}

	/**
	 * @param coord the coord to set
	 */
	public void setCoord(Coordinates coord) {
		this.coord = coord;
	}

	
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("IndexCoordinates [index=");
		builder.append(index);
		builder.append(", index orig=");
		builder.append(indexOriginalAtom);
		builder.append(", coord=");
		builder.append(coord);
		builder.append("]");
		return builder.toString();
	}

	
}
