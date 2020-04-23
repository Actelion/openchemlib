/*
 * Copyright (c) 2018.
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

package com.actelion.research.calc.statistics.median;

/**
 * ModelMedianInteger
 * Feb 9, 2012 MvK: Start implementation
 */
public class ModelMedianInteger {
	
	public int lowerQuartile;
	
	public int median;
	
	public int upperQuartile;
	
	public int id;
	
	public int size;
	
	public int range(){
		return upperQuartile-lowerQuartile;
	}
	
	public String toString() {
		
		StringBuilder sb = new StringBuilder(); 
		
		sb.append(lowerQuartile);
		sb.append("\t");
		sb.append(median);
		sb.append("\t");
		sb.append(upperQuartile);
		
		return sb.toString();
	}

}
