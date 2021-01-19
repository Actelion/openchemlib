/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.calc;

public interface SOMController {
	public int getInputVectorCount();		// may be -1 if input vectors are generated on the fly

	/**
	 * In VectorSOM the returned Object is changed by the method.
	 * --> getInputVector(int row) has to provide a deep copy of the object.
	 * @param row
	 * @return deep copy of the Object for VectorSOM
	 */
	public Object getInputVector(int row);	// row may not be used if input vectors are generated on the fly
	}
