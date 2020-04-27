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

import com.actelion.research.util.graph.complete.ICompleteGraph;

public interface IMolDistHist extends ICompleteGraph {
	
	public double getRelMaxDistInHist(int indexAt1, int indexAt2);
	
	public PPNode getNode(int i);
	
	public byte [] getDistHist(int indexAt1, int indexAt2, byte[] arr);
	
	public boolean isInevitablePharmacophorePoint(int indexNode);
	
	public int getNumInevitablePharmacophorePoints();

	
}
