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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ClusterNode implements Comparable<ClusterNode> {
	
	private int indexCenter;
	
	private List<Integer> liIndexMember;
	
	private int rmsd;
	
	public ClusterNode(int indexCenter) {
		this.indexCenter = indexCenter;
		
		liIndexMember = new ArrayList<Integer>();
		
		rmsd = 0;
	}
	
	public void add(int index){
		liIndexMember.add(index);
	}
	
	public int compareTo(ClusterNode cl) {
		return Double.compare(rmsd, cl.rmsd);
	}
	
	@Override
	public boolean equals(Object obj) {
		boolean bEqual = true;
		
		ClusterNode cl = (ClusterNode)obj;
		
		if(liIndexMember.size()!=cl.liIndexMember.size()){
			return false;
		}
		
		List<Integer> l1 = new ArrayList<Integer>(liIndexMember);
		l1.add(indexCenter);
		
		List<Integer> l2 = new ArrayList<Integer>(cl.liIndexMember);
		l2.add(cl.indexCenter);
		
		Collections.sort(l1);
		
		Collections.sort(l2);
		
		for (int i = 0; i < l1.size(); i++) {
			int v1 = l1.get(i);
			int v2 = l2.get(i);
			
			if(v1!=v2){
				bEqual=false;
				break;
			}
		}
		
		return bEqual;
	}

	public int getIndexCenter() {
		return indexCenter;
	}

	public int getRMSD() {
		return rmsd;
	}
	
	/**
	 * Sum of squared distances of the center to the other cluster members
	 * @param rmsd
	 */ 
	public void setRMSD(int rmsd) {
		this.rmsd = rmsd;
	}
	
	public static int getRMSD(ClusterNode cluster, DistHist mdh){
		int rmsd=0;
		
		List<Integer> liIndex = cluster.getClusterMember();
		
		for (int i = 0; i < liIndex.size(); i++) {
			int dist = mdh.getMinDist(cluster.getIndexCenter(), liIndex.get(i));
			rmsd+=dist*dist;
		}
		
		return rmsd;
	}

	public List<Integer> getClusterMember() {
		return liIndexMember;
	}
	
	public boolean isCluster(){
		if(liIndexMember.size()>0)
			return true;
		else
			return false;
	}

}
