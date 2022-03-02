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
