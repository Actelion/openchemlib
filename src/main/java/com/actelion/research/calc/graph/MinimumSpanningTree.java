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

package com.actelion.research.calc.graph;

import com.actelion.research.calc.Matrix;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * MinimumSpanningTree
 * Kruskal's algorithm
 * @author Modest von Korff
 * @version 1.0
 * Dec 17, 2012 MvK Start implementation
 */
public class MinimumSpanningTree {
	
	private Matrix maAdjacency;
	
	/**
	 * 
	 * @param maAdjacency adjacency matrix, only the upper triangle is used. 
	 * Non bonding fields have to be indicated as <code>NaN</code>.
	 */
	public MinimumSpanningTree(Matrix maAdjacency) {
		this.maAdjacency = maAdjacency;
	}
	
	/**
	 * 
	 * @return symmetric adjacency matrix. If an edge is a member of the tree the upper and 
	 * lower triangle contain the corresponding value from the input adjacency matrix. 
	 */
	public Matrix getMST(){
		Matrix maMST = new Matrix(maAdjacency.rows(), maAdjacency.cols());
		
		int [] arrNode = new int [maAdjacency.rows()];
		
		Point pStart = getStartEdge();
		
		
		arrNode[pStart.x]=1;
		arrNode[pStart.y]=1;
		
		List<Point> liMST = new ArrayList<Point>();
		
		liMST.add(pStart);
		
		int edgesMinSpanningTree = maAdjacency.rows()-1;
		
		int ccAddedEdges2MST = 0;
		
		while(ccAddedEdges2MST < edgesMinSpanningTree){
			
			double min = Double.MAX_VALUE;
			int rowMin=-1;
			int colMin=-1;
			
			for (int i = 0; i < maAdjacency.rows(); i++) {
				
				if(arrNode[i]==1) {
				
					for (int j = 0; j < maAdjacency.cols(); j++) {
						
						if(i==j){
							continue;
						}
						
						if(arrNode[j]==1) {
							continue;
						}
						
						if(!Double.isNaN(maAdjacency.get(i, j))){
								
							if(maAdjacency.get(i, j) < min){
								min = maAdjacency.get(i, j);
								
								rowMin = i;
								
								colMin = j;
								
							}
						}
					}
				}
			}
			
			// No neighbor node found 
			if(rowMin==-1){
				break;
			} else {
				arrNode[rowMin]=1;
				
				arrNode[colMin]=1;
				
				liMST.add(new Point(rowMin, colMin));

				ccAddedEdges2MST++;
			} 
			
		}
		
		maMST.set(Double.NaN);
		
		for (Point p : liMST) {
			
			maMST.set(p.x, p.y, maAdjacency.get(p.x, p.y));
			maMST.set(p.y, p.x, maAdjacency.get(p.x, p.y));
			
		}
		
		return maMST;
	}
	
	
	private Point getStartEdge(){
		
		double min = Double.MAX_VALUE;

		Point pMin = new Point();

		for (int i = 0; i < maAdjacency.rows(); i++) {

			for (int j = i + 1; j < maAdjacency.cols(); j++) {

				if (!Double.isNaN(maAdjacency.get(i, j))) {

					if (maAdjacency.get(i, j) < min) {
						min = maAdjacency.get(i, j);

						pMin.x = i;

						pMin.y = j;

					}
				}
			}
		}
		
		return pMin;
	}
}
