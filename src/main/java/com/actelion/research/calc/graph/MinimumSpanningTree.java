package com.actelion.research.calc.graph;

import java.awt.Point;
import java.util.ArrayList;
import java.util.List;

import com.actelion.research.calc.Matrix;

/**
 * MinimumSpanningTree
 * Kruskal's algorithm
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
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
