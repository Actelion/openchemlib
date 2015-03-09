package com.actelion.research.util.datamodel;

import java.awt.Point;

/**
 * 
 * 
 * PointDouble
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Sep 26, 2011 MvK: Start implementation
 */
public class PointDouble {
	
	public double x;
	public double y;
	
	public PointDouble() {
		
	}
	
	public PointDouble(double x, double y) {
		this.x=x;
		this.y=y;
	}
	
	public PointDouble(Point p) {
		this.x=p.x;
		this.y=p.y;
	}
	
	public void set(double d1, double d2) {
		this.x = d1;
		this.y = d2;
	}

	public void set(PointDouble td) {
		this.x = td.x;
		this.y = td.y;
	}
	
	public static double euclideanDist(PointDouble dt1, PointDouble dt2) {

		double dDist = Math.sqrt((dt1.x-dt2.x)*(dt1.x-dt2.x) + (dt1.y-dt2.y)*(dt1.y-dt2.y));

		return dDist;
	}

}