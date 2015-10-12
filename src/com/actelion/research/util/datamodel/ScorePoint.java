

package com.actelion.research.util.datamodel;

import java.awt.Point;


import com.actelion.research.util.Formatter;

/*
 * Project: DD_core
 * @(#)ScorePoint.java
 *
 * Copyright (c) 2003 - 2015
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * Author: MvK
 */
public class ScorePoint extends Point {

	private static final long serialVersionUID = 24052013;
	
	private double score;
	
	public ScorePoint() {
		super(-1,-1);
	}

	public ScorePoint(Point p) {
		super(p);
	}
	
	public ScorePoint(Point p, double value) {
		super(p);
		score = value;
	}
	
	public ScorePoint(int x, int y) {
		super(x,y);
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}
	
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("x " + x);
		sb.append(" y " + y);
		sb.append(" " + Formatter.format3(score));
		
		return sb.toString();
	}
	
}
