package com.actelion.research.util.datamodel;

import java.awt.Point;

import com.actelion.research.util.Formatter;

/**
 * 
 * ScorePoint
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * 10 Dec 2010 MvK: Start implementation
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
