package com.actelion.research.gui.generic;

public class GenericPoint {
	public double x,y;

	public GenericPoint() {}

	public GenericPoint(double x, double y) {
		this.x = x;
		this.y = y;
	}

	public double getX() {
		return x;
	}

	public double getY() {
		return y;
	}

	public double distance(GenericPoint p) {
		double dx = x - p.x;
		double dy = y - p.y;
		return Math.sqrt(dx*dx+dy*dy);
	}

	@Override
	public String toString() {
		return "x:"+x+" y:"+y;
	}
}
