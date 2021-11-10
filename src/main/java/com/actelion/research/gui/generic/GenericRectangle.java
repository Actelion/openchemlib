package com.actelion.research.gui.generic;

public class GenericRectangle implements GenericShape {
	private double x1,y1,x2,y2;

	public GenericRectangle(double x, double y, double w, double h) {
		x1 = x;
		y1 = y;
		x2 = x + w;
		y2 = y + h;
	}

	@Override
	public boolean contains(double x, double y) {
		return x>=x1 && x<=x2 && y>=y1 && y<=y2;
	}
}
