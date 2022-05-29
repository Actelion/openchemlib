package com.actelion.research.gui.generic;

public class GenericRectangle implements GenericShape {
	public double x,y,width,height;

	public GenericRectangle() {}

	public GenericRectangle(double x, double y, double w, double h) {
		this.x = x;
		this.y = y;
		this.width = w;
		this.height = h;
	}

	@Override
	public boolean contains(double x, double y) {
		return x>=this.x && x<=this.x+this.width && y>=this.y && y<=this.y+this.height;
	}

	public boolean contains(GenericRectangle r) {
		return contains(r.x, r.y) && contains(r.x+r.width, r.y+r.height);
	}

	public void set(double x, double y, double w, double h) {
		this.x = x;
		this.y = y;
		this.width = w;
		this.height = h;
		}

	public GenericRectangle union(GenericRectangle r) {
		double x = Math.min(this.x, r.x);
		double y = Math.min(this.y, r.y);
		double w = Math.max(this.x+this.width, r.x+r.width) - x;
		double h = Math.max(this.y+this.height, r.y+r.height) - y;
		return new GenericRectangle(x, y, w, h);
	}

	public boolean intersects(GenericRectangle r) {
		return (x < r.x+r.width) && (y < r.y+r.height) && (r.x < x+width) && (r.y < y+height);
	}

	public GenericRectangle intersection(GenericRectangle r) {
		if (!intersects(r))
			return null;

		double x = Math.max(this.x, r.x);
		double y = Math.max(this.y, r.y);
		double w = Math.min(this.x+this.width, r.x+r.width) - x;
		double h = Math.min(this.y+this.height, r.y+r.height) - y;
		return new GenericRectangle(x, y, w, h);
	}

	public double getX() {
		return x;
	}

	public double getY() {
		return y;
	}

	public double getWidth() {
		return width;
	}

	public double getHeight() {
		return height;
	}

	public double getCenterX() {
		return x + width / 2.0;
	}

	public double getCenterY() {
		return y + height / 2.0;
	}

	@Override
	public String toString() {
		return "x:"+x+" y:"+y+" w:"+width+" h:"+height;
	}
}
