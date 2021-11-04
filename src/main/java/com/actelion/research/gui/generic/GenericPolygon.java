package com.actelion.research.gui.generic;

import java.util.Arrays;

public class GenericPolygon implements GenericShape {
	private double[] mPoint;
	private int mIndex;

	public GenericPolygon(int size) {
		mPoint = new double[2*size];
		mIndex = 0; // is twice of used point count
	}

	public GenericPolygon() {
		mPoint = new double[128];
		mIndex = 0; // is twice of used point count
	}

	public void addPoint(double x, double y) {
		if (mIndex == mPoint.length)
			mPoint = Arrays.copyOf(mPoint, 2 * mIndex);

		mPoint[mIndex++] = x;
		mPoint[mIndex++] = y;
	}

	public void removeLastPoint() {
		if (mIndex >= 2)
			mIndex -= 2;
	}

	public void clear() {
		mIndex = 0;
	}

	public boolean contains(double x, double y) {
		boolean result = false;
//		for (int i=0, j=nvert-1; i<nvert; j=i++) {
//			if (((verty[i]>y) != (verty[j]>y))
//			 && (x < (vertx[j]-vertx[i]) * (y-verty[i]) / (verty[j]-verty[i]) + vertx[i]))
		int j = mIndex - 2;
		for (int i=0; i<mIndex; i+=2) {
			if (((mPoint[i+1]>y) != (mPoint[j+1]>y))
					&& (x < (mPoint[j]-mPoint[i]) * (y-mPoint[i+1]) / (mPoint[j+1]-mPoint[i+1]) + mPoint[i]))
				result = !result;
			j = i;
		}
		return result;
	}

	public int getSize() {
		return mIndex / 2;
	}

	public double getX(int i) {
		return mPoint[2*i];
	}

	public double getY(int i) {
		return mPoint[2*i+1];
	}

	public double[] getPoints() {
		return mPoint;
	}
}
