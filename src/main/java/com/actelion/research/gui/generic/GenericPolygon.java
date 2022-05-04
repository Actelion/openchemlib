package com.actelion.research.gui.generic;

import java.util.Arrays;

public class GenericPolygon implements GenericShape {
	private double[] mX,mY;
	private int mIndex;

	public GenericPolygon(int size) {
		mX = new double[size];
		mY = new double[size];
		mIndex = 0;
	}

	public GenericPolygon() {
		mX = new double[64];
		mY = new double[64];
		mIndex = 0;
	}

	public void addPoint(double x, double y) {
		if (mIndex == mX.length) {
			mX = Arrays.copyOf(mX, 2*mIndex);
			mY = Arrays.copyOf(mY, 2*mIndex);
			}

		mX[mIndex] = x;
		mY[mIndex] = y;
		mIndex++;
	}

	public void removeLastPoint() {
		if (mIndex > 0)
			mIndex--;
	}

	public void clear() {
		mIndex = 0;
	}

	public boolean contains(double x, double y) {
		boolean result = false;
//		for (int i=0, j=nvert-1; i<nvert; j=i++) {
//			if (((verty[i]>y) != (verty[j]>y))
//			 && (x < (vertx[j]-vertx[i]) * (y-verty[i]) / (verty[j]-verty[i]) + vertx[i]))
		int j = mIndex - 1;
		for (int i=0; i<mIndex; i++) {
			if (((mY[i]>y) != (mY[j]>y))
					&& (x < (mX[j]-mX[i]) * (y-mY[i]) / (mY[j]-mY[i]) + mX[i]))
				result = !result;
			j = i;
		}
		return result;
	}

	public int getSize() {
		return mIndex;
	}

	public double getX(int i) {
		return mX[i];
	}

	public double getY(int i) {
		return mY[i];
	}

	public double[] getX() {
		return mX;
	}

	public double[] getY() {
		return mY;
	}
}
