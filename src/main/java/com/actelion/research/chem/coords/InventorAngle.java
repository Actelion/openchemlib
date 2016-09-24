package com.actelion.research.chem.coords;

/**
 * Created by thomas on 9/23/16.
 */
public class InventorAngle {
	double mAngle;
	double mLength;

	protected static double getAngle(double x1, double y1, double x2, double y2) {
		double angle;
		double xdif = x2 - x1;
		double ydif = y2 - y1;

		if (ydif != 0) {
			angle = Math.atan(xdif/ydif);
			if (ydif < 0) {
				if (xdif < 0)
					angle -= Math.PI;
				else
					angle += Math.PI;
			}
		}
		else
			angle = (xdif >0) ? Math.PI/2 : -Math.PI/2;

		return angle;
	}

	protected InventorAngle(double angle, double length) {
		mAngle = angle;
		mLength = length;
	}

	protected InventorAngle(double x1, double y1, double x2, double y2) {
		mAngle = getAngle(x1, y1, x2, y2);
		double xdif = x2 - x1;
		double ydif = y2 - y1;
		mLength = Math.sqrt(xdif * xdif + ydif * ydif);
	}
}
