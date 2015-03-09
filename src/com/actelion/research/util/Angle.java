/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.util;

/**
 * This class provides various operations on angles,
 * while making sure that the values is always in the rangle
 * from -pi >= v > pi.
 */
public class Angle {
	public final double TWO_PI = Math.PI + Math.PI;          
	private double  mValue;

	public Angle() {
		}

	public Angle(Angle a) {
		this(a.mValue);
		}

	public Angle(double value) {
		mValue = value;
		normalize();
		}

	public void normalize() {
		while (mValue < -Math.PI) mValue += TWO_PI;
		while (mValue >= Math.PI) mValue -= TWO_PI;
		}

	public double toDegrees() {
		return 180.0 * mValue / Math.PI;
		}

	public double getValue() {
		return  mValue;
		}

	public void setValue(double value) {
		mValue = value;
		normalize();
		}

	public double cos() {
		return Math.cos(mValue);
		}

	public double sin() {
		return Math.sin(mValue);
		}

	public double tan() {
		return Math.tan(mValue);
		}

	public static Angle	arcsin(double x) {
		return new Angle(Math.asin(x));
		}

	public static Angle	arccos(double x) {
		return new Angle(Math.acos(x));
		}

	public static Angle	arctan(double y, double x) {
		return new Angle(Math.atan2(y, x));
		}

	/**
	 * Determines the mean angle between a1 and a2 on that
	 * side of the circle where a1 and a2 have the shorter connection.
	 * @param a1
	 * @param a2
	 * @return
	 */
	public static double mean(Angle a1, Angle a2) {
		double mean = (a1.mValue + a2.mValue) / 2.0;
		double dif = a2.mValue - a1.mValue;
		if (Math.abs(dif) > Math.PI) {
			if (mean < 0.0)
				mean += Math.PI;
			else
				mean -= Math.PI;
			}
		return mean;
		}

	/**
	 * Determines the normalized difference between a2 and a1,
	 * i.e. the difference on that side where a1 and a2 are closer.
	 * @param a2
	 * @param a1
	 * @return signed value between -pi and +pi
	 */
	public static double difference(double a2, double a1) {
		double a = a2 - a1;
		if (a >= Math.PI)
			a -= 2.0 * Math.PI;
		else if (a < -Math.PI)
			a += 2.0 * Math.PI;
		return a;
		}

	/**
	 * Determines the normalized difference between a2 and a1,
	 * i.e. the difference on that side where a1 and a2 are closer.
	 * @param a2
	 * @param a1
	 * @return signed value between -pi and +pi
	 */
	public static double difference(Angle a2, Angle a1) {
		return difference(a2.mValue, a1.mValue);
		}

	public void add(double value) {
		mValue += value;
		normalize();
		}

	public void add(Angle a) {
		add(a.mValue);
		}

	public void subtract(double value) {
		mValue -= value;
		normalize();
		}

	public void subtract(Angle a) {
		subtract(a.mValue);
		}

	/**
	 * Determines whether this angle is on the smaller value side of angle a,
	 * i.e. whether it is reached faster, when walking the circle towards lower numbers.
	 * @param a
	 * @return
	 */
	public boolean isSmallerThan(Angle a) {
		double dif = a.mValue - mValue;
		return ((dif > 0.0 && dif < Math.PI) || (dif < 0.0 && dif > Math.PI));
		}

	/**
	 * Determines whether this angle is on the greater value side of angle a,
	 * i.e. whether it is reached faster, when walking the circle towards higher numbers.
	 * @param a
	 * @return
	 */
	public boolean isGreaterThan(Angle a) {
		double dif = mValue - a.mValue;
		return ((dif > 0.0 && dif < Math.PI) || (dif < 0.0 && dif > Math.PI));
		}

	public String toString() {
//        return DoubleFormat.toString(toDegrees()) + " degrees";
		return toDegrees() + " degrees";
		}
	}
