/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
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
