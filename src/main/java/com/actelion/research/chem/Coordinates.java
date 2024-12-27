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
package com.actelion.research.chem;

import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;
import java.util.Random;

/**
 * Class to encapsulate 3D coordinates
 */
public final class Coordinates implements Serializable, Comparable<Coordinates> {

	public double x,y,z;

	public Coordinates() {}

	public Coordinates(Coordinates c) {
		this(c.x, c.y, c.z);
	}

	public Coordinates(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	/**
	 * Copies x,y,z from c to this
	 * @param c
	 * @return this after copying c
	 */
	public Coordinates set(Coordinates c) {
		set(c.x, c.y, c.z);
		return this;
	}

	public void set(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public double getLength() {
		return dist();
	}

	public double dist() {
		return Math.sqrt(distSq());
	}
	public double distSq() {
		return  x*x + y*y + z*z;
	}

	public double distanceSquared(Coordinates c) {
		return (c.x-x)*(c.x-x) + (c.y-y)*(c.y-y) + (c.z-z)*(c.z-z);
	}
	public double distSquareTo(Coordinates c) {
		return distanceSquared(c);
	}
	public double distance(Coordinates c) {
		return Math.sqrt(distanceSquared(c));
	}

	public double dot(Coordinates c) {
		return x*c.x + y*c.y + z*c.z;
	}

	public Coordinates cross(Coordinates c) {
		return new Coordinates(y*c.z-z*c.y, -(x*c.z-z*c.x), x*c.y-y*c.x);
	}

	/**
	 * Gets the angle formed between the 2 vectors ([0,PI])
	 * @param c
	 * @return angle in radian
	 */
	public double getAngle(Coordinates c) {
		double d1 = distSq();
		double d2 = c.distSq();
		if(d1==0 || d2==0) return 0;
		double d = dot(c) / Math.sqrt(d1*d2);
		if(d>=1) return 0;
		if(d<=-1) return Math.PI;
		return Math.acos(d);
	}

	/**
	 * Calculates the angle of the line from this location to c
	 * projected into the x/y plane. With Y facing upwards and X right,
	 * if the line points in Y direction, ten angle is 0.0 increasing
	 * in clockwise direction.
	 * @param c
	 * @return -PI < angle < PI
	 */
	public double getAngleXY(Coordinates c) {
		double dx = c.x - x;
		double dy = c.y - y;

		if (dy == 0.0)
			return (dx > 0.0) ? Math.PI/2.0 : -Math.PI/2.0;

		double angle = Math.atan(dx/dy);
		if (dy < 0.0)
			return (dx < 0.0) ? angle - Math.PI : angle + Math.PI;

		return angle;
	}

	public double getDihedral(Coordinates c2, Coordinates c3, Coordinates c4) {
		return getDihedral(this, c2, c3, c4);
	}

	public Coordinates subC(Coordinates c) {
		return new Coordinates(x-c.x, y-c.y, z-c.z);
	}

	public Coordinates addC(Coordinates c) {
		return new Coordinates(x+c.x, y+c.y, z+c.z);
	}

	public Coordinates scaleC(double scale) {
		return new Coordinates(x*scale, y*scale, z*scale);
	}

	/**
	 * @param c
	 * @return this after subtracting c
	 */
	public final Coordinates sub(Coordinates c) {
		x-=c.x;
		y-=c.y;
		z-=c.z;
		return this;
	}

	/**
	 * @param c
	 * @return this after subtracting c
	 */
	public Coordinates add(Coordinates c) {
		x+=c.x;
		y+=c.y;
		z+=c.z;
		return this;
	}

	public void add(double dx, double dy, double dz) {
		x += dx;
		y += dy;
		z += dz;
	}

	/**
	 * @param scale
	 * @return this after scaling this
	 */
	public Coordinates scale(double scale) {
		x*=scale;
		y*=scale;
		z*=scale;
		return this;
	}
	public final void negate() {
		x=-x;
		y=-y;
		z=-z;
	}

	/**
	 * @param m
	 * @return this after rotating it with rotation matrix m
	 */
	public Coordinates rotate(double[][] m) {
		double x0 = x;
		double y0 = y;
		double z0 = z;
		x = x0*m[0][0]+y0*m[1][0]+z0*m[2][0];
		y = x0*m[0][1]+y0*m[1][1]+z0*m[2][1];
		z = x0*m[0][2]+y0*m[1][2]+z0*m[2][2];
		return this;
	}

	/**
	 * @param m
	 * @return new Coordinates created from this point rotated by rotation matrix m
	 */
	public Coordinates rotateC(double[][] m) {
		return new Coordinates(x*m[0][0]+y*m[1][0]+z*m[2][0],
							   x*m[0][1]+y*m[1][1]+z*m[2][1],
							   x*m[0][2]+y*m[1][2]+z*m[2][2]);
	}

	public Coordinates rotate(Coordinates normal, double theta) {
		if(Math.abs(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z-1)>1E-6) throw new IllegalArgumentException("normal needs to a unit vector: "+normal);
		double x = normal.x;
		double y = normal.y;
		double z = normal.z;
		double c = Math.cos(theta);
		double s = Math.sin(theta);
		double t = 1-c;
		Coordinates opp = new Coordinates(
				(t*x*x+c)*this.x 	+ (t*x*y+s*z)*this.y + (t*x*z-s*y)*this.z,
				(t*x*y-s*z)*this.x	+ (t*y*y+c)*this.y 	+ (t*y*z+s*x)*this.z,
				(t*x*z+s*y)*this.x 	+ (t*z*y-s*x)*this.y + (t*z*z+c)*this.z
		);
		return opp;
	}


	/**
	 * @return new Coordinates with a copy of this scaled to length=1.0
	 */
	public Coordinates unitC() {
		double d = dist();
		if(d==0) {
			System.err.println("Cannot call unitC() on a null vector");
//			Thread.dumpStack();
			return new Coordinates(1,0,0);

		}
		return new Coordinates(x/d, y/d, z/d);
	}

	/**
	 * @return this after scaling it to length=1.0
	 */
	public Coordinates unit() {
		double d = dist();
		if(d==0) {
			System.out.println("Cannot call unit() on a null vector. Returned (1,0,0)");
//			Thread.dumpStack();
			x = 1;
			y = 0;
			z = 0;
//			throw new IllegalArgumentException("Cannot call unit() on a null vector.");
			return this;
		}
		x/=d;
		y/=d;
		z/=d;
		return this;
	}

	/**
	 * Calculates the center point between this and c and sets this to the center point.
	 * @param c
	 * @return this after updating it to the center position
	 */
	public Coordinates center(Coordinates c) {
		x = (x + c.x) / 2.0;
		y = (y + c.y) / 2.0;
		z = (z + c.z) / 2.0;
		return this;
	}

	/**
	 * Updates this to contains the center between c1 and c2.
	 * @param c1
	 * @param c2
	 */
	public void center(Coordinates c1, Coordinates c2) {
		x = (c1.x + c2.x) / 2.0;
		y = (c1.y + c2.y) / 2.0;
		z = (c1.z + c2.z) / 2.0;
	}

	/**
	 * Updates this to contain a point on the straight line through c1 and c2.
	 * @param c1
	 * @param c2
	 * @param f location on line 0.0 -> c1, 1.0 -> c2
	 * @return this after updating to be a point on the line
	 */
	public Coordinates between(Coordinates c1, Coordinates c2, double f) {
		x = c1.x + f * (c2.x - c1.x);
		y = c1.y + f * (c2.y - c1.y);
		z = c1.z + f * (c2.z - c1.z);
		return this;
	}

	public boolean insideBounds(Coordinates[] bounds) {
		return bounds!=null && bounds[0].x<=x && x<=bounds[1].x && bounds[0].y<=y && y<=bounds[1].y && bounds[0].z<=z && z<=bounds[1].z;
	}

	@Override
	public String toString() {
		DecimalFormat df = new DecimalFormat("0.00");
		return "[" + df.format(x) + ", " + df.format(y) + ", " + df.format(z) + "]";
	}

	public String toStringSpaceDelimited() {
		DecimalFormat df = new DecimalFormat("0.00");
		return df.format(x) + " " + df.format(y) + " " + df.format(z);
	}

	public String toStringSpaceDelimited(Locale locale) {
		NumberFormat nf = NumberFormat.getNumberInstance(locale);
		DecimalFormat df = (DecimalFormat)nf;
		df.applyPattern("0.00");
		return df.format(x) + " " + df.format(y) + " " + df.format(z);
	}

	@Override
	public boolean equals(Object o) {
		if(o==null || !(o instanceof Coordinates)) return false;
		Coordinates c = (Coordinates) o;
		return Math.abs(c.x-x) + Math.abs(c.y-y) + Math.abs(c.z-z) < 1E-6;
	}

	public boolean isNaN() {
		return Double.isNaN(x) || Double.isNaN(y) || Double.isNaN(z);
	}

	public Coordinates min(Coordinates c) {
		return new Coordinates(Math.min(x, c.x), Math.min(y, c.y), Math.min(z, c.z));
	}
	public Coordinates max(Coordinates c) {
		return new Coordinates(Math.max(x, c.x), Math.max(y, c.y), Math.max(z, c.z));
	}

	public double cosAngle(Coordinates c) {
		double d = dist() * c.dist();
		if(d<=0) return 0;
		return dot(c) / d;
	}

	/////////////////// UTILITITIES ///////////////////////////////////////
	public static Coordinates min(Coordinates[] c) {
		Coordinates min = new Coordinates(c[0]);
		for (int i = 1; i < c.length; i++) {
			min.x = Math.min(c[i].x, min.x);
			min.y = Math.min(c[i].y, min.y);
			min.z = Math.min(c[i].z, min.z);
		}
		return min;
	}

	public static Coordinates max(Coordinates[] c) {
		Coordinates max = new Coordinates(c[0]);
		for (int i = 1; i < c.length; i++) {
			max.x = Math.max(c[i].x, max.x);
			max.y = Math.max(c[i].y, max.y);
			max.z = Math.max(c[i].z, max.z);
		}
		return max;
	}

	public static Coordinates createBarycenter(Coordinates... coords) {
		if(coords==null) throw new IllegalArgumentException("The coordinates are null");
		Coordinates res = new Coordinates();
		for(int i=0; i<coords.length; i++) {
			res.x += coords[i].x;
			res.y += coords[i].y;
			res.z += coords[i].z;
		}
		res.x /= coords.length;
		res.y /= coords.length;
		res.z /= coords.length;
		return res;
	}

	/**
	 * Get the mirror image of p through the plane defined by c1, c2, c3
	 * @param p
	 * @param c1
	 * @param c2
	 * @param c3
	 * @return
	 */
	public static Coordinates getMirror(Coordinates p, Coordinates c1, Coordinates c2, Coordinates c3) {
		//define a unit normal vector to the plane
		Coordinates r31 = new Coordinates(c3);
		r31.sub(c1);
		Coordinates r21 = new Coordinates(c2);
		r21.sub(c1);
		Coordinates c = r31.cross(r21);
		if(c.distSq()<0.05) return new Coordinates(p);
		Coordinates n = c.unitC();

		Coordinates pc1 = new Coordinates(c1);
		pc1.sub(p);
		double l = pc1.dot(n);
		n.scale(2*l);
		Coordinates pp = new Coordinates(p);
		pp.add(n);
		return pp;
	}

/*	public static double getDihedral(Coordinates c1, Coordinates c2, Coordinates c3, Coordinates c4) {
		//Coordinates c1 = this;
		// Calculate the vectors between atoms
		Coordinates r12 = new Coordinates(c1);
		r12.sub(c2);
		Coordinates r23 = new Coordinates(c2);
		r23.sub(c3);
		Coordinates r34 = new Coordinates(c3);
		r34.sub(c4);

		//  Calculate the cross products
		Coordinates A = r12.cross(r23);
		Coordinates B = r23.cross(r34);
		Coordinates C = r23.cross(A);

		//  Calculate the distances
		double rA = A.dist();
		double rB = B.dist();
		double rC = C.dist();

		//  Calculate the sin and cos
		//  cos = A*B/(rA*rB)
		//  sin = C*B/(rC*rB)
		double cos_phi = A.dot(B) / (rA * rB);
		double sin_phi = C.dot(B) / (rC * rB);

		//  Get phi, assign the sign based on the sine value
		return -Math.atan2(sin_phi, cos_phi);
	}*/

	/**
	 * Calculates a signed torsion as an exterior spherical angle
	 * from a valid sequence of 4 points in space.
	 * Looking along the line from c2 to c3, the torsion angle is 0.0, if the
	 * projection of c2->c1 and c3->c4 point in the same direction.
	 * If the projection of vector c2-c1 is rotated in clockwise direction,
	 * the angle increases, i.e. has a positive value.
	 * http://en.wikipedia.org/wiki/Dihedral_angle
	 * @param c1
	 * @param c2
	 * @param c3
	 * @param c4
	 * @return torsion in the range: -pi <= torsion <= pi
	 */
	public static double getDihedral(Coordinates c1, Coordinates c2, Coordinates c3, Coordinates c4) {
		// changed from above, because it seems a little more efficient; TLS 2-Nov-2016
		Coordinates v1 = c2.subC(c1);
		Coordinates v2 = c3.subC(c2);
		Coordinates v3 = c4.subC(c3);

		Coordinates n1 = v1.cross(v2);
		Coordinates n2 = v2.cross(v3);

		return -Math.atan2(v2.getLength() * v1.dot(n2), n1.dot(n2));
	}

	//@Override Annotation incompatible with 1.5
	public int compareTo(Coordinates o) {
		if(x!=o.x) return x<o.x?-1:1;
		if(y!=o.y) return y<o.y?-1:1;
		if(z!=o.z) return z<o.z?-1:1;
		return 0;
	}

	public static Coordinates random() {
		Random random = new Random();
		return new Coordinates(random.nextDouble()*2-1, random.nextDouble()*2-1, random.nextDouble()*2-1);
	}
	public static double getRmsd(Coordinates[] c1, Coordinates[] c2) {
		return getRmsd(c1, c2,Math.min(c1.length, c2.length));
	}

	public static double getRmsd(Coordinates[] c1, Coordinates[] c2, int l) {
		double sum = 0;
		for (int i = 0; i < l; i++) {
			sum+= c1[i].distanceSquared(c2[i]);
		}
		return l>0? Math.sqrt(sum/l): 0;
	}
}