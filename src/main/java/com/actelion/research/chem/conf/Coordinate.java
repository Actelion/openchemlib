package com.actelion.research.chem.conf;

public class Coordinate {
	public double x,y,z;

	public Coordinate() {
		}

	public Coordinate(Coordinate c) {
		this.x = c.x;
		this.y = c.y;
		this.z = c.z;
		}

	public Coordinate(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
		}

	/**
	 * Copies all components from c to this.
	 * @param c
	 * @return this after after copying c
	 */
	public Coordinate copyFrom(Coordinate c) {
		this.x = c.x;
		this.y = c.y;
		this.z = c.z;
		return this;
		}

	public double getX() {
		return x;
		}

	public double getY() {
		return y;
	}

	public double getZ() {
		return z;
	}

	public void setZ(double z) {
		this.z = z;
		}

	public void setX(double x) {
		this.x = x;
		}

	public void setY(double y) {
		this.y = y;
		}

	/**
	 * @param dx
	 * @param dy
	 * @param dz
	 * @return this after translating in three dimensions
	 */
	public Coordinate translate(double dx, double dy, double dz) {
		x += dx;
		y += dy;
		z += dz;
		return this;
		}

	public void translateX(double dx) {
		x += dx;
		}

	public void translateY(double dy) {
		y += dy;
		}

	public void translateZ(double dz) {
		z += dz;
		}

	public double distanceTo(Coordinate c) {
		double dx = x - c.x;
		double dy = y - c.y;
		double dz = z - c.z;
		return Math.sqrt(dx*dx + dy*dy + dz*dz);
		}

	/**
	 * Adds the components of c to this and returns this.
	 * @param c vector to add
	 * @return this after addition
	 */
	public Coordinate add(Coordinate c) {
		x += c.x;
		y += c.y;
		z += c.z;
		return this;
		}

	/**
	 * Subtracts the components of c from this and returns this.
	 * @param c vector to subtract
	 * @return this after subtraction
	 */
	public Coordinate subtract(Coordinate c) {
		x -= c.x;
		y -= c.y;
		z -= c.z;
		return this;
		}

	/**
	 * Multiplies the components of this with factor f.
	 * @param f factor
	 * @return this after multiplication
	 */
	public Coordinate multiply(double f) {
		x *= f;
		y *= f;
		z *= f;
		return this;
		}

	/**
	 * Considering this and the passed Coordinate as vectors,
	 * this method returns the angle between them.
	 * @param v
	 * @return
	 */
	public double getAngle(Coordinate v) {
		return Math.acos(dotProduct(v)/(getLength()*v.getLength()));
		}

	/**
	 * Calculates the angle of the line from this location to c
	 * projected into the x/y plane. With Y facing upwards and X right,
	 * if the line points in Y direction, ten angle is 0.0 increasing
	 * in clockwise direction.
	 * @param c
	 * @return -PI < angle < PI
	 */
	public double getAngleXY(Coordinate c) {
		double dx = c.x - x;
		double dy = c.y - y;

		if (dy == 0.0)
			return (dx > 0.0) ? Math.PI/2.0 : -Math.PI/2.0;

		double angle = Math.atan(dx/dy);
		if (dy < 0.0)
			return (dx < 0.0) ? angle - Math.PI : angle + Math.PI;

		return angle;
		}

	/**
	 * @return distance between this coordinate and P(0,0,0)
	 */
	public double getLength() {
		return Math.sqrt(x*x + y*y + z*z);
	}

	/**
	 * Calculates the center point between this and c and sets this to the center point.
	 * @param c
	 * @return this after updating it to the center position
	 */
	public Coordinate center(Coordinate c) {
		x = (x + c.x) / 2.0;
		y = (y + c.y) / 2.0;
		z = (z + c.z) / 2.0;
		return this;
		}

	/**
	 * Updates this to contains the center between c1 and c2.
	 * @param c1
	 * @param c2
	 * @return this after calculating the center
	 */
	public Coordinate center(Coordinate c1, Coordinate c2) {
		x = (c1.x + c2.x) / 2.0;
		y = (c1.y + c2.y) / 2.0;
		z = (c1.z + c2.z) / 2.0;
		return this;
		}

	/**
	 * Updates this to contain a point on the straight line through c1 and c2.
	 * @param c1
	 * @param c2
	 * @param f location on line 0.0 -> c1, 1.0 -> c2
	 * @return this after updating to be a point on the line
	 */
	public Coordinate between(Coordinate c1, Coordinate c2, double f) {
		x = c1.x + f * (c2.x - c1.x);
		y = c1.y + f * (c2.y - c1.y);
		z = c1.z + f * (c2.z - c1.z);
		return this;
		}

	/**
	 * Calculates and returns the dot product between this and c.
	 * @param c
	 * @return
	 */
	public double dotProduct(Coordinate c) {
		return x*c.x + y*c.y + z*c.z;
	}

	/**
	 * Calculates and returns the cross product between this and c
	 * and returns it as a new vector.
	 * @param c
	 * @return
	 */
	public Coordinate crossProduct(Coordinate c) {
		return new Coordinate(y*c.z - z*c.y, z*c.x - x*c.z, x*c.y - y*c.x);
		}

	/**
	 * Converts this vector to the unit vector (length==1.0) in place.
	 * @return this vector after changing its length to 1.0
	 */
	public Coordinate normalize() {
		assert(x!=0 || y!=0 || z!=0);

		double d = Math.sqrt(x*x + y*y + z*z);

		x /= d;
		y /= d;
		z /= d;
		return this;
		}
	}
