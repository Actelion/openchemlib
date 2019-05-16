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

package com.actelion.research.chem.forcefield.mmff;

import java.lang.Math;
import java.lang.Double;

import com.actelion.research.chem.ExtendedMolecule;

/**
 * The Vector3 class provides common vector operations used throughout the
 * MMFF codebase. This is not supposed to be a comprehensive vector
 * library;  instead, only those functions which were needed have been
 * implemented.
 */
public class Vector3 {
    private final static double TOL = 0.01;

    public final double x;
    public final double y;
    public final double z;

    /**
     * Constructs a new vector with default coordinates.
     */
    public Vector3() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }

    /**
     * Constructs a new vector with given x, y and z coordinates.
     *  @param x X coordinate.
     *  @param y Y coordinate.
     *  @param z Z coordinate.
     */
    public Vector3(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    /**
     * Copy constructor, copies the x,y,z coordinates from another vector.
     */
    public Vector3(Vector3 that) {
        x = that.x;
        y = that.y;
        z = that.z;
    }

    /**
     * Constructs a new vector from the x, y, z coordinates of an atom in
     * a molecule.
     *  @param mol The molecule that the atom is in.
     *  @param atom The atom index.
     */
    public Vector3(MMFFMolecule mol, int atom) {
        this.x = mol.getAtomX(atom);
        this.y = mol.getAtomY(atom);
        this.z = mol.getAtomZ(atom);
    }

    /**
     * Constructs a new vector from the x, y, z coordinates of an atom in
     * a molecule.
     *  @param mol The molecule that the atom is in.
     *  @param atom The atom index.
     */
    public Vector3(ExtendedMolecule mol, int atom) {
        this.x = mol.getAtomX(atom);
        this.y = mol.getAtomY(atom);
        this.z = mol.getAtomZ(atom);
    }

    /**
     * Constructs a new vector from three consecutive doubles in a
     * positions array.
     *  @param pos The double array of position coordinates.
     *  @param offset Where in the array to fetch the values from.
     */
    public Vector3(double[] pos, int offset) {
        x = pos[3*offset    ];
        y = pos[3*offset + 1];
        z = pos[3*offset + 2];
    }

    /**
     * Constructs a new vector starting at atom1 position and ending at
     * atom2 position. The vector is constructed using point information
     * from a position array.
     */
    public Vector3(double[] pos, int atom1, int atom2) {
        this(new Vector3(pos, atom1).to(new Vector3(pos, atom2)));
    }

    /**
     * Constructs a vector that represents the line travelling from the
     * position of atom1 to the position of atom2. The resultant vector
     * points towards atom2 from atom1.
     *  @param mol The molecule that the atoms are in.
     *  @param atom1 The first atom and starting point of the line.
     *  @param atom2 The second atom and the end point of the line.
     */
    public Vector3(MMFFMolecule mol, int atom1, int atom2) {
        this(new Vector3(mol, atom1).to(new Vector3(mol, atom2)));
    }

    /**
     * Constructs a vector that represents the line travelling from the
     * position of atom1 to the position of atom2. The resultant vector
     * points towards atom2 from atom1.
     *  @param mol The molecule that the atoms are in.
     *  @param atom1 The first atom and starting point of the line.
     *  @param atom2 The second atom and the end point of the line.
     */
    public Vector3(ExtendedMolecule mol, int atom1, int atom2) {
        this(new Vector3(mol, atom1).to(new Vector3(mol, atom2)));
    }

    /**
     * Checks for equality between two vectors. Vectors are considered
     * equal if each component value is equal to its corresponding
     * component value. Doubleing point values are considered equal if
     * their difference is below some tolerance value (set in TOL).
     *  @param obj The vector to compare this vector with.
     *  @return True if they are equal, false otherwise.
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == this)
            return true;

        if (!(obj instanceof Vector3))
            return false;

        Vector3 that = (Vector3)obj;

        return Math.abs(x-that.x) < TOL
            && Math.abs(y-that.y) < TOL
            && Math.abs(z-that.z) < TOL;
    }

    /**
     * Returns the hashcode of this vector.
     *  @return The hashcode.
     */
    @Override
    public int hashCode() {
        return new Double(x).hashCode()
            + new Double(y).hashCode()
            + new Double(z).hashCode();
    }

    /**
     * Writes the contents of this vector to the specified positions in a
     * positions array.
     *  @param pos The positions array to write to.
     *  @param offset Where in the array to write to.
     */
    public void write(double[] pos, int offset) {
        pos[3*offset    ] = x;
        pos[3*offset + 1] = y;
        pos[3*offset + 2] = z;
    }

    /**
     * Returns this vector negated (signs swapped).
     *  @return The negated vector.
     */
    public Vector3 negate() {
        return new Vector3(-x, -y, -z);
    }

    /**
     * Returns a vector that points from this vector to that vector. The
     * resultant vector will have a length equal to the distance between
     * the two vectors and will have a direction pointing towards that
     * vector.
     */
    public Vector3 to(Vector3 that) {
        return new Vector3(that.x-x, that.y-y, that.z-z);
    }

    /**
     * Returns the length (or magnitude) of this vector. This is the
     * euclidean distance from this vector to the origin (A vector at
     * x=0, y=0, z=0).
     *  @return The length.
     */
    public double length() {
        return Math.sqrt(x*x+y*y+z*z);
    }

    /**
     * Normalises a vector. If the vector to be normalised has length
     * zero, then the zero vector is returned.
     *  @return This vector normalised.
     */
    public Vector3 normalise() {
        if (length() > 0.0)
            return new Vector3(x/length(), y/length(), z/length());
        return new Vector3(0.0, 0.0, 0.0);
    }

    /**
     * Adds two vectors together by component.
     *  @param that The other vector.
     *  @return The component sum of the two vectors.
     */
    public Vector3 add(Vector3 that) {
        return new Vector3(x+that.x, y+that.y, z+that.z);
    }

    /**
     * Subtracts that vector from this vector by component.
     *  @param that The other vector.
     *  @return The component difference of the two vectors.
     */
    public Vector3 sub(Vector3 that) {
        return new Vector3(x-that.x, y-that.y, z-that.z);
    }

    /**
     * Computes the euclidean distance between two Vector3s.
     *  @param that The other Vector3 to calculate the distance to.
     *  @return The distance between the two vectors.
     */
    public double distance(Vector3 that) {
        return Math.sqrt(
                  (x-that.x)*(x-that.x)
                + (y-that.y)*(y-that.y)
                + (z-that.z)*(z-that.z) );
    }

    /**
     * Computes the vector dot product between this vector and that
     * vector.
     *  @param that The other vector.
     *  @return The dot product of the two vectors.
     */
    public double dot(Vector3 that) {
        return x*that.x + y*that.y + z*that.z;
    }

    /**
     * Computes the vector cross product between this vector and that
     * vector.
     *  @param that The other vector.
     *  @return The cross product of the two vectors.
     */
    public Vector3 cross(Vector3 that) {
        return new Vector3(y*that.z - z*that.y,
                z*that.x - x*that.z,
                x*that.y - y*that.x);
    }

    /**
     * Returns the angle between this vector and that vector.
     *  @param that The other vector.
     *  @return The angle between the two vectors in radians.
     */
    public double angle(Vector3 that) {
        return Math.acos(dot(that) / (length() * that.length()));
    }

    /**
     * Returns the cosine of the angle between this vector and that
     * vector.
     *  @param that The other vector.
     *  @return The cosine of the angle between the two vectors in
     *     radians.
     */
    public double cosAngle(Vector3 that) {
        return dot(that) / (length() * that.length());
    }

    /**
     * Returns a string form of this vector.
     *  @return The string form of the vector.
     */
    public String toString() {
        return "("+String.format("%.3f", x)+","+String.format("%.3f", y)+","
            +String.format("%.3f", z)+")";
    }
}
