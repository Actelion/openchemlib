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

import java.util.*;

/**
 * Represents a polynomial spline function.
 */
public final class FastSpline {
   
	public final static class Polynome  {
		private final double coeffs[]; //x0 x1 x2 x3 
		
		public Polynome(double[] coeffs) {
			this.coeffs = coeffs;			
		}
		
		public final Polynome derivative() {
			return new Polynome(new double[] {coeffs[1], 2*coeffs[2], 3*coeffs[3], 0});
		}

		public final double value(double x) {
			return coeffs[0]+x*(coeffs[1]+x*(coeffs[2]+x*coeffs[3]));			
		}
		public final double[] getCoefficients() {
			return coeffs;
		}
	}
    
	
    /** Spline segment interval delimiters (knots).   Size is n+1 for n segments. */
    private final double knots[];

    /**
     * The polynomial functions that make up the spline.  The first element
     * determines the value of the spline over the first subinterval, the
     * second over the second, etc.   Spline function values are determined by
     * evaluating these functions at <code>(x - knot[i])</code> where i is the
     * knot segment to which x belongs.
     */
    private final Polynome polynomials[];
    
    /** 
     * Number of spline segments = number of polynomials
     *  = number of partition points - 1 
     */
    private final int n;
    

    /**
     * Construct a polynomial spline function with the given segment delimiters
     * and interpolating polynomials.
     */
    public FastSpline(double knots[], Polynome polynomials[]) {
        this.n = knots.length -1;
        this.knots = new double[n + 1];
        this.polynomials = new Polynome[n];
        
        System.arraycopy(knots, 0, this.knots, 0, n + 1);
        System.arraycopy(polynomials, 0, this.polynomials, 0, n);
    }

    /**
     * Compute the value for the function.
     */
    public final double value(double v) throws ArrayIndexOutOfBoundsException {    	
        int i = Arrays.binarySearch(knots, v);
        if (i < 0) i = -i - 2;
        if (i < 0) i = 0; 
        return polynomials[i].value(v - knots[i]);
    }
        
    /**
     * Returns the derivative of the polynomial spline function as a PolynomialSplineFunction
     */
    public final FastSpline derivative() {
    	Polynome derivativePolynomials[] = new Polynome[n];
        for (int i = 0; i < n; i++) derivativePolynomials[i] = polynomials[i].derivative();
        return new FastSpline(knots, derivativePolynomials);
    }

}
