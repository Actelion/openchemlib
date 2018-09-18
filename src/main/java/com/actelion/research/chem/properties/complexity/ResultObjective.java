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

package com.actelion.research.chem.properties.complexity;

public class ResultObjective {

	private static final double DEFAULT_FACTOR = -1;
	
	private double score;
	
	private double slope;
	
	private double slopeFragments;
	
	private double slopeR2;
	
	private int nRegressionPoints;
	
	private double factorIsomorphSymmetric;
	
	public ResultObjective() {
		factorIsomorphSymmetric = DEFAULT_FACTOR;
	}

	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @param score the score to set
	 */
	protected void setScore(double score) {
		this.score = score;
	}

	/**
	 * @return the slopeR2
	 */
	protected double getSlopeR2() {
		return slopeR2;
	}

	/**
	 * @param slopeR2 the slopeR2 to set
	 */
	protected void setSlopeR2(double slopeR2) {
		this.slopeR2 = slopeR2;
	}

	/**
	 * @return the slopeFragments
	 */
	public double getSlopeFragments() {
		return slopeFragments;
	}

	/**
	 * @param slopeFragments the slopeFragments to set
	 */
	public void setSlopeFragments(double slopeFragments) {
		this.slopeFragments = slopeFragments;
	}

	/**
	 * @return the slope
	 */
	public double getSlope() {
		return slope;
	}

	/**
	 * @param slope the slope to set
	 */
	protected void setSlope(double slope) {
		this.slope = slope;
	}

	/**
	 * @return the factorIsomorphSymmetric
	 */
	public double getFactorIsomorphSymmetric() {
		return factorIsomorphSymmetric;
	}

	/**
	 * @param factorIsomorphSymmetric the factorIsomorphSymmetric to set
	 */
	protected void setFactorIsomorphSymmetric(double factorIsomorphSymmetric) {
		this.factorIsomorphSymmetric = factorIsomorphSymmetric;
	}

	/**
	 * @return the nRegressionPoints
	 */
	public int getNumRegressionPoints() {
		return nRegressionPoints;
	}

	/**
	 * @param nRegressionPoints the nRegressionPoints to set
	 */
	protected void setNumRegressionPoints(int nRegressionPoints) {
		this.nRegressionPoints = nRegressionPoints;
	}

	
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("ResultObjective [score=");
		builder.append(score);
		builder.append(", slope=");
		builder.append(slope);
		builder.append(", slopeFrags=");
		builder.append(slopeFragments);
		builder.append(", points=");
		builder.append(nRegressionPoints);
		builder.append(", fac isomorph symm=");
		builder.append(factorIsomorphSymmetric);
		builder.append("]");
		return builder.toString();
	}
	
	

}
