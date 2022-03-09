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
 */

package com.actelion.research.calc.classification;

import com.actelion.research.util.Formatter;

import java.io.Serializable;
import java.util.List;

/**
 * PrecisionAndRecall
 * @author Modest von Korff
 * @version 1.0
 * Jun 18, 2015 MvK Start implementation
 */
public class PrecisionAndRecall implements Serializable {

	public static final String ATTR_TP = "TP";
	public static final String ATTR_TN = "TN";
	public static final String ATTR_FP = "FP";
	public static final String ATTR_FN = "FN";


	public int truePositive=0;
	public int trueNegative=0;
	public int falsePositive=0;
	public int falseNegative=0;
	
	/**
	 * @param truePositive
	 * @param trueNegative
	 * @param falsePositive
	 * @param falseNegative
	 */
	public PrecisionAndRecall(int truePositive, int trueNegative,
			int falsePositive, int falseNegative) {
		super();
		this.truePositive = truePositive;
		this.trueNegative = trueNegative;
		this.falsePositive = falsePositive;
		this.falseNegative = falseNegative;
	}
	/**
	 * 
	 */
	public PrecisionAndRecall() {
	}

	public void add(PrecisionAndRecall p){
		this.truePositive+=p.truePositive;
		this.trueNegative+=p.trueNegative;
		this.falsePositive+=p.falsePositive;
		this.falseNegative+=p.falseNegative;
	}

	public void parse2PrecisionAndRecall(String sVal){

		if(PrecisionAndRecall.ATTR_TP.equals(sVal)) {
			truePositive++;
		} else if(PrecisionAndRecall.ATTR_TN.equals(sVal)) {
			trueNegative++;
		} else if(PrecisionAndRecall.ATTR_FP.equals(sVal)) {
			falsePositive++;
		} else if(PrecisionAndRecall.ATTR_FN.equals(sVal)) {
			falseNegative++;
		} else {
			throw new RuntimeException("Parsing error for '" + sVal + "'.");
		}
	}



	/**
	 * @return the truePositive
	 */
	public int getTruePositive() {
		return truePositive;
	}

	/**
	 * @param truePositive the truePositive to set
	 */
	public void setTruePositive(int truePositive) {
		this.truePositive = truePositive;
	}

	/**
	 * @return the trueNegative
	 */
	public int getTrueNegative() {
		return trueNegative;
	}

	/**
	 * @param trueNegative the trueNegative to set
	 */
	public void setTrueNegative(int trueNegative) {
		this.trueNegative = trueNegative;
	}

	/**
	 * @return the falsePositive
	 */
	public int getFalsePositive() {
		return falsePositive;
	}

	/**
	 * @param falsePositive the falsePositive to set
	 */
	public void setFalsePositive(int falsePositive) {
		this.falsePositive = falsePositive;
	}

	/**
	 * @return the falseNegative
	 */
	public int getFalseNegative() {
		return falseNegative;
	}

	/**
	 * @param falseNegative the falseNegative to set
	 */
	public void setFalseNegative(int falseNegative) {
		this.falseNegative = falseNegative;
	}

	public int getSum(){
		return trueNegative+truePositive+falsePositive+falseNegative;
	}

	public double calculatePrecision(){
		
		int tp_fp = truePositive+falsePositive;
		
		if(tp_fp==0){
			return 0;
		}
		
		return truePositive/(double)(truePositive+falsePositive);
		
	}
	
	public double calculateRecall(){
				
		int tp_fn = truePositive+falseNegative;
		
		if(tp_fn==0){
			return 0;
		}
		
		return truePositive/(double)(truePositive+falseNegative);
		
	}
	public double calculateAccuracy(){

		double acc = (truePositive+trueNegative)/(double)getSum();

		return acc;

	}
	
	/**
	 * F1 score
	 * @return
	 */
	public double calculateHarmonicMean(){
		
		double p = calculatePrecision();
		
		double r = calculateRecall();
		
		if((p+r)==0){
			return 0;
		}
		
		double f1 = 2.0 * (p*r/(p+r));
		
		return f1;
		
	}

	public double calculateCohensKappa(){

		double all = getSum();

		double pYes = ((truePositive+falseNegative)/all) * ((truePositive+falsePositive)/all);

		double pNo = ((falsePositive+trueNegative)/all) * ((falseNegative+trueNegative)/all);

		double pe = pYes+pNo;

		double po = calculateAccuracy();

		double k = (po-pe)/(1.0-pe);

		return k;
	}

	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		
		StringBuilder sb = new StringBuilder();
		
		sb.append("Accuracy " + Formatter.format3(calculateAccuracy()));
		sb.append(", precision " + Formatter.format3(calculatePrecision()));
		sb.append(", recall " + Formatter.format3(calculateRecall()));
		sb.append(", F1 "+ Formatter.format3(calculateHarmonicMean()));
		sb.append(", Cohens' kappa " + Formatter.format3(calculateCohensKappa()));
		sb.append("\n");
		sb.append("True positive " + truePositive);
		sb.append("\n");
		sb.append("True negative " + trueNegative);
		sb.append("\n");
		sb.append("False positive " + falsePositive);
		sb.append("\n");
		sb.append("False negative " + falseNegative);
		
		return sb.toString();
	}


	public static double getHarmonicMean(List<PrecisionAndRecall> li){

		double sum=0;
		for (PrecisionAndRecall precisionAndRecall : li) {

			sum += precisionAndRecall.calculateHarmonicMean();
		}

		return sum / li.size();
	}

	public static void main(String[] args) {

		PrecisionAndRecall pr = new PrecisionAndRecall();

		pr.truePositive=20;
		pr.trueNegative=15;
		pr.falseNegative=5;
		pr.falsePositive=10;

		System.out.println("Cohen's kappa: " + pr.calculateCohensKappa());

	}

	
}
