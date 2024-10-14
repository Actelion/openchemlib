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

package com.actelion.research.chem.descriptor.flexophore.completegraphmatcher;

import com.actelion.research.chem.descriptor.flexophore.DistHist;
import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;

/**
 * 
 * 
 * HistogramMatchCalculator
 * @author Modest von Korff
 * Oct 2, 2012 MvK: Start implementation
 * May 15 2013 MvK: Heavy bug detected. Wrong similarity results. reset() added.
 * Mar 01 2016 MvK sliding filter added.
 * 2020 March, re-written
 */
public class HistogramMatchCalculator {

	/**
	 *
	 * @param dh1
	 * @param indexDistHist1At1
	 * @param indexDistHist1At2
	 * @param dh2
	 * @param indexDistHist2At1
	 * @param indexDistHist2At2
	 * @return integral of overlap
	 */
	public static double getFractionOverlapIntegral(DistHist dh1, int indexDistHist1At1, int indexDistHist1At2, DistHist dh2, int indexDistHist2At1, int indexDistHist2At2){
		int indexPostStartDistHist1 = dh1.getIndexPosStartForDistHist(indexDistHist1At1, indexDistHist1At2);
		int indexPostStartDistHist2 = dh2.getIndexPosStartForDistHist(indexDistHist2At1, indexDistHist2At2);
		int n = ConstantsFlexophoreGenerator.BINS_HISTOGRAM;
		double sumMin = 0;
		double sumMax = 0;

		for (int i = 0; i < n; i++) {
			int v1 = dh1.getValueAtAbsolutePosition(indexPostStartDistHist1+i);
			int v2 = dh2.getValueAtAbsolutePosition(indexPostStartDistHist2+i);
			sumMin += Math.min(v1, v2);
			sumMax += Math.max(v1, v2);
		}
		double score = sumMin / sumMax;
		return score;
	}
	public static double getFractionOverlappingBins(DistHist dh1, int indexDistHist1At1, int indexDistHist1At2, DistHist dh2, int indexDistHist2At1, int indexDistHist2At2){
		int indexPostStartDistHist1 = dh1.getIndexPosStartForDistHist(indexDistHist1At1, indexDistHist1At2);
		int indexPostStartDistHist2 = dh2.getIndexPosStartForDistHist(indexDistHist2At1, indexDistHist2At2);
		int n = ConstantsFlexophoreGenerator.BINS_HISTOGRAM;
		double sumMin = 0;
		double sumMax = 0;

		double sum1=0;
		double sum2=0;
		for (int i = 0; i < n; i++) {
			int v1 = (dh1.getValueAtAbsolutePosition(indexPostStartDistHist1+i)==0)?0:1;
			int v2 = (dh2.getValueAtAbsolutePosition(indexPostStartDistHist2+i)==0)?0:1;
			sum1 += v1;
			sum2 += v2;
			sumMin += Math.min(v1, v2);
			sumMax += Math.max(v1, v2);
		}
		double score = sumMin / Math.min(sum1, sum2);
		return score;
	}

	public static double getPercentageOverlap(DistHist dh1, int indexDistHist1At1, int indexDistHist1At2, DistHist dh2, int indexDistHist2At1, int indexDistHist2At2){
		int indexPostStartDistHist1 = dh1.getIndexPosStartForDistHist(indexDistHist1At1, indexDistHist1At2);
		int indexPostStartDistHist2 = dh2.getIndexPosStartForDistHist(indexDistHist2At1, indexDistHist2At2);
		int n = ConstantsFlexophoreGenerator.BINS_HISTOGRAM;
		double sumMin = 0;
		double sum1 = 0;
		double sum2 = 0;
		for (int i = 0; i < n; i++) {
			int v1 = dh1.getValueAtAbsolutePosition(indexPostStartDistHist1+i);
			int v2 = dh2.getValueAtAbsolutePosition(indexPostStartDistHist2+i);
			sumMin += Math.min(v1, v2);

			sum1 += v1;
			sum2 += v2;
		}
		double score = sumMin / Math.max(sum1, sum2);;
		return score;
	}



}
