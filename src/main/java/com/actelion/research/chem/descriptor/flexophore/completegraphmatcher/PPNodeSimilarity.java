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

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.statistics.median.MedianStatisticFunctions;
import com.actelion.research.chem.descriptor.flexophore.IPPNode;
import com.actelion.research.chem.descriptor.flexophore.PPNode;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.interactionstatistics.InteractionDistanceStatistics;
import com.actelion.research.chem.interactionstatistics.InteractionSimilarityTable;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.datamodel.table.TableModelString;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * 
 * PPNodeSimilarity
 * @author Modest von Korff
 * @version 1.0
 * Jan 7, 2013 MvK Start implementation
 * Dec 2020, MvK unskewed similarity calculation for similarity hard thresh. Node similarity is now independent from
 * query base order.
 */
public class PPNodeSimilarity implements IPPNodeSimilarity {

	public static final int SIMILARITY_MODE_SIMPLE = 0;

	/**
	 * Hard thresh means if the similarity of two atom types in a node to node comparison is below a threshold
	 * the complete node to node similarity becomes 0.
	 */

	public static final int SIMILARITY_MODE_HARD_THRESH = 1;
	public static final int SIMILARITY_MODE_HARD_THRESH_AVR = 2;
	public static final int SIMILARITY_MODE_HARD_THRESH_OPTIMISTIC = 3;
	public static final int SIMILARITY_MODE_CARBON = 4;


	/**
	 * Similarity value distribution of the interaction table (April 2020). First two rows bin borders. Third row counts.
	 * .900	.905	.910	.915	.920	.925	.930	.935	.940	.945	.950	.955	.960	.965	.970	.975	.980	.985	.990	.995
	 * .905	.910	.915	.920	.925	.930	.935	.940	.945	.950	.955	.960	.965	.970	.975	.980	.985	.990	.995	1.000
	 *  191	 239	 202	 280	 213	 207	 193	 165	 229	 221	 157	 139	 131	  79	  74	  23	  17	   5	   2	  56

	 A similarity thresh of 0.99 allows 58 pairwise interactions.
	 */
	public static final double THRESH_SIMILARITY_HARD_MATCH = 0.9;

	public static final double HARD_MATCH_OPTIMISTIC_PERCENTILE = 0.75;

	private static final double TINY = 0.0000001;
	private static final int SIZE_SIM_MATRIX = 20;


	private static final double THRESH_CARBON_INTERACTIONS = 0.6;

	private static PPNodeSimilarity INSTANCE = null;


	private InteractionSimilarityTable interactionSimilarityTable;

	private Matrix maSimilarity;

	private int similarityMode;

	private double threshSimilarityHardMatch;

	private boolean verbose;
	
	/**
	 * This constructor is used for parallel mode.
	 */
	public PPNodeSimilarity(int versionInteractionTable, int modePPNodeSimilarity){

		maSimilarity = new Matrix(SIZE_SIM_MATRIX, SIZE_SIM_MATRIX);

		interactionSimilarityTable = InteractionSimilarityTable.getInstance();

		similarityMode = modePPNodeSimilarity;

		threshSimilarityHardMatch = THRESH_SIMILARITY_HARD_MATCH;
	}

	public void setThreshSimilarityHardMatch(double threshSimilarityHardMatch) {
		this.threshSimilarityHardMatch = threshSimilarityHardMatch;
	}


	public String toStringParameter(){
		StringBuilder sb = new StringBuilder();

		sb.append("ObjectiveFlexophoreHardMatchUncovered, similarity mode ");
		switch (similarityMode){
			case SIMILARITY_MODE_SIMPLE:
				sb.append("simple");
				break;
			case SIMILARITY_MODE_HARD_THRESH:
				sb.append("hard thresh multiplicative");
				sb.append(", threshold=");
				sb.append(threshSimilarityHardMatch);
				break;
			case SIMILARITY_MODE_HARD_THRESH_AVR:
				sb.append("hard thresh average");
				sb.append(", threshold=");
				sb.append(threshSimilarityHardMatch);
				break;
			case SIMILARITY_MODE_HARD_THRESH_OPTIMISTIC:
				sb.append("hard thresh optimistic");
				sb.append(", percentile=");
				sb.append(HARD_MATCH_OPTIMISTIC_PERCENTILE);
				sb.append(", threshold=");
				sb.append(threshSimilarityHardMatch);
				break;
			case SIMILARITY_MODE_CARBON:
				sb.append("carbon");
				break;
		}

		return sb.toString();
	}

	/**
	 * Use this as constructor for serial mode.
	 * @return
	 */
	public static PPNodeSimilarity getInstance(int versionInteractionTable, int modePPNodeSimilarity){
		if(INSTANCE == null) {
			synchronized(PPNodeSimilarity.class) {
				INSTANCE = new PPNodeSimilarity(versionInteractionTable, modePPNodeSimilarity);
			}
		}
		return INSTANCE;
	}

	public void setVerbose(boolean v) {
		this.verbose = v;
	}

	public double getSimilarity(IPPNode query, IPPNode base) {

		double similarity = 0;

		switch (similarityMode){
			case SIMILARITY_MODE_SIMPLE:
				similarity = getSimilaritySimple((PPNode)query, (PPNode)base);
				break;
			case SIMILARITY_MODE_HARD_THRESH:
				similarity = getSimilarityHardMatchMultiplicative((PPNode)query, (PPNode)base);
				break;
			case SIMILARITY_MODE_HARD_THRESH_AVR:
				similarity = getSimilarityHardMatchAverage((PPNode)query, (PPNode)base);
				break;
			case SIMILARITY_MODE_HARD_THRESH_OPTIMISTIC:
				similarity = getSimilarityHardMatchOptimistic((PPNode)query, (PPNode)base);
				break;
			case SIMILARITY_MODE_CARBON:
				similarity = getSimilarityExtraCarbonConsideration((PPNode)query, (PPNode)base);
				break;
		}

		return similarity;
	}


		/**
         *
         * @param query
         * @param base
         * @return
         * @throws Exception
         */
	public double getSimilaritySimple(PPNode query, PPNode base) {

		maSimilarity.set(0);

		for (int i = 0; i < query.getInteractionTypeCount(); i++) {

			int interactionTypeQuery = query.getInteractionType(i);

			for (int j = 0; j < base.getInteractionTypeCount(); j++) {

				int interactionTypeBase = base.getInteractionType(j);

				double similarity = 1.0 - interactionSimilarityTable.getDistance(interactionTypeQuery, interactionTypeBase);

				maSimilarity.set(i,j,similarity);
			}
		}

		if(verbose) {
			System.out.println("PPNodeSimilarityMultiplicative");

			TableModelString tableModelString = new TableModelString(query.getInteractionTypeCount(), base.getInteractionTypeCount());

			for (int i = 0; i < query.getInteractionTypeCount(); i++) {
				int interactionType = query.getInteractionType(i);
				String s = InteractionAtomTypeCalculator.getString(interactionType);
				tableModelString.setRowName(i, s);
			}

			for (int i = 0; i < base.getInteractionTypeCount(); i++) {
				int interactionType = base.getInteractionType(i);
				String s = InteractionAtomTypeCalculator.getString(interactionType);
				tableModelString.setColName(i, s);
			}

			tableModelString.set(maSimilarity, 2);

			System.out.println(tableModelString.toString());

		}


		List<Double> liSimilarities = new ArrayList<>();

		if(base.getInteractionTypeCount() > query.getInteractionTypeCount()) {

			for (int col = 0; col < base.getInteractionTypeCount(); col++) {

				double maxSimInCol = 0;
				for (int row = 0; row < query.getInteractionTypeCount(); row++) {
					if(maSimilarity.get(row,col)>maxSimInCol){
						maxSimInCol = maSimilarity.get(row,col);
					}
				}

				liSimilarities.add(maxSimInCol);

				// System.out.println("Sim maxSimInCol " + Formatter.format2(maxSimInCol) + "\t" + InteractionAtomTypeCalculator.getString(interactionTypeBase) + "\t" + InteractionAtomTypeCalculator.getString(interactionTypeQuery));
			}
		} else {
			for (int row = 0; row < query.getInteractionTypeCount(); row++) {

				double maxSimInRow = 0;
				for (int col = 0; col < base.getInteractionTypeCount(); col++) {
					if(maSimilarity.get(row,col) > maxSimInRow){
						maxSimInRow = maSimilarity.get(row,col);
					}
				}
				// System.out.println("Sim maxSimInRow " + Formatter.format2(maxSimInRow) + "\t" + InteractionAtomTypeCalculator.getString(interactionTypeBase) + "\t" + InteractionAtomTypeCalculator.getString(interactionTypeQuery));

				liSimilarities.add(maxSimInRow);							}
		}

		double sim = 0;

		if(liSimilarities.size()>0) {
			sim = 1;
			for (Double simPart : liSimilarities) {
				sim *= simPart;
			}
		}


		if(verbose) {
			System.out.println("Sim " + Formatter.format2(sim));
			System.out.println();
		}

		return sim;
	}
	public double getSimilarityHardMatchMultiplicative(PPNode query, PPNode base) {

		IPPNode queryCmp = query;

		PPNode baseCmp = base;

		if(base.hasHeteroAtom() && query.hasHeteroAtom()){

			queryCmp = PPNode.getHeteroOnlyNode(query);

			baseCmp = PPNode.getHeteroOnlyNode(base);

		}

		List<Double> liSimilarities = getSimilarityList((PPNode)queryCmp, baseCmp);

		double sim = 0;

		if(liSimilarities.size()>0) {
			sim = 1;
			for (double simPart : liSimilarities) {
				if(simPart< threshSimilarityHardMatch){
					sim = 0;
					break;
				}
				sim *= simPart;
			}
		}

		if(verbose) {
			System.out.println("Sim " + Formatter.format2(sim));
			System.out.println();
		}

		return sim;
	}

	public List<Double> getSimilarityList(PPNode query, PPNode base) {

		maSimilarity.set(0);

		for (int i = 0; i < query.getInteractionTypeCount(); i++) {
			int interactionTypeQuery = query.getInteractionType(i);
			for (int j = 0; j < base.getInteractionTypeCount(); j++) {
				int interactionTypeBase = base.getInteractionType(j);
				try {
					double similarity = 1.0 - interactionSimilarityTable.getDistance(interactionTypeQuery, interactionTypeBase);
					maSimilarity.set(i,j,similarity);
				} catch (Exception e) {
					System.err.println("Error in PPNodeSimilarity");
					System.err.println("interactionTypeQuery " + interactionTypeQuery);
					System.err.println("interactionTypeBase " + interactionTypeBase);
					throw new RuntimeException(e);
				}
			}
		}

		if(verbose) {
			System.out.println("PPNodeSimilarity");
			TableModelString tableModelString = new TableModelString(query.getInteractionTypeCount(), base.getInteractionTypeCount());
			for (int i = 0; i < query.getInteractionTypeCount(); i++) {
				int interactionType = query.getInteractionType(i);
				String s = InteractionAtomTypeCalculator.getString(interactionType);
				tableModelString.setRowName(i, s);
			}
			for (int i = 0; i < base.getInteractionTypeCount(); i++) {
				int interactionType = base.getInteractionType(i);
				String s = InteractionAtomTypeCalculator.getString(interactionType);
				tableModelString.setColName(i, s);
			}
			tableModelString.set(maSimilarity, 2);
			System.out.println(tableModelString.toString());
		}

		double [] arrMaxSim = getTopValues(maSimilarity, query.getInteractionTypeCount(), base.getInteractionTypeCount(), threshSimilarityHardMatch);

		List<Double> liSimilarities = new ArrayList<>(arrMaxSim.length);
		for (double v : arrMaxSim) {
			liSimilarities.add(v);
		}

		return liSimilarities;
	}


	private static double [] getTopValues(Matrix maSimilarity, int rows, int cols, double thresh){

		double [] arrTopSim = new double[Math.max(rows, cols)];

		if(rows==1 && cols==1){
			arrTopSim[0]=maSimilarity.get(0,0);
			return arrTopSim;
		}

		if(cols > rows) {

			for (int col = 0; col < cols; col++) {

				double maxSimInCol = 0;
				for (int row = 0; row < rows; row++) {
					if(maSimilarity.get(row,col)>maxSimInCol){
						maxSimInCol = maSimilarity.get(row,col);
					}
				}
				arrTopSim[col]=maxSimInCol;
			}
		} else if(cols < rows){
			for (int row = 0; row < rows; row++) {
				double maxSimInRow = 0;
				for (int col = 0; col < cols; col++) {
					if(maSimilarity.get(row,col) > maxSimInRow){
						maxSimInRow = maSimilarity.get(row,col);
					}
				}
				arrTopSim[row]=maxSimInRow;
			}
		} else {
			// For rows=cols we search for minimum values in rows and in cols. We take the best result and
			// exclude results that contain a similarity value below the threshold.

			boolean invalidCol=false;
			double sumCol=0;

			double [] arrTopSimCol = new double[cols];

			for (int col = 0; col < cols; col++) {
				double maxSimInCol = 0;
				for (int row = 0; row < rows; row++) {
					double v = maSimilarity.get(row,col);
					if(v>maxSimInCol){
						maxSimInCol = v;
					}
				}
				arrTopSimCol[col]=maxSimInCol;
				sumCol += maxSimInCol;
				if(maxSimInCol<thresh){
					invalidCol=true;
				}
			}

			double [] arrTopSimRow = new double[rows];
			boolean invalidRow=false;
			double sumRow=0;
			for (int row = 0; row < rows; row++) {
				double maxSimInRow = 0;
				for (int col = 0; col < cols; col++) {

					double v = maSimilarity.get(row,col);
					if(v > maxSimInRow){
						maxSimInRow = v;
					}
				}
				arrTopSimRow[row]=maxSimInRow;
				sumRow += maxSimInRow;
				if(maxSimInRow<thresh){
					invalidRow=true;
				}
			}

			if(invalidCol && invalidRow){
				if(sumCol>sumRow){
					arrTopSim=arrTopSimCol;
				} else {
					arrTopSim=arrTopSimRow;
				}
			} else if(invalidCol){
				arrTopSim=arrTopSimRow;
			} else if(invalidRow){
				arrTopSim=arrTopSimCol;
			} else{
				if(sumCol>sumRow){
					arrTopSim=arrTopSimCol;
				} else {
					arrTopSim=arrTopSimRow;
				}
			}
		}

		return arrTopSim;
	}


	public double getSimilarityHardMatchAverage(PPNode query, PPNode base) {
		List<Double> liSimilarities = getSimilarityList(query, base);

		double sumSim = 0;

		if(liSimilarities.size()>0) {

			for (double simPart : liSimilarities) {
				if(simPart< threshSimilarityHardMatch){
					sumSim = 0;
					break;
				}
				sumSim += simPart;
			}
		}

		double sim = sumSim/liSimilarities.size();

		if(verbose) {
			System.out.println("Sim " + Formatter.format2(sim));
			System.out.println();
		}

		return sim;
	}

	public double getSimilarityHardMatchOptimistic(PPNode query, PPNode base) {
		List<Double> liSimilarities = getSimilarityList(query, base);

		double sim = 0;

		if(liSimilarities.size()>0) {
			if(liSimilarities.size()==1) {
				sim = liSimilarities.get(0);
			} else {
				Collections.sort(liSimilarities);

				sim = MedianStatisticFunctions.getPercentileFromSorted(liSimilarities, HARD_MATCH_OPTIMISTIC_PERCENTILE);
			}
		}

		if(verbose) {
			System.out.println("Sim " + Formatter.format2(sim));
			System.out.println();
		}

		return sim;
	}

	public double getSimilarityExtraCarbonConsideration(PPNode query, PPNode base) {

		maSimilarity.set(0);

		final double valNoInteraction = -1;

		boolean lowCarbonFractionQuery = (query.getFractionCarbonInteractions()< THRESH_CARBON_INTERACTIONS)?true:false;
		boolean lowCarbonFractionBase = (base.getFractionCarbonInteractions()< THRESH_CARBON_INTERACTIONS)?true:false;

//		for (int i = 0; i < query.getInteractionTypeCount(); i++) {
//			int interactionTypeQuery = query.getInteractionId(i);
//			if(InteractionAtomTypeCalculator.getAtomicNumber(interactionTypeQuery)==16){
//				System.out.println("Sulfur interaction");
//			}
//		}


		for (int i = 0; i < query.getInteractionTypeCount(); i++) {

			int interactionTypeQuery = query.getInteractionType(i);
			boolean carbonInteractionQuery = InteractionAtomTypeCalculator.isCarbonInteraction(interactionTypeQuery);
			if(carbonInteractionQuery){
				if(lowCarbonFractionQuery) {
					for (int j = 0; j < base.getInteractionTypeCount(); j++) {
						maSimilarity.set(i,j,valNoInteraction);
					}
					continue;
				}
			}

			for (int j = 0; j < base.getInteractionTypeCount(); j++) {

				int interactionTypeBase = base.getInteractionType(j);
				boolean carbonInteractionBase = InteractionAtomTypeCalculator.isCarbonInteraction(interactionTypeBase);

				if(carbonInteractionBase){
					if(lowCarbonFractionBase) {
						maSimilarity.set(i,j,valNoInteraction);
						continue;
					}
				}

				double similarity = valNoInteraction;
				if(carbonInteractionBase == carbonInteractionQuery){
					similarity = 1.0 - interactionSimilarityTable.getDistance(interactionTypeQuery, interactionTypeBase);
				}

				maSimilarity.set(i,j,similarity);
			}
		}

		if(verbose) {
			System.out.println("PPNodeSimilarityMultiplicative");

			TableModelString tableModelString = new TableModelString(query.getInteractionTypeCount(), base.getInteractionTypeCount());

			for (int i = 0; i < query.getInteractionTypeCount(); i++) {
				int interactionType = query.getInteractionType(i);
				String s = InteractionAtomTypeCalculator.getString(interactionType);
				tableModelString.setRowName(i, s);
			}

			for (int i = 0; i < base.getInteractionTypeCount(); i++) {
				int interactionType = base.getInteractionType(i);
				String s = InteractionAtomTypeCalculator.getString(interactionType);
				tableModelString.setColName(i, s);
			}

			tableModelString.set(maSimilarity, 2);

			System.out.println(tableModelString.toString());

		}


		List<Double> liSimilarities = new ArrayList<>();

		if(base.getInteractionTypeCount() > query.getInteractionTypeCount()) {

			for (int col = 0; col < base.getInteractionTypeCount(); col++) {

				int interactionTypeQuery=-1;
				int interactionTypeBase=-1;

				double maxSimInCol = valNoInteraction;
				for (int row = 0; row < query.getInteractionTypeCount(); row++) {
					if(maSimilarity.get(row,col)>maxSimInCol){
						maxSimInCol = maSimilarity.get(row,col);

						interactionTypeBase = base.getInteractionType(col);
						interactionTypeQuery = query.getInteractionType(row);
					}
				}


				if(Math.abs(maxSimInCol-valNoInteraction)<TINY){
					continue;
				}

				liSimilarities.add(maxSimInCol);

				// System.out.println("Sim maxSimInCol " + Formatter.format2(maxSimInCol) + "\t" + InteractionAtomTypeCalculator.getString(interactionTypeBase) + "\t" + InteractionAtomTypeCalculator.getString(interactionTypeQuery));
			}
		} else {
			for (int row = 0; row < query.getInteractionTypeCount(); row++) {

				int interactionTypeQuery=-1;
				int interactionTypeBase=-1;

				double maxSimInRow = valNoInteraction;
				for (int col = 0; col < base.getInteractionTypeCount(); col++) {
					if(maSimilarity.get(row,col) > maxSimInRow){
						maxSimInRow = maSimilarity.get(row,col);
						interactionTypeBase = base.getInteractionType(col);
						interactionTypeQuery = query.getInteractionType(row);

					}
				}
				// System.out.println("Sim maxSimInRow " + Formatter.format2(maxSimInRow) + "\t" + InteractionAtomTypeCalculator.getString(interactionTypeBase) + "\t" + InteractionAtomTypeCalculator.getString(interactionTypeQuery));


				if(Math.abs(maxSimInRow-valNoInteraction)<TINY){
					continue;
				}

				liSimilarities.add(maxSimInRow);							}
		}


		double sim = 0;

		if(liSimilarities.size()>0) {
			sim = 1;
			for (Double simPart : liSimilarities) {
				sim *= simPart;
			}
		}


		if(verbose) {
			System.out.println("Sim " + Formatter.format2(sim));
			System.out.println();
		}

		return sim;
	}


	public boolean isValidType(int type){
		boolean valid = true;

		try {
			int key = InteractionDistanceStatistics.getInstance().getKey(type);
			interactionSimilarityTable.getDistance(key, key);
		} catch (Exception e) {
			valid = false;
		}

		return valid;
	}
	
}
