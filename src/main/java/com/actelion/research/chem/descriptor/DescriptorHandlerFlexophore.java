/*
 * Copyright (c) 2020.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
 */

package com.actelion.research.chem.descriptor;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.descriptor.flexophore.*;
import com.actelion.research.chem.descriptor.flexophore.completegraphmatcher.ObjectiveFlexophoreHardMatchUncovered;
import com.actelion.research.chem.descriptor.flexophore.completegraphmatcher.PPNodeSimilarity;
import com.actelion.research.chem.descriptor.flexophore.generator.CreatorMolDistHistViz;
import com.actelion.research.util.CommandLineParser;
import com.actelion.research.util.graph.complete.CompleteGraphMatcher;

import java.util.Arrays;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 *
 * DescriptorHandlerFlexophore
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * 29 Jan 2009 MvK: Start implementation
 * 15 Oct 2012 MvK renamed DescriptorHandler3DMM2PPInteract-->DescriptorHandlerFlexophore
 * 19 Apr 2013 MvK major changes in Flexophore encoding decoding
 * 25 Apr 2013 MvK Flexophore version changed --> 3.0
 * 07 May 2013 MvK bug fixes in encoding. Flexophore version changed --> 3.1
 * 15 May 2013 MvK Bug fix for the objective function
 * 03 May 2016 MvK versioning for interaction tables from Joel introduced.
 * 10 Jun 2016 MvK DescriptorHandlerFlexophoreV4 --> DescriptorHandlerFlexophore, V4 becomes today the new Flexophore
 * 18 Jun 2016 MvK New ConformationGenerator from TS for Flexophore creation. V 4.1
 * 29 Jun 2016 MvK number of histogram bins and range increased in CGMult. V 4.2
 * 15 Jul 2016 MvK if generation of conformer failed a new seed is injected and the generation is tried again.
 * 11 Aug 2016 MvK number of bins increase from 50 to 80, histogram range increased from 25 to 40 Angstroem. --> V.4.3
 * 30 Jan 2017 MvK minor bug fix. Two constants for the number of conformations. --> V.4.4. Compatible with V.4.3
 * April 2020 Version 5.0, new interaction tables from Joel Wahl, algorithmic changes.
 * 20.05.2020 Changed definition of end-standing aliphatic groups. Center of gravity moved one bond to the outside.
 * Aromatic rings are not aliphatic any more. Methyl groups connected to a ring are not considered. The next O or N
 * must have a minimum distance of two bonds to the end standing atom.
 */
public class DescriptorHandlerFlexophore implements DescriptorHandler {

	public static final boolean DEBUG = false;

	public static final String SEP_PARAMETER = ",";

	public static final String PARA_TAG_PPNODE_SIMILARITY_THRESH = "ThreshNodeSimilarity";
	public static final String PARA_TAG_HISTOGRAM_SIMILARITY_THRESH = "ThreshHistogramSimilarity";
	public static final String PARA_TAG_MODE_MATCH = "PPNodeMatchMode";

	public static final String PARA_TAG_MODE_SINGLE_CONFORMATION_QUERY = "SingleConfQuery";

	protected static final int MAX_NUM_HEAVY_ATOMS = 70;

	// until 24.03.2020
	// private static final double CORRECTION_FACTOR = 0.40;

	private static final double CORRECTION_FACTOR = 0.4;

	private static final int MAX_TRIES_TO_GENERATE_CONFORMER = 25;

	private static final int MAX_TRIES_TO_GENERATE_CONFORMER_ONE_CONF = 11;

	protected static final int MIN_NUM_ATOMS = 6;

	// 250
	public static final int NUM_CONFORMATIONS = 250;

	public static final int MAX_NUM_SOLUTIONS = 1000;

	public static final MolDistHist FAILED_OBJECT = new MolDistHist();

	// Version 3.0 after definition of new interaction types by Joel Freyss.
	// 07.05.2013 Version 3.1 after bug fixes in encoding.
	// 17.09.2015 Version 3.2. Joel re-calculated interaction tables. Differences in atom types.
	// 10.06.2016 Version 4.0. Total re-implementation of the Flexophore. The pharmacophore point recognition is now
	// more generic.
	// April 2020 Version 5.0, new interaction tables from Joel Wahl
	// Version not defined. -1
	public static final int VERSION_INTERACTION_TABLES = -1;

	public static final int MODE_PPNODE_SIMILARITY_COMPARISON = PPNodeSimilarity.SIMILARITY_MODE_HARD_THRESH;

	// Production
	public static final double THRESH_SIMILARITY_COMPARISON_NODE = PPNodeSimilarity.THRESH_SIMILARITY_HARD_MATCH;

	// Test
	// public static final double THRESH_SIMILARITY_COMPARISON_NODE = 0.001;

	// Production
	public static final double THRESH_HISTOGRAM_SIMILARITY = ObjectiveFlexophoreHardMatchUncovered.THRESH_HISTOGRAM_SIMILARITY;

	// Test
	// public static final double THRESH_HISTOGRAM_SIMILARITY = 0.5;

	private static DescriptorHandlerFlexophore INSTANCE;


	private ConcurrentLinkedQueue<CompleteGraphMatcher<IMolDistHist>> queueCGM;

	private MolDistHistEncoder molDistHistEncoder;

	protected CreatorMolDistHistViz creatorMolDistHistViz;



	//
	// If you change this, do not forget to change the objective in CompleteGraphMatcher<IMolDistHist> getNewCompleteGraphMatcher().
	//
	private ObjectiveFlexophoreHardMatchUncovered objectiveCompleteGraphHard;

	protected Exception recentException;

	private int versionInteractionTable;
	private int modePPNodeSimilarityComparison;
	private double threshSimilarityHardMatch;
	private double threshHistogramSimilarity;

	private boolean singleConformationModeQuery;


	private ThreadMaster threadMaster;

	public DescriptorHandlerFlexophore(String parameter) {
		CommandLineParser cmd = new CommandLineParser(parameter, SEP_PARAMETER);
		int versionInteractionTable = VERSION_INTERACTION_TABLES;
		int modePPNodeSimilarityComparison = MODE_PPNODE_SIMILARITY_COMPARISON;
		double threshSimilarityHardMatch = THRESH_SIMILARITY_COMPARISON_NODE;
		double threshHistogramSimilarity = THRESH_HISTOGRAM_SIMILARITY;
		boolean singleConformationModeQuery = false;

		if(cmd.contains(PARA_TAG_PPNODE_SIMILARITY_THRESH)) {
			threshSimilarityHardMatch =  cmd.getAsDouble(PARA_TAG_PPNODE_SIMILARITY_THRESH);
		}
		if(cmd.contains(PARA_TAG_MODE_MATCH)) {
			modePPNodeSimilarityComparison =  cmd.getAsInt(PARA_TAG_MODE_MATCH);
		}
		if(cmd.contains(PARA_TAG_HISTOGRAM_SIMILARITY_THRESH)) {
			threshHistogramSimilarity =  cmd.getAsDouble(PARA_TAG_HISTOGRAM_SIMILARITY_THRESH);
		}

		if(cmd.contains(PARA_TAG_MODE_SINGLE_CONFORMATION_QUERY)) {
			singleConformationModeQuery =  cmd.getAsBoolean(PARA_TAG_MODE_SINGLE_CONFORMATION_QUERY);
		}

		init(versionInteractionTable, modePPNodeSimilarityComparison, threshSimilarityHardMatch, threshHistogramSimilarity, singleConformationModeQuery);
	}

	public DescriptorHandlerFlexophore() {
		init(VERSION_INTERACTION_TABLES,
				MODE_PPNODE_SIMILARITY_COMPARISON,
				THRESH_SIMILARITY_COMPARISON_NODE,
				THRESH_HISTOGRAM_SIMILARITY, false);
	}

	public DescriptorHandlerFlexophore(
			int versionInteractionTable,
			int modePPNodeSimilarityComparison,
			double threshSimilarityHardMatch,
			double threshHistogramSimilarity,
			boolean singleConformationModeQuery) {
		init(
				versionInteractionTable,
				modePPNodeSimilarityComparison,
				threshSimilarityHardMatch,
				threshHistogramSimilarity,
				singleConformationModeQuery);
	}

	private void init(
			int versionInteractionTable,
			int modePPNodeSimilarityComparison,
			double threshSimilarityHardMatch,
			double threshHistogramSimilarity,
			boolean singleConformationModeQuery){

		this.versionInteractionTable = versionInteractionTable;
		this.modePPNodeSimilarityComparison = modePPNodeSimilarityComparison;
		this.threshSimilarityHardMatch = threshSimilarityHardMatch;
		this.threshHistogramSimilarity = threshHistogramSimilarity;
		this.singleConformationModeQuery = singleConformationModeQuery;

		MolDistHistViz.createIndexTables();

		queueCGM = new ConcurrentLinkedQueue<>();

		queueCGM.add(getNewCompleteGraphMatcher());

		molDistHistEncoder = new MolDistHistEncoder();

		objectiveCompleteGraphHard = new ObjectiveFlexophoreHardMatchUncovered(
				versionInteractionTable, modePPNodeSimilarityComparison, threshSimilarityHardMatch, threshHistogramSimilarity);

		creatorMolDistHistViz = new CreatorMolDistHistViz();


	}

	public void setModeQuery(boolean modeQuery){
		objectiveCompleteGraphHard.setModeQuery(modeQuery);
	}

	public void setThreadMaster(ThreadMaster threadMaster) {
		this.threadMaster = threadMaster;
		creatorMolDistHistViz.setThreadMaster(threadMaster);
	}

	public boolean isSingleConformationModeQuery() {
		return singleConformationModeQuery;
	}

	public void setSingleConformationModeQuery(boolean singleConformationModeQuery) {
		this.singleConformationModeQuery = singleConformationModeQuery;
	}

	private CompleteGraphMatcher<IMolDistHist> getNewCompleteGraphMatcher(){

		ObjectiveFlexophoreHardMatchUncovered objective =
				new ObjectiveFlexophoreHardMatchUncovered(
						versionInteractionTable,
						modePPNodeSimilarityComparison,
						threshSimilarityHardMatch,
						threshHistogramSimilarity);

		CompleteGraphMatcher<IMolDistHist> cgMatcher = new CompleteGraphMatcher<>(objective);

		cgMatcher.setMaxNumSolutions(MAX_NUM_SOLUTIONS);

		return cgMatcher;
	}

	public ObjectiveFlexophoreHardMatchUncovered getObjectiveCompleteGraphHard() {
		return objectiveCompleteGraphHard;
	}

	public DescriptorInfo getInfo() {
		return DescriptorConstants.DESCRIPTOR_Flexophore;
	}

	public String getVersion() {
		return DescriptorConstants.DESCRIPTOR_Flexophore.version;
	}

	public String toStringParameter() {
		StringBuilder sb = new StringBuilder();

		sb.append("DescriptorHandlerFlexophoreV5");
		sb.append(" ");
		if(singleConformationModeQuery){
			sb.append("singleConformationModeQuery=true");
		}else{
			sb.append("singleConformationModeQuery=false");
		}
		sb.append(", ");
		sb.append(objectiveCompleteGraphHard.toStringParameter());

		return sb.toString();
	}

	public static DescriptorHandlerFlexophore getDefaultInstance(){
		if(INSTANCE==null){
			synchronized (DescriptorHandlerFlexophore.class) {
				INSTANCE = new DescriptorHandlerFlexophore();
			}
		}
		return INSTANCE;
	}

	public String encode(Object o) {

		if(calculationFailed(o)){
			return FAILED_STRING;
		}

		MolDistHist mdh = null;

		if(o instanceof MolDistHist){

			mdh = (MolDistHist)o;

		} else if(o instanceof MolDistHistViz){

			mdh = ((MolDistHistViz)o).getMolDistHist();

		} else {
			return FAILED_STRING;
		}

		return molDistHistEncoder.encode(mdh);

	}

	public MolDistHist decode(byte[] bytes) {
		try {
			return bytes == null || bytes.length == 0 ? null : Arrays.equals(bytes, FAILED_BYTES) ? FAILED_OBJECT : molDistHistEncoder.decode(bytes);
		} catch (RuntimeException e1) {
			return FAILED_OBJECT;
		}
	}

	public MolDistHist decode(String s) {
		try {
			return s == null || s.length() == 0 ? null
					: s.equals(FAILED_STRING) ? FAILED_OBJECT
					:                           molDistHistEncoder.decode(s);
		} catch (RuntimeException e1) {
			return FAILED_OBJECT;
		}
	}

	public MolDistHist createDescriptorSingleConf(StereoMolecule mol) {

		MolDistHistViz mdhv = creatorMolDistHistViz.createFromGivenConformation(mol);

		return mdhv.getMolDistHist();
	}

	public MolDistHist createDescriptorSingleConf(ConformerSet conformerSet) {

		if(conformerSet==null){
			return FAILED_OBJECT;
		}

		MolDistHistViz mdhv = createVisualDescriptorSingleConf(conformerSet);

		return mdhv.getMolDistHist();
	}

	public MolDistHistViz createVisualDescriptorSingleConf(ConformerSet conformerSet) {

		if(conformerSet==null){
			return new MolDistHistViz();
		}

		StereoMolecule mol = conformerSet.first().toMolecule();

		MolDistHistViz mdhv = creatorMolDistHistViz.createFromGivenConformation(mol);

		return mdhv;
	}

	public MolDistHist createDescriptor(Object mol) {

		StereoMolecule fragBiggest = (StereoMolecule)mol;

		fragBiggest.stripSmallFragments();

		fragBiggest.ensureHelperArrays(StereoMolecule.cHelperCIP);

		if(fragBiggest.getAtoms() < MIN_NUM_ATOMS){
			return FAILED_OBJECT;
		} else if(fragBiggest.getAtoms() > MAX_NUM_HEAVY_ATOMS){
			return FAILED_OBJECT;
		}

//        IDCodeParser parser = new IDCodeParser();
//		fragBiggest = parser.getCompactMolecule(new Canonizer(fragBiggest).getIDCode());

		MolDistHistViz mdhv = createVisualDescriptor(fragBiggest);
		MolDistHist mdh = (mdhv == null) ? null : mdhv.getMolDistHist();

		recentException = null;

		if(mdh == null) {
			mdh = FAILED_OBJECT;
		} else if (mdh.getNumPPNodes() > ObjectiveFlexophoreHardMatchUncovered.MAX_NUM_NODES_FLEXOPHORE) {

			String msg = "Flexophore exceeded maximum number of nodes.";

			recentException = new RuntimeException(msg);

			mdh = FAILED_OBJECT;;
		}

		return mdh;
	}

	public MolDistHistViz createVisualDescriptor(StereoMolecule fragBiggest) {

		MolDistHistViz mdhv = null;

		boolean conformationGenerationFailed = true;

		int ccFailed = 0;

		recentException = null;

		while (conformationGenerationFailed) {

			conformationGenerationFailed = true;

			try {

				mdhv = creatorMolDistHistViz.create(fragBiggest);

				if(creatorMolDistHistViz.isOnlyOneConformer() && (ccFailed < MAX_TRIES_TO_GENERATE_CONFORMER_ONE_CONF)){
					conformationGenerationFailed = true;
				} else {
					conformationGenerationFailed = false;
				}

			} catch (ExceptionConformationGenerationFailed e) {
				recentException = e;
			} catch (Exception e) {
				recentException = e;
				break;
			}

			if(threadMaster!=null && threadMaster.threadMustDie()){
				break;
			}

			if(conformationGenerationFailed) {
				if(DEBUG) {
					System.out.println("DescriptorHandlerFlexophore Inject new seed");
				}

				creatorMolDistHistViz.injectNewSeed();

				ccFailed++;
			}

			if(ccFailed==MAX_TRIES_TO_GENERATE_CONFORMER){

				try {
					if(DEBUG) {

						Canonizer can = new Canonizer(fragBiggest);

						String msg = "DescriptorHandlerFlexophore Impossible to generate conformer for\n" + can.getIDCode();

						System.err.println(msg);
					}

				} catch (Exception e) {
					e.printStackTrace();
					recentException = e;
				}
				break;
			}
		}

		return mdhv;
	}

	public Exception getRecentException() {
		return recentException;
	}


	public float getSimilarity(Object query, Object base) {

		float sc=0;

		if(base == null
				|| query == null
				|| ((IMolDistHist)base).getNumPPNodes() == 0
				|| ((IMolDistHist)query).getNumPPNodes() == 0) {
			sc = 0;

		} else {

			IMolDistHist mdhvBase = (IMolDistHist)base;

			IMolDistHist mdhvQuery = (IMolDistHist)query;

			if(mdhvBase.getNumPPNodes() > ObjectiveFlexophoreHardMatchUncovered.MAX_NUM_NODES_FLEXOPHORE){

				System.out.println("DescriptorHandlerFlexophore getSimilarity(...) mdhvBase.getNumPPNodes() " + mdhvBase.getNumPPNodes());

				return 0;
			} else if(mdhvQuery.getNumPPNodes() > ObjectiveFlexophoreHardMatchUncovered.MAX_NUM_NODES_FLEXOPHORE){
				System.out.println("DescriptorHandlerFlexophore getSimilarity(...) mdhvQuery.getNumPPNodes() " + mdhvQuery.getNumPPNodes());
				return 0;
			}

			sc = (float)getMinimumSimilarity(mdhvBase, mdhvQuery);

		}

		return normalizeValue(sc);

		// return sc;
	}

	private double getMinimumSimilarity(IMolDistHist mdhvBase, IMolDistHist mdhvQuery){

		double sc = 0;

		if(mdhvBase.getNumPPNodes() == mdhvQuery.getNumPPNodes()){
			double s1 = getSimilarity(mdhvBase, mdhvQuery);
			double s2 = getSimilarity(mdhvQuery, mdhvBase);

			sc = Math.max(s1, s2);
		} else {
			sc = getSimilarity(mdhvBase, mdhvQuery);
		}

		return sc;
	}


	private double getSimilarity(IMolDistHist iBase, IMolDistHist iQuery){

		MolDistHistViz mdhvBase = null;
		if(iBase instanceof MolDistHist){
			mdhvBase = new MolDistHistViz((MolDistHist)iBase);
		} else if(iBase instanceof MolDistHistViz){
			mdhvBase = new MolDistHistViz((MolDistHistViz)iBase);
		}

		MolDistHistViz mdhvQuery = null;
		if(iQuery instanceof MolDistHist){
			mdhvQuery = new MolDistHistViz((MolDistHist)iQuery);
		} else if(iQuery instanceof MolDistHistViz){
			mdhvQuery = new MolDistHistViz((MolDistHistViz)iQuery);
		}

		CompleteGraphMatcher<IMolDistHist> cgMatcher = queueCGM.poll();

		if(cgMatcher == null){
			cgMatcher = getNewCompleteGraphMatcher();
		}

		cgMatcher.set(mdhvBase, mdhvQuery);

		double sc = (float)cgMatcher.calculateSimilarity();

		queueCGM.add(cgMatcher);

		return sc;
	}

	public double getSimilarityNodes(PPNode query, PPNode base) {
		return objectiveCompleteGraphHard.getSimilarityNodes(query, base);
	}


	public float normalizeValue(double value) {
		return value <= 0.0f ? 0.0f
				: value >= 1.0f ? 1.0f
				: (float)(1.0-Math.pow(1-Math.pow(value, CORRECTION_FACTOR) ,1.0/CORRECTION_FACTOR));
	}

	public boolean calculationFailed(Object o) {

		if(o instanceof MolDistHist){
			return ((MolDistHist)o).getNumPPNodes() == 0;
		} else if(o instanceof MolDistHistViz){
			return ((MolDistHistViz)o).getNumPPNodes() == 0;
		}

		return true;

	}

	public DescriptorHandlerFlexophore getThreadSafeCopy() {

		DescriptorHandlerFlexophore dh = new DescriptorHandlerFlexophore(
				versionInteractionTable,
				modePPNodeSimilarityComparison,
				threshSimilarityHardMatch,
				threshHistogramSimilarity,
				singleConformationModeQuery);

		dh.setModeQuery(objectiveCompleteGraphHard.isModeQuery());

		return dh;
	}
	public void setObjectiveQueryBiased(boolean enable){
		throw new RuntimeException("setObjectiveQueryBiased(...) is not implemented");
	}

}
