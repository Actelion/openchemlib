/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * 3. Neither the name of the copyright holder nor the
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
 * @author Modest v. Korff
 */

package com.actelion.research.chem.descriptor;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.descriptor.flexophore.*;
import com.actelion.research.chem.descriptor.flexophore.completegraphmatcher.ObjectiveBlurFlexophoreHardMatchUncovered;
import com.actelion.research.chem.descriptor.flexophore.completegraphmatcher.PPNodeSimilarity;
import com.actelion.research.chem.descriptor.flexophore.generator.CreatorMolDistHistViz;
import com.actelion.research.util.ArrayUtils;
import com.actelion.research.util.CommandLineParser;
import com.actelion.research.util.graph.complete.CompleteGraphMatcher;
import com.actelion.research.util.graph.complete.SolutionCompleteGraph;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 *
 * DescriptorHandlerFlexophore
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
public class DescriptorHandlerFlexophore implements IDescriptorHandlerFlexophore {

	public static final boolean DEBUG = false;

	public static final String SEP_PARAMETER = ",";

	public static final String PARA_TAG_PPNODE_SIMILARITY_THRESH = "ThreshNodeSimilarity";
	public static final String PARA_TAG_HISTOGRAM_SIMILARITY_THRESH = "ThreshHistogramSimilarity";
	public static final String PARA_TAG_MODE_MATCH = "PPNodeMatchMode";

	public static final String PARA_TAG_MODE_SINGLE_CONFORMATION_QUERY = "SingleConfQuery";

	protected static final int MAX_NUM_HEAVY_ATOMS = 70;

	// until 24.03.2020
	// until 17.12.2020
	// private static final double CORRECTION_FACTOR = 0.40;

	private static final double CORRECTION_FACTOR = 0.4;

	private static final int MAX_TRIES_TO_GENERATE_CONFORMER = 25;

	private static final int MAX_TRIES_TO_GENERATE_CONFORMER_ONE_CONF = 11;

	protected static final int MIN_NUM_ATOMS = 6;

	/**
	 *
	 * 06.09.2024 Set to 200 to be aligned with PheSA
	 */
	public static final int NUM_CONFORMATIONS = 200;

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
	public static final double THRESH_HISTOGRAM_SIMILARITY = ObjectiveBlurFlexophoreHardMatchUncovered.THRESH_HISTOGRAM_SIMILARITY;

	// Test
	// public static final double THRESH_HISTOGRAM_SIMILARITY = 0.5;

	private static DescriptorHandlerFlexophore INSTANCE;


	private ConcurrentLinkedQueue<CompleteGraphMatcher<IMolDistHist>> queueCGM;

	private MolDistHistEncoder molDistHistEncoder;

	protected CreatorMolDistHistViz creatorMolDistHistViz;



	//
	// If you change this, do not forget to change the objective in CompleteGraphMatcher<IMolDistHist>
	// getNewCompleteGraphMatcher().
	//
	private ObjectiveBlurFlexophoreHardMatchUncovered objectiveCompleteGraphHard;

	protected Exception recentException;

	private int versionInteractionTable;
	private int modePPNodeSimilarityComparison;
	private double threshSimilarityHardMatch;
	private double threshHistogramSimilarity;

	private boolean singleConformationModeQuery;
	private boolean includeNodeAtoms;

	private SolutionCompleteGraph solution;

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

		objectiveCompleteGraphHard = new ObjectiveBlurFlexophoreHardMatchUncovered(
				versionInteractionTable, modePPNodeSimilarityComparison, threshSimilarityHardMatch, threshHistogramSimilarity);

		queueCGM.add(getNewCompleteGraphMatcher());

		molDistHistEncoder = new MolDistHistEncoder();


		creatorMolDistHistViz = new CreatorMolDistHistViz();
	}

	public void setIncludeNodeAtoms(boolean b) {
		includeNodeAtoms = b;
	}

	public void setModeQuery(boolean modeQuery){
		objectiveCompleteGraphHard.setModeQuery(modeQuery);
		queueCGM.clear();
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

	public CompleteGraphMatcher<IMolDistHist> getNewCompleteGraphMatcher(){

		ObjectiveBlurFlexophoreHardMatchUncovered objective =
				new ObjectiveBlurFlexophoreHardMatchUncovered(
						versionInteractionTable,
						modePPNodeSimilarityComparison,
						threshSimilarityHardMatch,
						threshHistogramSimilarity);


		objective.setModeQuery(objectiveCompleteGraphHard.isModeQuery());

		CompleteGraphMatcher<IMolDistHist> cgMatcher = new CompleteGraphMatcher<>(objective);

		cgMatcher.setMaxNumSolutions(MAX_NUM_SOLUTIONS);

		return cgMatcher;
	}

	public ObjectiveBlurFlexophoreHardMatchUncovered getObjectiveCompleteGraph() {
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

		int[] oldToNewAtom = fragBiggest.stripSmallFragments();

		fragBiggest.ensureHelperArrays(StereoMolecule.cHelperCIP);

		if(fragBiggest.getAtoms() < MIN_NUM_ATOMS){
			return FAILED_OBJECT;
		} else if(fragBiggest.getAtoms() > MAX_NUM_HEAVY_ATOMS){
			return FAILED_OBJECT;
		}

		MolDistHistViz mdhv = createVisualDescriptor(fragBiggest);
		MolDistHist mdh = (mdhv == null) ? null : mdhv.getMolDistHist();

		recentException = null;

		if(mdh == null) {
			mdh = FAILED_OBJECT;
		} else if (mdh.getNumPPNodes() > ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE) {
			String msg = "Flexophore exceeded maximum number of nodes.";
			recentException = new RuntimeException(msg);
			mdh = FAILED_OBJECT;
		}
		else if (includeNodeAtoms) {
			List<PPNodeViz> nodeList = mdhv.getNodes();
			int[][] nodeAtom = new int[nodeList.size()][];
			for (int i=0; i<nodeList.size(); i++) {
				nodeAtom[i] = ArrayUtils.toIntArray(nodeList.get(i).getListIndexOriginalAtoms());
				if (oldToNewAtom != null) {
					int[] newToOldAtom = new int[fragBiggest.getAtoms()];
					for (int j=0; j<oldToNewAtom.length; j++)
						if (oldToNewAtom[j] != -1)
							newToOldAtom[oldToNewAtom[j]] = j;
					for (int j=0; j<nodeAtom[i].length; j++)
						nodeAtom[i][j] = newToOldAtom[nodeAtom[i][j]];
				}
			}
			mdh.setNodeAtoms(nodeAtom);
		}

		return mdh;
	}

	/**
	 * This descriptor contains the molecule used for construction. The descriptor also contains information about
	 * corresponding atoms in the molecule.
	 * @param fragBiggest
	 * @return
	 */
	public MolDistHistViz createVisualDescriptor(StereoMolecule fragBiggest) {

		MolDistHistViz mdhv = null;
		recentException = null;
		try {
			mdhv = creatorMolDistHistViz.create(fragBiggest);
		} catch (Exception e) {
			try {
				if(DEBUG) {
					Canonizer can = new Canonizer(fragBiggest);
					String msg = "DescriptorHandlerFlexophore Impossible to generate conformer for\n" + can.getIDCode();
					System.err.println(msg);
				}
			} catch (Exception e2) {
				e2.printStackTrace();
				recentException = e;
			}
		}
		return mdhv;
	}

	public MolDistHistViz createVisualDescriptor(ConformerSet cs){
		return creatorMolDistHistViz.createFromConformerSet(cs);
	}

	public Exception getRecentException() {
		return recentException;
	}

	/**
	 * In query mode, all pharmacophore points in the query must be found in base.
	 * @param query
	 * @param base
	 * @return
	 */
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

			if(mdhvBase.getNumPPNodes() > ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE){
				System.out.println("DescriptorHandlerFlexophore getSimilarity(...) mdhvBase.getNumPPNodes() " + mdhvBase.getNumPPNodes());
				return 0;
			} else if(mdhvQuery.getNumPPNodes() > ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE){
				System.out.println("DescriptorHandlerFlexophore getSimilarity(...) mdhvQuery.getNumPPNodes() " + mdhvQuery.getNumPPNodes());
				return 0;
			}

			sc = (float) getSimilarity(mdhvBase, mdhvQuery);
		}

		return normalizeValue(sc);
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

		solution = cgMatcher.getBestMatchingSolution();

		queueCGM.add(cgMatcher);

		return sc;
	}

	public SolutionCompleteGraph getRecentSolution(){
		return solution;
	}

	/**
	 *
	 * @param mdhvBase
	 * @param mdhvQuery
	 * @return
	 */
	public ModelSolutionSimilarity getBestMatch(MolDistHistViz mdhvBase, MolDistHistViz mdhvQuery){

		CompleteGraphMatcher<IMolDistHist> cgMatcher = queueCGM.poll();

		if(cgMatcher == null){
			cgMatcher = getNewCompleteGraphMatcher();
		}

		cgMatcher.set(mdhvBase, mdhvQuery);
		cgMatcher.calculateSimilarity();

		SolutionCompleteGraph scgBest = cgMatcher.getBestMatchingSolution();

		int n = scgBest.getSizeHeap();

		if (n == 0) // added to prevent NegativeArraySizeException in next statement; TLS 10Jan2022
			return null;

		float [] arrSimNode = new float[scgBest.getNodesQuery()];

		Arrays.fill(arrSimNode, -1);

		ObjectiveBlurFlexophoreHardMatchUncovered objectiveCompleteGraphHard = (ObjectiveBlurFlexophoreHardMatchUncovered)cgMatcher.getObjectiveCompleteGraph();
		objectiveCompleteGraphHard.setMatchingInfoInQueryAndBase(scgBest);

		for (int indexHeap = 0; indexHeap < n; indexHeap++) {

			int indexQuery = scgBest.getIndexQueryFromHeap(indexHeap);

			arrSimNode[indexQuery] = mdhvQuery.getNode(indexQuery).getSimilarityMappingNodes();
		}

		ModelSolutionSimilarity mss = new ModelSolutionSimilarity(scgBest, arrSimNode);

		return mss;
	}

	public ModelSolutionSimilarity getBestMatch(MolDistHist mdhBase, MolDistHist mdhQuery){

		MolDistHistViz mdhvBase = new MolDistHistViz(mdhBase);
		MolDistHistViz mdhvQuery = new MolDistHistViz(mdhQuery);

		CompleteGraphMatcher<IMolDistHist> cgMatcher = queueCGM.poll();

		if(cgMatcher == null){
			cgMatcher = getNewCompleteGraphMatcher();
		}

		cgMatcher.set(mdhvBase, mdhvQuery);
		cgMatcher.calculateSimilarity();

		SolutionCompleteGraph scgBest = cgMatcher.getBestMatchingSolution();

		int n = scgBest.getSizeHeap();

		if (n == 0) // added to prevent NegativeArraySizeException in next statement; TLS 10Jan2022
			return null;

		float [] arrSimNode = new float[scgBest.getNodesQuery()];

		Arrays.fill(arrSimNode, -1);

		ObjectiveBlurFlexophoreHardMatchUncovered objectiveCompleteGraphHard = (ObjectiveBlurFlexophoreHardMatchUncovered)cgMatcher.getObjectiveCompleteGraph();
		objectiveCompleteGraphHard.setMatchingInfoInQueryAndBase(scgBest);

		for (int indexHeap = 0; indexHeap < n; indexHeap++) {

			int indexQuery = scgBest.getIndexQueryFromHeap(indexHeap);

			arrSimNode[indexQuery] = mdhvQuery.getNode(indexQuery).getSimilarityMappingNodes();
		}

		ModelSolutionSimilarity mss = new ModelSolutionSimilarity(scgBest, arrSimNode);

		return mss;
	}

	public static float normalizeValue(double value) {
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
		dh.setIncludeNodeAtoms(includeNodeAtoms);

		return dh;
	}
	public void setObjectiveQueryBiased(boolean enable){
		throw new RuntimeException("setObjectiveQueryBiased(...) is not implemented");
	}

}
