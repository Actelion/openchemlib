package com.actelion.research.chem.descriptor.flexophore.completegraphmatcher;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.graph.MinimumSpanningTree;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.chem.descriptor.flexophore.*;
import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.graph.complete.IObjectiveCompleteGraph;
import com.actelion.research.util.graph.complete.SolutionCompleteGraph;

import java.util.Arrays;
import java.util.List;

/**
 * 
 * 
 * ObjectiveBlurFlexophoreHardMatchUncovered
 * The weighting of the coverage is hard. Which means that uncovered nodes 
 * strongly change the final similarity score.
 * The distance histograms are heavily blurred. This increases the probability of matching.
 * look in <code>getScoreUncoveredNearestNodesBase(SolutionCompleteGraph solution)</code> 
 * and <code>getScoreUncoveredNearestNodesQuery(SolutionCompleteGraph solution)</code>.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * Oct 2, 2012 MvK: Start implementation
 * Mar 3. 2016 MvK: updates. Lowered thresh for histogram similarity.
 * Mar 31. 2020 MvK: fraction of carbon is considered in pharmacophore node similarity.
 */
public class ObjectiveBlurFlexophoreHardMatchUncovered implements IObjectiveCompleteGraph<IMolDistHist>{

	public static final String VERSION = "07.09.2022";
	public static final String INFO = "";

	// 0.15 best thresh tested on 31.03.2016
	public final static double THRESH_HISTOGRAM_SIMILARITY	= 0.15;

	public final static double THRESH_HISTOGRAM_SIMILARITY_OPTIMISTIC	= 0.0000001;


	// The thresh for the node similarity depends on the number of interaction types in the node.
	// 03.03.2016 Top result so far for 0.9
	// 13.04.2020 Maybe obsolete
	// ToDo
	// final static double THRESH_NODE_SIMILARITY_START = 0.5;

	// Changed to 0.9 21.08.2024 MvK
	final static double THRESH_NODE_SIMILARITY_START = 0.9;
	final static double OPTIMISTIC_HISTOGRAM_THRESH = 0.0;

	private static final float INIT_VAL = -1;

	private boolean modeQuery;

	private int marginQuery;

	private MolDistHistViz mdhvBase;
	private MolDistHistViz mdhvBaseBlurredHist;

	private MolDistHistViz mdhvQuery;
	private MolDistHistViz mdhvQueryBlurredHist;

	private int nodesBase;

	private int nodesQuery;


	private byte [] arrTmpHist;

	private boolean validHelpersQuery;

	private boolean validHelpersBase;

	private boolean resetSimilarityArrays;

	private double threshNodeMinSimilarityStart;

	private double threshHistogramSimilarity;

	private PPNodeSimilarity nodeSimilarity;

	private boolean optimisticHistogramSimilarity;

	private boolean excludeHistogramSimilarity;

	// For small Flexophore descriptor comparison.
	private boolean fragmentNodesMapping;

	boolean verbose;


	private float [][] arrSimilarityNodes;

	private float [][] arrSimilarityHistograms;

	private double [][] arrRelativeDistanceMatrixQuery;

	private double [][] arrRelativeDistanceMatrixBase;

	private double [] arrSimilarityTmp;

	private Matrix maHelperAdjacencyQuery;

	private Matrix maHelperAdjacencyBase;

	private double sumDistanceMinSpanTreeQuery;

	private double sumDistanceMinSpanTreeBase;

	private int numMandatoryPPPoints;

	private boolean [] arrMandatoryPPPoint;

	private double avrPairwiseMappingScaled;

	private double coverageQuery;

	private double coverageBase;

	private double similarity;

	private long deltaNanoQueryBlur;
	private long deltaNanoBaseBlur;
	private long deltaNanoSimilarity;


	private SlidingWindowDistHist slidingWindowDistHist;

	public ObjectiveBlurFlexophoreHardMatchUncovered(){
		this(DescriptorHandlerFlexophore.VERSION_INTERACTION_TABLES,
				DescriptorHandlerFlexophore.MODE_PPNODE_SIMILARITY_COMPARISON,
				DescriptorHandlerFlexophore.THRESH_SIMILARITY_COMPARISON_NODE,
				THRESH_HISTOGRAM_SIMILARITY);

	}

	public ObjectiveBlurFlexophoreHardMatchUncovered(
			int versionInteractionTable,
			int modePPNodeSimilarity,
			double threshSimilarityNodeHardMatch,
			double threshHistogramSimilarity) {
		
		arrTmpHist =  new byte[ConstantsFlexophoreGenerator.BINS_HISTOGRAM];
		nodeSimilarity = new PPNodeSimilarity(versionInteractionTable, modePPNodeSimilarity);
		nodeSimilarity.setThreshSimilarityHardMatch(threshSimilarityNodeHardMatch);
		this.threshHistogramSimilarity = threshHistogramSimilarity;
		threshNodeMinSimilarityStart = THRESH_NODE_SIMILARITY_START;
		slidingWindowDistHist = new SlidingWindowDistHist(ConstantsFlexophoreGenerator.FILTER);
		modeQuery = false;
		deltaNanoQueryBlur=0;
		deltaNanoBaseBlur=0;
		deltaNanoSimilarity=0;
		setFragmentNodesMapping(false);
		setOptimisticHistogramSimilarity(false);
		marginQuery = 0;
		excludeHistogramSimilarity=false;
		initSimilarityMatrices();

	}

	/**
	 * Allows mapping of small Flexophores, up to one pharmacophore node.
	 * @param fragmentNodesMapping
	 */
	public void setFragmentNodesMapping(boolean fragmentNodesMapping) {
		this.fragmentNodesMapping = fragmentNodesMapping;
	}

	/**
	 * Only used in mode query
	 * The query must hit with all pharmacophore nodes except margin. Margin gives the number of nodes that need
	 * @param marginQuery
	 */
	public void setMarginQuery(int marginQuery) {
		this.marginQuery = marginQuery;
	}

	/**
	 * An overlap of two compared histograms is scored as a full match.
	 * @param optimisticHistogramSimilarity
	 */
	public void setOptimisticHistogramSimilarity(boolean optimisticHistogramSimilarity) {
		this.optimisticHistogramSimilarity = optimisticHistogramSimilarity;
		if(optimisticHistogramSimilarity){
			threshHistogramSimilarity = THRESH_HISTOGRAM_SIMILARITY_OPTIMISTIC;
		}
	}

	public void setExcludeHistogramSimilarity(boolean excludeHistogramSimilarity) {
		this.excludeHistogramSimilarity = excludeHistogramSimilarity;
	}
	public boolean isExcludeHistogramSimilarity() {
		return excludeHistogramSimilarity;
	}

	private void initSimilarityMatrices(){
		
		arrSimilarityNodes = new float [ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE][];
		for (int i = 0; i < ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE; i++) {
			arrSimilarityNodes[i] = new float [ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE];
			Arrays.fill(arrSimilarityNodes[i], INIT_VAL);
		}

		int maxNumHistograms = ((ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE* ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE)- ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE)/2;
		
		arrSimilarityHistograms = new float [maxNumHistograms][];
		for (int i = 0; i < maxNumHistograms; i++) {
			arrSimilarityHistograms[i] = new float [maxNumHistograms];
			Arrays.fill(arrSimilarityHistograms[i], INIT_VAL);
		}

		arrSimilarityTmp=new double[ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE];
	}

	/**
	 *
	 * @param modeQuery if set true the query must hit with all pharmacophore nodes. No penalty terms apply for not
	 * matched base nodes.
	 */
	public void setModeQuery(boolean modeQuery) {
		this.modeQuery = modeQuery;
	}

	public boolean isModeQuery() {
		return modeQuery;
	}

	public String toStringParameter(){
		StringBuilder sb = new StringBuilder();

		sb.append("ObjectiveFlexophoreHardMatchUncovered Thresh histogram similarity ");
		sb.append(threshHistogramSimilarity);
		sb.append(", thresh node similarity start " + threshNodeMinSimilarityStart);

		sb.append(", ");
		sb.append(nodeSimilarity.toStringParameter());

		return sb.toString();
	}

	private void resetSimilarityMatrices(){
		
		for (int i = 0; i < nodesQuery; i++) {
			for (int j = 0; j < nodesBase; j++) {
				arrSimilarityNodes[i][j] = INIT_VAL;	
			}
		}

		int numHistogramsQuery = ((nodesQuery*nodesQuery)-nodesQuery)/2;
		
		int numHistogramsBase = ((nodesBase*nodesBase)-nodesBase)/2;
		
		for (int i = 0; i < numHistogramsQuery; i++) {
			for (int j = 0; j < numHistogramsBase; j++) {
				arrSimilarityHistograms[i][j] = INIT_VAL;	
			}
		}
		
		resetSimilarityArrays = false;
	}

	@Override
	public void setVerbose(boolean v) {
		verbose=v;
		nodeSimilarity.setVerbose(v);
	}

	/**
	 * If a single histogram is not matching the solution is invalid.
	 * 
	 * If at least one node is not matching the solution is invalid.
	 */
	public boolean isValidSolution(SolutionCompleteGraph solution) {
		
		boolean mapping = true;
		if(!validHelpersQuery){
			calculateHelpersQuery();
		}
		if(!validHelpersBase){
			calculateHelpersBase();
		}
		if(resetSimilarityArrays){
			resetSimilarityMatrices();
		}
		
		// 
		// Should contain at least one pppoint with a hetero atom.
		// 
		
		int heap = solution.getSizeHeap();
		//
		// Check for inevitable pharmacophore points.
		//
		if(numMandatoryPPPoints > 0) {
			int ccInevitablePPPointsInSolution = 0;
			for (int i = 0; i < heap; i++) {
				int indexNodeQuery = solution.getIndexQueryFromHeap(i);
				if(mdhvQueryBlurredHist.isMandatoryPharmacophorePoint(indexNodeQuery)){
					ccInevitablePPPointsInSolution++;
				}
			}
			int neededMinInevitablePPPoints = Math.min(heap, numMandatoryPPPoints);
			if(ccInevitablePPPointsInSolution < neededMinInevitablePPPoints){
				mapping = false;
			}
		}

		//
		// Check for one hetero atom in solution.
		// Not checked for fragment mapping!

		if (!fragmentNodesMapping && mapping) {
			boolean heteroNodeQuery = false;
			boolean heteroNodeBase = false;
			for (int i = 0; i < heap; i++) {
				int indexNodeQuery = solution.getIndexQueryFromHeap(i);
				PPNode nodeQuery = mdhvQueryBlurredHist.getNode(indexNodeQuery);
				if (nodeQuery.hasHeteroAtom()) {
					heteroNodeQuery = true;
				}
				int indexNodeBase = solution.getIndexCorrespondingBaseNode(indexNodeQuery);
				PPNode nodeBase = mdhvBaseBlurredHist.getNode(indexNodeBase);
				if (nodeBase.hasHeteroAtom()) {
					heteroNodeBase = true;
				}
				if (heteroNodeQuery && heteroNodeBase) {
					break;
				}
			}
			if (!heteroNodeQuery || !heteroNodeBase) {
				mapping = false;
			}
		}


		//
		// Check for matching nodes.
		//
		if(mapping){
			for (int i = 0; i < heap; i++) {
				int indexNodeQuery = solution.getIndexQueryFromHeap(i);
				int indexNodeBase = solution.getIndexCorrespondingBaseNode(indexNodeQuery);
				if(!areNodesMapping(indexNodeQuery, indexNodeBase)) {
					mapping = false;
					break;
				}
			}
		}
		
		//
		// Check for matching histograms.
		//
		if(mapping && !excludeHistogramSimilarity){
			outer:
			for (int i = 0; i < heap; i++) {

				int indexNode1Query = solution.getIndexQueryFromHeap(i);

				int indexNode1Base = solution.getIndexCorrespondingBaseNode(indexNode1Query);

				for (int j = i+1; j < heap; j++) {
					int indexNode2Query = solution.getIndexQueryFromHeap(j);

					int indexNode2Base = solution.getIndexCorrespondingBaseNode(indexNode2Query);

					if(!areHistogramsMapping(indexNode1Query, indexNode2Query, indexNode1Base, indexNode2Base)){
						mapping = false;
						break outer;
					} else {
						// System.out.println("Match");
					}
				}
			}
		}

		return mapping;
	}

	/**
	 * Dynamic calculation of similarity threshold. Depends on the number of interaction types in the nodes.
	 * @param indexNodeQuery
	 * @param indexNodeBase
	 * @return
	 */
	public boolean areNodesMapping(int indexNodeQuery, int indexNodeBase) {
		
		if(!validHelpersQuery){
			calculateHelpersQuery();
		}
		
		if(!validHelpersBase){
			calculateHelpersBase();
		}
		
		if(resetSimilarityArrays){
			resetSimilarityMatrices();
		}

		boolean match = true;

		double simNodes = getSimilarityNodes(indexNodeQuery, indexNodeBase);

		//
		// Dynamic calculation of similarity threshold.
		//
		PPNode ppNodeBase = mdhvBaseBlurredHist.getNode(indexNodeBase);

		PPNode ppNodeQuery = mdhvQueryBlurredHist.getNode(indexNodeQuery);

		int interactionTypeCountBase = ppNodeBase.getInteractionTypeCount();

		int interactionTypeCountQuery = ppNodeQuery.getInteractionTypeCount();

		double threshCalc = 0;

		if(interactionTypeCountBase > interactionTypeCountQuery) {
			threshCalc = Math.pow(threshNodeMinSimilarityStart, interactionTypeCountBase);
		} else {
			threshCalc = Math.pow(threshNodeMinSimilarityStart, interactionTypeCountQuery);
		}

		// System.out.println("threshCalc " + Formatter.format3(threshCalc));

		if(simNodes < threshCalc){
			match=false;
		}
		
		return match;
	}
	
	private boolean areHistogramsMapping(int indexNode1Query, int indexNode2Query, int indexNode1Base, int indexNode2Base) {
		
		boolean match = true;

		double simHistograms = getSimilarityHistogram(indexNode1Query, indexNode2Query, indexNode1Base, indexNode2Base);

		if(simHistograms < threshHistogramSimilarity){
			match=false;
		}
		
		return match;
	}

	/**
	 * Calculates the similarity for the pharmacophore nodes and the distance histograms.
	 * @param solution
	 * @return
	 */
	public float getSimilarity(SolutionCompleteGraph solution) {

		long t0 = System.nanoTime();

		if(!validHelpersQuery){
			calculateHelpersQuery();
		}
		
		if(!validHelpersBase){
			calculateHelpersBase();
		}

		if(resetSimilarityArrays){
			resetSimilarityMatrices();
		}
		
		int heap = solution.getSizeHeap();

		//
		// the query must hit with all pharmacophore nodes
		//
		if(modeQuery) {
			if (nodesQuery != heap) {
				similarity=0;
				return (float)similarity;
			}
		}

		if(fragmentNodesMapping && heap==1){
			int indexNodeQuery = solution.getIndexQueryFromHeap(0);
			int indexNodeBase = solution.getIndexCorrespondingBaseNode(indexNodeQuery);
			similarity = getSimilarityNodes(indexNodeQuery, indexNodeBase);
			return (float)similarity;
		}

		if(numMandatoryPPPoints>0) {
			int ccMandatoryPPPoints = 0;
			for (int i = 0; i < heap; i++) {
				int indexNode1Query = solution.getIndexQueryFromHeap(i);
				if (arrMandatoryPPPoint[indexNode1Query]) {
					ccMandatoryPPPoints++;
				}
			}
			if(numMandatoryPPPoints>ccMandatoryPPPoints){
				similarity=0;
				return (float)similarity;
			}
		}

		int cc=0;
		int nMappings = ((heap * heap)-heap) / 2;
		double [] arrMappingWeights = new double[nMappings];
		double [] arrSimilarityWeighted = new double[nMappings];

		for (int i = 0; i < heap; i++) {
			int indexNode1Query = solution.getIndexQueryFromHeap(i);
			int indexNode1Base = solution.getIndexCorrespondingBaseNode(indexNode1Query);
			
			for (int j = i+1; j < heap; j++) {
				int indexNode2Query = solution.getIndexQueryFromHeap(j);
				int indexNode2Base = solution.getIndexCorrespondingBaseNode(indexNode2Query);
				double scorePairwiseMapping = getScorePairwiseMapping(indexNode1Query, indexNode2Query, indexNode1Base, indexNode2Base);
				double w =
						mdhvQuery.getWeightPharmacophorePoint(indexNode1Query)
								* mdhvQuery.getWeightPharmacophorePoint(indexNode2Query);

				arrMappingWeights[cc]=w;
				arrSimilarityWeighted[cc++]=scorePairwiseMapping * w;
				if(verbose) {
					System.out.println("scorePairwiseMapping " + Formatter.format2(scorePairwiseMapping));
				}
			}
		}
		// double mappings = ((heap * heap)-heap) / 2.0;

		double sumMappingWeights = ArrayUtilsCalc.sum(arrMappingWeights);
		double sumSimilarityWeighted = ArrayUtilsCalc.sum(arrSimilarityWeighted);

		avrPairwiseMappingScaled = sumSimilarityWeighted/sumMappingWeights;
				
		coverageQuery = getRatioMinimumSpanningTreeQuery(solution);
		
		coverageBase = getRatioMinimumSpanningTreeBase(solution);
		
		double coverage = coverageQuery * coverageBase;

		// double ratioNodes = Math.min(nodesQuery, nodesBase) / (double)Math.max(nodesQuery, nodesBase);

		double ratioNodesMatchQuery = Math.min(nodesQuery, heap) / (double)Math.max(nodesQuery, heap);
		double ratioNodesMatchBase = Math.min(heap, nodesBase) / (double)Math.max(heap, nodesBase);

		if(modeQuery) {
			similarity = avrPairwiseMappingScaled * coverageQuery * coverageQuery * ratioNodesMatchQuery * ratioNodesMatchQuery;
		} else {
			similarity = avrPairwiseMappingScaled * coverage * ratioNodesMatchQuery * ratioNodesMatchBase;
		}

		if(verbose) {
			StringBuilder sb = new StringBuilder();
			sb.append("ObjectiveFlexophoreHardMatchUncovered");
			sb.append(" similarity");
			sb.append("\t");
			sb.append(Formatter.format2(similarity));
			sb.append("\t");
			sb.append("avrPairwiseMappingScaled");
			sb.append("\t");
			sb.append(Formatter.format2(avrPairwiseMappingScaled));
			sb.append("\t");
			sb.append("coverage");
			sb.append("\t");
			sb.append(Formatter.format2(coverage));
			sb.append("\t");
			sb.append("ratioNodesMatchQuery");
			sb.append("\t");
			sb.append(Formatter.format2(ratioNodesMatchQuery));
			sb.append("\t");
			sb.append("ratioNodesMatchBase");
			sb.append("\t");
			sb.append(Formatter.format2(ratioNodesMatchBase));

			System.out.println(sb.toString());
		}

		deltaNanoSimilarity += System.nanoTime()-t0;

		return (float)similarity;

		// For testing
		// return (float)1.0;

	}
	public float getSimilarityHistograms(SolutionCompleteGraph solution) {

		long t0 = System.nanoTime();

		if(!validHelpersQuery){
			calculateHelpersQuery();
		}

		if(!validHelpersBase){
			calculateHelpersBase();
		}

		if(resetSimilarityArrays){
			resetSimilarityMatrices();
		}

		int heap = solution.getSizeHeap();

		//
		// the query must hit with all pharmacophore nodes
		//
		if(modeQuery) {
			if (nodesQuery != heap) {
				similarity=0;
				return (float)similarity;
			}
		}

		double sumSimDistHist = 0;

		for (int i = 0; i < heap; i++) {

			int indexNode1Query = solution.getIndexQueryFromHeap(i);

			int indexNode1Base = solution.getIndexCorrespondingBaseNode(indexNode1Query);

			for (int j = i+1; j < heap; j++) {
				int indexNode2Query = solution.getIndexQueryFromHeap(j);

				int indexNode2Base = solution.getIndexCorrespondingBaseNode(indexNode2Query);

				double simDistHist = getSimilarityHistogram(indexNode1Query, indexNode2Query, indexNode1Base, indexNode2Base);

				sumSimDistHist += simDistHist;

				if(verbose) {
					System.out.println("scorePairwiseMapping " + Formatter.format2(simDistHist));
				}
			}
		}

		double mappings = ((heap * heap)-heap) / 2.0;

		float simDistHistAvr = (float)(sumSimDistHist/mappings);

		return simDistHistAvr;

	}
	public float getSimilarityNodes(SolutionCompleteGraph solution) {

		long t0 = System.nanoTime();

		if(resetSimilarityArrays){
			resetSimilarityMatrices();
		}

		int heap = solution.getSizeHeap();

		//
		// the query must hit with all pharmacophore nodes except margin. Margin gives the number of nodes that need
		// not to hit.
		//
		if(modeQuery) {
			if (Math.abs(nodesQuery-heap)>marginQuery) {
				similarity=0;
				return (float)similarity;
			}
		}

		if(fragmentNodesMapping && heap==1){
			int indexNodeQuery = solution.getIndexQueryFromHeap(0);
			int indexNodeBase = solution.getIndexCorrespondingBaseNode(indexNodeQuery);
			similarity = getSimilarityNodes(indexNodeQuery, indexNodeBase);
			return (float)similarity;
		}

		double sumSimilarityNodesWeighted = 0;

		double sumWeights = 0;
		for (int i = 0; i < heap; i++) {
			int indexNodeQuery = solution.getIndexQueryFromHeap(i);
			double w = mdhvQuery.getWeightPharmacophorePoint(indexNodeQuery);
			int indexNodeBase = solution.getIndexCorrespondingBaseNode(indexNodeQuery);
			double similarityNodePairWeighted = getSimilarityNodes(indexNodeQuery, indexNodeBase)*w;
			sumSimilarityNodesWeighted += similarityNodePairWeighted;
			sumWeights+=w;
		}

		avrPairwiseMappingScaled = sumSimilarityNodesWeighted / sumWeights;
		coverageQuery = 0;
		coverageBase = 0;
		double ratioNodesMatchQuery = Math.min(nodesQuery, heap) / (double)Math.max(nodesQuery, heap);
		double ratioNodesMatchBase = Math.min(heap, nodesBase) / (double)Math.max(heap, nodesBase);

		if(modeQuery) {
			similarity = avrPairwiseMappingScaled * ratioNodesMatchQuery * ratioNodesMatchQuery;
		} else {
			similarity = avrPairwiseMappingScaled * ratioNodesMatchQuery * ratioNodesMatchBase;
		}

		if(verbose) {
			StringBuilder sb = new StringBuilder();
			sb.append("ObjectiveFlexophoreHardMatchUncovered");
			sb.append(" similarity nodes");
			sb.append("\t");
			sb.append(Formatter.format2(similarity));
			sb.append("\t");
			sb.append("avrPairwiseMappingScaled");
			sb.append("\t");
			sb.append(Formatter.format2(avrPairwiseMappingScaled));
			sb.append("\t");
			sb.append("ratioNodesMatchQuery");
			sb.append("\t");
			sb.append(Formatter.format2(ratioNodesMatchQuery));
			sb.append("\t");
			sb.append("ratioNodesMatchBase");
			sb.append("\t");
			sb.append(Formatter.format2(ratioNodesMatchBase));

			System.out.println(sb.toString());
		}

		deltaNanoSimilarity += System.nanoTime()-t0;

		return (float)similarity;

		// For testing
		// return (float)1.0;

	}


	/**
	 *
	 * @param solution a valid solution for IMolDistHis base and IMolDistHis query. Query and base must be set before
	 *                 starting the similarity calculation.
	 * @param indexHeap heap index of the node for which the histogram similarities to all other nodes will be calculated.
	 *                  indexHeap is the index from the list of matching query and base nodes.
	 * @return histogram similarity for a single node. Calculated as average from the sum of histogram similarities.
	 */
	public float getSimilarityHistogramsForNode(SolutionCompleteGraph solution, int indexHeap) {
		int heap = solution.getSizeHeap();

		//
		// the query must hit with all pharmacophore nodes
		//
		if(modeQuery) {
			if (nodesQuery != heap) {
				return 0;
			}
		}

		double sumPairwiseMapping = 0;

		int indexNode1Query = solution.getIndexQueryFromHeap(indexHeap);

		int indexNode1Base = solution.getIndexCorrespondingBaseNode(indexNode1Query);

		for (int i = 0; i < heap; i++) {
			if(indexHeap==i)
				continue;

			int indexNode2Query = solution.getIndexQueryFromHeap(i);

			int indexNode2Base = solution.getIndexCorrespondingBaseNode(indexNode2Query);

			double simHists = getSimilarityHistogram(indexNode1Query, indexNode2Query, indexNode1Base, indexNode2Base);

			sumPairwiseMapping += simHists;

			if(verbose) {
				System.out.println("scorePairwiseMapping " + Formatter.format2(simHists));
			}
		}

		double mappings = heap-1;

		double avrPairwiseMapping = sumPairwiseMapping/mappings;

		return (float)avrPairwiseMapping;
	}

	public long getDeltaNanoQueryBlur() {
		return deltaNanoQueryBlur;
	}

	public long getDeltaNanoBaseBlur() {
		return deltaNanoBaseBlur;
	}

	public long getDeltaNanoSimilarity() {
		return deltaNanoSimilarity;
	}

	/**
	 * Sets the color information for the visualization of the Flexophore PPPoints.
	 * Call before visualization. Method sets identical info values for corresponding nodes.
	 * @param solution
	 */
	public void setMatchingInfoInQueryAndBase(SolutionCompleteGraph solution){

		mdhvQuery.resetInfoColor();
		
		mdhvBase.resetInfoColor();
		
		int heap = solution.getSizeHeap();
		
		for (int i = 0; i < heap; i++) {
			
			int indexNodeQuery = solution.getIndexQueryFromHeap(i);
			
			int indexNodeBase = solution.getIndexCorrespondingBaseNode(indexNodeQuery);
			
			double similarityMappingNodes = getSimilarityNodes(indexNodeQuery, indexNodeBase);

			mdhvQuery.setSimilarityMappingNodes(indexNodeQuery, (float)similarityMappingNodes);
						
			mdhvQuery.setMappingIndex(indexNodeQuery, i);
			
			mdhvBase.setMappingIndex(indexNodeBase, i);
			
			mdhvBase.setSimilarityMappingNodes(indexNodeBase, (float)similarityMappingNodes);
		}

	}

	@Override
	public boolean isModeFragment() {
		return fragmentNodesMapping;
	}


	public IMolDistHist getBase() {
		return mdhvBase;
	}

	public IMolDistHist getQuery() {
		return mdhvQuery;
	}

	public void setBase(IMolDistHist iMolDistHistBase) {

		if(iMolDistHistBase.getNumPPNodes()>=ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE){
			throw new RuntimeException("Number of base pharmacophore nodes (" +iMolDistHistBase.getNumPPNodes() + ") exceeds limit of " + ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE + ".");
		}

		long t0 = System.nanoTime();

		if(iMolDistHistBase instanceof MolDistHistViz) {
			mdhvBase = (MolDistHistViz) iMolDistHistBase;
			mdhvBaseBlurredHist = new MolDistHistViz((MolDistHistViz) iMolDistHistBase);
		} else if(iMolDistHistBase instanceof MolDistHist) {
			mdhvBase = new MolDistHistViz((MolDistHist) iMolDistHistBase);
			mdhvBaseBlurredHist = new MolDistHistViz((MolDistHist) iMolDistHistBase);
		}


		if((slidingWindowDistHist!=null) && !fragmentNodesMapping)
			slidingWindowDistHist.apply(mdhvBaseBlurredHist);

		nodesBase = iMolDistHistBase.getNumPPNodes();
		
		validHelpersBase = false;
		
		resetSimilarityArrays = true;

		if(!checkAtomTypes(mdhvBase)) {
			throw new RuntimeException("Base contains Invalid atom type for similarity calculation " + mdhvBase.getMolDistHist().toString() + ".");
		}

		deltaNanoBaseBlur += System.nanoTime()-t0;
	}

	public void setSlidingWindowDistHistNull() {
		this.slidingWindowDistHist = null;
	}

	public void setQuery(IMolDistHist iMolDistHistQuery) {

		if(iMolDistHistQuery.getNumPPNodes()>=ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE){
			throw new RuntimeException("Number of query pharmacophore nodes (" +iMolDistHistQuery.getNumPPNodes() + ") exceeds limit of " + ConstantsFlexophore.MAX_NUM_NODES_FLEXOPHORE + ".");
		}

		long t0 = System.nanoTime();

		if(iMolDistHistQuery instanceof MolDistHistViz) {
			mdhvQuery = (MolDistHistViz) iMolDistHistQuery;
			mdhvQueryBlurredHist = new MolDistHistViz((MolDistHistViz) iMolDistHistQuery);
		} else if(iMolDistHistQuery instanceof MolDistHist) {
			mdhvQuery = new MolDistHistViz((MolDistHist) iMolDistHistQuery);
			mdhvQueryBlurredHist = new MolDistHistViz((MolDistHist) iMolDistHistQuery);
		}

		if((slidingWindowDistHist!=null) && !fragmentNodesMapping)
			slidingWindowDistHist.apply(mdhvQueryBlurredHist);

		nodesQuery = iMolDistHistQuery.getNumPPNodes();

		arrMandatoryPPPoint = new boolean[nodesQuery];
		for (int i = 0; i < arrMandatoryPPPoint.length; i++) {
			arrMandatoryPPPoint[i]=iMolDistHistQuery.isMandatoryPharmacophorePoint(i);
		}
		numMandatoryPPPoints = iMolDistHistQuery.getNumMandatoryPharmacophorePoints();
		
		validHelpersQuery = false;
		
		resetSimilarityArrays = true;

		if(!checkAtomTypes(mdhvQuery)) {
			throw new RuntimeException("Base contains Invalid atom type for similarity calculation " + mdhvQuery.getMolDistHist().toStringNodes() + ".");
		}

		deltaNanoQueryBlur += System.nanoTime()-t0;
	}

	private boolean checkAtomTypes(MolDistHistViz mdhv) {

		boolean valid = true;
		List<PPNodeViz> liNode = mdhv.getNodes();

		stop:
		for (PPNodeViz ppNodeViz : liNode) {
			int c = ppNodeViz.getInteractionTypeCount();
			for (int i = 0; i < c; i++) {
				int type = ppNodeViz.getInteractionType(i);
				if(!nodeSimilarity.isValidType(type)){
					valid = false;
					break stop;
				}
			}
		}

		return valid;
	}

	
	private void calculateHelpersQuery(){
		
		arrRelativeDistanceMatrixQuery = calculateRelativeDistanceMatrix(mdhvQueryBlurredHist);
		
		maHelperAdjacencyQuery = new Matrix(arrRelativeDistanceMatrixQuery);
		
		MinimumSpanningTree mst = new MinimumSpanningTree(maHelperAdjacencyQuery);
		
		Matrix maMST = mst.getMST();
		
		sumDistanceMinSpanTreeQuery = maMST.getSumUpperTriangle();
		
		validHelpersQuery = true;
	}


	private void calculateHelpersBase(){
		
		arrRelativeDistanceMatrixBase = calculateRelativeDistanceMatrix(mdhvBaseBlurredHist);

		maHelperAdjacencyBase = new Matrix(arrRelativeDistanceMatrixBase);
		
		
		MinimumSpanningTree mst = new MinimumSpanningTree(maHelperAdjacencyBase);
		
		Matrix maMST = mst.getMST();
		
		sumDistanceMinSpanTreeBase = maMST.getSumUpperTriangle();
		
		validHelpersBase = true;
	}
	
	private double getRatioMinimumSpanningTreeQuery(SolutionCompleteGraph solution) {
		
		double ratioCovered2Total = 0;
		
		int heap = solution.getSizeHeap();

		maHelperAdjacencyQuery.set(Double.NaN);
		
		for (int i = 0; i < heap; i++) {
			
			int indexNode1 = solution.getIndexQueryFromHeap(i);

			for (int j = i+1; j < heap; j++) {
				
				int indexNode2 = solution.getIndexQueryFromHeap(j);
				
				maHelperAdjacencyQuery.set(indexNode1, indexNode2, arrRelativeDistanceMatrixQuery[indexNode1][indexNode2]);
				maHelperAdjacencyQuery.set(indexNode2, indexNode1, arrRelativeDistanceMatrixQuery[indexNode1][indexNode2]);
				
			}
		}
		
		MinimumSpanningTree mstSolution = new MinimumSpanningTree(maHelperAdjacencyQuery);
		
		Matrix maMST = mstSolution.getMST();
		
		double sum = maMST.getSumUpperTriangle();
				
		// ratioCovered2Total = (sum*sum)/(sumDistanceMinSpanTreeQuery*sumDistanceMinSpanTreeQuery);
		
		double sumMSTSquared = sum*sum;
		
		double sumMSTSQuerySquared = sumDistanceMinSpanTreeQuery*sumDistanceMinSpanTreeQuery;
		
		ratioCovered2Total = Math.min(sumMSTSquared, sumMSTSQuerySquared)/Math.max(sumMSTSquared, sumMSTSQuerySquared);
	
		return ratioCovered2Total;
	}
	
	private double getRatioMinimumSpanningTreeBase(SolutionCompleteGraph solution) {
		
		double ratioCovered2Total = 0;
		
		int heap = solution.getSizeHeap();

		maHelperAdjacencyBase.set(Double.NaN);
		
		for (int i = 0; i < heap; i++) {
			
			int indexNode1 = solution.getIndexBaseFromHeap(i);

			for (int j = i+1; j < heap; j++) {
				
				int indexNode2 = solution.getIndexBaseFromHeap(j);
				
				maHelperAdjacencyBase.set(indexNode1, indexNode2, arrRelativeDistanceMatrixBase[indexNode1][indexNode2]);
				maHelperAdjacencyBase.set(indexNode2, indexNode1, arrRelativeDistanceMatrixBase[indexNode1][indexNode2]);
				
			}
		}
		
		MinimumSpanningTree mstSolution = new MinimumSpanningTree(maHelperAdjacencyBase);
		
		Matrix maMST = mstSolution.getMST();
		
		double sum = maMST.getSumUpperTriangle();
		
		double sumMSTSquared = sum*sum;
		
		double sumMSTSBaseSquared = sumDistanceMinSpanTreeBase*sumDistanceMinSpanTreeBase;
		
		ratioCovered2Total = Math.min(sumMSTSquared, sumMSTSBaseSquared)/Math.max(sumMSTSquared, sumMSTSBaseSquared);
		// ratioCovered2Total = (sum)/(sumDistanceMinSpanTreeBase);
		
		return ratioCovered2Total;
	}
	
	/**
	 * Calculate a distance matrix from the center of gravity distance bins. 
	 * The values are standardized by dividing them by the highest center of gravity value.  
	 * @param mdh
	 * @return
	 */
	private double [][] calculateRelativeDistanceMatrix(IMolDistHist mdh){
		
		int nodes = mdh.getNumPPNodes();
		
		double [][] arrDist = new double [nodes][nodes];
		
		int maxMedianDistanceBin = 0; 
		
		for (int i = 0; i < arrDist.length; i++) {
			
			for (int j = i+1; j < arrDist.length; j++) {
				
				int medianDistanceBin = getCenterOfGravityDistanceBin(mdh, i, j);
				
				arrDist[i][j] = medianDistanceBin;
				arrDist[j][i] = medianDistanceBin;
				
				if(medianDistanceBin > maxMedianDistanceBin){
					maxMedianDistanceBin = medianDistanceBin;
				}
			}
		}
		
		for (int i = 0; i < arrDist.length; i++) {
			
			for (int j = i+1; j < arrDist.length; j++) {
				
				arrDist[i][j] = arrDist[i][j] / maxMedianDistanceBin;
				
				arrDist[j][i] = arrDist[j][i] / maxMedianDistanceBin;
				
			}
		}
		
		return arrDist;
	}

	/**
	 * 
	 * @param mdh
	 * @param indexNode1
	 * @param indexNode2
	 * @return the index of the distance bin with the center of gravity for the histogram values.
	 */
	private int getCenterOfGravityDistanceBin(IMolDistHist mdh, int indexNode1, int indexNode2) {
		
		byte [] arr = mdh.getDistHist(indexNode1, indexNode2, arrTmpHist);

		double sum=0;
		
		for (int i = 0; i < arr.length; i++) {
				sum += arr[i]; 
		}
		
		double center = sum / 2.0;
		
		sum=0;
		
		int bin = -1;
		for (int i = arr.length-1; i >= 0; i--) {
			
			sum += arr[i]; 
			
			if(sum >= center) {
				bin=i;
				break;
			}
		}
		
		return bin;
	}

	/**
	 * Compares nodes and histograms
	 * @param indexNode1Query
	 * @param indexNode2Query
	 * @param indexNode1Base
	 * @param indexNode2Base
	 * @return
	 */
	private double getScorePairwiseMapping(int indexNode1Query, int indexNode2Query, int indexNode1Base, int indexNode2Base) {
		double score = 0;

		double simNodePair1 = getSimilarityNodes(indexNode1Query, indexNode1Base);
		
		double simNodePair2 = getSimilarityNodes(indexNode2Query, indexNode2Base);
		
		double simHists = getSimilarityHistogram(indexNode1Query, indexNode2Query, indexNode1Base, indexNode2Base);

		if(optimisticHistogramSimilarity) {
			if (simHists >= OPTIMISTIC_HISTOGRAM_THRESH) {
				simHists = 1.0;
			}
		}

		if(verbose){
			System.out.println("simHists " + Formatter.format2(simHists));
		}
		score = simNodePair1 * simNodePair1 * simNodePair2 * simNodePair2 * simHists * simHists;
		return score;
	}

	public double getSimilarityNodes(int indexNodeQuery, int indexNodeBase) {
		if(arrSimilarityNodes[indexNodeQuery][indexNodeBase] < 0 || verbose){
			float similarity = (float)nodeSimilarity.getSimilarity(mdhvQueryBlurredHist.getNode(indexNodeQuery), mdhvBaseBlurredHist.getNode(indexNodeBase));
			arrSimilarityNodes[indexNodeQuery][indexNodeBase]=similarity;
		} 
		return arrSimilarityNodes[indexNodeQuery][indexNodeBase];
	}
	
	public float getSimilarityHistogram(int indexNode1Query, int indexNode2Query, int indexNode1Base, int indexNode2Base) {

		int indexHistogramQuery = DistHist.getIndex(indexNode1Query, indexNode2Query, nodesQuery);
		int indexHistogramBase = DistHist.getIndex(indexNode1Base, indexNode2Base, nodesBase);
		if(arrSimilarityHistograms[indexHistogramQuery][indexHistogramBase] < 0){
			float similarityHistogram = 0;
			similarityHistogram =
					(float)HistogramMatchCalculator.getFractionOverlappingBins(
							mdhvQueryBlurredHist, indexNode1Query, indexNode2Query, mdhvBaseBlurredHist, indexNode1Base, indexNode2Base);
			arrSimilarityHistograms[indexHistogramQuery][indexHistogramBase]=similarityHistogram;
		}

		return arrSimilarityHistograms[indexHistogramQuery][indexHistogramBase];
	}
	
	public double getSimilarityNodes(IPPNode query, IPPNode base) {
		return nodeSimilarity.getSimilarity(query, base);
	}

	/**
	 * @return
	 */
	public String toStringRecentSimilarityResults() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("ObjectiveFlexophoreHardMatchUncovered toStringRecentSimilarityResults()");
		sb.append("avr pairwise mapping " + Formatter.format3(avrPairwiseMappingScaled) + "\n");
		sb.append("coverage query " + Formatter.format3(coverageQuery) + "\n");
		sb.append("coverage base " + Formatter.format3(coverageBase) + "\n");
		sb.append("similarity " + Formatter.format3(similarity));
		
		return sb.toString();
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		
		Matrix maSimNodes = new Matrix(arrSimilarityNodes);
		
		int rowEnd = -1;
		for (int i = 0; i < maSimNodes.rows(); i++) {
			if(maSimNodes.get(i, 0)<0){
				rowEnd = i;
				break;
			}
		}
		
		int colEnd = -1;
		for (int i = 0; i < maSimNodes.cols(); i++) {
			if(maSimNodes.get(0, i)<0){
				colEnd = i;
				break;
			}
		}
		
		for (int i = 0; i < rowEnd; i++) {
			for (int j = 0; j < colEnd; j++) {
				sb.append(Formatter.format2(maSimNodes.get(i, j)));
				sb.append("  ");
			}
			sb.append("\n");
		}
		
		return sb.toString();
	}
	
}
