package com.actelion.research.chem.descriptor.flexophore.completegraphmatcher;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.graph.MinimumSpanningTree;
import com.actelion.research.chem.descriptor.DescriptorHandlerFlexophore;
import com.actelion.research.chem.descriptor.flexophore.*;
import com.actelion.research.chem.descriptor.flexophore.generator.ConstantsFlexophoreGenerator;
import com.actelion.research.util.Formatter;
import com.actelion.research.util.graph.complete.IObjectiveCompleteGraph;
import com.actelion.research.util.graph.complete.SolutionCompleteGraph;

import java.util.Arrays;

/**
 * 
 * 
 * ObjectiveFlexophoreHardMatchUncovered
 * The weighting of the coverage is hard. Which means that uncovered nodes 
 * strongly change the final similarity score.
 * look in <code>getScoreUncoveredNearestNodesBase(SolutionCompleteGraph solution)</code> 
 * and <code>getScoreUncoveredNearestNodesQuery(SolutionCompleteGraph solution)</code>.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * Oct 2, 2012 MvK: Start implementation
 * Mar 3. 2016 MvK: updates. Lowered thresh for histogram similarity.
 * Mar 31. 2020 MvK: fraction of carbon is considered in pharmacophore node similarity.
 */
public class ObjectiveFlexophoreHardMatchUncovered implements IObjectiveCompleteGraph<IMolDistHist>{

	public static final String VERSION = "02.04.2020 08:00";
	public static final String INFO = "";

	public static final int MAX_NUM_NODES_FLEXOPHORE = 64;

	// 0.15 best thresh tested on 31.03.2016
	// private final static double THRESH_HISTOGRAM_SIMILARITY	= 0.15;
	public final static double THRESH_HISTOGRAM_SIMILARITY	= 0.1;

	// The thresh for the node similarity depends on the number of interaction types in the node.
	// 03.03.2016 Top result so far for 0.9
	// 13.04.2020 Maybe obsolete
	// ToDo
	final static double THRESH_NODE_SIMILARITY_START = 0.5;



	private static final float INIT_VAL = -1;

	
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

	// private IPPNodeSimilarity nodeSimilarity;
	private PPNodeSimilarity nodeSimilarity;

	boolean verbose;

	// private ScaleClasses scaleClassesSimilarityNodes;
	
	// private ScaleClasses scaleClassesFinalSimilarity;

	private float [][] arrSimilarityNodes;
	
	private float [][] arrSimilarityHistograms;
	
	private double [][] arrRelativeDistanceMatrixQuery;
	
	private double [][] arrRelativeDistanceMatrixBase;
	
	private Matrix maHelperAdjacencyQuery;
	
	private Matrix maHelperAdjacencyBase;
	
	private double sumDistanceMinSpanTreeQuery;
	
	private double sumDistanceMinSpanTreeBase;

	private int numInevitablePPPoints;
	
	private double avrPairwiseMappingScaled;
	
	private double coverageQuery;
	
	private double coverageBase;
	
	private double similarity;

	private SlidingWindowDistHist slidingWindowDistHist;

	public ObjectiveFlexophoreHardMatchUncovered(){
		this(DescriptorHandlerFlexophore.VERSION_INTERACTION_TABLES,
				DescriptorHandlerFlexophore.MODE_PPNODE_SIMILARITY_COMPARISON,
				DescriptorHandlerFlexophore.THRESH_SIMILARITY_COMPARISON_NODE, THRESH_HISTOGRAM_SIMILARITY);

	}



	public ObjectiveFlexophoreHardMatchUncovered(
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

		initSimilarityMatrices();
	}
	
	private void initSimilarityMatrices(){
		
		arrSimilarityNodes = new float [MAX_NUM_NODES_FLEXOPHORE][];
		for (int i = 0; i < MAX_NUM_NODES_FLEXOPHORE; i++) {
			arrSimilarityNodes[i] = new float [MAX_NUM_NODES_FLEXOPHORE];
			Arrays.fill(arrSimilarityNodes[i], INIT_VAL);
			
		}

		int maxNumHistograms = ((MAX_NUM_NODES_FLEXOPHORE*MAX_NUM_NODES_FLEXOPHORE)-MAX_NUM_NODES_FLEXOPHORE)/2;
		
		arrSimilarityHistograms = new float [maxNumHistograms][];
		for (int i = 0; i < maxNumHistograms; i++) {
			arrSimilarityHistograms[i] = new float [maxNumHistograms];
			Arrays.fill(arrSimilarityHistograms[i], INIT_VAL);
		}
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
		if(numInevitablePPPoints > 0) {
			
			int ccInevitablePPPointsInSolution = 0;
			
			for (int i = 0; i < heap; i++) {
				int indexNodeQuery = solution.getIndexQueryFromHeap(i);
				
				if(mdhvQueryBlurredHist.isInevitablePharmacophorePoint(indexNodeQuery)){
					ccInevitablePPPointsInSolution++;
				}
			}
			
			int neededMinInevitablePPPoints = Math.min(heap, numInevitablePPPoints);
			
			if(ccInevitablePPPointsInSolution < neededMinInevitablePPPoints){
				mapping = false;
			}
			
		}

		//
		// Check for one hetero atom in solution
		//
		if(mapping){
			
			boolean heteroNodeQuery = false;
			
			boolean heteroNodeBase = false;

			for (int i = 0; i < heap; i++) {
				
				int indexNodeQuery = solution.getIndexQueryFromHeap(i);
				PPNode nodeQuery = mdhvQueryBlurredHist.getNode(indexNodeQuery);
				if(nodeQuery.hasHeteroAtom()){
					heteroNodeQuery = true;
				}
				
				int indexNodeBase = solution.getIndexCorrespondingBaseNode(indexNodeQuery);
				
				PPNode nodeBase = mdhvBaseBlurredHist.getNode(indexNodeBase);
				if(nodeBase.hasHeteroAtom()){
					heteroNodeBase = true;
				}

				if(heteroNodeQuery && heteroNodeBase){
					break;
				}

			}
			
			if(!heteroNodeQuery || !heteroNodeBase) {
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
		if(mapping){
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

	public float getSimilarity(SolutionCompleteGraph solution) {

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

		double sumPairwiseMapping = 0;

		// double productPairwiseMapping = 0;

		for (int i = 0; i < heap; i++) {
			
			int indexNode1Query = solution.getIndexQueryFromHeap(i);
			
			int indexNode1Base = solution.getIndexCorrespondingBaseNode(indexNode1Query);
			
			for (int j = i+1; j < heap; j++) {
				int indexNode2Query = solution.getIndexQueryFromHeap(j);
				
				int indexNode2Base = solution.getIndexCorrespondingBaseNode(indexNode2Query);

				double scorePairwiseMapping = getScorePairwiseMapping(indexNode1Query, indexNode2Query, indexNode1Base, indexNode2Base);

				sumPairwiseMapping += scorePairwiseMapping;

				if(verbose) {
					System.out.println("scorePairwiseMapping " + Formatter.format2(scorePairwiseMapping));
				}
			}
		}
		
		double mappings = ((heap * heap)-heap) / 2.0;
	
		avrPairwiseMappingScaled = sumPairwiseMapping/mappings;
				
		coverageQuery = getRatioMinimumSpanningTreeQuery(solution);
		
		coverageBase = getRatioMinimumSpanningTreeBase(solution);
		
		double coverage = coverageQuery * coverageBase;

		// double ratioNodes = Math.min(nodesQuery, nodesBase) / (double)Math.max(nodesQuery, nodesBase);

		double ratioNodesMatchQuery = Math.min(nodesQuery, heap) / (double)Math.max(nodesQuery, heap);
		double ratioNodesMatchBase = Math.min(heap, nodesBase) / (double)Math.max(heap, nodesBase);

		similarity = avrPairwiseMappingScaled * coverage * ratioNodesMatchQuery * ratioNodesMatchBase;

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

		return (float)similarity;

		// For testing
		// return (float)1.0;

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
	
	
	public IMolDistHist getBase() {
		return mdhvBase;
	}

	public IMolDistHist getQuery() {
		return mdhvQuery;
	}

	public void setBase(IMolDistHist iMolDistHistBase) {

		mdhvBase = (MolDistHistViz)iMolDistHistBase;
		mdhvBaseBlurredHist = new MolDistHistViz((MolDistHistViz)iMolDistHistBase);

		slidingWindowDistHist.apply(mdhvBaseBlurredHist);

		nodesBase = iMolDistHistBase.getNumPPNodes();
		
		validHelpersBase = false;
		
		resetSimilarityArrays = true;
	}
	
	public void setQuery(IMolDistHist iMolDistHistQuery) {
		
		mdhvQuery = (MolDistHistViz) iMolDistHistQuery;
		mdhvQueryBlurredHist = new MolDistHistViz((MolDistHistViz)iMolDistHistQuery);

		slidingWindowDistHist.apply(mdhvQueryBlurredHist);

		nodesQuery = iMolDistHistQuery.getNumPPNodes();
		
		numInevitablePPPoints = iMolDistHistQuery.getNumInevitablePharmacophorePoints();
		
		validHelpersQuery = false;
		
		resetSimilarityArrays = true;
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
		
		if(simHists==0){

			System.out.println("ObjectiveFlexophoreHardMatchUncovered getScorePairwiseMapping(int indexNode1Query, int indexNode2Query, int indexNode1Base, int indexNode2Base)");

			System.out.println("Sim hists = 0");
		}

		if(verbose){
			System.out.println("simHists " + Formatter.format2(simHists));
		}

		// score = simNodePair1 * simNodePair1 * simNodePair2 * simNodePair2 * simHists * simHists * simHists;

		score = simNodePair1 * simNodePair1 * simNodePair2 * simNodePair2 * simHists * simHists;


		return score;
	}

	private float getSimilarityNodes(int indexNodeQuery, int indexNodeBase) {
		
		if(arrSimilarityNodes[indexNodeQuery][indexNodeBase] < 0 || verbose){
			
			float similarity = (float)nodeSimilarity.getSimilarity(mdhvQueryBlurredHist.getNode(indexNodeQuery), mdhvBaseBlurredHist.getNode(indexNodeBase));
			
			arrSimilarityNodes[indexNodeQuery][indexNodeBase]=similarity;
		} 
		
		return arrSimilarityNodes[indexNodeQuery][indexNodeBase];
	}
	
	private float getSimilarityHistogram(int indexNode1Query, int indexNode2Query, int indexNode1Base, int indexNode2Base) {

		int indexHistogramQuery = DistHist.getIndex(indexNode1Query, indexNode2Query, nodesQuery);

		int indexHistogramBase = DistHist.getIndex(indexNode1Base, indexNode2Base, nodesBase);

		if(arrSimilarityHistograms[indexHistogramQuery][indexHistogramBase] < 0){

			float similarityHistogram =
					(float)HistogramMatchCalculator.getSimilarity((MolDistHistViz) mdhvQueryBlurredHist, indexNode1Query, indexNode2Query, (MolDistHistViz) mdhvBaseBlurredHist, indexNode1Base, indexNode2Base);

			arrSimilarityHistograms[indexHistogramQuery][indexHistogramBase]=similarityHistogram;
		}

		return arrSimilarityHistograms[indexHistogramQuery][indexHistogramBase];
	}
	
	public double getSimilarityNodes(PPNode query, PPNode base) {
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
