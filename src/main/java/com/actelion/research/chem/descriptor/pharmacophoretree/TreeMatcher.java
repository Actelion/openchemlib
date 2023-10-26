package com.actelion.research.chem.descriptor.pharmacophoretree;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;
import java.util.stream.Collectors;

import com.actelion.research.chem.descriptor.pharmacophoretree.TreeMatcher.FeatureMatch;


/**
 * Finds the optimal matching of nodes between two PharmacophoreTrees employing a dynamic programing scheme
 * termed as "match-search" algorithm in the original publication (DOI:10.1023/a:1008068904628).
 * Given a set of initial splits (a split is created by cutting an edge each in both trees), both pair of subtrees
 * resulting from the split are used as an input for the recursive extensionMatch method. Starting from the head 
 * nodes of the subtrees, extension matches up to a certain size containing the head node and connecting nodes are enumerated.
 * The best-scoring extension match are further considered. The resulting cut subtrees from the extension cuts form possible matches
 * and are also used as an input for recursive calls of the match-search routine. The highest-scoring assignments of
 * the subtrees are determined using a Hungarian Algorithm. If a match between two subtrees fulfills size and balance criteria, the recursion is stopped. 
 * The assignments of the nodes in the matching and the score of a subtree match are stored in TreeMatching objects. These
 * objects are stored in a dynamic programming matrix to facilitate memoization of the results.
 * @author joel
 *
 */


public class TreeMatcher {
	
	
	public static final int EXTENSION_MATCHES = 3; //number of explicitly considered extension matches at every recursion step
	public static final double ALPHA = 0.8; //weighting of source-tree match vs extension-tree match, takes values from 0 to 1
	public static final double NULL_MATCH_SCALING = 0.5; 
	public static final double SIMILARITY_SCALING_SPLIT_SCORE = 0.6;
	public static final double MATCH_BALANCE = 2.0; //named beta in the original publication
	public static final double MATCH_SIZE_LIMIT = 3.0;
	public static final int MATCH_NODE_NR_LIMIT = 2;
	public static final int EXTENSION_MATCH_NODE_NR_LIMIT = 3;
	public static final int INITIAL_SPLITS = 5;
	public static final double SIZE_RATIO = 2; //if ratio of sizes (nr of atoms) of two matches differs by more than that, the similarity is zero 
	private TreeMatching[][] dpMatchMatrix;
	private PharmacophoreTree queryTree;
	private PharmacophoreTree baseTree;
	private List<PharmacophoreNode> queryNodes;
	private List<PharmacophoreNode> baseNodes;

	public TreeMatcher(PharmacophoreTree queryTree, PharmacophoreTree baseTree) {
		
		this.queryTree = queryTree;
		this.baseTree = baseTree;
		queryNodes = queryTree.getNodes();
		baseNodes = baseTree.getNodes();
		dpMatchMatrix = new TreeMatching[2*queryTree.getEdges().size()][2*baseTree.getEdges().size()];
		
	}
	/**
	 * finds a set of balanced, high-scoring initial splits that are then used as an input for the extension match
	 * algorithm
	 * @return
	 */
	public TreeMatching matchSearch() {
		// search for initial splits
		double bestScore = 0.0;
		TreeMatching bestMatch = new TreeMatching();
		int[][] splits = findInitialSplits();
		for(int[] split : splits) {
			int index1 = split[0];
			int cut1 = PharmacophoreTree.CUT_LEFT;
			
			int index2 = split[1]/2;
			int cut2 = split[1]%2 == 0 ? PharmacophoreTree.CUT_LEFT : PharmacophoreTree.CUT_RIGHT;
		
			List<Integer> sourceTreeEdges1 = new ArrayList<Integer>();
			List<Integer> targetTreeEdges1 = new ArrayList<Integer>();
			List<Integer> sourceTreeEdgeParents1 = new ArrayList<Integer>();
			List<Integer> targetTreeEdgeParents1 = new ArrayList<Integer>();
			int [] headNodes1 = queryTree.initialCut(cut1, index1,sourceTreeEdges1, sourceTreeEdgeParents1, targetTreeEdges1, targetTreeEdgeParents1);
			List<Integer> sourceTreeEdges2 = new ArrayList<Integer>();
			List<Integer> targetTreeEdges2 = new ArrayList<Integer>();
			List<Integer> sourceTreeEdgeParents2 = new ArrayList<Integer>();
			List<Integer> targetTreeEdgeParents2 = new ArrayList<Integer>();
			int [] headNodes2 = baseTree.initialCut(cut2, index2 , sourceTreeEdges2, sourceTreeEdgeParents2, targetTreeEdges2, targetTreeEdgeParents2);		


			//compare source trees
			
			TreeMatching match1 = extensionMatch(headNodes1[0], headNodes2[0],index1,index2,cut1,cut2,
					sourceTreeEdges1,sourceTreeEdges2,sourceTreeEdgeParents1, sourceTreeEdgeParents2);


			
			//compare target trees
			
			TreeMatching match2 = extensionMatch(headNodes1[1], headNodes2[1],index1,index2,cut1*-1,cut2*-1,
					targetTreeEdges1,targetTreeEdges2,targetTreeEdgeParents1, targetTreeEdgeParents2);

			match1.addMatching(match2);
			match1.calculate();
			if(match1.sim>bestScore) {
				bestScore = match1.sim;
				bestMatch = match1;
			}

			

		}
		bestMatch.calculate();

		return bestMatch;
		}
		
		public int[][] findInitialSplits() {
			int[] cuts = {PharmacophoreTree.CUT_LEFT,PharmacophoreTree.CUT_RIGHT};
			double[][] initialSplitScores = new double[queryTree.getEdges().size()][2*baseTree.getEdges().size()];
			for(int i=0;i<queryTree.getEdges().size();i++) {
					List<Integer> sourceTreeEdges1 = new ArrayList<Integer>();
					List<Integer> targetTreeEdges1 = new ArrayList<Integer>();
					List<Integer> sourceTreeEdgeParents1 = new ArrayList<Integer>();
					List<Integer> targetTreeEdgeParents1 = new ArrayList<Integer>();
					int [] headNodes1 = queryTree.initialCut(cuts[0],i,sourceTreeEdges1, sourceTreeEdgeParents1, targetTreeEdges1, targetTreeEdgeParents1);
			
					for(int j=0;j<baseTree.getEdges().size();j++) {
						for(int cut2: cuts) {
							List<Integer> sourceTreeEdges2 = new ArrayList<Integer>();
							List<Integer> targetTreeEdges2 = new ArrayList<Integer>();
							List<Integer> sourceTreeEdgeParents2 = new ArrayList<Integer>();
							List<Integer> targetTreeEdgeParents2 = new ArrayList<Integer>();
							int [] headNodes2 = baseTree.initialCut(cut2, j , sourceTreeEdges2, sourceTreeEdgeParents2, targetTreeEdges2, targetTreeEdgeParents2);		
							Collection<Integer> nodesQuerySource = queryTree.getNodesFromEdges(sourceTreeEdges1);
							nodesQuerySource.add(headNodes1[0]);
							Collection<Integer> nodesQueryTarget = queryTree.getNodesFromEdges(targetTreeEdges1);
							nodesQueryTarget.add(headNodes1[1]);
							
							Collection<Integer> nodesBaseSource = baseTree.getNodesFromEdges(sourceTreeEdges2);
							nodesBaseSource.add(headNodes2[0]);
							Collection<Integer> nodesBaseTarget = baseTree.getNodesFromEdges(targetTreeEdges2);
							nodesBaseTarget.add(headNodes2[1]);
							int index1 = i;
							int index2 = cut2 == PharmacophoreTree.CUT_LEFT ? j*2 : j*2+1;

							initialSplitScores[index1][index2] = getSplitScore(queryTree,baseTree,nodesQuerySource,nodesBaseSource,
									nodesQueryTarget,nodesBaseTarget);

						}
					}
					
			}
			int[][] bestSplits = new int[INITIAL_SPLITS][2];
			double[] bestScores = new double[INITIAL_SPLITS];
			TreeUtils.retrieveHighestValuesFrom2DArray(initialSplitScores, bestScores, bestSplits);
			return bestSplits;
			
		}
		
		/**
		 * recursive extension match procedure, taking two subtrees as an input
		 * @param headNode1
		 * @param headNode2
		 * @param cutEdge1
		 * @param cutEdge2
		 * @param cutDir1
		 * @param cutDir2
		 * @param subTreeEdgeIndeces1
		 * @param subTreeEdgeIndeces2
		 * @param subTreeEdgeParentIndeces1
		 * @param subTreeEdgeParentIndeces2
		 * @return
		 */
		private TreeMatching extensionMatch(int headNode1, int headNode2, int cutEdge1, int cutEdge2, int cutDir1, int cutDir2,List<Integer> subTreeEdgeIndeces1, 
				List<Integer> subTreeEdgeIndeces2,List<Integer> subTreeEdgeParentIndeces1, List<Integer> subTreeEdgeParentIndeces2) {
			TreeMatching matching;

			
			int index1 = cutDir1 == PharmacophoreTree.CUT_LEFT ? cutEdge1*2 : cutEdge1*2+1;
			int index2 = cutDir2 == PharmacophoreTree.CUT_LEFT ? cutEdge2*2 : cutEdge2*2+1;
			Collection<Integer> nodes1 = queryTree.getNodesFromEdges(subTreeEdgeIndeces1);
			nodes1.add(headNode1);
			Collection<Integer> nodes2 = baseTree.getNodesFromEdges(subTreeEdgeIndeces2);
			nodes2.add(headNode2);
			if(dpMatchMatrix[index1][index2]!= null)
				matching = dpMatchMatrix[index1][index2];

			else {

				List<FeatureMatch> matches = assessMatch(nodes1,nodes2); 
				
				if(matches!=null) {
					matching = new TreeMatching();
					for(FeatureMatch fmatch : matches)
						matching.addFeatureMatch(fmatch);
					matching.calculate();
					dpMatchMatrix[index1][index2] = matching;
					} // go into recursion
				else {

					List<int[]> cuts1 = queryTree.getExtensionCuts(subTreeEdgeIndeces1,subTreeEdgeParentIndeces1);
					List<int[]> cuts2 = baseTree.getExtensionCuts(subTreeEdgeIndeces2,subTreeEdgeParentIndeces2);
					double[][] scores = new double[cuts1.size()][cuts2.size()];
					for(int i=0;i<cuts1.size();i++) {
						int[] cut1 = cuts1.get(i);
						Set<Integer> extensionNodes1 = new HashSet<Integer>();
						Set<Integer> sourceNodes1 = new HashSet<Integer>();
						queryTree.enumerateExtensionCutFast(headNode1,cut1, subTreeEdgeIndeces1,
							extensionNodes1, sourceNodes1);
						for(int j=0;j<cuts2.size();j++) {
							int[] cut2 = cuts2.get(j);
							Set<Integer> extensionNodes2 = new HashSet<Integer>();
							Set<Integer> sourceNodes2 = new HashSet<Integer>();
							baseTree.enumerateExtensionCutFast(headNode2,cut2, subTreeEdgeIndeces2,
									extensionNodes2, sourceNodes2);
							scores[i][j] = scoreExtensionMatch(queryTree,baseTree,extensionNodes1, extensionNodes2,
										sourceNodes1, sourceNodes2);

						}
					}
					int[][] bestCuts = new int[cuts1.size()*cuts2.size()][2];
					double[] bestScores = new double[cuts1.size()*cuts2.size()];
					TreeUtils.retrieveHighestValuesFrom2DArray(scores, bestScores, bestCuts);
					double bestScore = -Double.MAX_VALUE;
					TreeMatching bestMatching = null;
					int counter = 0;
					for(int[] cut:bestCuts) {
						if(counter>EXTENSION_MATCHES)
							break;
	
						int[] cut1 = cuts1.get(cut[0]);
						int[] cut2 = cuts2.get(cut[1]);
						List<List<Integer>> sourceTreeEdgeIndeces1 = new ArrayList<List<Integer>>();
						List<List<Integer>> sourceTreeEdgeParentIndeces1 = new ArrayList<List<Integer>>();
						List<Integer> sourceTreeHeadNodes1 = new ArrayList<Integer>();
						Set<Integer> extensionNodes1 = new HashSet<Integer>();
						List<Integer> cutEdges1 = new ArrayList<Integer>();
						List<Integer> cutDirs1 = new ArrayList<Integer>();
						queryTree.enumerateExtensionCutFull(headNode1,cut1, subTreeEdgeIndeces1,
								subTreeEdgeParentIndeces1, sourceTreeEdgeIndeces1,sourceTreeEdgeParentIndeces1, 
								sourceTreeHeadNodes1,extensionNodes1, cutEdges1, cutDirs1);
						List<List<Integer>> sourceTreeEdgeIndeces2 = new ArrayList<List<Integer>>();
						List<List<Integer>> sourceTreeEdgeParentIndeces2 = new ArrayList<List<Integer>>();
						List<Integer> sourceTreeHeadNodes2 = new ArrayList<Integer>();
						Set<Integer> extensionNodes2 = new HashSet<Integer>();
						List<Integer> cutEdges2 = new ArrayList<Integer>();
						List<Integer> cutDirs2 = new ArrayList<Integer>();
						baseTree.enumerateExtensionCutFull(headNode2,cut2, subTreeEdgeIndeces2,
								subTreeEdgeParentIndeces2, sourceTreeEdgeIndeces2,sourceTreeEdgeParentIndeces2, 
								sourceTreeHeadNodes2,extensionNodes2, cutEdges2, cutDirs2);
						FeatureMatch extensionMatch = assessExtensionMatch(extensionNodes1,extensionNodes2);

						if(extensionMatch==null)
							continue;
						counter++;
						TreeMatching[][] sourceTreeMatches = new TreeMatching[sourceTreeHeadNodes1.size()][sourceTreeHeadNodes2.size()];
						double[][] sourceTreeScores = new double[sourceTreeHeadNodes1.size()][sourceTreeHeadNodes2.size()];

							
						
						for(int i=0;i<sourceTreeHeadNodes1.size();i++) {
							for(int j=0;j<sourceTreeHeadNodes2.size();j++) {
								TreeMatching m = extensionMatch(sourceTreeHeadNodes1.get(i), sourceTreeHeadNodes2.get(j),
										cutEdges1.get(i), cutEdges2.get(j), cutDirs1.get(i), cutDirs2.get(j),
										sourceTreeEdgeIndeces1.get(i), sourceTreeEdgeIndeces2.get(j), sourceTreeEdgeParentIndeces1.get(i),
										sourceTreeEdgeParentIndeces2.get(j));
								sourceTreeMatches[i][j] = m;
								sourceTreeScores[i][j] = m.sim;
							}						
						}
						int[][] assignment = new int[0][0];	
						boolean transpose = false;
						if  (sourceTreeScores.length > 0 &&  sourceTreeScores[0].length>0) {
							if (sourceTreeScores.length > sourceTreeScores[0].length)
							{	//Cols must be >= Rows.
								sourceTreeScores = HungarianAlgorithm.transpose(sourceTreeScores);
								transpose = true;
							}

							if(sourceTreeScores.length>0 && sourceTreeScores[0].length>0);
								assignment = HungarianAlgorithm.hgAlgorithm(sourceTreeScores, "max");		
						
								if(transpose) {
									sourceTreeScores = HungarianAlgorithm.transpose(sourceTreeScores);
									for(int a=0;a<assignment.length;a++) {
										int[] pair = assignment[a];
										int ele = pair[0];
										pair[0] = pair[1];
										pair[1] = ele;
									}
								}
						}
							

						// find null matches:
						Set<Integer> matchedSourceTrees1 = new HashSet<Integer>();
						Set<Integer> matchedSourceTrees2 = new HashSet<Integer>();
						

						TreeMatching extensionMatching = new TreeMatching();
						extensionMatching.addFeatureMatch(extensionMatch);
						for(int i=0;i<assignment.length;i++) {
							matchedSourceTrees1.add(assignment[i][0]);
							matchedSourceTrees2.add(assignment[i][1]);
							extensionMatching.addMatching(sourceTreeMatches[assignment[i][0]][assignment[i][1]]);
						}

						for(int i=0;i<sourceTreeHeadNodes1.size();i++) {
							if(!matchedSourceTrees1.contains(i)) { // null match
								//System.out.println("null match1");
								FeatureMatch nullMatch = getMatch(sourceTreeHeadNodes1.get(i),sourceTreeEdgeIndeces1.get(i), 
										-1,new ArrayList<Integer>());
								extensionMatching.addFeatureMatch(nullMatch);
							}
						}
						for(int j=0;j<sourceTreeHeadNodes2.size();j++) {
							if(!matchedSourceTrees2.contains(j)) { // null match
								//System.out.println("null match2");
								FeatureMatch nullMatch = getMatch(-1,new ArrayList<Integer>(), 
										sourceTreeHeadNodes2.get(j),sourceTreeEdgeIndeces2.get(j));
								extensionMatching.addFeatureMatch(nullMatch);
							}
						}
						
						extensionMatching.calculate();
						double extensionScore = extensionMatching.sim;
						if(extensionScore>=bestScore) { 
							bestScore = extensionScore;
							bestMatching = extensionMatching;
						}
						
					}

					matching = bestMatching;	


				}
			}
			return matching;
			
		
		}
		
		/**
		 * accept match if:
		 * match is a nullMatch 
		 * least one of the subtrees has a size of less than 3 atoms AND trees are balanced or
		 * at least one of the subtrees contains only one node AND trees are balanced
		 * if the number of nodes and size criterion is fulfilled, but the trees are not balanced:
		 * null-matches are formed!
		 * @return
		 */
		private List<FeatureMatch> assessMatch(Collection<Integer> nodes1,
				Collection<Integer> nodes2) {
			List<FeatureMatch> matches = null;
			double size1 = getSizeOfNodeCollection(nodes1,queryTree);
			double size2 = getSizeOfNodeCollection(nodes2,baseTree);
			boolean balanced = isMatchBalanced(size1,size2);
			//if(nodes1.size()==0 || nodes2.size()==0) {
			//	matches = new ArrayList<FeatureMatch>();
			//	matches.add(getMatch(nodes1, nodes2));
			//}
			
			 if ((size1<MATCH_SIZE_LIMIT || size2<MATCH_SIZE_LIMIT) ||
					(nodes1.size()<MATCH_NODE_NR_LIMIT || nodes2.size()<MATCH_NODE_NR_LIMIT)) {
				if(balanced) {
					matches = new ArrayList<FeatureMatch>();
					matches.add(getMatch(nodes1, nodes2));
				}
				else {
					matches = new ArrayList<FeatureMatch>();
					matches.add(getMatch(nodes1,new ArrayList<Integer>()));
					matches.add(getMatch(new ArrayList<Integer>(),nodes2));
				}
					
			}
					
			return matches;
		}
		
		private FeatureMatch assessExtensionMatch(Collection<Integer> nodes1,
				Collection<Integer> nodes2) {
			FeatureMatch match = null;
			double size1 = getSizeOfNodeCollection(nodes1,queryTree);
			double size2 = getSizeOfNodeCollection(nodes2,baseTree);
			if(nodes1.size()!=0 && nodes2.size()!=0) {
				if ((size1<MATCH_SIZE_LIMIT || size2<MATCH_SIZE_LIMIT) ||
					(nodes1.size()<MATCH_NODE_NR_LIMIT || nodes2.size()<MATCH_NODE_NR_LIMIT)) {
				 		
				 		match = getMatch(nodes1, nodes2);
				
				}
			}
		

			return match;
		}
		
		
		private FeatureMatch getMatch(int headNode1,List<Integer> subTreeEdgeIndeces1, 
				int headNode2,List<Integer> subTreeEdgeIndeces2) {

			FeatureMatch m = null;
			int[][] match = new int[2][];
			if(headNode1==-1) {
				Collection<Integer> nodes2 = baseTree.getNodesFromEdges(subTreeEdgeIndeces2);
				nodes2.add(headNode2);
				match[0] = new int[0];
				match[1] = nodes2.stream().mapToInt(x -> x).toArray();
				m = new FeatureMatch(match);
				m.calculate(queryNodes,baseNodes);
			}
			else if(headNode2==-1) {
				Collection<Integer> nodes1 = queryTree.getNodesFromEdges(subTreeEdgeIndeces1);
				nodes1.add(headNode1);
				match[1] = new int[0];
				match[0] = nodes1.stream().mapToInt(x -> x).toArray();
				m = new FeatureMatch(match);
				m.calculate(queryNodes,baseNodes);
			}
			else {	
				Collection<Integer> nodes1 = queryTree.getNodesFromEdges(subTreeEdgeIndeces1);
				nodes1.add(headNode1);
				Collection<Integer> nodes2 = baseTree.getNodesFromEdges(subTreeEdgeIndeces2);
				nodes2.add(headNode2);
				match[0] = nodes1.stream().mapToInt(x -> x).toArray();
				match[1] = nodes2.stream().mapToInt(x -> x).toArray();
				m = new FeatureMatch(match);
				m.calculate(queryNodes,baseNodes);
			}

			

			return m;
		}
			
		
		private FeatureMatch getMatch(Collection<Integer> nodes1,Collection<Integer> nodes2) {
			FeatureMatch m = null;
			int[][] match = new int[2][];
			match[0] = nodes1.stream().mapToInt(x -> x).toArray();
			match[1] = nodes2.stream().mapToInt(x -> x).toArray();
			m = new FeatureMatch(match);
			m.calculate(queryNodes,baseNodes);
				
			return m;
				
		}
		

		
		private double scoreExtensionMatch(PharmacophoreTree pTree1, PharmacophoreTree pTree2, Set<Integer> extensionNodes1, 
				Set<Integer> extensionNodes2, Set<Integer> sourceNodes1, Set<Integer> sourceNodes2) {
			double extensionScore = PharmacophoreNode.getSimilarity(extensionNodes1, extensionNodes2, queryNodes, baseNodes);
			double sourceScore = PharmacophoreNode.getSimilarity(sourceNodes1, sourceNodes2, queryNodes, baseNodes);

			return ALPHA*extensionScore+(1-ALPHA)*sourceScore;
		}
		
		private boolean isMatchBalanced(double size1,double size2) {
			boolean isBalanced = true;

			double ratio = size1/size2;
			if(ratio > MATCH_BALANCE || ratio < 1.0/MATCH_BALANCE)
				isBalanced = false;
			return isBalanced;

		}
		
		private double getSizeOfNodeCollection(Collection<Integer> nodes, PharmacophoreTree pTree) {
			double size = 0;
			List<PharmacophoreNode> n = pTree.getNodes(nodes);
			for(PharmacophoreNode node : n)
				size += node.getSize();
			
			return size;
			
		}
		
		private double getCutBalance(Collection<Integer> nodes1, Collection<Integer> nodes2) {
			
			double bal = 1.0;
			int abs = Math.abs(nodes1.size()-nodes2.size());
			if(abs>2) 
				bal = 1.0-((abs-2.0)/(nodes1.size()+nodes2.size()-2.0));
			return bal;
				
		}
		
		private double getSplitScore(PharmacophoreTree pTree1, PharmacophoreTree pTree2, Collection<Integer> a1, Collection<Integer> b1,
				Collection<Integer> a2,Collection<Integer> b2) {
			double score = 0.0;
			int na1 = a1.size();
			int na2 = a2.size();
			int nb1 = b1.size();
			int nb2 = b2.size();
			double balance = 0.0;
			if(na1+na2<nb1+nb2)
				balance = getCutBalance(a1,a2);
			else if(na1+na2==nb1+nb2)
				balance = 0.5*getCutBalance(a1,a2)+0.5*getCutBalance(b1,b2);
			else 
				balance = getCutBalance(b1,b2);

			FeatureMatch match1 = getMatch(a1,b1);
			FeatureMatch match2 = getMatch(a2,b2);
			TreeMatching matching = new TreeMatching();
			matching.addFeatureMatch(match1);
			matching.addFeatureMatch(match2);
			matching.calculate();
			score = (1-SIMILARITY_SCALING_SPLIT_SCORE)*(matching.sim)+SIMILARITY_SCALING_SPLIT_SCORE*balance;
			return score;
			
			
		}
		
		
		
		public static class TreeMatching {
			private List<FeatureMatch> matches;
			private double sim;
			private double size1;
			private double size2;
			
			public TreeMatching() {
				matches = new ArrayList<FeatureMatch>();
			}
			
			public void addFeatureMatch(FeatureMatch match) {
				matches.add(match);
			}
			
			public void addMatching(TreeMatching matching) {

				matches.addAll(matching.matches);

			
			}
			
			public void calculate() {
				sim = 0.0;
				size1 = 0.0;
				size2 = 0.0;
				for(FeatureMatch match : matches) {
					sim += match.size*match.sim;
					size1+=match.sizes[0];
					size2+=match.sizes[1];
				}
				// Formula for similarity taken from Langer and Hoffmann: Pharmacophores and Pharmacophore Searches, p. 86
				// replacing the formula from the original publication
				double nom = 0.5*sim;
				double denom = ((NULL_MATCH_SCALING*Math.max(size1, size2)+(1.0-NULL_MATCH_SCALING)*Math.min(size1, size2)));
				sim = nom/denom;
				
			}
			
			public List<FeatureMatch> getMatches() {
				return matches;
			}
			
			public double getSim() {
				return sim; 
			}
			
			public double getSize1() {
				return size1;
			}
			
			public double getSize2() {
				return size2;
			}
			
		}
		/**
		 * TODO: don't add null-matches!
		 * @author joel
		 *
		 */
		public static class FeatureMatch {
			private double sim;
			private double size;
			double[] sizes = new double[2];
			private int[][] match;
	
			public FeatureMatch(int[][] match) {
				this.match = match;
			}
		

			public void calculate(List<PharmacophoreNode> treeNodes1, List<PharmacophoreNode> treeNodes2) {
				sim = 0.0;
				sizes[0] = 0.0;
				sizes[1] = 0.0;
				size = 0.0;
	
				List<PharmacophoreNode> nodes1 = new ArrayList<PharmacophoreNode>();
				List<PharmacophoreNode> nodes2 = new ArrayList<PharmacophoreNode>();
				if(match[0].length!=0) {
					for(Integer node: match[0]) {
						nodes1.add(treeNodes1.get(node));
					}
				}
				if(match[1].length!=0) {
					for(Integer node: match[1]) {
						nodes2.add(treeNodes2.get(node));
					}
				}
				if(nodes1.size()==0 || nodes2.size()==0 ) 
					sim = 0.0;
				
				else 
					sim = PharmacophoreNode.getSimilarity(nodes1, nodes2);
				if(nodes1.size()>0) {
					for(PharmacophoreNode pn : nodes1) {
						sizes[0]+=pn.getSize();
					}
				}
				
				if(nodes2.size()>0) {
					for(PharmacophoreNode pn : nodes2) {
						sizes[1]+=pn.getSize();
					}
				}
				size = sizes[0]+sizes[1];
			}
			
			public double[] getSizes() {
				return sizes;
			}
			
			public double getSim() {
				return sim;
			}
			
			public int[][] getMatch() {
				return match;
			}
			
			public void setSizes(double[] sizes) {
				this.sizes = sizes;
				size = sizes[0] + sizes[1];
			}
			
			public void setSim(double sim) {
				this.sim = sim;
			}
			
		
		
		}
	
		
}
