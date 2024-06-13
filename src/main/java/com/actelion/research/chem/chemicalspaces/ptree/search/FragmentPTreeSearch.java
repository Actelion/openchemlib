package com.actelion.research.chem.chemicalspaces.ptree.search;
import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.chemicalspaces.ptree.PharmTreeSynthonReactionHelper;
import com.actelion.research.chem.chemicalspaces.ptree.synthon.PharmTreeSynthon;
import com.actelion.research.chem.chemicalspaces.ptree.synthon.PharmTreeSynthonLibrary;
import com.actelion.research.chem.chemicalspaces.synthon.SynthonReactor;
import com.actelion.research.chem.descriptor.DescriptorHandlerSkeletonSpheres;
import com.actelion.research.chem.descriptor.pharmacophoretree.DescriptorHandlerPTree;
import com.actelion.research.chem.descriptor.pharmacophoretree.HungarianAlgorithm;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreNode;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreTree;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreTreeGenerator;
import com.actelion.research.chem.descriptor.pharmacophoretree.TreeMatcher;
import com.actelion.research.chem.descriptor.pharmacophoretree.TreeUtils;


import com.actelion.research.chem.descriptor.pharmacophoretree.TreeMatcher.FeatureMatch;
import com.actelion.research.chem.descriptor.pharmacophoretree.TreeMatcher.TreeMatching;


/**
 * Based on: https://doi.org/10.1023/A:1011144622059
 * Given a query molecule on one hand and lists of building blocks that can react according to given rules defined 
 * by a rxnHelper, the FragmentPTreeSearch searches the virtual space defined by the fragments and the reactions for compounds
 * that have similar features compared to the query. The first step is the construction of the potential solution space,
 * consisting of an edge-link table, which stores for every edge and direction in the query tree a list of matched fragments
 * together with the score of the matchings, matched to the subtree of the query tree created by cutting this edge.
 * Compatible fragments from the high-scoring solutions are reacted (product enumeration) and the similarity to the query is
 * calculated. Solutions with similarities higher than a defined cutoff are stored and returned.
 * Two Modes: Scaffold Hopping returns results with high similarity regarding PharmTree, Hit Expansion aims to find
 * hits with high chemical similarity (SkelSpheres) 
 * In order to make the search more efficient, initially only a chemically diverse set of fragments is matched (cluster centroids).
 * Only if the similarity of the cluster centroid is below a given threshold, the remaining members of the cluster are matched.
 *   
 * @author Joel Wahl
 *
 */

public class FragmentPTreeSearch {
	

	
	// synthon id and linker id are merged into one integer, in order to have a fast lookup for compatible synthons
	private static final int LINKER_ID_MASK = 7;
	private static final int SYNTHON_ID_MASK = 56;
	private static final int SYNTHON_SHIFT = 3;


	private static final double SUBTREE_MATCHING_BUFFER = 0.3;
	private static final int N_BEST_SOLUTIONS = 2000;
	
	private PharmacophoreTree queryTree;
	private PharmTreeSynthonLibrary synthonLib;
	private PharmTreeSynthonReactionHelper rxnHelper;
	private SearchResult[][] edgeLinkTable;
	private double pTreeSimilarityCutoff;
	private Map<Integer,Map<PharmTreeSynthon,int[]>> linkerToSynthons; //maps the linkerID+synthonID to the corresponding trees, the index of the linker node and the index of the edge attached to the linker node are stored

	private int[] cuts;
	private int nBestSolutions;

	
	public FragmentPTreeSearch(StereoMolecule queryMol,PharmacophoreTree queryTree, PharmTreeSynthonLibrary synthonLib,
			double pTreeSimilarityCutoff) {
		this.queryTree = queryTree;
		this.rxnHelper = synthonLib.getReactionHelper();
		this.pTreeSimilarityCutoff = pTreeSimilarityCutoff;
		this.synthonLib = synthonLib;
		cuts = new int[] {PharmacophoreTree.CUT_LEFT,PharmacophoreTree.CUT_RIGHT};
		nBestSolutions = N_BEST_SOLUTIONS;
		linkerToSynthons = new HashMap<Integer,Map<PharmTreeSynthon,int[]>>();

	}
	
	private void processFragments() {

		Set<Integer> allLinkers = new HashSet<Integer>();
		int highestLinkerID = 0;
		rxnHelper.getReactantsWithLinkers().values().stream().forEach(e -> allLinkers.addAll(e));
		for(int linkerID : allLinkers)
			if(linkerID>highestLinkerID)
				highestLinkerID = linkerID;
		List<List<PharmTreeSynthon>> allSynthons = synthonLib.getSynthons();
		edgeLinkTable = new SearchResult[highestLinkerID*allSynthons.size()][queryTree.getEdges().size()*2];
		for(int i=0;i<allSynthons.size();i++) {
			List<PharmTreeSynthon> synthons = allSynthons.get(i);
			for(int j=0;j<synthons.size();j++) {
				PharmTreeSynthon synthon = synthons.get(j);
				PharmacophoreTree fragmentTree = synthon.getPharmTree();
				for(int n=0;n<fragmentTree.getNodes().size();n++) {
					PharmacophoreNode node = fragmentTree.getNodes().get(n);
					if(node.isLinkNode()) {
						int linkerID = node.getFunctionalities()[0];
						int linkerEdge = -1;
						for(int e=0;e<fragmentTree.getEdges().size();e++) {
							int[] edge = fragmentTree.getEdges().get(e);
							if(edge[0]==n || edge[1]==n) {
								linkerEdge = e;
								break;
							}
						}
						int id = linkerID + (i<<SYNTHON_SHIFT);
						linkerToSynthons.putIfAbsent(id, new HashMap<PharmTreeSynthon,int[]>());
						Map<PharmTreeSynthon,int[]> map = linkerToSynthons.get(id);
						map.put(synthon,new int[] {n,linkerEdge});
						
					}
				}
				
			}
		}

		
	}
	
	
	public Map<String,Double> search() {
		LinkedHashMap<String,Double> hits = new LinkedHashMap<String,Double>();
		processFragments();
	
		
		for(int i=0;i<queryTree.getEdges().size();i++) {
			List<Integer> querySourceTreeEdgeIndeces = new ArrayList<Integer>();
			List<Integer> queryTargetTreeEdgeIndeces = new ArrayList<Integer>();
			List<Integer> querySourceTreeEdgeParentIndeces = new ArrayList<Integer>();
			List<Integer> queryTargetTreeEdgeParentIndeces = new ArrayList<Integer>();
			int [] headNodes = queryTree.initialCut(cuts[0],i,querySourceTreeEdgeIndeces, querySourceTreeEdgeParentIndeces, 
					queryTargetTreeEdgeIndeces , queryTargetTreeEdgeParentIndeces);
			for(int l : linkerToSynthons.keySet()) {
				int synthonID = (l & SYNTHON_ID_MASK)>>SYNTHON_SHIFT ;
				int linkerID = l & LINKER_ID_MASK;
				for(int cutDirIndex=0;cutDirIndex<cuts.length;cutDirIndex++) {
					SearchResult result = edgeLinkTable[(linkerID-1)*synthonLib.getSynthons().size()+synthonID][2*i+cutDirIndex];
					if(result==null) {
						List<Integer> querySubTreeEdgeIndeces;
						List<Integer> querySubTreeEdgeParentIndeces;
						int queryCutDir; 
						int querySubTreeHeadNode;
						//match to source side of cut query tree
						if(cuts[cutDirIndex] == PharmacophoreTree.CUT_LEFT) {
							querySubTreeEdgeIndeces = querySourceTreeEdgeIndeces;
							querySubTreeEdgeParentIndeces = querySourceTreeEdgeParentIndeces;
							queryCutDir = PharmacophoreTree.CUT_LEFT;
							querySubTreeHeadNode = headNodes[0];
						}

						else {
							//match to target side of cut query tree
							querySubTreeEdgeIndeces = queryTargetTreeEdgeIndeces;
							querySubTreeEdgeParentIndeces = queryTargetTreeEdgeParentIndeces;
							queryCutDir = PharmacophoreTree.CUT_RIGHT;
							querySubTreeHeadNode = headNodes[1];
						}
						result = new SearchResult();
						Map<PharmTreeSynthon,int[]> treeToLinkerHead = linkerToSynthons.get(l);
						for(PharmTreeSynthon synthon: treeToLinkerHead.keySet()) {
							FragmentMatching matching = matchFragmentSubtree(synthon, treeToLinkerHead,synthonID, linkerID, queryCutDir,
									querySubTreeHeadNode, i,  querySubTreeEdgeIndeces,
									querySubTreeEdgeParentIndeces);
							if(matching.sim > (pTreeSimilarityCutoff-SUBTREE_MATCHING_BUFFER)) {
								//similarity of cluster centroid below threshold
								result.addResult(matching);		

								 
							}
						
						}
						if(result.getResults().size()> nBestSolutions) { //prune matchings, store only best solutionsJu
							List<FragmentMatching> prunedMatchings = result.getResults().subList(0, nBestSolutions);
							result.setResult(prunedMatchings);
						}


						edgeLinkTable[(linkerID-1)*synthonLib.getSynthons().size()+synthonID][2*i+cutDirIndex] = result;
					}
					
			}
			
		
		}
		}

		getHits(hits);
		hits = hits.entrySet().stream().sorted(Map.Entry.<String,Double>comparingByValue().reversed()).collect(Collectors.toMap(Map.Entry::getKey,Map.Entry::getValue,
				(e1, e2) -> e1, LinkedHashMap::new));


		return hits;
		
		
		
	}
	
	private FragmentMatching matchFragmentSubtree(PharmTreeSynthon synthon, Map<PharmTreeSynthon,int[]> treeToLinkerHead,
			int synthonID, int linkerID, int queryCutDir, int queryTreeHeadNode, int queryCutEdge, List<Integer> querySubTreeEdgeIndeces,
			List<Integer> querySubTreeEdgeParentIndeces){
		int fragTreeHeadNode = treeToLinkerHead.get(synthon)[0];
		List<Integer> fragmentTreeEdgeIndeces = new ArrayList<Integer>();
		List<Integer> fragmentTreeEdgeParentIndeces = new ArrayList<Integer>();
		PharmacophoreTree fragmentTree = synthon.getPharmTree();
		fragTreeHeadNode = processFragmentTree(fragTreeHeadNode,fragmentTree,fragmentTreeEdgeIndeces,fragmentTreeEdgeParentIndeces);
		FragmentMatchSearch matchSearch;
		//match to source side of cut query tree
		matchSearch = new FragmentMatchSearch(this, queryTree, synthon, synthonID, linkerID, queryTreeHeadNode, fragTreeHeadNode, queryCutEdge,
					treeToLinkerHead.get(synthon)[1],queryCutDir,querySubTreeEdgeIndeces, fragmentTreeEdgeIndeces,
					querySubTreeEdgeParentIndeces, fragmentTreeEdgeParentIndeces);

		return matchSearch.matchSearch();
	}
	
	public void getHits(LinkedHashMap<String,Double> hits) {

		final double buffer = 0.05; 
		for(int i=0;i<edgeLinkTable.length;i++) {
			for(int j=0;j<edgeLinkTable[0].length;j++) {
				SearchResult result = edgeLinkTable[i][j];
				if(result==null)
					continue;

				int edge = j/2;
				int cutDirIndex = j%2; 
				List<List<FragmentMatching>> solutions = constructSolutions(result);
				int linkerID = (i/synthonLib.getSynthons().size()) + 1;
				int synthonID = i%synthonLib.getSynthons().size();
				// find compatible matches 
				int compatibleCutDirIndex = cutDirIndex == 0 ? 1 : 0;

				for(int k=0;k<synthonLib.getSynthons().size();k++) {
					if(k==synthonID)
						continue;
					int id = (linkerID-1)*synthonLib.getSynthons().size()+k;
						//same linker ID, different synthonID
					
					SearchResult compatibleResult = edgeLinkTable[id][2*edge+compatibleCutDirIndex];
					if(compatibleResult==null)
						continue;
					List<List<FragmentMatching>> compatibleSolutions = constructSolutions(compatibleResult);
					for(List<FragmentMatching> solution1 : solutions) {
						double bestScore = 0.0;
						for(List<FragmentMatching> solution2 : compatibleSolutions) {
							List<FragmentMatching> combinedSolution = new ArrayList<FragmentMatching>();	
							combinedSolution.addAll(solution1);
							combinedSolution.addAll(solution2);
							double sim = getTotalSimilarity(combinedSolution);

							if(sim<(pTreeSimilarityCutoff-buffer))
								break;
							else {
								List<StereoMolecule> reactants = combinedSolution.stream().map(r -> r.synthon.getStructure()).collect(Collectors.toList());
								if(reactants.size()!=synthonLib.getSynthons().size())
									continue;
								StereoMolecule product = SynthonReactor.react(reactants);

								boolean accept = true;
	

									if(accept) {	
										StringBuilder resultString = new StringBuilder();
										resultString.append(product.getIDCode());

										resultString.append("____");
										combinedSolution.stream().forEach(r -> {
											resultString.append(r.synthon.getId());
											resultString.append("____");
										});
										resultString.append(synthonLib.getReactionID());

										String rs = resultString.toString();
										if(hits.containsKey(rs)) {
											double oldSim = hits.get(rs);
											if(oldSim<sim)
												hits.put(rs, sim);
										}
										else {
											hits.put(rs, sim);
										}
										if(sim>bestScore) {
											bestScore = sim;
										}

								}
								}
								}

					}
				}
				
			}
		}

	}
	
	
	
	private List<List<FragmentMatching>> constructSolutions(SearchResult result) {
		List<List<FragmentMatching>> allSolutions = new ArrayList<List<FragmentMatching>>();
		List<FragmentMatching> matchings = result.getResults();
		for(FragmentMatching fm : matchings) {
			List<List<FragmentMatching>> solutionSet = new ArrayList<List<FragmentMatching>>();
			List<FragmentMatching> solution = new ArrayList<FragmentMatching>();
			solution.add(fm);
			solutionSet.add(solution);
			if(fm.getFurtherMatches()!=null) {
				Map<Integer,List<FragmentMatching>> furtherMatches = fm.furtherMatchings;
				for(int key : furtherMatches.keySet()) {
					List<List<FragmentMatching>> toBeDeleted = new ArrayList<List<FragmentMatching>>();
					List<List<FragmentMatching>> toBeAdded = new ArrayList<List<FragmentMatching>>();
					for(List<FragmentMatching> oneSolution : solutionSet) {
						toBeDeleted.add(oneSolution);
						for(FragmentMatching ffm : furtherMatches.get(key)) {
							List<FragmentMatching> newSolution = new ArrayList<FragmentMatching>(oneSolution);
							newSolution.add(ffm);
							toBeAdded.add(newSolution);
						}
					}
					solutionSet.removeAll(toBeDeleted);
					solutionSet.addAll(toBeAdded);
						
				}
			}

			allSolutions.addAll(solutionSet);
		}
		return allSolutions;
	}
	

	
	
		
	/**
	 *  
	 * @param headNode
	 * @param treeEdges
	 * @param treeEdgeParents
	 */
	private int processFragmentTree(int linkerNode, PharmacophoreTree fragmentTree, List<Integer> treeEdges,List<Integer> treeEdgeParents) {
		int cutEdge = -1;
		int headNode = -1;
		for(int e=0;e<fragmentTree.getEdges().size();e++) {
			int[] edge = fragmentTree.getEdges().get(e);
			if(edge[0]==linkerNode) {
				cutEdge = e;
				headNode = edge[1];
				break;
			}
			if(edge[1]==linkerNode) {
				cutEdge = e;
				headNode = edge[0];
				break;
			}
		}
		
		fragmentTree.treeWalkBFS(headNode, cutEdge, treeEdges, treeEdgeParents);
		return headNode;
		}
	
	private List<FragmentMatching> matchCompatibleFragments(PharmacophoreTree queryTree, int querySubTreeHeadNode, 
			int querySubtreeCutEdge, int querySubtreeCutEdgeDir,int linkerID, int fragmentID) {
		List<Integer> querySubTreeEdgeIndeces = new ArrayList<Integer>();
		List<Integer> querySubTreeEdgeParentIndeces = new ArrayList<Integer>();
		queryTree.treeWalkBFS(querySubTreeHeadNode, querySubtreeCutEdge, querySubTreeEdgeIndeces, querySubTreeEdgeParentIndeces);
		List<FragmentMatching> allCompatibleMatchings = new ArrayList<FragmentMatching>();
		int cutDirIndex = querySubtreeCutEdgeDir == cuts[0] ? 0: 1;
		for(int i=0;i<synthonLib.getSynthons().size();i++) {
			List<FragmentMatching> matchings = new ArrayList<FragmentMatching>();
			if(fragmentID==i)  //fragments with same ID are not compatible
				continue;
			int index = (i<<SYNTHON_SHIFT) + linkerID;
			Map<PharmTreeSynthon,int[]> compatibleTrees = linkerToSynthons.get(index);
			if(compatibleTrees==null)
				continue;
			SearchResult result = edgeLinkTable[(linkerID-1)*synthonLib.getSynthons().size()+i][2*querySubtreeCutEdge+cutDirIndex];
			if(result!=null) 
				matchings = result.getResults();
			
			else {
				for(PharmTreeSynthon synthon : compatibleTrees.keySet()) {
					PharmacophoreTree fragTree = synthon.getPharmTree();
					List<Integer> fragmentTreeEdgeIndeces = new ArrayList<Integer>();
					List<Integer> fragmentTreeEdgeParentIndeces = new ArrayList<Integer>();
					int[] res = compatibleTrees.get(synthon);
					int fragmentHeadLinkerNode = res[0];
					int fragmentLinkerEdge = res[1];
					int[] linkerEdge = fragTree.getEdges().get(fragmentLinkerEdge);
					int fragmentHeadNode = linkerEdge[0] ==  fragmentHeadLinkerNode ? linkerEdge[1] : linkerEdge [0];
					fragTree.treeWalkBFS(fragmentHeadNode, fragmentLinkerEdge, 
							fragmentTreeEdgeIndeces, fragmentTreeEdgeParentIndeces);
					FragmentMatchSearch fms = new FragmentMatchSearch(this,queryTree, synthon,i, linkerID, querySubTreeHeadNode,fragmentHeadNode, 
								querySubtreeCutEdge, fragmentLinkerEdge,querySubtreeCutEdgeDir,querySubTreeEdgeIndeces,fragmentTreeEdgeIndeces,
								querySubTreeEdgeParentIndeces,fragmentTreeEdgeParentIndeces);
					FragmentMatching matching = fms.matchSearch();
					matching.calculate();
					matchings.add(matching);
				}
			    matchings.sort((e1,e2) -> {
				return Double.compare(e2.sim, e1.sim);}); // reverse order;
				
			    if(matchings.size()> nBestSolutions) 
			    	matchings = matchings.subList(0,nBestSolutions);
			    SearchResult sr = new SearchResult();
			    sr.setResult(matchings);
			    edgeLinkTable[(linkerID-1)*synthonLib.getSynthons().size()+i][2*querySubtreeCutEdge+cutDirIndex] = sr;
			}
			allCompatibleMatchings.addAll(matchings);
			}
		

		return allCompatibleMatchings;
	
	}
	
	private double getTotalSimilarity(List<FragmentMatching> matchings) { //calculate total similarity from a list of matchings
		double sim = 0.0;
		double size1 = 0.0;
		double size2 = 0.0;
		for(FragmentMatching matching : matchings) {
			for(FeatureMatch match : matching.getTreeMatching().getMatches()) {
				double[] sizes = match.getSizes();
				double s = match.getSim();
				sim +=(sizes[0]+sizes[1])*s;
				size1+=sizes[0];
				size2+=sizes[1];
		}
		}
		
		return 0.5*sim/
				((TreeMatcher.NULL_MATCH_SCALING*Math.max(size1, size2)+(1.0-TreeMatcher.NULL_MATCH_SCALING)*Math.min(size1, size2)));
		
	}



public static class FragmentMatchSearch {
	
	
	public static final int EXTENSION_MATCHES = 3; //number of explicitly considered extension matches at every recursion step
	public static final double ALPHA = 0.8; //weighting of source-tree match vs extension-tree match, takes values from 0 to 1
	public static final double NULL_MATCH_SCALING = 0.3; 
	public static final double SIMILARITY_SCALING_SPLIT_SCORE = 0.6;
	public static final double MATCH_BALANCE = 2.0; //named beta in the original publication
	public static final double MATCH_SIZE_LIMIT = 3.0;
	public static final int MATCH_NODE_NR_LIMIT = 2;
	public static final int EXTENSION_MATCH_NODE_NR_LIMIT = 3;
	
	private PharmacophoreTree queryTree;
	private PharmacophoreTree fragmentTree;
	private PharmTreeSynthon synthon;
	private int fragmentTreeSynthonID;
	private int fragmentTreeLinkerID;
	private int queryTreeHeadNode;
	private int cutEdgeQueryTree;
	private int cutEdgeFragmentTree;
	private int cutDirQueryTree;
	private int cutDirFragmentTree;
	private List<Integer> querySubTreeEdgeIndeces;
	private List<Integer> querySubTreeEdgeParentIndeces;
	private int fragTreeHeadNode;
	private List<Integer> fragTreeEdgeIndeces;
	private List<Integer> fragTreeEdgeParentIndeces;
	private TreeMatching[][] dpMatchMatrix;
	private FragmentPTreeSearch pTreeSearch;
	private List<PharmacophoreNode> queryNodes;
	private List<PharmacophoreNode> fragmentNodes;
	
	
	public FragmentMatchSearch(FragmentPTreeSearch pTreeSearch, PharmacophoreTree queryTree, PharmTreeSynthon synthon, int fragmentTreeSynthonID, int fragmentTreeLinkerID, int queryTreeHeadNode, int fragTreeHeadNode, int cutEdgeQueryTree, 
			int cutEdgeFragmentTree, int cutDirQueryTree, List<Integer> querySubTreeEdgeIndeces, List<Integer> fragTreeEdgeIndeces,
			List<Integer> querySubTreeEdgeParentIndeces, List<Integer> fragTreeEdgeParentIndeces) {
		this.queryTree = queryTree;
		this.synthon = synthon;
		this.queryTreeHeadNode = queryTreeHeadNode;
		this.cutEdgeQueryTree = cutEdgeQueryTree;
		this.cutEdgeFragmentTree = cutEdgeFragmentTree;
		this.cutDirQueryTree = cutDirQueryTree;
		this.querySubTreeEdgeIndeces = querySubTreeEdgeIndeces;
		this.querySubTreeEdgeParentIndeces = querySubTreeEdgeParentIndeces;
		this.fragTreeHeadNode = fragTreeHeadNode;
		this.fragTreeEdgeIndeces = fragTreeEdgeIndeces;
		this.fragTreeEdgeParentIndeces = fragTreeEdgeParentIndeces ;
		this.fragmentTreeSynthonID = fragmentTreeSynthonID;
		this.fragmentTreeLinkerID = fragmentTreeLinkerID;
		this.pTreeSearch = pTreeSearch;
		this.queryNodes = queryTree.getNodes();
		this.fragmentTree = synthon.getPharmTree();
		this.fragmentNodes = fragmentTree.getNodes();

		if(fragmentTree.getEdges().get(cutEdgeFragmentTree)[1] == fragTreeHeadNode)
			cutDirFragmentTree = PharmacophoreTree.CUT_RIGHT;
		else 
			cutDirFragmentTree = PharmacophoreTree.CUT_LEFT;
			
		
		dpMatchMatrix = new TreeMatching[2*queryTree.getEdges().size()][2*fragmentTree.getEdges().size()];
	}
	
	public FragmentMatching matchSearch() {
		 return recMatchSearch(queryTreeHeadNode,fragTreeHeadNode,cutEdgeQueryTree,cutEdgeFragmentTree,cutDirQueryTree,cutDirFragmentTree,
				 querySubTreeEdgeIndeces,fragTreeEdgeIndeces,querySubTreeEdgeParentIndeces,fragTreeEdgeParentIndeces);
	}
	
	
	
	
	/**
	 * /recursive part of match-search algorithm, as described in: 
	 * https://doi.org/10.1023/A:1008068904628
	 * The algorithm uses a dynamic-programming approach, whereby the results of matching subtrees are stored in a matrix
	 * and can be reused for increased performance
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
	private FragmentMatching recMatchSearch(int headNode1, int headNode2, int cutEdge1, int cutEdge2, int cutDir1, int cutDir2,List<Integer> subTreeEdgeIndeces1, 
			List<Integer> subTreeEdgeIndeces2,List<Integer> subTreeEdgeParentIndeces1, List<Integer> subTreeEdgeParentIndeces2) {
		TreeMatching treeMatching = new TreeMatching();
		FragmentMatching fragmentMatching = new FragmentMatching(synthon, fragmentTreeSynthonID, fragmentTreeLinkerID );
		int index1 = cutDir1 == PharmacophoreTree.CUT_LEFT ? cutEdge1*2 : cutEdge1*2+1;
		int index2 = cutDir2 == PharmacophoreTree.CUT_LEFT ? cutEdge2*2 : cutEdge2*2+1;
		Set<Integer> nodes1 = queryTree.getNodesFromEdges(subTreeEdgeIndeces1);
		nodes1.add(headNode1);
		Set<Integer> nodes2 = fragmentTree.getNodesFromEdges(subTreeEdgeIndeces2);
		nodes2.add(headNode2);
		if(fragmentTree.getNodes().get(headNode2).isLinkNode() && nodes2.size()==1) { //node is link node -> needs a new call for a match search
			// look for fragments that have compatible links!
			int linkerID = fragmentTree.getNodes().get(headNode2).getFunctionalities()[0];
			List<FragmentMatching> matchings = pTreeSearch.matchCompatibleFragments(queryTree, headNode1, 
					cutEdge1, cutDir1,linkerID,fragmentTreeSynthonID);
			fragmentMatching.setTreeMatching(treeMatching);
			fragmentMatching.addFurtherMatches(linkerID, matchings);
			
		}
		else {
			
			if(dpMatchMatrix[index1][index2]!= null) {
				//result found in dynamic-programing matrix
				treeMatching = dpMatchMatrix[index1][index2];
				fragmentMatching.setTreeMatching(treeMatching);
			}
			else {
				//check if match fulfills criteria, if not, create extension match
				List<FeatureMatch> matches = assessMatch(nodes1,nodes2); 
				if(matches!=null) {
					treeMatching = new TreeMatching();
					for(FeatureMatch fmatch : matches)
						treeMatching.addFeatureMatch(fmatch);
					treeMatching.calculate();
					dpMatchMatrix[index1][index2] = treeMatching;
					fragmentMatching.setTreeMatching(treeMatching);
					} 
				else { // create extension match
					List<int[]> cuts1 = queryTree.getExtensionCuts(subTreeEdgeIndeces1,subTreeEdgeParentIndeces1);
					List<int[]> cuts2 = fragmentTree.getExtensionCuts(subTreeEdgeIndeces2,subTreeEdgeParentIndeces2);
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
							fragmentTree.enumerateExtensionCutFast(headNode2,cut2, subTreeEdgeIndeces2,
									extensionNodes2, sourceNodes2);
							scores[i][j] = scoreExtensionMatch(queryTree,fragmentTree,extensionNodes1, extensionNodes2,
										sourceNodes1, sourceNodes2);
	
						}
					}
					int[][] bestCuts = new int[cuts1.size()*cuts2.size()][2];
					double[] bestScores = new double[cuts1.size()*cuts2.size()];
					TreeUtils.retrieveHighestValuesFrom2DArray(scores, bestScores, bestCuts);

					double bestScore = -Double.MAX_VALUE;
					FragmentMatching bestMatching = null;
					int counter = 0;
					// fully enumerate best extension cuts and find the one with the best match
					for(int[] cut:bestCuts) {
						if(counter>EXTENSION_MATCHES)
							break;
				
						if(cut[0]==-1 || cut[1]==-1)
							continue;
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
						fragmentTree.enumerateExtensionCutFull(headNode2,cut2, subTreeEdgeIndeces2,
								subTreeEdgeParentIndeces2, sourceTreeEdgeIndeces2,sourceTreeEdgeParentIndeces2, 
								sourceTreeHeadNodes2,extensionNodes2, cutEdges2, cutDirs2);
						FeatureMatch extensionMatch = assessExtensionMatch(extensionNodes1,extensionNodes2);
						if(extensionMatch==null)
							continue;
						counter++;
						FragmentMatching[][] sourceTreeMatches = new FragmentMatching[sourceTreeHeadNodes1.size()][sourceTreeHeadNodes2.size()];
						double[][] sourceTreeScores = new double[sourceTreeHeadNodes1.size()][sourceTreeHeadNodes2.size()];
						for(int i=0;i<sourceTreeHeadNodes1.size();i++) {
							for(int j=0;j<sourceTreeHeadNodes2.size();j++) {
								FragmentMatching m = recMatchSearch(sourceTreeHeadNodes1.get(i), sourceTreeHeadNodes2.get(j),
										cutEdges1.get(i), cutEdges2.get(j), cutDirs1.get(i), cutDirs2.get(j),
										sourceTreeEdgeIndeces1.get(i), sourceTreeEdgeIndeces2.get(j), sourceTreeEdgeParentIndeces1.get(i),
										sourceTreeEdgeParentIndeces2.get(j));
								sourceTreeMatches[i][j] = m;
								sourceTreeScores[i][j] = m.sim;						
							}						
						}
						int[][] assignment = new int[0][0];	
						boolean transpose = false;
						//cut subtrees created by the extension cuts of both subtrees should be assigned by a
						// bipartite graph matching, we used the hungarian algorithm here
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
						
	
						FragmentMatching extensionMatching = new FragmentMatching(synthon,fragmentTreeSynthonID,fragmentTreeLinkerID);
						TreeMatching extensionTreeMatching = new TreeMatching();
						extensionMatching.setTreeMatching(extensionTreeMatching);
						extensionTreeMatching.addFeatureMatch(extensionMatch);
						for(int i=0;i<assignment.length;i++) {
							matchedSourceTrees1.add(assignment[i][0]);
							matchedSourceTrees2.add(assignment[i][1]);
							extensionMatching.addFragmentMatch(sourceTreeMatches[assignment[i][0]][assignment[i][1]]);
						}
						for(int i=0;i<sourceTreeHeadNodes1.size();i++) {
							if(!matchedSourceTrees1.contains(i)) { // null match
								FeatureMatch nullMatch = getMatch(sourceTreeHeadNodes1.get(i),sourceTreeEdgeIndeces1.get(i), 
										-1,new ArrayList<Integer>());
								extensionTreeMatching.addFeatureMatch(nullMatch);
							}
						}
						for(int j=0;j<sourceTreeHeadNodes2.size();j++) {
							if(!matchedSourceTrees2.contains(j)) { // null match
								FeatureMatch nullMatch = getMatch(-1,new ArrayList<Integer>(), 
										sourceTreeHeadNodes2.get(j),sourceTreeEdgeIndeces2.get(j));
								extensionTreeMatching.addFeatureMatch(nullMatch);
							}
						}
						extensionTreeMatching.calculate();
						extensionMatching.calculate();
						double extensionScore = extensionMatching.sim;
						if(extensionScore>=bestScore) { 
							bestScore = extensionScore;
							bestMatching = extensionMatching;
						}
						
					}

							
					fragmentMatching = bestMatching;	
	
				}
			}
		}
		fragmentMatching.calculate();

		return fragmentMatching;
		
	
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
	private List<FeatureMatch> assessMatch(Set<Integer> nodes1,
			Set<Integer> nodes2) {

		List<FeatureMatch> matches = null;

		double size1 = getSizeOfNodeSet(nodes1,queryTree);

		double size2 = getSizeOfNodeSet(nodes2,fragmentTree);


		boolean balanced = isMatchBalanced(size1,size2);
		boolean containsLinkNodes = false; //if one of the two collections contain a link nodes, they cannot be matched directly
		for(int n : nodes2) {
			if(fragmentTree.getNodes().get(n).isLinkNode()) {
				containsLinkNodes = true;
				break;
			}
		}
		if(!containsLinkNodes) {
			
			if ((size1<MATCH_SIZE_LIMIT || size2<MATCH_SIZE_LIMIT) ||
					(nodes1.size()<MATCH_NODE_NR_LIMIT || nodes2.size()<MATCH_NODE_NR_LIMIT)) {
				if(balanced) {
					matches = new ArrayList<FeatureMatch>();
					matches.add(getMatch(nodes1, nodes2));
				}
				else {
					matches = new ArrayList<FeatureMatch>();
					matches.add(getMatch(nodes1,new HashSet<Integer>()));
					matches.add(getMatch(new HashSet<Integer>(),nodes2));
				}
					
			}
		}

		return matches;
	}
	

	
	private FeatureMatch assessExtensionMatch( Set<Integer> nodes1,
			Set<Integer> nodes2) {
		FeatureMatch match = null;
		double size1, size2;


		size1 = getSizeOfNodeSet(nodes1,queryTree);

		size2 = getSizeOfNodeSet(nodes2,fragmentTree);


		boolean containsLinkNodes = false; //if one of the two collections contain a link nodes, they cannot be matched directly
		for(int n : nodes2) {
			if(fragmentTree.getNodes().get(n).isLinkNode()) {
				containsLinkNodes = true;
				break;
			}
		}
		if(!containsLinkNodes) {
		if(nodes1.size()!=0 && nodes2.size()!=0) {
			if ((size1<MATCH_SIZE_LIMIT || size2<MATCH_SIZE_LIMIT) ||
				(nodes1.size()<MATCH_NODE_NR_LIMIT || nodes2.size()<MATCH_NODE_NR_LIMIT)) {
			 		
			 		match = getMatch(nodes1, nodes2);
			
			}
		}
		}
		return match;
	}
	
	
	private FeatureMatch getMatch(int headNode1,List<Integer> subTreeEdgeIndeces1, 
			int headNode2,List<Integer> subTreeEdgeIndeces2) {
		
		FeatureMatch m = null;
		int[][] match = new int[2][];
		if(headNode1==-1) {
			Set<Integer> nodes2 = fragmentTree.getNodesFromEdges(subTreeEdgeIndeces2);
			nodes2.add(headNode2);
			match[0] = new int[0];
			match[1] = nodes2.stream().mapToInt(x -> x).toArray();
			m = new FeatureMatch(match);

			m.calculate(queryNodes,fragmentNodes);
		}

		else if(headNode2==-1) {
			Set<Integer> nodes1 = queryTree.getNodesFromEdges(subTreeEdgeIndeces1);
			nodes1.add(headNode1);
			match[1] = new int[0];
			match[0] = nodes1.stream().mapToInt(x -> x).toArray();
			m = new FeatureMatch(match);
			m.calculate(queryNodes,fragmentNodes);

		}
		else {	
			Set<Integer> nodes1 = queryTree.getNodesFromEdges(subTreeEdgeIndeces1);
			nodes1.add(headNode1);
			Set<Integer> nodes2 = fragmentTree.getNodesFromEdges(subTreeEdgeIndeces2);
			nodes2.add(headNode2);
			match[0] = nodes1.stream().mapToInt(x -> x).toArray();
			match[1] = nodes2.stream().mapToInt(x -> x).toArray();
			m = new FeatureMatch(match);
			m.calculate(queryNodes,fragmentNodes);
			
		}
		

		return m;
	}
		
	
	private FeatureMatch getMatch(Set<Integer> nodes1,Set<Integer> nodes2) {

		FeatureMatch m = null;
		int[][] match = new int[2][];
		match[0] = nodes1.stream().mapToInt(x -> x).toArray();
		match[1] = nodes2.stream().mapToInt(x -> x).toArray();
		m = new FeatureMatch(match);
		m.calculate(queryNodes,fragmentNodes);

			
		return m;
			
	}
	

	
	private double scoreExtensionMatch(PharmacophoreTree pTree1, PharmacophoreTree pTree2, Set<Integer> extensionNodes1, 
			Set<Integer> extensionNodes2, Set<Integer> sourceNodes1, Set<Integer> sourceNodes2) {
		double extensionScore = 0.0;
		double sourceScore = 0.0;

		extensionScore = PharmacophoreNode.getSimilarity(extensionNodes1, extensionNodes2, queryNodes, fragmentNodes);
		sourceScore = PharmacophoreNode.getSimilarity(sourceNodes1, sourceNodes2, queryNodes, fragmentNodes);


		return ALPHA*extensionScore+(1-ALPHA)*sourceScore;
	}
	
	private static boolean isMatchBalanced(double size1,double size2) {
		boolean isBalanced = true;

		double ratio = size1/size2;
		if(ratio > MATCH_BALANCE || ratio < 1.0/MATCH_BALANCE)
			isBalanced = false;
		return isBalanced;

	}
	
	public static double getSizeOfNodeSet(Set<Integer> nodes, PharmacophoreTree pTree) {
		double size = 0;
		List<PharmacophoreNode> n = pTree.getNodes(nodes);
		for(PharmacophoreNode node : n)
			size += node.getSize();
		
		return size;
		
	}
	
	public static int[] getFunctionalitiesOfNodeSet(Set<Integer> nodes, PharmacophoreTree pTree) {
		int[] functionalities = new int[PharmacophoreNode.FUNCTIONALITY_WEIGHTS.length];
		List<PharmacophoreNode> n = pTree.getNodes(nodes);
		for(PharmacophoreNode node : n) {
			int[] functionalities2 = node.getFunctionalities();
			for(int i=0;i<functionalities.length;i++)
				functionalities[i]+=functionalities2[i];
		}
		
		return functionalities;
		
	}

}
	
	public static class FragmentMatching {
		private PharmTreeSynthon synthon;
		private int fragmentSynthonID;
		private int fragmentLinkerID;
		private TreeMatching treeMatching;
		private Map<Integer,List<FragmentMatching>> furtherMatchings;
		private double sim;
		private double size1;
		private double size2;
		
		public FragmentMatching(PharmTreeSynthon synthon, int synthonID, int linkerID) {
			this.synthon = synthon;
			this.fragmentSynthonID = synthonID;
			this.fragmentLinkerID = linkerID;
			furtherMatchings = new HashMap<Integer,List<FragmentMatching>>();
		}
		
		public void addFragmentMatch(FragmentMatching fragmentMatching) {
			if(fragmentMatching.treeMatching != null)
				treeMatching.addMatching(fragmentMatching.treeMatching);
			if(fragmentMatching.furtherMatchings!=null) {
				fragmentMatching.furtherMatchings.forEach((key,value) -> furtherMatchings.merge(key, value, (v1,v2) ->  {
					v1.addAll(v2);
					return v1;}));
			}
		}
		
		
		public void setTreeMatching(TreeMatching treeMatching) {
			this.treeMatching = treeMatching;
		}
		
		public TreeMatching getTreeMatching() {
			return treeMatching;
		}
		
		public void addFurtherMatch(int linkerID, FragmentMatching fragmentMatching) {
			furtherMatchings.putIfAbsent(linkerID, new ArrayList<FragmentMatching>());
			furtherMatchings.get(linkerID).add(fragmentMatching);
		}
		
		public void addFurtherMatches(int linkerID, List<FragmentMatching> fragmentMatchings) {
			furtherMatchings.putIfAbsent(linkerID, new ArrayList<FragmentMatching>());
			furtherMatchings.get(linkerID).addAll(fragmentMatchings);
		}
		
		public Map<Integer,List<FragmentMatching>> getFurtherMatches() {
			return furtherMatchings;
		}
		
		public int getFragmentLinkerID() {
			return fragmentLinkerID;
		}
		
		public int getFragmentSynthonID() {
			return fragmentSynthonID;
		}
		
		public PharmTreeSynthon getFragmentPTree() {
			return this.synthon;
		}
		
		public void calculate() {
			if(furtherMatchings.keySet().size()==0) {
				sim = treeMatching.getSim();
				size1 = treeMatching.getSize1();
				size2 = treeMatching.getSize2();
			}
			else {
				sim = 0.0;
				size1 = 0.0;
				size2 = 0.0;
				for(FeatureMatch match : treeMatching.getMatches()) {
					double[] sizes = match.getSizes();
					double s = match.getSim();
					sim +=(sizes[0]+sizes[1])*s;
					size1+=sizes[0];
					size2+=sizes[1];
				}
	
				for(int linkerID : furtherMatchings.keySet()) {
					if(furtherMatchings.get(linkerID).size()==0)
						continue;
					FragmentMatching furtherMatch = furtherMatchings.get(linkerID).get(0); //calculate from highest scored additional match
					for(FeatureMatch match : furtherMatch.treeMatching.getMatches()) {
						double[] sizes = match.getSizes();
						double s = match.getSim();
						sim +=(sizes[0]+sizes[1])*s;
						size1+=sizes[0];
						size2+=sizes[1];
					}
				}
				sim = 0.5*sim/
						((TreeMatcher.NULL_MATCH_SCALING*Math.max(size1, size2)+(1.0-TreeMatcher.NULL_MATCH_SCALING)*Math.min(size1, size2)));
				if(size1==0 && size2==0)
					sim = 0.0;
			}
			
		}
		


		
		
	}
	
	public static class SearchResult {
		
		List<FragmentMatching> results;
		
		private SearchResult() {
			results = new ArrayList<FragmentMatching>();
		}
		
		public List<FragmentMatching> getResults() {
			return results;
		}
		
		public void setResult(List<FragmentMatching> results) {
			this.results = results;
		}
		
		public FragmentMatching getResult(int index) {
			return results.get(index);
		}
		
		public void addResult(List<FragmentMatching> results) {
			this.results.addAll(results);
			this.results.sort((c1,c2) -> {
				return Double.compare(c2.sim,c1.sim);} //reverse order
				);
		}
		
		public void addResult(FragmentMatching result) {
			results.add(result);
			results.sort((c1,c2) -> {
				return Double.compare(c2.sim,c1.sim);} //reverse order
				);
		}
	}
	
		
		
		
	}
	
	


