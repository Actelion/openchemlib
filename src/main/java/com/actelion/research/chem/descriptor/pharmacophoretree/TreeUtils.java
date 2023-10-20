package com.actelion.research.chem.descriptor.pharmacophoretree;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class TreeUtils {
	
	
	/*
	 * insert an element into an sorted array of length N. Only insert if index is smaller than N
	 */
	
	private static void binaryInsert(double[] arr, int[][] indexPairs, double val, int[] indexPair) {
		int n = arr.length;
		int index = Arrays.binarySearch(arr, val);
		if(index<0) 
			index=-index-1;		
		else 
			index = index +1;
		if(index<n) {
			double currentVal = arr[index];
			int[] currentIndexPair = indexPairs[index];
			arr[index] = val;
			indexPairs[index] = indexPair;
			for(int i=index+1;i<n;i++) {
				double previousVal = arr[i];
				int[] previousIndexPair = indexPairs[i];
				arr[i] = currentVal;
				indexPairs[i] = currentIndexPair;
				currentVal = previousVal;
				currentIndexPair = previousIndexPair;
			}
			
	}

		
	}
	
	public static void retrieveHighestValuesFrom2DArray(double[][] arr,double[] val, int[][] indeces) {
		Arrays.fill(val, 1);
		for(int[] indecesRow:indeces)
			Arrays.fill(indecesRow, -1);
		for(int i=0;i<arr.length;i++) {
			for(int j=0;j<arr[0].length;j++) {
				binaryInsert(val,indeces,-arr[i][j],new int[] {i,j});
			}
		}
	}	
	
	public static Map<Integer,List<Integer>> getAdjacencyList(int n, List<int[]> edges) {
		
		Map<Integer,List<Integer>> adjacencyList = new HashMap<Integer,List<Integer>>();
		for(int i=0;i<n;i++)
			adjacencyList.putIfAbsent(i, new ArrayList<Integer>());
		
		for(int[] edge : edges) {
			int n1 = edge[0];
			int n2 = edge[1];
			adjacencyList.get(n1).add(n2);
			adjacencyList.get(n2).add(n1);
		}
		return adjacencyList;
		
	}
	
	public static Map<Integer,Map<Integer,Integer>> getAdjacencyListWithBondOrders(int n, List<PharmacophoreTree.BiGramInt> edges) {
		Map<Integer, Map<Integer,Integer>> adjacencyList = new HashMap<Integer,Map<Integer,Integer>>();
		for(int i=0;i<n;i++)
			adjacencyList.putIfAbsent(i, new HashMap<Integer, Integer>());
		
		for(PharmacophoreTree.BiGramInt edge : edges) {
			int n1 = edge.edge[0];
			int n2 = edge.edge[1];
			adjacencyList.get(n1).put(n2,edge.order);
			adjacencyList.get(n2).put(n1, edge.order);
		}
		return adjacencyList;
	}

	

	
	

}
